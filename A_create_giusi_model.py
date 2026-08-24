from os.path import exists
from cobrak.thermokinetic_data_retrieval import get_database_dG0s_for_cobrak_model
from copy import deepcopy
from cobrak.dataclasses import Model, Reaction, EnzymeReactionData, Enzyme, Metabolite
from cobrak.io import json_load, json_write, load_annotated_sbml_model_as_cobrak_model
from math import log
from statistics import median

model: Model = json_load("./iCH360_cobrak.json", Model)
median_kcat = median(reac.enzyme_reaction_data.k_cat for reac in model.reactions.values() if reac.enzyme_reaction_data)

iml = load_annotated_sbml_model_as_cobrak_model("iML1515.xml")
json_write("A_iml_cobrak.json", iml)

# Common manual kcat changes (used for published iCH analyses)
model.reactions["NADTRHD_fw"].enzyme_reaction_data.k_cat = 706056.6952036729 / 10
model.reactions["ME1_fw"].enzyme_reaction_data.k_cat =  1181225.0119759631 / 100


# Add Giusi's SACA pathway :-)
# Second, methanol uptake via diffusion :D
# →meoh_e
model.reactions["EX_meoh_e"] = deepcopy(iml.reactions["MEOHtrpp_ENZ_s0001_FWD"])
model.reactions["EX_meoh_e"].stoichiometries = { "meoh_e" : 1.0 }
model.reactions["EX_meoh_e"].max_flux = 0.0
# meoh_e → meoh_p
model.reactions["MEOHtex"] = deepcopy(iml.reactions["MEOHtex_ENZ_b0929_FWD"])
# meoh_p → meoh_c
model.reactions["MEOHtrpp"] = deepcopy(iml.reactions["MEOHtrpp_ENZ_s0001_FWD"])
model.reactions["MEOHtrpp"].enzyme_reaction_data = None

# ALCD1 (MDH): 2 MeOH + 2 NAD+ → 2 FALD + 2 NADH + 2H+;
# http://bigg.ucsd.edu/models/iAF987/reactions/ALCD1, EC 1.1.1.244
# meoh_c + nad_c → fald_c + h_c + nadh_c
model.reactions["MDHx"] = Reaction( # !"x" as MDH already exists
    name="Alcohol dehydrogenase (methanol)",
    stoichiometries={
        "meoh_c": -1.0, # !All stoichiometries set to 1
        "nad_c": -1.0,
        "fald_c": 1.0,
        "h_c": 1.0,
        "nadh_c": 1.0,
    },
    min_flux=0.0,
    max_flux=1000.0,
    dG0=0.0,
    dG0_uncertainty=0.0,
    enzyme_reaction_data=EnzymeReactionData(
        identifiers=["mdhct41"],
        k_cat=0.20 * 3_600, # !We use the mean :3
        k_ms={
            "meoh_c": 21.6 / 1000,
            "nad_c": 0.093 / 1000, # !Can't find the value in the publication xD, only BRENDA/SABIO-RK entry for NAD⁺ and C. necator; see https://www.brenda-enzymes.org/literature.php?e=1.1.1.244&r=741715 or https://pubmed.ncbi.nlm.nih.gov/26846745
            "fald_c": 3e-4, # !"Filled" value
            "nadh_c": 3e-4, # "Filled" value
        },
    ),
    annotation={
        "ec-code": "1.1.1.244",
        "uniprot": "F8GNE5",
    },
)
model.enzymes["mdhct41"] = Enzyme(molecular_weight=40.689) # !might be a bit different as engineered variant, see https://www.uniprot.org/uniprotkb/F8GNE5/entry
model.metabolites["fald_c"] = deepcopy(iml.metabolites["fald_c"])
model.metabolites["meoh_e"] = deepcopy(iml.metabolites["meoh_e"])
model.metabolites["meoh_p"] = deepcopy(iml.metabolites["meoh_p"])
model.metabolites["meoh_c"] = deepcopy(iml.metabolites["meoh_c"])

# GALS: 2 FALD → 1 GALD;
# No BiGG equivalent!
# 2 fald_c → 1 gcald_c
# kcat = 1.58 s⁻¹, Km = 165 mM, MW = 56.32 kDa
model.reactions["GALS"] = Reaction(
    name="Glycolaldehyde synthase",
    stoichiometries={
        "fald_c": -2.0,
        "gcald_c": 1.0,
    },
    min_flux=0.0,
    max_flux=1000.0,
    dG0=0.0,
    dG0_uncertainty=0.0,
    enzyme_reaction_data=EnzymeReactionData(
        identifiers=["gals_enzmye"],
        k_cat=1.58 * 3600,
        k_ms={
            "fald_c": 165 / 1000,
            "gcald_c": 3e-4, # "Filled" value
        },
    ),
)
model.enzymes["gals_enzmye"] = Enzyme(molecular_weight=56.32)
model.metabolites["gcald_c"] = deepcopy(iml.metabolites["gcald_c"])

# ACPS: 1 GALD + 1 Pi → 1 Acetyl Phosphate;
# No BiGG equivalent! But another ACPS reaction exists, so we add an "x" :-)
# 1 gcald_c + 1 pi_c → actp_c
# kcat = 0.16 s⁻¹, Km = 51 mM, MW = 87.76 kDa
model.reactions["ACPSx"] = Reaction(
    name="Acetyl-phosphate synthase",
    stoichiometries={
        "gcald_c": -1.0,
        "pi_c": -1.0,
        "actp_c": 1.0,
        "h2o_c": 1.0, # !Balanced with H₂O
    },
    min_flux=0.0,
    max_flux=1000.0,
    dG0=0.0,
    dG0_uncertainty=0.0,
    enzyme_reaction_data=EnzymeReactionData(
        identifiers=["acpsx_enzyme"],
        k_cat=0.16 * 3600,
        k_ms={
            "gcald_c": 51 / 1000,
            "pi_c": 3e-4, # !"Filled" value; is ignored anyway :3
            "actp_c": 3e-4, # "Filled" value
        },
    ),
)
model.enzymes["acpsx_enzyme"] = Enzyme(molecular_weight=87.76)

# PTA: 1 Acetyl Phosphate + CoA → 1 Acetyl-CoA
# already exists as PTAr_bw :-)

if not exists("create_giusi_model_dG0s.json"):
    dG0s, dG0_uncertainties = get_database_dG0s_for_cobrak_model(
        cobrak_model=model,
    )
    json_write("create_giusi_model_dG0s.json", dG0s)
    json_write("create_giusi_model_dG0_uncertainties.json", dG0_uncertainties)
else:
    dG0s = json_load("create_giusi_model_dG0s.json", dict[str, float])
model.reactions["MDHx"].dG0 = dG0s["MDHx"]
model.reactions["GALS"].dG0 = dG0s["GALS"]
model.reactions["ACPSx"].dG0 = dG0s["ACPSx"]
model.reactions["MEOHtrpp"].dG0 = dG0s["MEOHtrpp"]
model.reactions["MEOHtex"].dG0 = dG0s["MEOHtex"]

# Set M9 concentrations, see "M9_medium_giusi.xls"
data = """
#na1_e	0,101 !not in model
pi_e	0,07
nh4_e	0,0200
#cl_e	0,021253922
#k_e	0,0200
#mg2_e	0,002
#ca2_e	0,0001
so4_e	0,002
#zn2_e	0,0000062
#cu2_e	0,00000076
#cobalt2_e	0,00000042
#mn2_e	0,000000081
#fe2_e	0,000013
#fe3_e	0,000013
"""
for dataline in data.split("\n"):
    if "\t" not in dataline:
        continue
    if dataline.startswith("#"):
        continue
    met_id = dataline.split("\t")[0]
    conc = float(dataline.split("\t")[1].replace(",", "."))
    model.metabolites[met_id].log_min_conc = log(conc)
    model.metabolites[met_id].log_min_conc = log(conc)
    model.metabolites[met_id.replace("_e", "_c")].log_min_conc = min(model.metabolites[met_id.replace("_e", "_c")].log_min_conc, log(conc / 10))
    model.metabolites[met_id.replace("_e", "_c")].log_min_conc = min(model.metabolites[met_id.replace("_e", "_c")].log_min_conc, log(conc / 10))
    model.metabolites[met_id.replace("_e", "_p")].log_min_conc = min(model.metabolites[met_id.replace("_e", "_p")].log_min_conc, log(conc / 10))
    model.metabolites[met_id.replace("_e", "_p")].log_min_conc = min(model.metabolites[met_id.replace("_e", "_p")].log_min_conc, log(conc / 10))

# GlyoxAux strain
# the addition of the genes MeFtfL-MeFch-MeMtdA to produce methylene-THF
# from formate and of the gene bhcA to convert glyoxylate to glycine.
# →Enzymes required for production of methylene-THF from formate:
# ⋅formate-tetrahydrofolate ligase (ftfL, UniProt Q83WS0),
# ⋅5,10-methenyl-tetrahydrofolate cyclohydrolase (fchA, UniProt Q49135)
# ⋅5,10-methylene-tetrahydrofolate cyclohydrolase (mtdA, UniProt P55818).
# MeFtfL
"""
2. Glyoxylate to Glycine (bhcA)The bhcA gene functions as a transaminase. In BiGG models, this specific activity is typically represented by IDs that reflect the amino donor (usually Alanine or Glutamate).Reaction TypeBiGG Reaction IDBiGG EquationAlanine-Glyoxylate TransaminaseAGATala-L_c + glx_c <=> gly_c + pyr_cGlutamate-Glyoxylate TransaminaseGATHglu-L_c + glx_c <=> gly_c + akg_c
"""
# ⋅formate-tetrahydrofolate ligase (ftfL, UniProt Q83WS0), EC 3.5.4.9
model.reactions["FTHFLi"] = deepcopy(iml.reactions["FTHFLi"])
# MW: # https://www.uniprot.org/uniprotkb/Q83WS0/entry
model.enzymes["FTHFLi_enzyme"] = Enzyme(molecular_weight=59.391, name="Q83WS0")
# kcat (!other strain): https://www.sciencedirect.com/science/article/pii/S0006291X20311451?via%3Dihub#sec3
model.reactions["FTHFLi"].enzyme_reaction_data = EnzymeReactionData(
    k_cat=52.17 * 3_600, # !lowest given kcat for substrates
    identifiers=["FTHFLi_enzyme"],
)
# kms: https://pmc.ncbi.nlm.nih.gov/articles/PMC296244/
model.reactions["FTHFLi"].enzyme_reaction_data.k_ms["atp_c"] = 21 / 1_000
model.reactions["FTHFLi"].enzyme_reaction_data.k_ms["for_c"] = 22 /1_000
model.reactions["FTHFLi"].enzyme_reaction_data.k_ms["thf_c"] = 0.8 / 1_000
model.reactions["FTHFLi"].enzyme_reaction_data.k_ms["10fthf_c"] = 3e-4
model.reactions["FTHFLi"].enzyme_reaction_data.k_ms["adp_c"] = 3e-4
model.reactions["FTHFLi"].enzyme_reaction_data.k_ms["adp_c"] = 3e-4
model.reactions["FTHFLi"].enzyme_reaction_data.k_ms["pi_c"] = 3e-4

# ⋅5,10-methenyl-tetrahydrofolate cyclohydrolase (fchA, UniProt Q49135)
# (6R)-5,10-methenyltetrahydrofolate + H2O = (6R)-10-formyltetrahydrofolate + H+
# MTHFC already exists in model!?
model.reactions["MTHFC_bw"].enzyme_reaction_data.identifiers = ["MTHFCenzyme"]
# MW: https://www.uniprot.org/uniprotkb/Q49135/entry
model.enzymes["MTHFCenzyme"] = Enzyme(molecular_weight=21.723, name="Q49135")
# kms, kcats not measured for this enzyme in the organism, we use what was found in databases!

# ⋅5,10-methylene-tetrahydrofolate cyclohydrolase (mtdA, UniProt P55818).
# MTHFD already exists in model!?
model.reactions["MTHFD_bw"].enzyme_reaction_data.identifiers = ["MTHFDenzyme"]
# MW: https://www.uniprot.org/uniprotkb/P55818/entry
model.enzymes["MTHFDenzyme"] = Enzyme(molecular_weight=29.736, name="P55818")
# kms and kcat: https://pmc.ncbi.nlm.nih.gov/articles/PMC107583/ for the wrong direction!
# model.reactions["MTHFD_bw"].enzyme_reaction_data.k_ms["methf_c"] = 20 / 1_000_000
# model.reactions["MTHFD_bw"].enzyme_reaction_data.k_ms["nadph_c"] = 22 /1_000
model.reactions["MTHFD_bw"].enzyme_reaction_data.k_ms["nadp_c"] = 3e-4
model.reactions["MTHFD_bw"].enzyme_reaction_data.k_ms["mlthf_c"] = 3e-4

"""
2. Glyoxylate to Glycine (bhcA)The bhcA gene functions as a transaminase. In BiGG models, this specific activity is typically represented by IDs that reflect the amino donor (usually Alanine or Glutamate).Reaction TypeBiGG Reaction IDBiGG EquationAlanine-Glyoxylate TransaminaseAGATala-L_c + glx_c <=> gly_c + pyr_cGlutamate-Glyoxylate TransaminaseGATHglu-L_c + glx_c <=> gly_c + akg_c
oxaloacetate + glycine = glyoxylate + L-aspartate
"""
# i) Formaldehyde → Glycoaldehyde: Done with GALS
# 2 fald_c → gcald_c

# ii) Glycoaldehyde → Glycolate: With aldA (GCALDD_ENZ_b1415 in iML)
# gcalc_c + h2o_c + nad_c → glyclt_c + 2 h_c + nadh_c
model.reactions["GCALDD_ENZ_b1415"] = deepcopy(iml.reactions["GCALDD_ENZ_b1415"])
model.reactions["GCALDD_ENZ_b1415"].enzyme_reaction_data.k_cat=18.3 * 3_600 # https://www.brenda-enzymes.org/literature.php?e=1.2.1.22&r=670833
model.reactions["GCALDD_ENZ_b1415"].enzyme_reaction_data.k_ms = {
        "gcalc_c": 0.14 / 1_000, # https://www.brenda-enzymes.org/literature.php?e=1.2.1.22&r=670833
        "nad_c": 0.047 / 1_000, # https://www.brenda-enzymes.org/literature.php?e=1.2.1.22&r=695015
        "glyclt_c": 3e-4, # Standard product km value
        "nadh_c": 3e-4, # Standard product km value
}
model.enzymes["b1415"] = deepcopy(iml.enzymes["b1415"])
model.enzymes["b1415"].molecular_weight = 52.273 # https://www.uniprot.org/uniprotkb/P25553/entry
model.metabolites["glyclt_c"] = deepcopy(iml.metabolites["glyclt_c"])

# iii) Glycolate → Glyoxylate: With glcDEF (GLYCTO2_ENZ_b2979_and_b4468_and_b4467, GLYCTO3_ENZ_b2979_and_b4468_and_b4467, GLYCTO4_ENZ_b2979_and_b4468_and_b4467)
# glyclt_c + q8_c (or mqn8_c) → glx_c + q8h2_c (or mql8_c)
model.reactions["GLYCTO2_ENZ_b2979_and_b4468_and_b4467"] = deepcopy(iml.reactions["GLYCTO2_ENZ_b2979_and_b4468_and_b4467"])
model.reactions["GLYCTO2_ENZ_b2979_and_b4468_and_b4467"].enzyme_reaction_data.k_cat = 1.93 * 3_600 # https://www.brenda-enzymes.org/literature.php?e=1.1.99.14&r=741190
model.reactions["GLYCTO2_ENZ_b2979_and_b4468_and_b4467"].enzyme_reaction_data.k_ms = {
    "glyclt_c": 0.04 / 1_000, # https://www.brenda-enzymes.org/literature.php?e=1.1.99.14&r=348285
    "q8_c": 0.00016, # Standard value
    "glx_c": 0.0003,
    "q8h2_c": 0.0003,
}
model.reactions["GLYCTO3_ENZ_b2979_and_b4468_and_b4467"] = deepcopy(iml.reactions["GLYCTO3_ENZ_b2979_and_b4468_and_b4467"])
model.reactions["GLYCTO3_ENZ_b2979_and_b4468_and_b4467"].enzyme_reaction_data.k_cat = 1.93 * 3_600 # https://www.brenda-enzymes.org/literature.php?e=1.1.99.14&r=741190
model.reactions["GLYCTO3_ENZ_b2979_and_b4468_and_b4467"].enzyme_reaction_data.k_ms = {
    "glyclt_c": 0.04 / 1_000, # https://www.brenda-enzymes.org/literature.php?e=1.1.99.14&r=348285
    "q8_c": 0.00016, # Standard value
    "glx_c": 0.0003,
    "q8h2_c": 0.0003,
}
# model.reactions["GLYCTO4_ENZ_b2979_and_b4468_and_b4467"] = deepcopy(iml.reactions["GLYCTO4_ENZ_b2979_and_b4468_and_b4467"]), not included as 2dmmq8_c is missing in model

model.enzymes["b2979"] = deepcopy(iml.enzymes["b2979"])
model.enzymes["b2979"].molecular_weight = 53.812 # https://www.uniprot.org/uniprotkb/P0AEP9/entry
model.enzymes["b4468"] = deepcopy(iml.enzymes["b4468"])
model.enzymes["b4468"].molecular_weight = 38.361 # https://www.uniprot.org/uniprotkb/P52073/entry
model.enzymes["b4467"] = deepcopy(iml.enzymes["b4467"])
model.enzymes["b4467"].molecular_weight = 45.11 # https://www.uniprot.org/uniprotkb/P52073/entry

# See https://www.nature.com/articles/s41467-025-57407-3
model.reactions["GLXCL_ENZ_b0507"] = deepcopy(iml.reactions["GLXCL_ENZ_b0507"])
model.enzymes["b0507"] = deepcopy(iml.enzymes["b0507"])
model.enzymes["b0507"].molecular_weight = 64.732
model.reactions["GLXCL_ENZ_b0507"].enzyme_reaction_data = EnzymeReactionData(
    identifiers=["b0507"],
    k_cat=18.9 * 3_600, # https://www.brenda-enzymes.org/literature.php?e=4.1.1.47&r=705846
    k_ms={
        "glx_c": 0.9 / 1_000, # https://www.brenda-enzymes.org/literature.php?e=4.1.1.47&r=705846
        # "h_c" is ignored anyway
        "2h3oppan_c": 3e-4,
        "co2_c": 3e-4,
    },
)


model.reactions["TRSARr_ENZ_b0509_FWD"] = deepcopy(iml.reactions["TRSARr_ENZ_b0509_FWD"])
model.reactions["TRSARr_ENZ_b0509_FWD"].enzyme_reaction_data = EnzymeReactionData(
    k_cat=23 * 3_600, # SABIO 67732
    k_ms={
        "2h3oppan_c": 0.00019 / 1_000, # https://www.brenda-enzymes.org/literature.php?e=1.1.1.60&r=726270
        # "h_c" ignored
        "nadh_c": 0.0143 / 1_000, # https://www.brenda-enzymes.org/literature.php?e=1.1.1.60&r=726270
        "glyc__R_c": 3e-4,
        "nad_c": 3e-4,
    },
    identifiers=["b0509"],
)

model.reactions["TRSARr_ENZ_b0509_REV"] = deepcopy(iml.reactions["TRSARr_ENZ_b0509_REV"])
model.reactions["TRSARr_ENZ_b0509_REV"].enzyme_reaction_data = EnzymeReactionData(
    k_cat=median_kcat, # No value given! Using median
    k_ms={
        "glyc__R_c": 0.28 / 1_000, # https://www.brenda-enzymes.org/literature.php?e=1.1.1.60&r=287167
        "nad_c": 0.04 / 1_000, # https://www.brenda-enzymes.org/literature.php?e=1.1.1.60&r=287164
        "2h3oppan_c": 3e-4,
        # "h_c" is ignored
        "nadh_c": 3e-4,
    },
    identifiers=["b0509"],
)

model.enzymes["b0509"] = deepcopy(iml.enzymes["b0509"])
model.enzymes["b0509"].molecular_weight = 30.801
model.metabolites["2h3oppan_c"] = deepcopy(iml.metabolites["2h3oppan_c"])

dG0s, dG0_uncertainties = get_database_dG0s_for_cobrak_model(
    cobrak_model=model,
    inclusion_prefixes=["TRSARr_ENZ_b0509_FWD", "TRSARr_ENZ_b0509_REV", "GLXCL_ENZ_b0507", "GCALDD_ENZ_b1415", "GLYCTO2_ENZ_b2979_and_b4468_and_b4467", "GLYCTO3_ENZ_b2979_and_b4468_and_b4467"],
)
model.reactions["TRSARr_ENZ_b0509_FWD"].dG0 = dG0s["TRSARr_ENZ_b0509_FWD"]
model.reactions["TRSARr_ENZ_b0509_REV"].dG0 = dG0s["TRSARr_ENZ_b0509_REV"]
model.reactions["GLXCL_ENZ_b0507"].dG0 = dG0s["GLXCL_ENZ_b0507"]
model.reactions["GCALDD_ENZ_b1415"].dG0 = dG0s["GCALDD_ENZ_b1415"]
model.reactions["GLYCTO2_ENZ_b2979_and_b4468_and_b4467"].dG0 = dG0s["GLYCTO2_ENZ_b2979_and_b4468_and_b4467"]
model.reactions["GLYCTO3_ENZ_b2979_and_b4468_and_b4467"].dG0 = dG0s["GLYCTO3_ENZ_b2979_and_b4468_and_b4467"]


json_write("iCH360_cobrak_giusi.json", model)

del model.reactions["ACPSx"]

# iv) Glyoxylate → Glycine
model.reactions["bhca_reaction"] = Reaction(
    stoichiometries={
        "glx_c": -1.0,
        "asp__L_c": -1.0,
        "oaa_c": 1.0,
        "gly_c": 1.0,
    },
    min_flux=0.0,
    max_flux=1000.0,
    enzyme_reaction_data=EnzymeReactionData(
        k_cat=median_kcat, # !no real value existing :O Using median...
        k_ms={ # !no real values for this direction D:
            "oaa_c": 0.00016,
            "gly_c": 0.00016,
            "glx_c": 0.0003,
            "asp__L_c": 0.0003,
        },
        identifiers=["bhca_enzyme"],
    ),
)
dG0s, dG0_uncertainties = get_database_dG0s_for_cobrak_model(
    cobrak_model=model,
    inclusion_prefixes=["bhca_reaction", ],
)
model.reactions["bhca_reaction"].dG0 = dG0s["bhca_reaction"]
model.enzymes["bhca_enzyme"] = Enzyme(molecular_weight=42.507) # from here: https://www.uniprot.org/uniprotkb/A1B8Z3/entry

json_write("iCH360_cobrak_giusi_glyoxaux.json", model)
