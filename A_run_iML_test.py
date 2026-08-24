from cobrak.standard_solvers import GUROBI
from cobrak.utilities import delete_orphaned_metabolites_and_enzymes, create_cnapy_scenario_out_of_optimization_dict
from cobrak.constants import OBJECTIVE_VAR_NAME
from cobrak.lps import perform_lp_optimization
from cobrak.dataclasses import Model, Reaction, EnzymeReactionData, Enzyme, ExtraLinearConstraint, Metabolite
from cobrak.io import json_write, load_annotated_sbml_model_as_cobrak_model, save_cobrak_model_as_annotated_sbml_model

model: Model = load_annotated_sbml_model_as_cobrak_model("iML1515.xml")
json_write("A_iml_cobrak.json", model)

# Add Giusi's SACA pathway :-)

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


result = perform_lp_optimization(
    model,
    objective_target="BIOMASS_Ec_iML1515_core_75p37M",
    objective_sense=+1,
    solver=GUROBI,
)
print("µ under Glc [1/h]:", result[OBJECTIVE_VAR_NAME])

old_glc_flux = model.reactions["EX_glc__D_e_REV"].max_flux
model.reactions["EX_glc__D_e_REV"].max_flux = 0.0
model.reactions["EX_meoh_e"].max_flux = 15.0
model.reactions["EX_meoh_e"].stoichiometries = {
    "meoh_e": +1.0,
}

result = perform_lp_optimization(
    model,
    objective_target="BIOMASS_Ec_iML1515_core_75p37M",
    objective_sense=+1,
    solver=GUROBI,
)
print("µ under MeOH [1/h]:", result[OBJECTIVE_VAR_NAME])
json_write("iML_cobrak_giusi.json", model)

model.reactions["EX_glc__D_e_REV"].max_flux = old_glc_flux
model.reactions["EX_glc__D_e_REV"].max_flux = 0.0

save_cobrak_model_as_annotated_sbml_model(
    cobrak_model=model,
    filepath="iML1515_plus_SACA.xml",
    combine_base_reactions=True,
)

# model.reactions["EX_glc__D_e_REV"].max_flux = 0.0
model.reactions["ATPM"].min_flux = 0.0
for reac_id, reaction in model.reactions.items():
    if reaction.stoichiometries.get("accoa_c", 0.0) > 0.0:
        reaction.stoichiometries["accoaX_c"] = reaction.stoichiometries["accoa_c"]
model.metabolites["accoaX_c"] = Metabolite()
model.reactions["EX_accoaX"] = Reaction(
    name="Acetyl-phosphate synthase",
    stoichiometries={
        "accoaX_c": -1.0,
    },
    min_flux=0.0,
    max_flux=1000.0,
    dG0=0.0,
    dG0_uncertainty=0.0,
    enzyme_reaction_data=None,
)
model.reactions["EX_akg_e"].stoichiometries = {
    "akg_e": +1.0
}
c2aux_dels = [
    "PFL_ENZ_b0902_and_b0903", # pflB
    "PFL_ENZ_b0902_and_b0903_and_b2579", # pflB
    "OBTFL_ENZ_b0902_and_b0903_and_b2579", # pflB
    "OBTFL_ENZ_b0902_and_b0903", # pflB
    "POX_ENZ_b0871", # poxB
    "PDH_ENZ_b0115_and_b0116_and_b0114", # aceE
    "CS_ENZ_b0720", # gltA
    "MCITS_ENZ_b0333", # prpC
]
for c2aux_del in c2aux_dels:
    del model.reactions[c2aux_del]
model = delete_orphaned_metabolites_and_enzymes(model)
model.reactions["EX_glc__D_e_REV"].max_flux = 1000.0
model.extra_linear_constraints = [
    ExtraLinearConstraint(
        stoichiometries={"BIOMASS_Ec_iML1515_core_75p37M": 1.0},
        lower_value=1.0,
    )
]
result = perform_lp_optimization(
    model,
    objective_target="EX_accoaX",
    objective_sense=-1,
    solver=GUROBI,
)
print("AcCoA usage (AcCoA MBD) in mmol/(gDW⋅h):", result[OBJECTIVE_VAR_NAME], "@ µ of", result["BIOMASS_Ec_iML1515_core_75p37M"])

model.reactions["ACS_ENZ_b4069"].max_flux = 0.0
model.reactions["PAI2T_ENZ_b1517"].max_flux = 0.0
model.reactions["POR5_ENZ_b1378_and_b0684_FWD"].max_flux = 0.0
model.reactions["ME1_ENZ_b1479"].max_flux = 0.0
model.reactions["ME2_ENZ_b2463"].max_flux = 0.0
model.reactions["PPCK_ENZ_b3403"].max_flux = 0.0
model.reactions["GLYK_ENZ_b3926"].max_flux = 0.0
model.reactions["POR5_ENZ_b1378_and_b2895_FWD"].max_flux = 0.0

# alle außer PDH
result = perform_lp_optimization(
    model,
    objective_target="EX_meoh_e",
    objective_sense=-1,
    solver=GUROBI,
    with_loop_constraints=True,
)
print("AcCoA producing reactions of last result:")
for reac_id, reaction in model.reactions.items():
    if result.get(reac_id, 0.0) <= 0.0:
        continue
    if reaction.stoichiometries.get("accoa_c", 0.0) <= 0.0:
        continue
    print("*", reac_id, result[reac_id])
