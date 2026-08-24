# See https://gitlab.com/elad.noor/glyoxylate-auxotrophy
from cobrak.standard_solvers import GUROBI
from cobrak.constants import OBJECTIVE_VAR_NAME
from cobrak.lps import perform_lp_optimization
from cobrak.dataclasses import Model
from cobrak.io import json_load

giusi_model: Model = json_load("iCH360_cobrak_giusi_glyoxaux.json", Model)

# 30mM Glycerol + 5mM Formate + 500mM methanol

del giusi_model.reactions["EX_glc__D_e_fw"]
for wt in (True, False):
    giusi_model.reactions["ICL_fw"].max_flux = 1000.0 if wt else 0.0
    giusi_model.reactions["ICL_bw"].max_flux = 1000.0 if wt else 0.0
    giusi_model.reactions["MALS_fw"].max_flux = 1000.0 if wt else 0.0
    giusi_model.reactions["GHMT2r_fw"].max_flux = 1000.0 if wt else 0.0
    giusi_model.reactions["GHMT2r_fw"].max_flux = 1000.0 if wt else 0.0

    print("===WILD-TYPE===" if wt else "===MUTANT===")
    for test_reaction in ("Biomass_fw", "GCALDD_ENZ_b1415", "GLYCTO2_ENZ_b2979_and_b4468_and_b4467", "GLYCTO3_ENZ_b2979_and_b4468_and_b4467", "MTHFC_bw", "MTHFD_bw", "bhca_reaction",):
        giusi_model.reactions["EX_glyc_e_bw"].max_flux = 1000.0
        giusi_model.reactions["EX_for_e_bw"].max_flux = 0.0

        giusi_model.reactions["EX_meoh_e"].max_flux = 0.0
        result = perform_lp_optimization(
            cobrak_model=giusi_model,
            objective_target=test_reaction,
            objective_sense=+1,
            with_enzyme_constraints=True,
            with_thermodynamic_constraints=True,
            solver=GUROBI,
        )
        print(f"No MeOH → max({test_reaction}):", round(result[OBJECTIVE_VAR_NAME], 2), "@", round(result["EX_glc__D_e_bw"], 2), "glc", round(result["EX_meoh_e"], 2), "meoh uptake")

        giusi_model.reactions["EX_meoh_e"].max_flux = 1000.0
        result = perform_lp_optimization(
            cobrak_model=giusi_model,
            objective_target=test_reaction,
            objective_sense=+1,
            with_enzyme_constraints=True,
            with_thermodynamic_constraints=True,
            solver=GUROBI,
        )
        print(f"MeOH → max({test_reaction}):", round(result[OBJECTIVE_VAR_NAME], 2), "@", round(result["EX_glc__D_e_bw"], 2), "glc", round(result["EX_meoh_e"], 2), "meoh uptake")
