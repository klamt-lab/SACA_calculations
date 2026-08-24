from cobrak.standard_solvers import GUROBI
from cobrak.constants import OBJECTIVE_VAR_NAME
from cobrak.lps import perform_lp_optimization
from cobrak.dataclasses import Model
from cobrak.io import json_load

original_model: Model = json_load("iCH360_cobrak.json", Model)
giusi_model: Model = json_load("iCH360_cobrak_giusi.json", Model)

growth_rates = []
for model in (original_model, giusi_model):
    result = perform_lp_optimization(
        cobrak_model=model,
        objective_target="Biomass_fw",
        objective_sense=+1,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=True,
        solver=GUROBI,
    )
    growth_rates.append(result[OBJECTIVE_VAR_NAME])

del giusi_model.reactions["EX_glc__D_e_fw"]
for wt in (True, False):
    # !POX is not in iCH360
    giusi_model.reactions["PDH_fw"].max_flux = 1000.0 if wt else 0.0
    giusi_model.reactions["PFL_fw"].max_flux = 1000.0 if wt else 0.0
    giusi_model.reactions["CS_fw"].max_flux = 1000.0 if wt else 0.0

    print("===WILD-TYPE===" if wt else "===MUTANT===")
    for test_reaction in ("MDHx", "GALS", "ACPSx", "PTAr_bw", "Biomass_fw",):
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
