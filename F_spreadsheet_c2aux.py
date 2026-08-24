import json
from cobrak.constants import OBJECTIVE_VAR_NAME
from cobrak.io import json_load, save_cobrak_model_as_annotated_sbml_model, get_files, ensure_folder_existence
from cobrak.dataclasses import Model
from cobrak.spreadsheet_functionality import create_cobrak_spreadsheet, VariabilityDataset, OptimizationDataset
from cobrak.utilities import create_cnapy_scenario_out_of_optimization_dict


def get_best_result(path: str) -> dict[str, float]:
    if not path.startswith("./"):
        path = "./" + path
    if not path.endswith("/"):
        path = path + "/"
    final_best_result_filenames = [
        filename for filename in get_files(path=path)
        if filename.startswith("final_best_result_") and filename.endswith(".json") and "time" not in filename
        # if filename.startswith("best_evolution_result_") and filename.endswith(".json.BIN") and "time" not in filename
    ]
    results: list[dict[str, float]] = [
        json_load(f"{path}{final_best_result_filename}")
        for final_best_result_filename in final_best_result_filenames
    ]
    # if "13_C2AUX_GLC_MEOH_AKGkmkcatACPS" in path:
    #     print(final_best_result_filenames)
    #     print(results)
    #     input(".")
    try:
        max_obj = max(result[OBJECTIVE_VAR_NAME] for result in results)
    except ValueError:
        print(f"ERROR in {path}")
        raise ValueError
    min_obj = min(result[OBJECTIVE_VAR_NAME] for result in results)
    print(f"{path.replace('/', '').replace('.', '')}: obj∈[{round(min_obj, 6)};{round(max_obj, 6)}] Δ: {round(max_obj - min_obj, 6)}")
    return [result for result in results if result[OBJECTIVE_VAR_NAME] == max_obj][0]


model: Model = json_load("13_C2AUX_GLC_MEOH_AKGnochanges/used_cobrak_model__1_giusi1000.json", Model)
save_cobrak_model_as_annotated_sbml_model(
    cobrak_model=model,
    filepath="F_ich_giusi.xml",
    combine_base_reactions=True,
    add_enzyme_constraints=False,
)
results: dict[str, dict[str, float]] = {
    "c2aux_glc_meoh": get_best_result("13_C2AUX_GLC_MEOH_AKGnochanges"),
    "c2aux_glc_meoh_GALSkcat": get_best_result("13_C2AUX_GLC_MEOH_AKGkcatGALS"),
    "c2aux_glc_meoh_GALSkm": get_best_result("13_C2AUX_GLC_MEOH_AKGkmGALS"),
    "c2aux_glc_meoh_ACPSkcat": get_best_result("13_C2AUX_GLC_MEOH_AKGkcatACPS"),
    "c2aux_glc_meoh_ACPSkm": get_best_result("13_C2AUX_GLC_MEOH_AKGkmACPS"),
    "c2aux_glc_meoh_MDHkcat": get_best_result("13_C2AUX_GLC_MEOH_AKGkcatMDH"),
    "c2aux_glc_meoh_MDHkm": get_best_result("13_C2AUX_GLC_MEOH_AKGkmMDH"),
    "c2aux_glc_meoh_ALLkmkcat": get_best_result("13_C2AUX_GLC_MEOH_AKGkmkcatALL"),
}
create_cobrak_spreadsheet(
    "F_spreadsheet_c2aux.xlsx",
    model,
    variability_datasets={
            "glc_meoh": VariabilityDataset(json_load(
                "13_C2AUX_GLC_MEOH_AKGnochanges/variability_dict__1_giusi1000.json"
            ),
            True),
        },
    optimization_datasets={
            key: OptimizationDataset(value, with_df=True, with_vplus=True, with_kappa=True, with_gamma=True, with_iota=False, with_alpha=False, with_efficiency_coefficient=True)
            for key, value in results.items()
        },
)

ensure_folder_existence("cnapy_scenarios")
for key, value in results.items():
    create_cnapy_scenario_out_of_optimization_dict(
        f"cnapy_scenarios/FZ_{key}.scen",
        cobrak_model=model,
        optimization_dict=value,
    )
