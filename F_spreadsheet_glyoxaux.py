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
    ]
    results: list[dict[str, float]] = [
        json_load(f"{path}{final_best_result_filename}")
        for final_best_result_filename in final_best_result_filenames
    ]
    max_obj = max(result[OBJECTIVE_VAR_NAME] for result in results)
    min_obj = min(result[OBJECTIVE_VAR_NAME] for result in results)
    print(f"{path.replace('/', '').replace('.', '')}: obj∈[{round(min_obj, 6)};{round(max_obj, 6)}] Δ: {round(max_obj - min_obj, 6)}")
    return [result for result in results if result[OBJECTIVE_VAR_NAME] == max_obj][0]


model: Model = json_load("13_GLYOXAUX_GLYC_MEOHnochanges/used_cobrak_model__1_giusi1000.json", Model)
save_cobrak_model_as_annotated_sbml_model(
    cobrak_model=model,
    filepath="F_ich_giusi.xml",
    combine_base_reactions=True,
    add_enzyme_constraints=False,
)
results: dict[str, dict[str, float]] = {
    "glyoxaux_nodel_glyc_meoh": get_best_result("13_GLYOXAUX_NODEL_GLYC_MEOHnochanges"),
    "glyoxaux_glyc_meoh": get_best_result("13_GLYOXAUX_GLYC_MEOHnochanges"),
    "glyoxaux_glyc_meohGALSkcat": get_best_result("13_GLYOXAUX_GLYC_MEOHkcatGALS"),
    "glyoxaux_glyc_meohGALSkm": get_best_result("13_GLYOXAUX_GLYC_MEOHkmGALS"),
    # "glyoxaux_glyc_meohACPSkcat": get_best_result("13_GLYOXAUX_GLYC_MEOHkcatACPS"),
    # "glyoxaux_glyc_meohACPSkm": get_best_result("13_GLYOXAUX_GLYC_MEOHkmACPS"),
    "glyoxaux_glyc_meohMDHkcat": get_best_result("13_GLYOXAUX_GLYC_MEOHkcatMDH"),
    "glyoxaux_glyc_meohMDHkm": get_best_result("13_GLYOXAUX_GLYC_MEOHkmMDH"),
    "glyoxaux_glyc_meohALLkmkcat": get_best_result("13_GLYOXAUX_GLYC_MEOHkmkcatALL"),
}
create_cobrak_spreadsheet(
    "F_spreadsheet_glyoxaux.xlsx",
    model,
    variability_datasets={
            "ecTFVA_wt_glc_meoh": VariabilityDataset(json_load(
                "13_GLYOXAUX_GLYC_MEOHnochanges/variability_dict__1_giusi1000.json"
            ), True),
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
