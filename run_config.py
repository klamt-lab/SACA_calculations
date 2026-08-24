from pydantic import Field
from dataclasses import dataclass


@dataclass
class RunConfig:  # noqa: D101
    # Model changes
    manually_changed_kms: dict[str, dict[str, float]]
    manually_changed_kcats: dict[str, float]
    manually_changed_dG0s: dict[str, float]
    # Folder settings
    results_folder: str
    # ecTFVA settings
    ectfva_active_reacs: list[str]
    # Evolutionary algorithm settings
    round_num: int
    objective_target: str | dict[str, float]
    objective_sense: int
    set_bounds: dict[str, tuple[float, float]]
    working_results: list[dict[str, float]]
    changed_flux_bounds: dict[str, tuple[float, float]] = Field(default_factory=list)
    sampling_rounds_per_metaround: int = 2
    sampling_wished_num_feasible_starts: int = 5
    sampling_max_metarounds: int = 3
    evolution_num_gens: int = 150
    pop_size: int = 32
    protein_pool: float | None = None
    log_min_max_conc_overwrites: dict[str, tuple[float, float]] = Field(default_factory=dict)
    max_conc_sum: float | None = None
    nameaddition: str | None = None
    kicked_reacs: list[str] = Field(default_factory=list)
    do_parameter_variation: bool = False
    varied_reacs: list[str] = Field(default_factory=list)
    max_km_variation: float | None = None
    max_kcat_variation: float | None = None
    max_ki_variation: float | None = None
    max_ka_variation: float | None = None
    max_dG0_variation: float | None = None
    with_kappa: bool = True
    with_gamma: bool = True
    with_iota: bool = False
    with_alpha: bool = False
    change_known_values: bool = True
    change_unknown_values: bool = True
    use_shuffling_instead_of_uniform_random: bool = False
    use_shuffling_with_putting_back: bool = False
    free_upper_unfixed_concentrations: bool = False
    json_path_model_to_merge: str = ""
    shuffle_using_distribution_of_values_with_reference: bool = True
    modelpath: str = ""
