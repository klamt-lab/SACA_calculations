"""Creates JSONs that describe the main paper calculations.

Use with argument 'local' to store these JSONs in the subfolder "main_paper_calculations". Without this
argument, this script tries to start these runs using SLURM as on HPC clusters.
With a JSON's path as argument, you can then run 'examples/iCH360/H_run_calculations.py' (to be found from
COBRA-k's repository).
"""
from cobrak.dataclasses import Model

from math import log
import contextlib
import os
import subprocess
import tempfile
from sys import argv
from run_config import RunConfig
from cobrak.io import ensure_folder_existence, json_write, json_load

if argv[-1] == "local":
    print("LOCAL MODE")
    print("In the 'main_paper_calculations_jsons' subfolder,")
    print("the .jsons for all main paper")
    print("are created. From the COBRA-k main folder,")
    print("run examples/iCH360/H_run_calculations.py")
    print("with a JSON path to recacalculate,")
    print(
        "e.g. uv run run_calculations.py ./given_folder/XXX.json"
    )

def create_and_submit_slurm_job(json_path: str) -> None:  # noqa: D103
    jobname = "PSBR"
    # Define the SLURM script content
    slurm_script_content = f"""#!/bin/bash

#SBATCH -J {jobname} # Your job name
#SBATCH -e {jobname}%j.err
#SBATCH -o {jobname}%j.out
#SBATCH --time=0-10:45:00 # Maximum expected runtime.
#SBATCH --nodes=1                 # Request 1 full node
#SBATCH --ntasks=1                 # Allocate 1 task
#SBATCH --cpus-per-task=128         # Number CPUs per task
#SBATCH --partition general           # Choose Partition (Queue)
#SBATCH --mail-type=FAIL,END       # An email is sent on begin, end, and failure of the job
#SBATCH --mail-user=bekiaris@mpi-magdeburg.mpg.de # E-Mail for notification

export OMP_NUM_THREADS=1

bash run_cobrak_templates_py_giusi.sh ./run_calculations.py {json_path}
"""

    # Create a temporary file to store the SLURM script
    with tempfile.NamedTemporaryFile(
        encoding="utf-8", mode="w", delete=False, suffix=".slurm", dir=os.getcwd()
    ) as temp_file:
        temp_file.write(slurm_script_content)
        temp_file_path = temp_file.name

    # Make the SLURM script executable
    os.chmod(temp_file_path, 0o755)

    # Submit the SLURM job using sbatch without blocking the Python script
    subprocess.Popen(["sbatch", temp_file_path])

    print(f"SLURM job script created and submitted: {temp_file_path}")


if __name__ == "__main__":
    if argv[-1] != "local":
        # POINT A: Change "pbekiaris" to your HPC username if you run this script
        # on an HPC cluster
        os.chdir("/u/pbekiaris/")

    # km and kcat changes as in original COBRA-k paper
    MANUALLY_CHANGED_KMS = {}
    MANUALLY_CHANGED_KCATS = {
        "NADTRHD_fw": 706056.6952036729 / 10,
        "ME1_fw": 1181225.0119759631 / 100,
    }

    run_configs: list[RunConfig] = []

    carbon_exchanges = (
        "EX_fru_e_bw",
        "EX_xyl__D_e_bw",
        "EX_for_e_bw",
        "EX_acald_e_bw",
        "EX_lac__D_e_bw",
        "EX_glyc__R_e_bw",
        "EX_mal__L_e_bw",
        "EX_etoh_e_bw",
        "EX_glcn_e_bw",
        "EX_gln__L_e_bw",
        "EX_glyc_e_bw",
        "EX_ac_e_bw",
        "EX_akg_e_bw",
        "EX_pyr_e_bw",
        "EX_succ_e_bw",
        "EX_rib__D_e_bw",
        "EX_fum_e_bw",
        "EX_glu__L_e_bw",
        "EX_fru_e_bw",
        "EX_meoh_e",
        "EX_glc__D_e_bw",
    )

    # POINT B: Change the values here to change internal or external concentrations in all calculations
    # Concentrations
    glc_concs = {
        "glc__D_e": (log(30/1000), log(30/1000)),
    }
    meoh_concs = {
        "meoh_e": (log(500/1000), log(500/1000)),
        "meoh_p": (log(1e-6), log(500/1000)),
        "meoh_c": (log(1e-6), log(500/1000)),
    }
    succ_concs = {
        "succ_e": (log(45 / 1000), log(45 / 1000)),
    }
    akg_concs = {
        "akg_e": (log(4/1000), log(4/1000)),
    }
    glyc_concs = {
        "glyc_e": (log(60/1000), log(60/1000)),
    }
    glyc_form_met_concs = {
        "glyc_e": (log(30/1000), log(30/1000)),
        "for_e": (log(5/1000), log(5/1000)),
    } | meoh_concs

    # POINT C: Change the values here to change the reaction knock-out sets
    # Deletion sets
    no_meoh = [
        "EX_meoh_e"
    ]
    c2aux_dels = [
        "PFL_fw", # pflB
        # POX (poxB) is not in model
        "PDH_fw", # aceE
        "CS_fw", # gltA
        # MITS (prpC) is not in model
    ]
    glyoxaux_dels = [
        # GLXCL Δgcl
        "ICL_fw", "ICL_bw", "MALS_fw", # ΔaceBAK (B: ICL; A: MALS, !aceK not in model)
        # ΔglcGB, (G! not in model, B is MALS)
        "GHMT2r_fw", "GHMT2r_bw", # ΔglyA (ALATA_L not in model, ALATA_D2 not in model, THRA2 not in model, THFAT not in model)
        # ΔltaE (4HTHRA not in model, THRA not in model, THRA2 not in model, ALATA_D2 not in model, ALATA_L2 not in model)
        # ΔghrB (2DGULRx not in model, DKGLCNR2y and DKGLCNR2x not in model, 2DGULRY not in model, GLYCLTDx not in model, HPYRRy not in model, HPYRRx not in model, GLYCLTDy not in model, 2DGULRGx and 2DGULRGy not in model)
        # ΔghrA (GLYCLTDx not in model, HPYRRy and HPYRRx not in model, GLYCLTDy not in model)
        # Δkbl (GLYAT not in model)
        # glcDEF is inside :-)
    ]
    glyoxaux_concs = {
        "glyc_e": (log(30 / 1000), log(30/1000)),
        "for_e": (log(5 / 1000), log(5/1000)),
        "meoh_e": (log(500/ 1000), log(500/1000)),
    }

    model: Model = json_load("./iCH360_cobrak_giusi.json" if "local" in argv else "./GIUSI-k/iCH360_cobrak_giusi.json", Model)

    # POINT E: Change the number to change the prefix of newly created folders
    num = 14
    """
    * succinate (45mM) as a second carbon source to replace glucose
    * glycerol (60mM) as a second carbon source to replace glucose.
    * growth on only methanol with the wt
    * wt on 30mM Glycerol + 5mM Formate + 500mM methanol for the second E. coli selection strain in the strain list (E. coli GlyoxAux).
      This strain contains several deletions listed in the excel file and the addition of the genes MeFtfL-MeFch-MeMtdA to produce methylene-THF
      from formate and of the gene bhcA to convert glyoxylate to glycine.
    """
    wtpath = "iCH360_cobrak_giusi.json"
    glyoxauxpath = "iCH360_cobrak_giusi_glyoxaux.json"
    for round_num in (1, 2, 3, 4, 5, 6,):
        for kineticteststate in (
            # POINT F: Outcomment every km- or kcat-change that
            # you *don't* want to test
            "nochanges",
            "kmkcatMDH",
            "kmMDH",
            "kcatMDH",
            "kmkcatGALS",
            "kmGALS",
            "kcatGALS",
            "kmkcatACPS",
            "kmACPS",
            "kcatACPS",
            "kmkcatALL",
            "kmALL",
            "kcatALL"
        ):
            extrakmchanges = {
                "MDHx" if ("ALL" in kineticteststate or "MDH" in kineticteststate) else "": {
                    substrate: model.reactions["MDHx"].enzyme_reaction_data.k_ms[substrate] / 10
                    for substrate in ["meoh_c", "nad_c"]
                },
                "GALS" if ("ALL" in kineticteststate or "GALS" in kineticteststate) else "": {
                    substrate: model.reactions["GALS"].enzyme_reaction_data.k_ms[substrate] / 10
                    for substrate in ["fald_c",]
                },
                "ACPSx" if ("ALL" in kineticteststate or "ACPS" in kineticteststate) else "": {
                    substrate: model.reactions["ACPSx"].enzyme_reaction_data.k_ms[substrate] / 10
                    for substrate in ["gcald_c",]
                }, # not in glyox!
            } if "km" in kineticteststate else {}
            extrakcatchanges = {
                reac_id: model.reactions[reac_id if reac_id else "MDHx"].enzyme_reaction_data.k_cat * 10
                for reac_id in [
                    "MDHx" if ("ALL" in kineticteststate or "MDH" in kineticteststate) else "",
                    "GALS" if ("ALL" in kineticteststate or "GALS" in kineticteststate) else "",
                    "ACPSx" if ("ALL" in kineticteststate or "ACPS" in kineticteststate) else "", # not in glyox!
                ]
            } if "kcat" in kineticteststate else {}

            for (modelpath, scenario_name, kicked_reacs, active_carbon_exchanges, log_min_max_conc_overwrites) in (
                # POINT G:
                # Outcomment all runs that you *don't* want to run; by standard, all
                # calculations shown in the publication are run
                # WT
                # →Methanol
                [wtpath, f"{num}_WT_MEOH", [], ["EX_meoh_e"], meoh_concs,],
                # →Methanol & Glucose
                [wtpath, f"{num}_WT_GLC", [], ["EX_glc__D_e_bw"], glc_concs,],
                [wtpath, f"{num}_WT_GLC_MEOH", [], ["EX_glc__D_e_bw", "EX_meoh_e"], glc_concs | meoh_concs,],
                # →Methanol & Succinate
                ## [wtpath, f"{num}_WT_SUCC", [], ["EX_succ_e_bw"], succ_concs,],
                ## [wtpath, f"{num}_WT_SUCC_MEOH", [], ["EX_succ_e_bw", "EX_meoh_e"], succ_concs | meoh_concs,],
                # →Methanol & Glycerol
                ## [wtpath, f"{num}_WT_GLYC", [], ["EX_glyc_e_bw"], glyc_concs,],
                ## [wtpath, f"{num}_WT_GLYC_MEOH", [], ["EX_glyc_e_bw", "EX_meoh_e"], glyc_concs | meoh_concs,],
                # Mutant C2Aux
                # →Methanol & αkg
                [wtpath, f"{num}_C2AUX_MEOH_AKG", c2aux_dels, ["EX_meoh_e", "EX_akg_e_bw"], meoh_concs | akg_concs,],
                # →Methanol & Glucose & αkg
                [wtpath, f"{num}_C2AUX_GLC_AKG", c2aux_dels, ["EX_glc__D_e_bw", "EX_akg_e_bw"], glc_concs | akg_concs,],
                [wtpath, f"{num}_C2AUX_GLC_MEOH_AKG", c2aux_dels, ["EX_glc__D_e_bw", "EX_meoh_e", "EX_akg_e_bw"], glc_concs | meoh_concs | akg_concs,],
                # →Methanol & Glycerol & αkg
                ## [wtpath, f"{num}_C2AUX_GLYC_AKG", c2aux_dels, ["EX_glyc_e_bw", "EX_akg_e_bw"], glyc_concs | akg_concs,],
                ## [wtpath, f"{num}_C2AUX_GLYC_MEOH_AKG", c2aux_dels, ["EX_glyc_e_bw", "EX_meoh_e", "EX_akg_e_bw"], glyc_concs | meoh_concs | akg_concs,],
                # →Methanol & Succinate & αkg
                ## [wtpath, f"{num}_C2AUX_SUCC_AKG", c2aux_dels, ["EX_succ_e_bw", "EX_akg_e_bw"], succ_concs | akg_concs,],
                ## [wtpath, f"{num}_C2AUX_SUCC_MEOH_AKG", c2aux_dels, ["EX_succ_e_bw", "EX_meoh_e", "EX_akg_e_bw"], succ_concs | meoh_concs | akg_concs,],
                # Mutant GlyoxAux
                [glyoxauxpath, f"{num}_GLYOXAUX_NODEL_GLYC_MEOH", [], ["EX_glyc_e_bw", "EX_meoh_e",], glyoxaux_concs], # No formiate!
                [glyoxauxpath, f"{num}_GLYOXAUX_GLYC_MEOH", glyoxaux_dels, ["EX_glyc_e_bw", "EX_meoh_e",], glyoxaux_concs], # No formiate!
            ):
                if modelpath == glyoxauxpath and "ACPS" in kineticteststate:
                    continue
                if modelpath == glyoxauxpath:
                    if "ACPSx" in extrakmchanges:
                        del extrakmchanges["ACPSx"]
                    if "ACPSx" in extrakcatchanges:
                        del extrakcatchanges["ACPSx"]

                kicked_carbon_exchanges: list[str] = []
                active_carbon_exchange_flux_bounds: dict[str, tuple[float, float]] = {}
                for carbon_exchange in carbon_exchanges:
                    if carbon_exchange not in active_carbon_exchanges:
                        kicked_carbon_exchanges.append(carbon_exchange)
                    else:
                        active_carbon_exchange_flux_bounds[carbon_exchange] = (0.0, 1000.0)
                run_configs.extend(
                    [
                        RunConfig(
                            manually_changed_kms=MANUALLY_CHANGED_KMS | extrakmchanges,
                            manually_changed_kcats=MANUALLY_CHANGED_KCATS | extrakcatchanges,
                            manually_changed_dG0s={},
                            results_folder=f"/{scenario_name}{kineticteststate}/",
                            ectfva_active_reacs=["Biomass_fw"],
                            round_num=round_num,
                            objective_target={
                                "Biomass_fw": 1.0,
                                # "prot_pool_delivery": -0.01,
                            },
                            objective_sense=+1,
                            evolution_num_gens=500,
                            pop_size=128,
                            working_results=[],
                            changed_flux_bounds=active_carbon_exchange_flux_bounds,
                            set_bounds={},
                            protein_pool=0.224,
                            max_conc_sum=0.4,
                            nameaddition=f"giusi{maxglcvalue}",
                            kicked_reacs=kicked_reacs + kicked_carbon_exchanges,
                            log_min_max_conc_overwrites=log_min_max_conc_overwrites,
                            varied_reacs=[],
                            modelpath=modelpath,
                        )
                        for maxglcvalue in (
                            1000,
                        )
                    ]
                )

    if argv[-1] == "local":
        ensure_folder_existence("./main_paper_calculations_jsons")
        for run_config in run_configs:
            json_path = (
                "./main_paper_calculations_jsons/"
                + f"run_config_{run_config.round_num}_{run_config.nameaddition}_{run_config.results_folder.replace('/', '')}.json"
            )
            json_write(json_path, run_config)
    else:
        for run_config in run_configs:
            with contextlib.suppress(FileExistsError):
                os.makedirs("./GIUSI-k" + run_config.results_folder)
            json_path = (
                "./GIUSI-k"
                + run_config.results_folder
                + f"run_config_{run_config.round_num}_{run_config.nameaddition}.json"
            )
            json_write(json_path, run_config)
            create_and_submit_slurm_job(
                json_path="."
                + run_config.results_folder
                + f"run_config_{run_config.round_num}_{run_config.nameaddition}.json",
            )
