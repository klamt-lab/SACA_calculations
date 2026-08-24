# Code repository for "Growth-coupled designs expose proteome burden as a key constraint on in vivo implementation of the synthetic acetyl-coenzyme A pathway" ([→Literature](#literature))

## Content

In this repository, one can find the following Python script and further data:

**Folders**: These contain the results from the COBRA-k calculations for this publication. They have the following format: "NUMBER_STRAIN_SUBSTRATES_CHANGES", whereby "NUMBER" is an arbitrary number, strain the wildtype "WT" or C2Aux or GlyoxAux, the substrates are shortened and CHANGES stands for the implemented km and kcat changes, whereby "nochanges" or changes to "ALL" SACA pathway reactions are possible. The final COBRA-k results starts with "final_best_result_" and are in the standard COBRA-k JSON result format with fluxes, kinetic values, concentrations and more. The files starting with "run_config_" show the relevant settings of each run.

**Scripts**: The "*.py" files are the main scripts, whose usage are explained in the installation section below. With "A_create_giusi_model.py", the model iCH360-COBRA-k is enriched with the SACA pathway and made ready for the GlyoxAux and C2Aux strains (the strain models are iCH360_cobrak_giusi_glyoxaux.json and iCH360_cobrak_giusi_c2aux.json, repspectively). "A_run_iML_test.py" runs the MBD calculations in iML1515. The .py files starting with "B_"  test these strain models for stoichiometric capabilities. "D_cobrak_run_definitions.py" (see its usage below) runs the COBRA-k evolutionary algorithm calculations. The scripts starting with "F_" create COBRA-k result overview spreadsheets for the strains.

**CNApy scenarios and models**: The files starting with "FZ_" in the `cnapy_scenarios` folder are [CNApy](https://github.com/cnapy-org/CNApy) scenarios that can be directly visualized in the iCH360-COBRA-k+SACA CNApy model `F_ich360_favoino_et_al.cna`. For iML calculations, the CNApy model `iML_with_SACA.cna` is ready with two `iML_with_SACA_...` CNApy scenarios which include the C2Aux strain deletions or additionally also all knock-outs needed to calculate the MeOH-dependent MBD value by then minimizing the MeOH consumption in CNApy with this scenario.

**Spreadsheets** The spreadsheets created with the "F_..." scripts ca be found as .xlsx files in the main folder. Also, for the determination of the M9 medium concentrations, one can also wfind the spreadsheet `M9_medium.xls`.

## Installation

The code is written in Python and primarily uses the [COBRA-k package](https://github.com/klamt-lab/COBRA-k). To run the scripts, it is recommended (if you haven't already) to install [uv](https://github.com/astral-sh/uv) on your device.

Now, to run all .py scripts starting with "A_", "B_" and "F_" (see above about their content), simply run `uv run ...` with ... as the name of the .py script.

To run the actual COBRA-k evolutionary algorithm calculations, "D_cobrak_run_definitions.py" is the first important script, as it creates JSON files with configurations for each defined run (a run is a calculation, e.g. a single run of iCH360-COBRA-k as wild-type with glucose as substrate).
For settings look at the comments starting with `# POINT` to determine for which kind of runs you wish to create the JSON files. The respective JSON files are stored in a newly created folder called `main_paper_calculations_jsons`.

Now, if you run bare `uv run D_cobrak_run_definitions.py` the JSONs are called by SLURM scripts (scripts used by HPC computer clusters) and then called via `sbatch`, whereby you have to set your HPC username where the comment `# POINT A` lies in the script. To just create the scripts without calling sbatch (e.g. if you work on a normal computer), run `uv run D_cobrak_run_definitions.py local`. In the case you start with `local`, you have to run the subscript `run_calculations.py` as follows: `uv run run_calculations.py JSON_NAME local`, whereby `local` again indicates that you run it on your computer, and the `JSON_NAME` is the path to the run config JSON in `main_paper_calculations_jsons` that you want to calculate.

**Important note regarding external solvers**: You need MILP and NLP solvers either on your device or on the HPC cluster to run `run_calculations.py`. By default, Gurobi is used on local devices with IPOPT MA57 as NLP solver; on HPC clusters, CPLEX is used. Please set the respective variables (which are COBRA-k standard Solver instances such as `IPOPT_MA57` or `CPLEX`) to the solvers that you use in `run_calculations.py`.

## Literature

This repository is part of the following publiation:

* Giusi Favoino, Reka Madarasz, Pavlos Stephanos Bekiaris, Shannara Kayleigh Taylor Parkins, Ava Jahad Mohammadi, Alberto de Maria, Carl Felix Freiesleben, Steffen Klamt, and Pablo I. Nikel, "Growth-coupled designs expose proteome burden as a key constraint on in vivo implementation of the synthetic acetyl-coenzyme A pathway", *in submission*

iML1515's publication is:
* Monk, J. M., Lloyd, C. J., Brunk, E., Mih, N., Sastry, A., King, Z., ... & Palsson, B. O. (2017). i ML1515, a knowledgebase that computes Escherichia coli traits. *Nature biotechnology*, 35(10), 904-908.

iCH360's publication is:

* Corrao, M., He, H., Liebermeister, W., Noor, E., & Bar-Even, A. (2025). A compact model of Escherichia coli core and biosynthetic metabolism. *PLoS computational biology*, 21(10), e1013564.

iCH360-COBRA-k's publication is:

* Bekiaris, P. S., & Klamt, S. (2026). COBRA-k: A powerful framework bridging constraint-based and kinetic metabolic modeling. *Science Advances*, 12(4), eaeb3022.
