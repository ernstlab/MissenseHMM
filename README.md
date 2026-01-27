# MissenseHMM
MissenseHMM: state-based annotations for missense variants through joint modeling of pathogenicity scores.
The pre-generated annotations can be found in the `annotations_v1.0` folder. 

### Requirements
First, download dbNSFP data from https://www.dbnsfp.org/. Modify the `IN_LOC` and `VERSON` variables in `config.yaml` accordingly.

Next, download the ChromHMM software from https://compbio.mit.edu/ChromHMM/. Modify the `CHROMHMM_LOC` variable in `config.yaml` to point to the path to the `ChromHMM.jar` file.

Make sure you have Java >= `1.8.0_281` installed and the Python packages listed in `requirements.txt`. 

### Training and applying the model
Use `./submit` to run the full Snakemake pipeline. 
