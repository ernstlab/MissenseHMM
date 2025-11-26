import os

configfile: "config.yaml"

Z = config["Z"]
CHROMS = config["CHROMS"]
IN_LOC = config["IN_LOC"]
SCRIPT_LOC = config["SCRIPT_LOC"]
OUT_LOC = config["OUT_LOC"]
VERSION = config["VERSION"]
ENV = config["ENV"]
NUM_STATES = config["NUM_STATES"]
CHROMHMM_LOC = config["CHROMHMM_LOC"]


os.makedirs(OUT_LOC, exist_ok=True)
# os.makedirs(OUT_LOC + "/score_logs", exist_ok=True)


rule all:
    input:
        expand("{output}/dbnsfp_chr{chrm}.txt", chrm=CHROMS, output=OUT_LOC),
        expand("{output}/score_logs/dbnsfp_{z}.log", z=range(0, Z), output=OUT_LOC),
        OUT_LOC + "/thresholds.txt", 
        expand("{output}/samples/log_chr{chrm}.log", chrm=CHROMS, output=OUT_LOC),
        expand("{output}/model.log", output=OUT_LOC),
        OUT_LOC + "/max_states/all.bed"


rule create_datafile:
    input:
        "%s/dbNSFP%s/dbNSFP%s_variant.chr{chrm}.gz" % (IN_LOC, VERSION, VERSION)
    output:
        "%s/dbnsfp_chr{chrm}.txt" % OUT_LOC
    params:
        script_loc = SCRIPT_LOC,
        version = VERSION,
        output_loc = OUT_LOC
    threads: 1
    shell:
        """
        {params.script_loc}/to_datafile.py -c {wildcards.chrm} -v {params.version} -g -l --output {params.output_loc}
        """


rule threshold:
    input:
        expand("{output}/dbnsfp_chr{chrm}.txt", chrm=CHROMS, output=OUT_LOC)
    output:
        "%s/score_logs/dbnsfp_{z}.log" % OUT_LOC
    threads: 1
    params:
        script_loc = SCRIPT_LOC,
        output_loc = OUT_LOC
    shell:
        """
        mkdir -p {params.output_loc}/score_logs
        {params.script_loc}/threshold.py {wildcards.z} {params.output_loc} {params.output_loc}/score_logs
        """


rule generate_threshold_file:
    input:
        expand("{output}/score_logs/dbnsfp_{z}.log", z=range(0, Z), output=OUT_LOC)
    output:
        "%s/thresholds.txt" % OUT_LOC
    threads: 1
    params:
        script_loc = SCRIPT_LOC,
        output_loc = OUT_LOC
    shell:
        """
        {params.script_loc}/generate_thresholds.py {params.output_loc}/score_logs {params.output_loc}
        """


rule generate_binary:
    input:
        "%s/thresholds.txt" % OUT_LOC,
        expand("{output}/dbnsfp_chr{chrm}.txt", chrm=CHROMS, output=OUT_LOC)
    output:
        "%s/samples/log_chr{chrm}.log" % OUT_LOC
    threads: 1
    params:
        script_loc = SCRIPT_LOC,
        output_loc = OUT_LOC
    shell:
        """
        mkdir -p {params.output_loc}/samples
        {params.script_loc}/to_binary.py -c {wildcards.chrm} -i {params.output_loc} -o {params.output_loc}/samples -t {params.output_loc}/thresholds.txt
        echo DONE >  {params.output_loc}/samples/log_chr{wildcards.chrm}.log
        """


rule train_hmm:
    input: 
        expand("{output}/samples/log_chr{chrm}.log", chrm=CHROMS, output=OUT_LOC)
    output:
        "%s/model.log" % OUT_LOC
    threads: 4
    params:
        chromhmm_loc = CHROMHMM_LOC,
        output_loc = OUT_LOC,
        env = ENV,
        num_states = NUM_STATES
    shell:
        """
        if [ "{params.env}" = "cluster" ]; then
            . /u/local/Modules/default/init/modules.sh; module load java
        fi
        
        java -jar {params.chromhmm_loc}/ChromHMM.jar LearnModel -nobed -nobrowser -noenrich -noimage -pseudo -b 1 -n 2000 -d "-1" -lowmem -p 4 {params.output_loc}/samples {params.output_loc}/model_{params.num_states} {params.num_states} hg38
        
        echo DONE > {params.output_loc}/model.log
        """


rule apply_model:
    input:
        "%s/model.log" % OUT_LOC,
        "%s/thresholds.txt" % OUT_LOC,
        "%s/dbnsfp_chr{chrm}.txt" % OUT_LOC
    output:
        "%s/max_states/max_states_chr{chrm}.bed" % OUT_LOC
    threads: 1
    params:
        script_loc = SCRIPT_LOC,
        output_loc = OUT_LOC,
        num_states = NUM_STATES
    shell:
        """
        mkdir -p {params.output_loc}/max_states
        {params.script_loc}/max_pos_state.py -c {wildcards.chrm} -m {params.output_loc}/model_{params.num_states}/model_{params.num_states}.txt -o {params.output_loc}/max_states -t {params.output_loc}/thresholds.txt
        """
        
        
rule combine_states:
    input:
        expand("{output}/max_states/max_states_chr{chrm}.bed", chrm=CHROMS, output=OUT_LOC)
    output:
        "%s/max_states/all.bed" % OUT_LOC
    threads: 1
    params:
        output_loc = OUT_LOC
    shell:
        """
        cat {params.output_loc}/max_states/max_states_chr*.bed > {params.output_loc}/max_states/all.bed
        """