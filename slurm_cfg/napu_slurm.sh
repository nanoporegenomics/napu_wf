#!/bin/bash

#This is a sample submission script for Napu pipeline on biowulf. At minimum, you need to set four parameters:
#CROMWELL_CFG: Cromwell configuration file.
#WORKFLOW: Path to workflow to run
#WORKDIR: execution directory, where cromwell outputs will also be written. Must have enough space and write access
#WF_INPUTS: json file with input parameters (paths to read files etc).
#Modify and submit as: sbatch --time 48:0:0 -c 4 --mem 8g napu_slurm.sh
#Then you can do: tail -f slurm-XX to monitor (XX is sbatch returned job id)

CROMWELL_CFG=PATH_TO_CFG
WORKFLOW=PATH_TO_cardEndToEndVcf.wdl
WORKDIR=PATH_TO_ECEX_DIR
WF_INPUTS=PATH_TO_INPUT_JSON

module load cromwell singularity
ulimit -u 10240 -n 16384
cd ${WORKDIR}

java -Dconfig.file=${NAPU_CFG} -jar ${CROMWELL_JAR} run -i ${WF_INPUTS} ${WORKFLOW}
