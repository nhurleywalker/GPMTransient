#! /bin/bash -l

echo "Loading the transient data reduction pipeline"

if [[ -z $GXBASE ]]
then 
    echo "The GXBASE variable is not available, implying the GLEAM-X pipeline is not available. Exiting. "
    return 1
fi

if [[ -z $GPMBASE ]]
then 
    echo "The GPMBASE variable is not available, implying the Galactic plane monitoring pipeline is not available. Exiting. "
    return 1
fi

# Who is running the pipeline, used below for base install
GPTUSER=$(whoami)
export GPTUSER  

# The location of where the GP Transient pipeline is installed
export GPTBASE="/software/projects/pawsey0272/nhurleywalker/GPMTransient"

# The location where slurm log files will be saved for GPT tasks
export GPTLOG="${GPTBASE}/logs"

# The locaiton where generated scripts submitted to slurm will be placed for GPT tasks
export GPTSCRIPT="${GPTBASE}/scripts"

# These are to see the mwa-asvo copy queue business in obs_mantra to use garrawarla 
# as the longer running copy queue
export GXCOPYA='pawsey0272'             # Account to submit obs_manta.sh job under, if time accounting is being performed by SLURM.
                            # Leave this empty if the job is to be submitted as the user and there is no time accounting.
export GXCOPYQ='copy'      # A required parameter directing the job to a particular queue on $GXCOPYM. Set as just the queue name, e.g. 'copyq'
export GXCOPYM='setonix'       # A required parameter directing the job to be submitted to a particular machine. Set as just the machine name, e.g. 'zeus'

# Making sure the path for tasks are available 
export PATH="${PATH}:${GPTBASE}/bin:${GPTBASE}/dynspec:${GPTBASE}"

# Setting up the path to ensure the base GPT directory available for python scripts
export SINGULARITY_BINDPATH="${SINGULARITY_BINDPATH},${GPTBASE}"

# Necessary from singularity v3.8 onward so that PATH is sent to the container properly
export SINGULARITYENV_APPEND_PATH=$PATH

# Manual setup of the python path, meant to be used _within_ the container context
# as at the moment the GPT python code is not a proper module nor built into the 
# container
export PYTHONPATH="$PYTHONPATH:${GPTBASE}:${GPTBASE}/dynspec"

# Co-ordinates of the transient, used in a bunch of other tasks
export COORDS="18h39m02s -10d31m49.5s"

# Creates directories as needed below for the mandatory paths if they do not exist
if [[ ! -d "${GPTLOG}" ]]
then
    mkdir -p "${GPTLOG}"
fi

if [[ ! -d "${GPTSCRIPT}" ]]
then
    mkdir -p "${GPTSCRIPT}"
fi
