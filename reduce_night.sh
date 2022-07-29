#!/bin/bash -l
#SBATCH --account=mwasci
#SBATCH --partition=workq
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --output=/software/projects/pawsey0272/nhurleywalker/GPMTransient/logs/reduce_D0038.o%A
#SBATCH --error=/software/projects/pawsey0272/nhurleywalker/GPMTransient/logs/reduce_D0038.e%A

# Manual setup: user must add their own output and error log destinations
# Otherwise they will go to /home/$USER 

# We might need to source to GLEAM-X and the GPM profiles here if this is
# a script kicked off from the slurm cron job. Not sure if bash profiles
# would be executed by the slurm magic. These are to ensure the 
# environment in the slurm cron job is set up correctly. 
GXPROFILE=/software/projects/pawsey0272/nhurleywalker/GLEAM-X-pipeline/GLEAM-X-pipeline-setonix.profile
GPMPROFILE=/software/projects/pawsey0272/nhurleywalker/MWA-Galactic-Plane-Monitoring/GP-Monitor-setonix.profile
GPTPROFILE=/software/projects/pawsey0272/nhurleywalker/GPMTransient/GPM-Transient-setonix.profile

for var in GXPROFILE GPMPROFILE GPTPROFILE
do
    if [[ -f "${!var}" ]]
    then 
        source "${!var}"
    else 
        echo "Unable to locate the ${var} profile. "
        echo "${var} currently set to ${!var} "
        return 1 
    fi
done

# Manual setup
datadir="${GXSCRATCH}"

# Common singularity command to run the python code
SINGCMD="singularity exec ${GXCONTAINER} "

chans="69 93 121 145 157 169"
for chan in $chans
do
    obsid=$(${SINGCMD} "${GPMBASE}/gp_monitor_lookup.py" --cal --cent-chan $chan --project D0038 --startdate '2022-07-25 00:00:00' --stopdate '2022-07-26 00:00:00')
    if [[ $obsid == "" ]]
    then
        obsid=$(${SINGCMD} "${GPMBASE}/gp_monitor_lookup.py" --cal --cent-chan $chan --project D0038 --calsrc 3C444)
    fi

    if [[ $obsid != "" ]]
    then
        project=D0038_$(${SINGCMD} "${GPTBASE}/determine_night.py" --obsid $obsid)

        # Note: output format is: D0038_YYYY-MM-DD
        pdir="$datadir/$project"

        if [[ ! -d "${pdir}" ]]
        then
            mkdir "${pdir}"
        fi

        echo "${obsid}" > "${pdir}/cal_${chan}.txt"

        dep=($(obs_manta.sh -p "${project}" -o "${pdir}/cal_${chan}.txt"))
        depend=${dep[3]}
        dep=($(obs_autoflag.sh -d ${depend} -p "${project}" "${obsid}"))
        depend=${dep[3]}
        dep=($(obs_autocal.sh -d ${depend} -p "${project}" -f 0.5 "${obsid}"))
        depend=${dep[3]}
        obs_followim.sh -d ${depend} -p "${project}" "${obsid}"
        
    fi
done
