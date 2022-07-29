#! /bin/bash

usage()
{
echo "obs_followim.sh [-d dep] [-p project] [-a account] [-z] [-t] obsnum
  -d dep     : job number for dependency (afterok)
  -p project : project, (must be specified, no default)
  -t         : test. Don't submit job, just make the batch file
               and then return the submission command
  obsnum     : the obsid of the calibrator, the data around which will be processed." 1>&2;
exit 1;
}

pipeuser="${GXUSER}"

#initial variables
dep=
tst=
debug=
# parse args and set options
while getopts ':tzd:a:p:' OPTION
do
    case "$OPTION" in
	d)
	    dep=${OPTARG}
	    ;;
    a)
        account=${OPTARG}
        ;;
    p)
        project=${OPTARG}
        ;;
    z)
        debug=1
        ;;
	t)
	    tst=1
	    ;;
	? | : | h)
	    usage
	    ;;
  esac
done
# set the obsid to be the first non option
shift  "$(($OPTIND -1))"
obsnum=$1

queue="-p ${GXSTANDARDQ}"
base="${GXSCRATCH}/$project"

# if obsid is empty then just print help

if [[ -z ${obsnum} ]] || [[ -z $project ]] || [[ ! -d ${base} ]]
then
    usage
fi

if [[ ! -z ${dep} ]]
then
    depend="--dependency=afterok:${dep}"
fi

if [[ ! -z ${GXACCOUNT} ]]
then
    account="--account=${GXACCOUNT}"
fi

# start the real program

output="${GPTLOG}/followim_${obsnum}.o%A"
error="${GPTLOG}/followim_${obsnum}.e%A"
script="${GPTSCRIPT}/followim_${obsnum}.sh"

cat "${GPTBASE}/templates/followim.tmpl" | sed -e "s:CALID:${obsnum}:g" \
                                 -e "s:BASEDIR:${base}:g" \
                                 -e "s:DEBUG:${debug}:g" \
                                 -e "s:SCRIPT:${script}:g" \
                                 -e "s:ERROR:${error}:g" \
                                 -e "s:OUTPUT:${output}:g" \
                                 -e "s:QUEUE:${queue}:g" \
                                 -e "s:ACCOUNT:${account}:g" \
                                 -e "s:PIPEUSER:${pipeuser}:g" > "${script}"

chmod 755 "${script}"

# sbatch submissions need to start with a shebang
#echo '#!/bin/bash' > "${script}.sbatch"
#echo "singularity run ${GXCONTAINER} ${script}" >> "${script}.sbatch"

sub="sbatch --begin=now+1minutes --export=ALL  --time=01:00:00 --mem=3G -M ${GXCOMPUTER} --output=${output} --error=${error}"
sub="${sub} --cores=1 --ntasks-per-node=1 ${account} ${depend} ${queue} ${script}"
if [[ ! -z ${tst} ]]
then
    echo "script is ${script}"
    echo "submit via:"
    echo "${sub}"
    exit 0
fi
    
# submit job
jobids=($(${sub}))
jobid=${jobids[3]}

output=${output//%A/"${jobid}"}
error=${error//%A/"${jobid}"}

echo "Submitted ${script} as ${jobid} . Follow progress here:"
echo "${output}"
echo "${error}"
