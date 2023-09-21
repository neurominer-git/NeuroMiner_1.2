#!/bin/bash 
echo
echo '****************************************'
echo '*** NeuroMiner                       ***'
echo '*** SGE joblist manager:             ***'
echo '*** Interpret model predictions      ***'
echo '*** (c) 2022 N. Koutsouleris         ***'
echo '****************************************'
echo '          NM VERSION 1.1  		      '
echo '****************************************'
echo   

# compiled with matlab R2022a so MCR main is v912. Needs to change if different MCR is used.
export LD_LIBRARY_PATH=/opt/matlab/v912/runtime/glnxa64:/opt/matlab/v912/bin/glnxa64:/opt/matlab/v912/sys/os/glnxa64:/opt/matlab/v912/sys/opengl/lib/glnxa64
# Preload glibc_shim in case of RHEL7 variants
export LD_PRELOAD=/opt/matlab/R2022a/bin/glnxa64/glibc-2.17_shim.so

export JOB_DIR=$PWD
export NEUROMINER=/opt/NM/NeuroMinerMCCMain_1.1_v912/for_testing
export ACTION=mli
read -e -p 'Path to NM structure: ' datpath
if [ ! -f $datpath ]; then
 	echo $datpath' not found.'
 	exit 
fi

read -e -p 'Path to job directory ['$JOB_DIR']: ' tJOB_DIR
if [ "$tJOB_DIR" != '' ]; then
	if [ -d $tJOB_DIR ]; then 
		export JOB_DIR=$tJOB_DIR
	else
		echo $tJOB_DIR' not found.'
		exit
	fi
fi 

read -e -p 'Change path to compiled NM directory ['$NEUROMINER']: ' tNEUROMINER
if [ "$tNEUROMINER" != '' ]; then     
  if [ -d $tNEUROMINER ]; then  
    export NEUROMINER=$tNEUROMINER
  else
    echo $tNEUROMINER' not found.'
    exit
  fi    
fi

read -p 'Provide your email address: ' EMAIL
echo '-----------------------'
echo 'PATH definitions:'
echo 'LOG directory: '$JOB_DIR
echo 'NeuroMiner directory: '$NEUROMINER
echo '-----------------------'
read -p 'Index to analysis container (NM.analysis{<index>}): ' analind
if [ "$analind" = '' ] ; then
	echo 'An analysis index is mandatory! Exiting program.'
	exit   
fi
read -p 'Index to independent data container (NM.OOCV{<index>} [0 => MLI will run in the discovery data, default]): ' oocvind
if [ "$oocvind" = '' ] ; then
	oocvind=0
fi
export optmodelspath=NaN 
export optparamspath=NaN 
read -p 'Save optimized preprocessing parameters and models to disk for future use [ 1 = yes, 2 = no ]: ' saveparam
if [ "$saveparam" = '2' ] ; then
  read -p 'Load optimized preprocessing parameters and models from disk [ 1 = yes, 2 = no ]: ' loadparam
  if [ "$loadparam" = '1' ] ; then
    read -e -p 'Path to OptPreprocParam master file: ' optparamspath
    if [ ! -f $optparamspath ] ; then
	    echo $optparamspath' not found.'
	    exit
    fi
    read -e -p 'Path to OptModels master file: ' optmodelspath
    if [ ! -f $optmodelspath ] ; then
	    echo $optmodelspath' not found.'
	    exit
    fi
  fi
else
  loadparam=2
fi
read -p 'Recompute interpretations based on predictions [ 1 = yes, 2 = no, default ]: ' reestimateflag
if [ "$reestimateflag" = '' ] ; then
	reestimateflag=2
fi

read -p 'CV2 grid start row: ' CV2x1
read -p 'CV2 grid end row: ' CV2x2
read -p 'CV2 grid start column: ' CV2y1
read -p 'CV2 grid end column: ' CV2y2
read -p 'No. of SGE jobs: ' numCPU
read -p 'Server to use [any=1, psy0cf20=2, mitnvp1=3]: ' sind
if [ "$sind" = '1' ]; then
        SERVER_ID='all.q'
        echo "Please estimate RAM accurately"
elif [ "$sind" = '2' ]; then
        SERVER_ID='psy0cf20'
        echo "Please estimate RAM accurately"
elif [ "$sind" = '3' ]; then
        SERVER_ID='mitnvp1-2'
        echo "Please estimate RAM accurately"
else
        echo "Enter a number between 1-3"
fi
read -p 'Enter "Max vmem" in GB from email sent from a test run of a single fold: ' GB
HALFGB=$(($GB / 2))
vGB=$(($GB + $HALFGB))'G'
mGB=$GB'G'

read -p 'Use OpenMP [yes = 1 | no = 0]: ' pLibsvm
if [ "$pLibsvm" = '1' ] ; then
	read -p 'Specify number of CPUs assigned to MATLAB job [4, 8, 16, 32]: ' pnum
	PMODE='#\$-pe omp '$pnum
else
	PMODE=''
	pnum=1
fi

read -p 'Submit jobs immediately (y): ' todo
for curCPU in $(seq $((numCPU)))
do
SD='_CPU'$curCPU
ParamFile=$JOB_DIR/Param_NM_$ACTION$SD
SGEFile=$JOB_DIR/NM_$ACTION$SD
echo 'Generate parameter file: NM_'$ACTION$SD' => '$ParamFile
# Generate parameter file
cat > $ParamFile <<EOF
$NEUROMINER
$datpath
$analind
$oocvind
$saveparam
$loadparam
$reestimateflag
$optparamspath
$optmodelspath
$curCPU
$numCPU
$CV2x1
$CV2x2
$CV2y1
$CV2y2
EOF
cat > $SGEFile <<EOF
#!/bin/bash
#\$-o $JOB_DIR/\$JOB_IDnm$ACTION$SD -j y
#\$-N nm$ACTION$SD
#\$-S /bin/bash
#\$-M $EMAIL
#\$-l mem_total=$mGB
#\$-l h_vmem=$vGB
#\$-q $SERVER_ID
$PMODE
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export LD_PRELOAD=$LD_PRELOAD
export OMP_NUM_THREADS=$pnum
cd $NEUROMINER 
./NeuroMinerMCCMain $ACTION $ParamFile
EOF
chmod u+x $SGEFile
datum=`date +"%Y%m%d"`
if [ "$todo" = 'y' -o "$todo" = 'Y' ] ; then
qsub $SGEFile >> NeuroMiner_MLI_$datum.log
fi
done
