#!/bin/sh

#Simple bash script to make fits with combine

if [ $# -ne 3 ]
then
	echo "Usage: $0 [datacard] [angular_cuts_string] [seed]"
	exit 10
fi

#path to the text datacard to be used
DATACARD=$1
#name of the output root file
CUTS=$2  #thmu%1.1f_the%1.0f
#seed of the toy
#if SEED_TOY < 0, then pseudodata will not be smeared in HistosForCombine2D_mappedIntoTH1.C
SEED_TOY=$3

#*************   TO CHANGE WITH YOUR FINAL DIR   *************
#FINAL_DIR=/afs/cern.ch/user/r/rpilato/CMSSW_10_2_13/src/MUonE/FINAL_THESIS/PostGraziano_2022-10-13/ULTIMO_FIT/fit2_K_MS1_Intr1_Ebeam20MeV/
FINAL_DIR=/bundle/scratch/30day/cezhang/MUonE/Combine


#list of nuisance parameters (see the datacard for the names)
#NUISANCE_PARS="provalnN_error,MultipleScattering,SingleHitRes,Ebeam"
NUISANCE_PARS="provalnN_error,MultipleScattering"


#echo "Setting the CMSSW environment..."
#CURRENT_DIR=$PWD
#ulimit -s unlimited
#set -e
#cd /afs/cern.ch/user/r/rpilato/CMSSW_10_2_13/src
#export SCRAM_ARCH=slc7_amd64_gcc700
#source /cvmfs/cms.cern.ch/cmsset_default.sh
#eval `scramv1 runtime -sh`

#cd $CURRENT_DIR
#pwd


#create directory to store the file with the best fit values created by combine
mkdir OutputCombine
#create directory to store the file with the covariance matrix created by combine
mkdir HesseMatrix

FIRST_TOY=$SEED_TOY
LAST_TOY=$SEED_TOY
NTEMPLATES=40
### fits with combine
for (( NTOY=$FIRST_TOY; NTOY<=$LAST_TOY; NTOY++ ))
do

	#1) take the histograms and make the toys
	root -l -b -q HistosForCombine2D_mappedIntoTH1.C'('${NTOY}', "'${CUTS}'")'

	#2) make the fit with combine for this toy
	echo "Submitting jobs to make template fits with combine..."
	for (( iK=0; iK<=$NTEMPLATES; iK++ ))
	do
		echo "iK = $iK"

        	#perform the fit with combine
        	combine -M MultiDimFit --datacard $DATACARD -v 1 \
        	--text2workspace --X-allow-no-background -n ${CUTS}_iK${iK}_TOY${NTOY} --saveNLL --keyword-value CUTS=${CUTS} \
        	--keyword-value iK=$iK --keyword-value NTOY=$NTOY --stepSize 0.00001 --cminDefaultMinimizerStrategy=0 \
        	--X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --rMin 1 --rMax 1 --robustHesse=1 --robustHesseSave 1 --saveFitResult --trackParameters $NUISANCE_PARS --trackErrors $NUISANCE_PARS 2> stderr_${CUTS}.txt 

        	echo ""
        	
        	#move the output file with best fit results in the directory created above
        	mv higgsCombine${CUTS}_iK${iK}_TOY${NTOY}.MultiDimFit.mH120.CUTS${CUTS}.iK${iK}.NTOY${NTOY}.root ./OutputCombine/higgsCombine2D_${CUTS}_iK${iK}_TOY${NTOY}.root
        	#move the output file with the covariance matrix in the directory created above
        	mv robustHesse${CUTS}_iK${iK}_TOY${NTOY}.root ./HesseMatrix/robustHesse${CUTS}_iK${iK}_TOY${NTOY}.root
			#remove the additional output file created by combine
        	rm multidimfit${CUTS}_iK${iK}_TOY${NTOY}.root # best fit in the nui para, but already avai in 72
	done

	#rm stderr_${CUTS}.txt


	#3) interpolate to calculate the best fit value of K and the nuisance parameters. Calculate the errors including the covariance matrix
	root -l -b -q getFitParameters_withCorrelations.C'('${NTOY}', "'${CUTS}'")'

	rm toy_data_${NTOY}.root

done




echo "Fit has been performed"
