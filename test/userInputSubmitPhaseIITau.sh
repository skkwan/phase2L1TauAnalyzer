#!/bin/sh                                                                                                                              
#SingleTau.txt
#voms-proxy-init --voms cms --valid 100:00                                                                                              

cat SUB-AnalyzerRates.py > SUBPhase2.py

rm -rf    /nfs_scratch/ojalvo/MinBias-0PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/ojalvo/MinBias-0PU-$1-SUBPhase2/dags/daginputs

rm -rf    /nfs_scratch/ojalvo/SingleNeutrino-140PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/ojalvo/SingleNeutrino-140PU-$1-SUBPhase2/dags/daginputs

rm -rf    /nfs_scratch/ojalvo/SingleNeutrino-200PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/ojalvo/SingleNeutrino-200PU-$1-SUBPhase2/dags/daginputs

#rm -rf    /nfs_scratch/ojalvo/QCD-0PU-$1-SUBPhase2/
#mkdir -p /nfs_scratch/ojalvo/QCD-0PU-$1-SUBPhase2/dags/daginputs

#rm -rf    /nfs_scratch/ojalvo/QCD-200PU-$1-SUBPhase2/
#mkdir -p /nfs_scratch/ojalvo/QCD-200PU-$1-SUBPhase2/dags/daginputs


farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/MinBias-93X-0PU.txt \
--submit-dir=/nfs_scratch/ojalvo/Minbias-0PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/MinBias-0PU-$1-SUBPhase2/dags/dag \
MinBias-0PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2


#farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/QCD-93X-200PU.txt \
#--submit-dir=/nfs_scratch/ojalvo/QCD-200PU-$1-SUBPhase2/submit \
#--output-dag-file=/nfs_scratch/ojalvo/QCD-200PU-$1-SUBPhase2/dags/dag \
#QCD-200PU-$1  \
#$CMSSW_BASE  \
#$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

#farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/QCD-93X-0PU.txt \
#--submit-dir=/nfs_scratch/ojalvo/QCD-0PU-$1-SUBPhase2/submit \
#--output-dag-file=/nfs_scratch/ojalvo/QCD-0PU-$1-SUBPhase2/dags/dag \
#QCD-0PU-$1  \
#$CMSSW_BASE  \
#$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleNeutrino-93X-140PU.txt \
--submit-dir=/nfs_scratch/ojalvo/SingleNeutrino-140PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/SingleNeutrino-140PU-$1-SUBPhase2/dags/dag \
SingleNeutrino-140PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleNeutrino-93X-200PU.txt \
--submit-dir=/nfs_scratch/ojalvo/SingleNeutrino-200PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/SingleNeutrino-200PU-$1-SUBPhase2/dags/dag \
SingleNeutrino-200PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2


rm SUBPhase2.py

