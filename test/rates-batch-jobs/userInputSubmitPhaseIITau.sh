#!/bin/sh                                                                                                                              
#SingleTau.txt
#voms-proxy-init --voms cms --valid 100:00                                                                                              

cat test-Analyzer-rates.py > SUBPhase2.py
cat submit.py >> SUBPhase2.py

rm -rf    /nfs_scratch/skkwan/QCD-200PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/skkwan/QCD-200PU-$1-SUBPhase2/dags/daginputs

#rm -rf    /nfs_scratch/skkwan/SingleNeutrino-140PU-$1-SUBPhase2/
#mkdir -p /nfs_scratch/skkwan/SingleNeutrino-140PU-$1-SUBPhase2/dags/daginputs

#rm -rf    /nfs_scratch/skkwan/SingleNeutrino-200PU-$1-SUBPhase2/
#mkdir -p /nfs_scratch/skkwan/SingleNeutrino-200PU-$1-SUBPhase2/dags/daginputs

#rm -rf    /nfs_scratch/skkwan/QCD-0PU-$1-SUBPhase2/
#mkdir -p /nfs_scratch/skkwan/QCD-0PU-$1-SUBPhase2/dags/daginputs

#rm -rf    /nfs_scratch/skkwan/QCD-200PU-$1-SUBPhase2/
#mkdir -p /nfs_scratch/skkwan/QCD-200PU-$1-SUBPhase2/dags/daginputs


farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/QCD-103X-200PU.txt \
--submit-dir=/nfs_scratch/skkwan/QCD-200PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/skkwan/QCD-200PU-$1-SUBPhase2/dags/dag \
QCD-200PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/rates-batch-jobs/SUBPhase2.py $2


#farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/QCD-93X-200PU.txt \
#--submit-dir=/nfs_scratch/skkwan/QCD-200PU-$1-SUBPhase2/submit \
#--output-dag-file=/nfs_scratch/skkwan/QCD-200PU-$1-SUBPhase2/dags/dag \
#QCD-200PU-$1  \
#$CMSSW_BASE  \
#$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

#farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/QCD-93X-0PU.txt \
#--submit-dir=/nfs_scratch/skkwan/QCD-0PU-$1-SUBPhase2/submit \
#--output-dag-file=/nfs_scratch/skkwan/QCD-0PU-$1-SUBPhase2/dags/dag \
#QCD-0PU-$1  \
#$CMSSW_BASE  \
#$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

#farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleNeutrino-93X-140PU.txt \
#--submit-dir=/nfs_scratch/skkwan/SingleNeutrino-140PU-$1-SUBPhase2/submit \
#--output-dag-file=/nfs_scratch/skkwan/SingleNeutrino-140PU-$1-SUBPhase2/dags/dag \
#SingleNeutrino-140PU-$1  \
#$CMSSW_BASE  \
#$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

#farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleNeutrino-93X-200PU.txt \
#--submit-dir=/nfs_scratch/skkwan/SingleNeutrino-200PU-$1-SUBPhase2/submit \
#--output-dag-file=/nfs_scratch/skkwan/SingleNeutrino-200PU-$1-SUBPhase2/dags/dag \
#SingleNeutrino-200PU-$1  \
#$CMSSW_BASE  \
#$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2


#rm SUBPhase2.py
