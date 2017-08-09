#!/bin/sh                                                                                                                              
#SingleTau.txt
#voms-proxy-init --voms cms --valid 100:00                                                                                              

cat test-Analyzer.py > SUBPhase2.py
cat submit.py >> SUBPhase2.py

rm -r    /nfs_scratch/ojalvo/SingleTau-1Prong-0PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/ojalvo/SingleTau-1Prong-0PU-$1-SUBPhase2/dags/daginputs

rm -r    /nfs_scratch/ojalvo/SingleTau-1Prong-140PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/ojalvo/SingleTau-1Prong-140PU-$1-SUBPhase2/dags/daginputs

rm -r    /nfs_scratch/ojalvo/SingleTau-1Prong-200PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/ojalvo/SingleTau-1Prong-200PU-$1-SUBPhase2/dags/daginputs

rm -r    /nfs_scratch/ojalvo/SingleTau-3Prong-0PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/ojalvo/SingleTau-3Prong-0PU-$1-SUBPhase2/dags/daginputs

rm -r    /nfs_scratch/ojalvo/SingleTau-3Prong-140PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/ojalvo/SingleTau-3Prong-140PU-$1-SUBPhase2/dags/daginputs

rm -r    /nfs_scratch/ojalvo/SingleTau-3Prong-200PU-$1-SUBPhase2/
mkdir -p /nfs_scratch/ojalvo/SingleTau-1Prong-200PU-$1-SUBPhase2/dags/daginputs


farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleTauOneProng-92X-0PU.txt \
--submit-dir=/nfs_scratch/ojalvo/SingleTau-1Prong-0PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/SingleTau-1Prong-0PU-$1-SUBPhase2/dags/dag \
SingleTau-1Prong-0PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleTauOneProng-92X-140PU.txt \
--submit-dir=/nfs_scratch/ojalvo/SingleTau-1Prong-140PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/SingleTau-1Prong-140PU-$1-SUBPhase2/dags/dag \
SingleTau-1Prong-140PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleTauOneProng-92X-200PU.txt \
--submit-dir=/nfs_scratch/ojalvo/SingleTau-1Prong-200PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/SingleTau-1Prong-200PU-$1-SUBPhase2/dags/dag \
SingleTau-1Prong-200PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2


farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleTauThreeProng-92X-0PU.txt \
--submit-dir=/nfs_scratch/ojalvo/SingleTau-3Prong-0PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/SingleTau-3Prong-0PU-$1-SUBPhase2/dags/dag \
SingleTau-3Prong-0PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleTauThreeProng-92X-140PU.txt \
--submit-dir=/nfs_scratch/ojalvo/SingleTau-3Prong-140PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/SingleTau-3Prong-140PU-$1-SUBPhase2/dags/dag \
SingleTau-3Prong-140PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

farmoutAnalysisJobs --assume-input-files-exist  --input-file-list=samples/SingleTauThreeProng-9220X-0PU.txt \
--submit-dir=/nfs_scratch/ojalvo/SingleTau-3Prong-200PU-$1-SUBPhase2/submit \
--output-dag-file=/nfs_scratch/ojalvo/SingleTau-3Prong-200PU-$1-SUBPhase2/dags/dag \
SingleTau-3Prong-200PU-$1  \
$CMSSW_BASE  \
$CMSSW_BASE/src/L1Trigger/phase2L1TauAnalyzer/test/SUBPhase2.py $2

rm SUBPhase2.py

