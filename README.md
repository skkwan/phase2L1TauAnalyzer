# phase2L1TauAnalyzer
Phase 2 L1 Tau Analyzer
Using the primary recipe from here:
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#CMSSW_10_6_1_patch2
```
cmsrel CMSSW_10_6_1_patch2
cd CMSSW_10_6_1_patch2/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline phase2-l1t-integration-CMSSW_10_6_1_patch2
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v2.22.1-CMSSW_10_6_1_patch2



git cms-addpkg L1Trigger/L1TCommon

cd L1Trigger
# to just clone the repo... though it is better if you fork it and clone from your own area!
git clone git@github.com:isobelojalvo/phase2L1TauAnalyzer.git

cd ../

USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable" scram b -j 8

cd L1Trigger/phase2L1TauAnalyzer/test

cmsRun test-Analyzer-reprocess.py
```
