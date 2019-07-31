# phase2L1TauAnalyzer
Phase2 L1Tau Analyzer
Using the primary recipe from here:
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#Phase_2_L1T_Development_and_MC_R
```
cmsrel CMSSW_10_6_0_pre4
cd CMSSW_10_6_0_pre4/src
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline phase2-l1t-integration-CMSSW_10_6_0_pre4
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v2.17.26.1



git cms-addpkg L1Trigger/L1TCommon

cd L1Trigger
# to just clone the repo... though it is better if you fork it and clone from your own area!
git clone git@github.com:isobelojalvo/phase2L1TauAnalyzer.git

cd ../

USER_CXXFLAGS="-lTMVA -Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable" scram b -j 8

cd L1Trigger/phase2L1TauAnalyzer/test

cmsRun test-Analyzer.py
```
