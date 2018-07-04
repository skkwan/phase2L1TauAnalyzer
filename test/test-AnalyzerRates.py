# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --python_filename=rerun_step2_L1_onMCL1_FEVTHLTDEBUG.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 1 --era Phase2_timing --eventcontent FEVTDEBUGHLT --filein file:/afs/cern.ch/user/r/rekovic/release/CMSSW_9_3_2/src/step2_DIGI_PU200_10ev.root --conditions 93X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D17 --fileout file:step2_ZEE_PU200_1ev_rerun-L1-L1Ntuple.root --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleEMU
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('L1',eras.Phase2_trigger)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/00157B11-405C-E811-89CA-0CC47AFB81B4.root',
        #'/store/mc/PhaseIISpring17D/SingleNeutrino/GEN-SIM-DIGI-RAW/PU200_90X_upgrade2023_realistic_v9-v1/70002/2E87B534-B526-E711-98D7-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/026A6839-352D-E811-B52C-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/02940DEE-332D-E811-8199-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/02C09B32-352D-E811-9AE0-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/04DEC942-3B2D-E811-A916-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/082D030B-3C2D-E811-9296-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0880AAC8-332D-E811-A5DD-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0A163304-3C2D-E811-A4FC-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0A8C41DD-332D-E811-B680-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0AAC1F0B-352D-E811-BD97-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0AC603F2-3B2D-E811-B00C-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0AD75FF3-3B2D-E811-9225-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0C9C6D39-352D-E811-BD40-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0C9E7DEE-3B2D-E811-B76A-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0E72DCF5-3B2D-E811-AC8F-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/0E8E0440-352D-E811-830D-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/10BEC4C7-332D-E811-8C9C-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/14650440-352D-E811-A8D2-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/14EF783D-2F2D-E811-B4F0-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/163965F8-3B2D-E811-A5FE-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/165F73CC-332D-E811-82F6-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/16CB6667-2F2D-E811-A502-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/16FFAB47-352D-E811-A331-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/1A051D2C-352D-E811-9098-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/1AE6D567-3C2D-E811-BCA4-0242AC130002.root',
        '/store/relval/CMSSW_9_3_7/RelValQCD_Pt-15To7000_Flat_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/1C3E06C8-332D-E811-B4F4-0242AC130002.root'
        #'/store/relval/CMSSW_9_3_7/RelValMinBias_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/040DB7E8-832C-E811-A87F-0CC47A4C8F08.root',
        ),
        #'file:pickevents-1prongPi0-fail-raw.root'),
        #'file:pickevents-mini-1prong-fail-raw.root'),
        #'file:pickevents-1prongPi0-fail-mini.root'),
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/MINIAODSIM/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/C033607A-8E2C-E811-B4EF-0CC47A78A478.root'),
    secondaryFileNames = cms.untracked.vstring(
        #'file:pickevents-1prongPi0-fail-raw.root'
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/02306B8E-6D2C-E811-9625-0025905B858C.root',
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/163A5C8E-6D2C-E811-B990-0025905B85B2.root',
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/1A7E6D8F-6D2C-E811-8DFB-0025905B85BE.root',
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/5864268A-6D2C-E811-90F5-0025905A60EE.root',
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/627E087F-6D2C-E811-8F81-0CC47A4C8E34.root',
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/64AE2B7F-6D2C-E811-800E-0CC47A4D76A0.root',
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/820C0278-6D2C-E811-B95B-0CC47A4C8F0C.root',
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/BE6A3C80-6D2C-E811-9143-0CC47A4D7668.root',
        #'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/DA7C408A-6D2C-E811-9BFF-0025905B858A.root'
        ),
    inputCommands = cms.untracked.vstring("keep *", 
        "drop l1tHGCalTowerMapBXVector_hgcalTriggerPrimitiveDigiProducer_towerMap_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT",
        "drop l1tEMTFHit2016s_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT",
        "drop l1tEMTFHitExtras_simEmtfDigis_CSC_HLT",
        "drop l1tEMTFHitExtras_simEmtfDigis_RPC_HLT",
        "drop l1tEMTFTrackExtras_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016s_simEmtfDigis__HLT")                            
     #skipEvents = cms.untracked.uint32(80)
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:test_reprocess.root'),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')


process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives)

process.load('SimCalorimetry.EcalEBTrigPrimProducers.ecalEBTriggerPrimitiveDigis_cff')
process.EcalEBtp_step = cms.Path(process.simEcalEBTriggerPrimitiveDigis)

#process.TTClusterAssociatorFromPixelDigis.digiSimLinks          = cms.InputTag( "simSiPixelDigis","Tracker" )
process.L1TrackTrigger_step = cms.Path(process.L1TrackletTracksWithAssociators)

process.VertexProducer.l1TracksInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
#process.schedule = cms.Schedule(process.EcalEBtp_step,process.hgcl1tpg_step,process.L1TrackTrigger_step,process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
#process.schedule = cms.Schedule(process.L1TrackTrigger_step,process.L1simulation_step,process.endjob_step) 
#,process.FEVTDEBUGHLToutput_step)

#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
#from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAWEMUGEN_MC 

#call to customisation function L1NtupleEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
#process = L1NtupleRAWEMUGEN_MC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion



############################################################
# L1 clusters
############################################################

process.load("L1Trigger.phase2Demonstrator.L1CaloClusterProducer_cff")
process.L1CaloClusterProducer.debug = cms.untracked.bool(False)
process.L1Clusters = cms.Path(process.L1CaloClusterProducer)

############################################################
# L1 pf object
############################################################

process.load("L1Trigger.phase2Demonstrator.L1PFProducer_cff")
process.L1PFProducer.debug = cms.untracked.bool(False)
process.L1PFObjects = cms.Path(process.L1PFProducer)

############################################################
# L1 Tau object
############################################################
process.load("L1Trigger.phase2Demonstrator.L1PFTauProducer_cff")
process.L1PFTauProducer.debug = cms.untracked.bool(False)
process.L1PFTaus = cms.Path(process.L1PFTauProducer)


# L1 Tau Analyzer
process.load("L1Trigger.phase2L1TauAnalyzer.phase2L1TauAnalyzerRates_cfi")
process.analyzer = cms.Path(process.L1TauAnalyzerRates)

# output module
#process.out = cms.OutputModule( "PoolOutputModule",
#                                fileName = cms.untracked.string("CaloClusters.root"),
#                                fastCloning = cms.untracked.bool( False),
#                                outputCommands = cms.untracked.vstring(#'drop *',
#                                                                       'keep *_*_L1Phase2CaloClusters_*', 
#                                                                       )
#)
#process.FEVToutput_step = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("analyzer_rates.root"), 
   closeFileFast = cms.untracked.bool(True)
)

# Schedule definition
#process.schedule = cms.Schedule(process.EcalEBtp_step,process.hgcl1tpg_step,process.L1TrackTrigger_step,process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
process.schedule = cms.Schedule(process.L1Clusters,process.L1TrackTrigger_step,process.L1simulation_step,process.L1PFObjects,process.L1PFTaus,process.analyzer,process.endjob_step) 

#,process.FEVTDEBUGHLToutput_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
#from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAWEMUGEN_MC 

#call to customisation function L1NtupleEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
#process = L1NtupleRAWEMUGEN_MC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
