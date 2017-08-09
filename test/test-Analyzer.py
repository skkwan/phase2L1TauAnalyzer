# define basic process
import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1CaloClusters")
 
# import standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D4_cff')
process.load('Configuration.Geometry.GeometryExtended2023D13Reco_cff') ## this needs to match the geometry you are running on
process.load('Configuration.Geometry.GeometryExtended2023D13_cff')     ## this needs to match the geometry you are running on

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


# input
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
Source_Files = cms.untracked.vstring(
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SinglePion_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/00469E3E-7D26-E711-BBE0-5065F37D9132.root'
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SinglePion0_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/0024CEE8-5829-E711-A2CA-5065F382C231.root'
    #'file:singleE-gen-sim-raw.root'
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleE_FlatPt-8to100/GEN-SIM-DIGI-RAW/PU200_90X_upgrade2023_realistic_v9-v1/120000/002A4121-132C-E711-87AD-008CFAFBF618.root'
    #'file:/hdfs/store/mc/PhaseIISpring17D/GluGluHToZZTo4L_M125_14TeV_powheg2_JHUgenV702_pythia8/GEN-SIM-DIGI-RAW/PU200_100M_90X_upgrade2023_realistic_v9-v1/120000/8882318A-8E44-E711-AA96-0242AC130002.root'
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/TauThreeProngsEnriched/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/042B46D1-8E55-E711-A6C7-0CC47AA53D5A.root'
    '/store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/PU140_MB100M_90X_upgrade2023_realistic_v9-v1/120000/0052D3DF-4C56-E711-84A9-0025907DE22C.root'

    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/02E3D76F-8A55-E711-B9F0-0025901FB188.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/044D566A-C157-E711-BEDE-0CC47A0AD6AA.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/069E68E8-C457-E711-8E15-0CC47A0AD742.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/08AF49F4-8A55-E711-A870-0CC47A0AD74E.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/0AFC8973-C557-E711-8578-0025907859B8.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/0E4E7912-7857-E711-A03C-0CC47AA53DBE.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/0E4FFCB6-8B57-E711-A904-0CC47A57CB8E.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/0E7A1838-5557-E711-8B7B-0CC47A57CB8E.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/0ECD0F32-B157-E711-B00C-0025901FB188.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/109938E4-7B57-E711-92C4-0CC47A57CB8E.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/120AE875-B357-E711-AAC3-0025907D24F0.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/129004EE-C557-E711-B4A3-00259075D62E.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/18113155-C457-E711-812B-0CC47A57CD6A.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/188CA032-C057-E711-913E-0CC47A57CD6A.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/18C21A5B-C457-E711-927F-002590FD5A4C.root',
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/1AB4EFC0-B357-E711-8479-0CC47A57CB8E.root'
    )
process.source = cms.Source("PoolSource", fileNames = Source_Files)

############################################################
# L1 tracking
############################################################

# remake stubs 
# ===> IMPORTANT !!! stub window tuning as is by default in CMSSW is incorrect !!! <===
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)

process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")

# run only the tracking (no MC truth associators)
process.TTTracks = cms.Path(process.L1TrackletTracks)

# run the tracking AND MC truth associators)
process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators)

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
process.load("L1Trigger.phase2L1TauAnalyzer.phase2L1TauAnalyzer_cfi")
process.analyzer = cms.Path(process.L1TauAnalyzer)

# output module
#process.out = cms.OutputModule( "PoolOutputModule",
#                                fileName = cms.untracked.string("CaloClusters.root"),
#                                fastCloning = cms.untracked.bool( False ),
#                                outputCommands = cms.untracked.vstring(#'drop *',
#                                                                       'keep *_*_L1Phase2CaloClusters_*', 
#                                                                       )
#)
#process.FEVToutput_step = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("analyzer-1Prong.root"), 
   closeFileFast = cms.untracked.bool(True)
)


process.schedule = cms.Schedule(process.L1Clusters,process.TTClusterStub,process.TTTracksWithTruth,process.L1PFObjects,process.L1PFTaus,process.analyzer)


