import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1EventDisplayPhase2")

 
# import standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D4_cff')

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

Source_Files = cms.untracked.vstring(
    #'file:event206074.root'
    #'file:1PiPMEvent.root'
    #'file:1pi0Event.root'
    #'file:4Pi0events.root'
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SinglePion0_FlatPt-8to100/GEN-SIM-DIGI-RAW/PU140_90X_upgrade2023_realistic_v9-v1/70000/0001E4B8-CD28-E711-B90A-E0071B7AA780.root'
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SinglePion0_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/0024CEE8-5829-E711-A2CA-5065F382C231.root'
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SinglePion_FlatPt-8to100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/00469E3E-7D26-E711-BBE0-5065F37D9132.root'
    #'file:/hdfs/store/mc/PhaseIISpring17D/GluGluHToZZTo4L_M125_14TeV_powheg2_JHUgenV702_pythia8/GEN-SIM-DIGI-RAW/PU200_100M_90X_upgrade2023_realistic_v9-v1/120000/8882318A-8E44-E711-AA96-0242AC130002.root'
    #'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/TauThreeProngsEnriched/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/042B46D1-8E55-E711-A6C7-0CC47AA53D5A.root'
    'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/02E3D76F-8A55-E711-B9F0-0025901FB188.root',
    'root://cms-xrd-global.cern.ch//store/mc/PhaseIISpring17D/SingleTauOneProngFlatPt10To100/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/120000/044D566A-C157-E711-BEDE-0CC47A0AD6AA.root'
    )
process.source = cms.Source("PoolSource", fileNames = Source_Files)

############################################################
# L1 tracking
############################################################

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
process.L1PFObjects = cms.Path(process.L1PFProducer)
process.L1PFProducer.debug = cms.untracked.bool(True)

############################################################
# L1 Tau object
############################################################
process.load("L1Trigger.phase2Demonstrator.L1PFTauProducer_cff")
process.L1PFTaus = cms.Path(process.L1PFTauProducer)
process.L1PFTauProducer.debug = cms.untracked.bool(True)


############################################################
# L1 Tau Analyzer
############################################################
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

# ----------------------------------------------------------------------------------------------
# 
# Analyzer starts here
process.tauAnalyzer = cms.EDAnalyzer('P2L1TEventDisplayGenerator',
                                     debug = cms.untracked.bool(True),
                                     L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),            
                                     L1TrackPrimaryVertexTag = cms.InputTag("L1TkPrimaryVertex"),
                                     ecalTPGsBarrel  = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
                                     hcalDigis       = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT"),
                                     L1Clusters      = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                                     l1TauObjects    = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                     vertices = cms.InputTag("vertices"),
                                     genParticles = cms.InputTag("genParticles"),
                                     L1CrystalClustersInputTag = cms.InputTag("L1EGammaCrystalsProducer","EGCrystalCluster"),
                                     genMatchDeltaRcut = cms.untracked.double(0.25),
                                     recoPtCut = cms.double(5),
                                     folderName = cms.untracked.string("dummy")
                                     )


process.panalyzer = cms.Path(process.tauAnalyzer)

process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("eventDisplay-1Prong.root"), 
   closeFileFast = cms.untracked.bool(True)
)


#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("fullDump-eventDisplay-1prong.root"),
#    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
#    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
#)

#process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.L1Clusters,process.TTTracksWithTruth,process.L1PFObjects,process.L1PFTaus,process.panalyzer)
