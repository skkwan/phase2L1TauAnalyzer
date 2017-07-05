import FWCore.ParameterSet.Config as cms

L1TauAnalyzer = cms.EDAnalyzer('phase2L1TauAnalyzer',
                               L1Clusters      = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                               l1PFObjects     = cms.InputTag("L1PFProducer","L1PFObjects"),
                               genParticles    = cms.InputTag("genParticles", "", "HLT"),
)
