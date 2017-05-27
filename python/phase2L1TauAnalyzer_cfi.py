import FWCore.ParameterSet.Config as cms

L1TauAnalyzer = cms.EDAnalyzer('phase2L1TauAnalyzer',
                               L1Clusters      = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                               genParticles    = cms.InputTag("genParticles", "", "HLT"),
)
