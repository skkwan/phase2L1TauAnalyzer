import FWCore.ParameterSet.Config as cms

L1TauAnalyzer = cms.EDAnalyzer('phase2L1TauAnalyzer',
                               L1Clusters       = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                               l1PFObjects      = cms.InputTag("L1PFProducer","L1PFObjects"),
                               l1TauObjects     = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                               L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),            
                               genParticles     = cms.InputTag("genParticles", "", "HLT"),
                               packedCandidates = cms.InputTag("packedPFCandidates","","RECO"),
                               ecalTPGsBarrel = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
                               miniTaus         = cms.InputTag("slimmedTaus", "", "RECO")
)
