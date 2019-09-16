import FWCore.ParameterSet.Config as cms

L1TauAnalyzer = cms.EDAnalyzer('phase2L1TauAnalyzer',
                               L1PFObjects      = cms.InputTag("l1pfCandidates","PF"),
                               l1TauObjects     = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                               L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),            
                               genParticles     = cms.InputTag("genParticles", "", "HLT"),
#                               packedCandidates = cms.InputTag("packedPFCandidates","","RECO"),
                               packedCandidates = cms.InputTag("packedPFCandidates","","PAT"),
                               ecalTPGsBarrel = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
#                               miniTaus         = cms.InputTag("slimmedTaus", "", "RECO"),
                               miniTaus         = cms.InputTag("slimmedTaus", "", "PAT"),
                               L1VertexInputTag = cms.InputTag("L1TkPrimaryVertex")
)
