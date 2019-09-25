import FWCore.ParameterSet.Config as cms

L1TauAnalyzerRates = cms.EDAnalyzer('phase2L1TauAnalyzerRates',
                                    L1Clusters      = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                                    L1PFObjects      = cms.InputTag("l1pfCandidates","PF"),
                                    l1TauObjects    = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                    L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                    ecalTPGsBarrel = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
                                    genParticles    = cms.InputTag("genParticles", "", "HLT"),
                                    L1VertexInputTag = cms.InputTag("L1TkPrimaryVertex")
)
