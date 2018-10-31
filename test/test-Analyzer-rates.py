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
    input = cms.untracked.int32(10000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:pickevents-1prong-pf-fail-raw.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/00036F95-8259-E811-A6EB-0025905D1E0A.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/009A1EDA-045B-E811-BD51-0025905C543A.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/00FF7D28-5D59-E811-BEEE-0025905C53F0.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0251273F-5359-E811-B536-0025905C5502.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0256C839-4959-E811-8B3D-0CC47AFB7FFC.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/026B3B35-2A59-E811-88D9-0025904C68DA.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/02CED5E8-8259-E811-8390-0025905C5488.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/045A6158-0459-E811-B0CE-0025905C53AA.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/047EC7E7-2B59-E811-94EF-0CC47AF9B1AE.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/04A6F11D-8A59-E811-A069-0CC47AF9B2CE.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/04AD1EE2-9659-E811-A6C3-0025905C2C84.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/04BE2323-4059-E811-95A3-0025905C3D3E.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/04F1C569-AF5A-E811-ABA9-0025905C5476.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0615710B-6B59-E811-8B2C-0025905C4262.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0620E820-7E59-E811-9454-0CC47AFB7D08.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/06340963-2A59-E811-80D3-0025904CF710.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/06771039-5059-E811-8AF8-0025905C3D3E.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/06895E4B-5659-E811-9F5F-0025905C3DD8.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/06C5259E-2859-E811-87D8-0025905C431C.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/085BAAD6-1E59-E811-915F-0025905C2CA6.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/08725B36-FC58-E811-96AC-0025905C4264.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0881F78F-8359-E811-AA56-0CC47AFB80B0.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/08C0A083-6659-E811-9420-0025905C9740.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/08EDE4AA-4659-E811-A74B-0CC47AFB80E4.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0A107647-2B59-E811-AF09-0025905C2CD2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0A2611C7-EF59-E811-93C5-0CC47AF9B2FE.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0A2B07E5-6D59-E811-9208-0025905C3E36.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0A37A541-7159-E811-9CB0-0CC47AF9B20A.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0A492C94-8059-E811-8F29-0025905C5502.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0A710EBA-6E59-E811-A160-0025905C2CBE.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0A7180BA-2459-E811-B25B-0025905C9724.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0C0D5448-E859-E811-A6C4-0CC47AFB7DA8.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0C3D629C-3859-E811-9929-0025905C94D2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0C3F8479-5A59-E811-BA58-0025905C3D98.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0C7A217D-8859-E811-BD6B-0CC47AF9B2C2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0C9F064B-AA5A-E811-9026-0CC47AFB7D88.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0E0A8BF2-8559-E811-886D-0025905C2CD2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0E24BA15-5659-E811-AEBC-0025904C6568.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU140_93X_upgrade2023_realistic_v5-v1/80000/0E3CF352-3B59-E811-BB02-0CC47AFB7D88.root'

        ),
    secondaryFileNames = cms.untracked.vstring(
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

process.L1TrackTrigger_step = cms.Path(process.L1TrackletTracksWithAssociators)

process.VertexProducer.l1TracksInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

############################################################
# L1 pf object
############################################################

process.load("L1Trigger.Phase2L1ParticleFlow.pfTracksFromL1Tracks_cfi")

from L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff import *
process.l1pf = cms.Path(process.pfTracksFromL1Tracks+process.l1ParticleFlow)

############################################################
# L1 Tau object
############################################################

process.load("L1Trigger.Phase2L1Taus.L1PFTauProducer_cff")
process.L1PFTauProducer.L1PFObjects = cms.InputTag("l1pfProducer","PF")
process.L1PFTauProducer.L1Neutrals = cms.InputTag("l1pfProducer")
process.L1PFTauProducer.L1Clusters = cms.InputTag("l1pfProducer","PF")
process.L1PFTaus = cms.Path(process.L1PFTauProducer)

# L1 Tau Analyzer
process.load("L1Trigger.phase2L1TauAnalyzer.phase2L1TauAnalyzerRates_cfi")
process.analyzer = cms.Path(process.L1TauAnalyzerRates)


process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("analyzer.root"), 
   closeFileFast = cms.untracked.bool(True)
)

# Schedule definition
process.schedule = cms.Schedule(process.EcalEBtp_step,process.L1TrackTrigger_step,process.L1simulation_step,process.l1pf,process.L1PFTaus,process.analyzer,process.endjob_step) 

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
