# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: repr --processName=REPR --python_filename=reprocess_test_10_5_0_pre1.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 2 --era Phase2 --eventcontent FEVTDEBUGHLT --filein root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext4-v1/280000/FF5C31D5-D96E-5E48-B97F-61A0E00DF5C4.root --conditions 103X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D28 --fileout file:step2_2ev_reprocess_slim.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('REPR', eras.Phase2C8_trigger)
#process = cms.Process('REPR',eras.Phase2C4_trigger)
#process = cms.Process('REPR',eras.Phase2C4_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource")

process.source.fileNames = cms.untracked.vstring(
 "root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRSpring19MiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8_TuneCP5/MINIAODSIM/NoPU_106X_upgrade2023_realistic_v3-v2/130000/66C6C21F-A481-0D4C-9205-0332759441F2.root")


process.source.secondaryFileNames = cms.untracked.vstring(
 "root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRSpring19DR/GluGluHToTauTau_M125_14TeV_powheg_pythia8_TuneCP5/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v2/130000/2B1379C8-2D5B-C044-85F6-AD3788F02B34.root","root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRSpring19D\
R/GluGluHToTauTau_M125_14TeV_powheg_pythia8_TuneCP5/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v2/130000/EA19C627-D81C-8744-B1AE-7097DD56B1E7.root")


process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:22","1:5")

process.options = cms.untracked.PSet(

)
 
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('repr nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step2_2ev_reprocess_slim.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

############################################################
# L1 Tau object
############################################################

process.load("L1Trigger.Phase2L1Taus.L1PFTauProducer_cff")
process.L1PFTauProducer.L1PFObjects = cms.InputTag("l1pfCandidates","PF")
process.L1PFTauProducer.L1Neutrals = cms.InputTag("l1pfCandidates")
process.L1PFTauProducer.L1Clusters = cms.InputTag("l1pfCandidates","PF")
process.L1PFTauProducer.min_pi0pt = cms.double(3)
process.L1PFTaus = cms.Path(process.L1PFTauProducer)

# L1 Tau Analyzer
process.load("L1Trigger.phase2L1TauAnalyzer.phase2L1TauAnalyzer_cfi")

process.analyzer = cms.Path(process.L1TauAnalyzer)

process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("analyzer.root"), 
   closeFileFast = cms.untracked.bool(True)
)

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("dump.root"),
#    outputCommands = cms.untracked.vstring('keep *') 
#)

#process.e = cms.EndPath(process.out)


# Schedule definition
#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
process.schedule = cms.Schedule(process.L1simulation_step,process.L1PFTaus,process.analyzer,process.endjob_step)#,process.e)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

from L1Trigger.Configuration.customiseUtils import L1TrackTriggerTracklet
process = L1TrackTriggerTracklet(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
