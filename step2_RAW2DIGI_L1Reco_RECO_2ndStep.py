# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --conditions 103X_dataRun2_Prompt_v2 -s RAW2DIGI,L1Reco,RECO --process reRECO -n 30 --data --era Run2_2018_pp_on_AA --eventcontent AOD --runUnscheduled --scenario pp --datatier AOD --repacked --filein file:out.root --fileout file:step2.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('APPROXIMATED',eras.Run2_2018_pp_on_AA)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_DataMapper_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.Timing = cms.Service("Timing",
#  summaryOnly = cms.untracked.bool(False),
#  useJobReport = cms.untracked.bool(True)
#)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step2.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:30'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

#RecoHiTrackerRECO = cms.PSet(
#    outputCommands = cms.untracked.vstring(
#	'keep *_*siStrip*_*_*', 
#	'keep *_*SiStrip*_*_*',
#        'keep *_*RawData*_*_*'
#    )
#)
#
#process.AODEventContent.outputCommands.extend(RecoHiTrackerRECO.outputCommands)

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AOD'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('file:step2_2nd_ModifiedRecHits.root'),
    outputCommands = process.AODEventContent.outputCommands
   # outputCommands = cms.untracked.vstring(
   #	'keep *', 
    #)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)

#cluster summary producer needs to ignore strips
#process.clusterSummaryProducer.doStrips= False



process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODoutput_step = cms.EndPath(process.AODoutput)

process.load('RecoLocalTracker.SiStripClusterizer.SiStripApprox2Clusters_cfi')
process.SiStripConvertedClusters = process.SiStripApprox2Clusters.clone(
  inputApproxClusters = 'SiStripClusters2ApproxClusters' 
)
process.ConvertedClusterPath = cms.Path(process.SiStripConvertedClusters)

process.MeasurementTrackerEventPreSplitting.stripClusterProducer = "SiStripConvertedClusters"
process.MeasurementTrackerEvent.stripClusterProducer = "SiStripConvertedClusters"

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.approxClusters,process.newClusters,process.endjob_step,process.AODoutput_step)
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.ApproximateRecHitsPath,process.reconstruction_step,process.endjob_step,process.AODoutput_step)
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODoutput_step)
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.ConvertedClusterPath,process.endjob_step,process.AODoutput_step)
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.ConvertedClusterPath,process.reconstruction_step,process.endjob_step,process.AODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


from Configuration.Applications.ConfigBuilder import MassReplaceInputTag
#MassReplaceInputTag(process, new="rawDataMapperByLabel", old="rawDataCollector")
MassReplaceInputTag(process, new="rawDataRepacker", old="rawDataCollector")
#MassReplaceInputTag(process, new="ApproximateRecHits", old="siStripMatchedRecHits")
MassReplaceInputTag(process, new="SiStripConvertedClusters", old="siStripClusters")

#do not add changes to your config after this point (unless you know what you are doing)
#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
#process = MassReplaceInputTag(process,"siStripClusters","SiStripApproxClusters2Clusters")
#process = MassReplaceInputTag(process,"siStripMatchedRecHits","ApproximateRecHits")
