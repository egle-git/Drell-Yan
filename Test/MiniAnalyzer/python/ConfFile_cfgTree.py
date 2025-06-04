import FWCore.ParameterSet.Config as cms

def read_filelist(file):
    with open(file, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

filelist = read_filelist("rootfiles/rootfiles_test_real.txt")

process = cms.Process("Test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(filelist)
)

process.hltHighLevel = cms.EDFilter("HLTHighLevel",
   TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
   HLTPaths = cms.vstring('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ*'),       # provide list of HLT paths (or patterns) you want
   eventSetupPathsKey = cms.string(''),      # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
   andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
   throw = cms.bool(True)    # throw exception on unknown path names
)

process.demo = cms.EDAnalyzer('MiniAnalyzerTree',
                              muons = cms.InputTag("slimmedMuons"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.p = cms.Path(process.hltHighLevel+process.demo)
