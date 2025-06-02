import FWCore.ParameterSet.Config as cms

def read_filelist(file):
    with open(file, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

filelist = read_filelist("rootfiles/rootfilessimtchanantitop.txt")

process = cms.Process("Test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(filelist)
)


process.hltHighLevel = cms.EDFilter("HLTHighLevel",
   TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
   HLTPaths = cms.vstring('HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ*'),
   eventSetupPathsKey = cms.string(''),
   andOr = cms.bool(True),
   throw = cms.bool(True)
)

process.demo = cms.EDAnalyzer('MiniAnalyzerSim',
                              muons = cms.InputTag("slimmedMuons"),
                              GenParticle = cms.untracked.InputTag("prunedGenParticles"),
                              GenEventInfo = cms.untracked.InputTag("generator"),
                              mcProcess = cms.string("tchanantitop"),
)


process.p = cms.Path(process.hltHighLevel+process.demo)
