import FWCore.ParameterSet.Config as cms

def read_filelist(file):
    with open(file, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

filelist = read_filelist("rootfiles/rootfilessim1.txt") #rootfilessim1.txt or rootfilessim2.txt

process = cms.Process("Test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(filelist)
)


process.demo = cms.EDAnalyzer('MiniAnalyzerSimTrue',
                              GenParticle = cms.untracked.InputTag("prunedGenParticles"),
                              GenEventInfo = cms.untracked.InputTag("generator")
)


process.p = cms.Path(process.demo)
