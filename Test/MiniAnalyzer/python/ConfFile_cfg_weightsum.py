import FWCore.ParameterSet.Config as cms

def read_filelist(file):
    with open(file, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

filelist = read_filelist("rootfiles/rootfiles_test_ttnew.txt")
# rootfilessim1.txt  -  sim1
# rootfilessim2.txt  -  sim2
# rootfilessimtt.txt  -  tt
# rootfilessimww.txt  -  ww
# rootfilessimwz.txt  -  wz
# rootfilessimzz.txt  -  zz
# rootfilessimtwtop.txt  -  twtop
# rootfilessimtwantitop.txt  -  twantitop
# rootfilessimtchantop.txt  -  tchantop
# rootfilessimtchanantitop.txt  -  tchanantitop
process = cms.Process("Test")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(filelist)
)


process.demo = cms.EDAnalyzer('MiniAnalyzer_weightsum',
                              muons = cms.InputTag("slimmedMuons"),
                              GenEventInfo = cms.untracked.InputTag("generator")
)


process.p = cms.Path(process.demo)
