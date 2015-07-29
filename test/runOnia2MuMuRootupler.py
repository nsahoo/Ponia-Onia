import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v0')

process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('/store/user/nsahoo/Charmonium/Onia2MuMuPAT-Run2015B-Charmonium-v3/150729_042156/0000/Onia2MuMuPAT-Run2015B_105.root')
    fileNames = cms.untracked.vstring('/store/user/nsahoo/MuOnia/Onia2MuMuPAT-Run2015B-MuOnia-v3/150729_044939/0000/Onia2MuMuPAT-Run2015B_112.root')
#    fileNames = cms.untracked.vstring('/store/user/nsahoo/MuOnia/Onia2MuMuPAT-Run2015B-MuOnia-v1/150727_073035/0000/Onia2MuMuPAT-Run2015B_102.root')
#    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/n/nsahoo/Run2-TRG-study/CMSSW_7_4_7_patch2/src/HeavyFlavorAnalysis/Skimming/test/Onia2MuMuPAT-Run2015B.root')

)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Onia2MuMuRootuple-Run2015B-MuOnia.root'),
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('Ponia.Onia.Onia2MuMuRootupler_cfi')
process.p = cms.Path(process.rootuple)

process.rootuple.isMC = cms.bool(False)                 # is mc?
process.rootuple.onia_mass_cuts = cms.vdouble(0.,300.)    # you may need to adjust this
