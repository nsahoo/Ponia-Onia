import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v0')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_10_1_Z7q.root')
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Onia2MuMuRootuple-Run2015B-MuOnia-v2.root'),
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('Ponia.Onia.Onia2MuMuRootupler_cfi')
process.p = cms.Path(process.rootuple)

process.rootuple.isMC = cms.bool(False)                 # is mc?
process.rootuple.onia_mass_cuts = cms.vdouble(0.,300.)    # you may need to adjust this
