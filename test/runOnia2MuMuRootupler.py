import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v0')

process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_10_1_Z7q.root')
    fileNames = cms.untracked.vstring(
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_10_1_DQL.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_11_1_Y3Q.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_12_1_ymd.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_13_1_Qh0.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_14_1_exW.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_15_1_0db.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_16_1_S1i.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_17_1_XLy.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_18_1_0ti.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_19_1_Xz0.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_20_1_vHN.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_21_1_hfj.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_22_1_r6g.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_24_1_zCB.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_25_1_VVK.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_2_1_5wQ.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_3_1_zzD.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_4_1_d3O.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_5_1_2mI.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_6_1_ozN.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_7_1_STR.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_8_1_CCA.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v3/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_9_1_aEb.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_10_1_Z7q.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_11_1_ybZ.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_12_1_BYl.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_13_1_EhN.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_1_1_Vbv.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_2_1_YZi.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_3_1_VTv.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_4_1_1bJ.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_5_1_7CE.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_6_1_rMh.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_7_1_CUt.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_8_1_3T3.root',
'/store/user/zhenhu/MuOnia/Onia2MuMuPAT-Run2015B-MuOina-v2/3c0bc5c741de78e17f27570e8d4bbe40/Onia2MuMuPAT_9_1_6Pu.root')

)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Onia2MuMuRootuple-Run2015B-MuOnia-v2.root'),
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('Ponia.Onia.Onia2MuMuRootupler_cfi')
process.p = cms.Path(process.rootuple)

process.rootuple.isMC = cms.bool(False)                 # is mc?
process.rootuple.onia_mass_cuts = cms.vdouble(0.,300.)    # you may need to adjust this
