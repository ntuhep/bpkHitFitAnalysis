import FWCore.ParameterSet.Config as cms

defaultBTagSFUtilParameters = cms.PSet(
    MaxNJets      = cms.untracked.int32(8),
    BTagEffSF     = cms.untracked.double(1.),
    BTagAlgorithm = cms.untracked.string('CSV'),
    BTagCuts      = cms.untracked.vdouble(0.679,0.244),
    BTagEffs      = cms.untracked.vdouble(0.70,0.85),
    Debug         = cms.untracked.bool(False)
    )
