import FWCore.ParameterSet.Config as cms

defaultJetMetSystematicsParameters = cms.PSet(
    Type             = cms.untracked.string('None'),
    Scale            = cms.untracked.double(0.),
    MinDRLepJet      = cms.untracked.double(0.3),
    JECfile          = cms.untracked.string('GR_R_42_V19_AK5PF_Uncertainty.txt'),
    nEtaBins         = cms.untracked.int32(28),
    nPtBins          = cms.untracked.int32(39),
	#https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    JEREtaMax        = cms.untracked.vdouble(0.5,1.1,1.7,2.3,5.0),
    JERNominal       = cms.untracked.vdouble(1.052,1.057,1.096,1.134,1.288),
	JERSigmaSym      = cms.untracked.vdouble(0.012,0.012,0.017,0.035,0.127),
    JERSigmaNeg      = cms.untracked.vdouble(0.061,0.055,0.062,0.085,0.153),
    JERSigmaPos      = cms.untracked.vdouble(0.062,0.056,0.063,0.087,0.155),
	JERMinGenJetPt   = cms.untracked.double(0.0),
    swUncertainty    = cms.untracked.double(0.015),
    puUncertainty    = cms.untracked.double(0.2),
    jetArea          = cms.untracked.double(0.8),
    averagePU        = cms.untracked.double(2.2),
    bjetUncertainty  = cms.untracked.double(0.02),
    bjetUncertainty2 = cms.untracked.double(0.03)
    )
