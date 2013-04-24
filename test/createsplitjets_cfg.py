import FWCore.ParameterSet.Config as cms

process = cms.Process("TPANA")

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

options.register('inFile',
				 'file:/afs/cern.ch/work/g/grundler/private/bprimeKit/tpWb_ntupl532p4_500-1.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input file name")

options.register('outFile',
				 'splitjets.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name")

options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )#don't change, this doesn't set max events

process.source = cms.Source("EmptySource")

process.demo = cms.EDAnalyzer(
    'createSplitJets',
    MaxEvents         = cms.untracked.int32(-1),#set to -1 for all events
    InputFile         = cms.untracked.string(options.inFile),
    OutputFile        = cms.untracked.string(options.outFile),
	InputStandardJets = cms.untracked.string("PFJetInfo"),
	InputWJets        = cms.untracked.string("WJetInfo"),
	OutputJets        = cms.untracked.string("SplitJetInfo"),
	MinWMass          = cms.untracked.double(60.),
	MaxWMass          = cms.untracked.double(100.),
	MaxMassDrop       = cms.untracked.double(-1.),
	MaxDRMatch        = cms.untracked.double(0.4),
    BTagCut           = cms.untracked.double(0.679), #CSVM=0.679
    Debug             = cms.untracked.bool(False)
)

from pprint import pprint
pprint (vars(process.demo))

process.p = cms.Path(process.demo)
