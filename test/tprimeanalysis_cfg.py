import FWCore.ParameterSet.Config as cms

process = cms.Process("TPANA")

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

options.register('inFile',
				 'results.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input file name")

options.register('outFile',
				 'tprime.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name")

options.parseArguments()
process.Timing = cms.Service("Timing")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.load("RecoBTag.PerformanceDB.BTagPerformanceDB2013")
process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB2013")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )#don't change, this doesn't set max events

process.source = cms.Source("EmptySource")

from MyAna.bprimeKit.HitFitParameters_cfi import *
from MyAna.bprimeKit.EventParameters_cfi import *
from MyAna.bpkHitFitAnalysis.JetMetSystematicsParameters_cfi import *
from MyAna.bpkHitFitAnalysis.BTagSFUtilParameters_cfi import *

#these constants should probably be moved somewhere else
mW = 80.4
mZ = 91.2
mH = 125.9

debug = cms.untracked.bool(False)

useSplitJets = False
decayToGeneric = False # T->Wq as opposed to T->Wb. Also for T->Zq
hadMass = mW # set to the mass of whatever boson we want for the hadronic decay


jetBranches = 'PFJetInfo'
if useSplitJets:
    jetBranches = 'SplitJetInfo'

ObjectSelection = defaultObjectParameters.clone(
    Debug = debug,
    UsePFIsolation = cms.untracked.bool(True)
    )

EventSelection = defaultEventParameters.clone(
    Debug = debug,
    MinNJets  = cms.untracked.int32(4),
    #MinJetPTs = cms.untracked.vdouble(100.,75.),
    #MinMET    = cms.untracked.double(20.),
    CutLevels = cms.untracked.vstring('Initial','Vertex','Lepton','MuonVeto','ElectronVeto','MinJets'),
    #CutLevels = cms.untracked.vstring('Initial','Trigger','Vertex','Lepton','MuonVeto','ElectronVeto','MinJets','TightJetPt','MinMET'),
    ObjectParameters=ObjectSelection.clone()
    )

JetMetSystematics = defaultJetMetSystematicsParameters.clone(
	JECfile = cms.untracked.string('Summer12_V2_DATA_AK5PF_UncertaintySources.txt'),
	#JECfile = cms.untracked.string('/afs/cern.ch/user/t/twang/public/GR_P_V40_AN2_Uncertainty_AK5PF.txt'),
    Type = cms.untracked.string('None'),
    Scale = cms.untracked.double(0.),
    Debug = debug
    )

BTagSFUtil = defaultBTagSFUtilParameters.clone(
    BTagAlgorithm = cms.untracked.string('CSV'),
    BTagCuts      = cms.untracked.vdouble(0.679,0.244),
    BTagEffs      = cms.untracked.vdouble(0.70,0.85),
    Debug         = debug
    )   


configFileHitFit = 'MyAna/bpkHitFit/data/setting/RunHitFitConfiguration.txt'
if decayToGeneric:
   configFileHitFit = 'MyAna/bpkHitFit/data/setting/RunHitFitConfiguration_bMassless.txt'

HitFit = defaultHitFitParameters.clone(
    Debug = debug,
	Default = cms.untracked.FileInPath(configFileHitFit),
    JetCorrectionLevel = cms.untracked.string('L3'),
	HadWMass = cms.untracked.double(hadMass),
    NuSolution = cms.untracked.int32(2),
    #TopMass = cms.untracked.double(172.9),
    MaxNJet = cms.untracked.uint32(5)
    )

process.demo = cms.EDAnalyzer(
    'bpkHitFitAnalysis',
    InputFile        = cms.untracked.string(options.inFile),
    MaxEvents        = cms.untracked.int32(-1),#set to -1 for all events
    OutputFile       = cms.untracked.string(options.outFile),
	NAutoSave        = cms.untracked.int32(0),#How often to AutoSave (0 leaves default)
    Debug            = debug,
    Channels         = cms.untracked.vint32(11,13),#11 for electron, 13 for muon, can do both
    LeptonCollection = cms.untracked.string('PFLepInfo'),
    JetCollection    = cms.untracked.string(jetBranches),
    CutFlow              = cms.untracked.bool(False),#provide histogram showing how many events passed each selection level
    Skim                 = cms.untracked.bool(False),
    SelectionParameters  = EventSelection.clone(),
    RunHitFit            = cms.untracked.bool(True),#run HitFit
    PriorityBTags        = cms.untracked.int32(0),# change order according to btag
    HitFitParameters     = HitFit.clone(),
    DoJetMetSystematics  = cms.untracked.bool(False),# applied systematic
    JetMetSystematicsParameters = JetMetSystematics.clone(),
    BTagSFUtilParameters = BTagSFUtil.clone(),
	BTagUtilitySigma     = cms.untracked.double(0),# determined the sigma of btag(btag uncertainty)
    ModifyBTags          = cms.untracked.bool(False),# applied btag scaling factor
    BTagAlgorithm        = cms.untracked.string('CSV'), #CSV or TCHE
    BTagCut              = cms.untracked.double(0.679), #CSVM=0.679, TCHEM=3.3
    StripBranches        = cms.untracked.vstring('PairInfo*','JetInfo*','LepInfo*')
)

##Output file
#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string(options.outFile)
#)

from pprint import pprint
pprint (vars(process.demo))

process.p = cms.Path(process.demo)
