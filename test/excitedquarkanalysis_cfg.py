import FWCore.ParameterSet.Config as cms

process = cms.Process("TPANA")

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

options.register('inFile',
                 #'rfio:/castor/cern.ch/user/y/ymtzeng/CASTOR_FOLDER/INPUT_FOLDER/ROOTFILE',
				 #'root://eoscms//eos/cms/store/user/twang/2012ExcitedQuark/ReducedFiles/HitFit_Muon/TGluon450GeV.root',
				 'file:/afs/cern.ch/work/g/grundler/private/bprimeKit/TT2tgtg_M_800_results_10.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input file name")

options.register('outFile',
				 'excited.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name")

options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("CondCore.DBCommon.CondDBCommon_cfi")
#Datameasurments from Summer11
#process.load("RecoBTag.PerformanceDB.BTagPerformanceDB1107")
#process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")
#process.load("RecoBTag.PerformanceDB.BTagPerformanceDB100426")
#process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB100426")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )#don't change, this doesn't set max events

process.source = cms.Source("EmptySource")

from bpkFrameWork.bpkHitFit.HitFitParameters_cfi import *
from bpkFrameWork.bprimeKit.EventParameters_cfi import *
from bpkFrameWork.bpkHitFitAnalysis.JetMetSystematicsParameters_cfi import *
from bpkFrameWork.bpkHitFitAnalysis.BTagSFUtilParameters_cfi import *

debug = cms.untracked.bool(True)

usePFLeptons = True
usePFJets = True

lepBranches = 'LepInfo'
if usePFLeptons:
    lepBranches = 'PFLepInfo'

jetBranches = 'JetInfo'
if usePFJets:
    jetBranches = 'PFJetInfo'


ObjectSelection = defaultObjectParameters.clone(
    Debug = debug,
    UsePFIsolation = cms.untracked.bool(usePFLeptons)
    )

EventSelection = defaultEventParameters.clone(
    Debug = debug,
    MinNJets  = cms.untracked.int32(6),
    #MinJetPTs = cms.untracked.vdouble(100.,75.),
    #MinMET    = cms.untracked.double(20.),
    CutLevels = cms.untracked.vstring('Initial','Vertex','Lepton','MuonVeto','ElectronVeto','MinJets'),
	#CutLevels = cms.untracked.vstring('Initial','Vertex','Lepton'),
    #CutLevels = cms.untracked.vstring('Initial','Trigger','Vertex','Lepton','MuonVeto','ElectronVeto','MinJets'),
    #CutLevels = cms.untracked.vstring('Initial','Trigger','Vertex','Lepton','MuonVeto','ElectronVeto','MinJets','TightJetPt','MinMET'),
    ObjectParameters=ObjectSelection.clone()
    )

JetMetSystematics = defaultJetMetSystematicsParameters.clone(
	#JECfile = cms.untracked.string('Summer12_V2_DATA_AK5PF_UncertaintySources.txt'),
	JECfile = cms.untracked.string('/afs/cern.ch/user/t/twang/public/GR_P_V40_AN2_Uncertainty_AK5PF.txt'),
    Type = cms.untracked.string('jer'),
    Scale = cms.untracked.double(0.),
    Debug = debug
    )

BTagSFUtil = defaultBTagSFUtilParameters.clone(
    BTagAlgorithm = cms.untracked.string('CSV'),
    BTagCuts      = cms.untracked.vdouble(0.679,0.244),
    BTagEffs      = cms.untracked.vdouble(0.70,0.85),
    Debug         = debug
    )   

HitFit = defaultHitFitParameters.clone(
    Debug = debug,
	FitType = cms.untracked.int32(1),#TSTAR fit
	SaveOnlyBest = cms.untracked.bool(True),#only take best permutation
	Default = cms.untracked.FileInPath("bpkFrameWork/bpkHitFit/data/setting/RunHitFitConfiguration_topgluon.txt"),
    JetCorrectionLevel = cms.untracked.string('L3'),
    NuSolution = cms.untracked.int32(2),
	RequireMatchedBtag = cms.untracked.bool(True),
    TopMass = cms.untracked.double(172.9),
    MaxNJet = cms.untracked.uint32(6)
    )

process.demo = cms.EDAnalyzer(
    'bpkHitFitAnalysis',
    InputFile        = cms.untracked.string(options.inFile),
    MaxEvents        = cms.untracked.int32(-1),#set to -1 for all events
    #MaxEvents        = cms.untracked.int32(100),#set to -1 for all events
    OutputFile       = cms.untracked.string(options.outFile),
    Debug            = debug,
    Channels         = cms.untracked.vint32(13),#11 for electron, 13 for muon, can do both
    LeptonCollection = cms.untracked.string(lepBranches),
    JetCollection    = cms.untracked.string(jetBranches),
    CutFlow              = cms.untracked.bool(True),#provide histogram showing how many events passed each selection level
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
