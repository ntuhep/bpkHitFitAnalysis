#!/bin/bash

scramv1 project CMSSW CMSSW_5_3_11
export SCRAM_ARCH=slc5_amd64_gcc462
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATRecipes#CMSSW_5_2_X_pro_2012
cd CMSSW_5_3_11/src/
#cmsenv
eval `scramv1 runtime -sh`
addpkg DataFormats/PatCandidates V06-05-06-07
addpkg DataFormats/StdDictionaries V00-02-14
addpkg FWCore/GuiBrowsers V00-00-70

addpkg CommonTools/ParticleFlow V00-03-16
addpkg RecoParticleFlow/PFProducer V15-02-09  
addpkg CommonTools/RecoUtils V00-00-09
addpkg RecoMET/METFilters V00-00-08
addpkg RecoMET/METAnalyzers V00-00-08
addpkg CommonTools/RecoAlgos V00-03-23 
#cvs co -d EGamma/EGammaAnalysisTools -r V00-00-08 UserCode/EGamma/EGammaAnalysisTools
addpkg JetMETCorrections/Type1MET V04-06-09
addpkg PhysicsTools/PatAlgos V08-09-52
addpkg PhysicsTools/PatUtils V03-09-22
#addpkg TopQuarkAnalysis/TopPairBSM tlbsm_53x_v1_004
addpkg TopQuarkAnalysis/TopPairBSM                      V04-02-09
addpkg RecoBTag/SecondaryVertex V01-10-02
addpkg RecoVertex/AdaptiveVertexFinder V02-02-06

cvs co -d EGamma/EGammaAnalysisTools -r V00-00-31 UserCode/EGamma/EGammaAnalysisTools
#cvs co -d EGamma/EGammaAnalysisTools -r V00-00-18 UserCode/EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools/data
cat download.url | xargs wget
cd -

cvs co -r HEAD UserCode/sixie/Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h

#HITFIT
cvs co TopQuarkAnalysis/TopHitFit
#git clone git@github.com:ntuhep/bpkHitFit.git MyAna/bpkHitFit
git clone https://github.com/ntuhep/bpkHitFit.git MyAna/bpkHitFit

#bprimeKit
#git clone git@github.com:ntuhep/bprimeKit.git MyAna/bprimeKit
git clone https://github.com/ntuhep/bprimeKit.git MyAna/bprimeKit

cvs co -r HEAD UserCode/NTUHEP/jetTools.py.txt 
# doubleSecondaryVertexHighEffBJetTags_cfi.py  inclusiveSecondaryVertexFinderTagInfos_cfi.py
cvs co -r HEAD UserCode/NTUHEP/doubleSecondaryVertexHighEffBJetTags_cfi.py
cvs co -r HEAD UserCode/NTUHEP/inclusiveSecondaryVertexFinderTagInfos_cfi.py

cp UserCode/NTUHEP/jetTools.py.txt PhysicsTools/PatAlgos/python/tools/jetTools.py 
cp UserCode/NTUHEP/doubleSecondaryVertexHighEffBJetTags_cfi.py RecoBTag/SecondaryVertex/python/
cp UserCode/NTUHEP/inclusiveSecondaryVertexFinderTagInfos_cfi.py RecoBTag/SecondaryVertex/python/

## for coherent noise
cvs co -r V01-00-11-01 DPGAnalysis/Skims
cvs co -r V00-10-10-06 DPGAnalysis/SiStripTools
cvs co -r V00-00-08 DataFormats/TrackerCommon
cvs co -r V01-09-05 RecoLocalTracker/SubCollectionProducers

## for tobtec fakes
cvs co -d KStenson/TrackingFilters UserCode/KStenson/TrackingFilters
cp KStenson/TrackingFilters/plugins/TobTecFakesFilter.cc RecoMET/METFilters/plugins/
cp KStenson/TrackingFilters/python/tobtecfakesfilter_cfi.py RecoMET/METFilters/python
rm -r KStenson/TrackingFilters

## gluon tag
cvs co -r v1-2-3 -d QuarkGluonTagger/EightTeV UserCode/tomc/QuarkGluonTagger/EightTeV

#git clone git@github.com:ntuhep/bpkHitFitAnalysis.git MyAna/bpkHitFitAnalysis
git clone https://github.com/ntuhep/bpkHitFitAnalysis.git MyAna/bpkHitFitAnalysis
sed -i 's/<flags/#<flags/g' MyAna/bprimeKit/BuildFile.xml

scramv1 b -j 8

