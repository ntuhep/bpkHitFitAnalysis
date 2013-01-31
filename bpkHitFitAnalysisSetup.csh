#!/bin/tcsh

setenv SCRAM_ARCH slc5_amd64_gcc462
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATRecipes#CMSSW_5_2_X_pro_2012
cmsrel CMSSW_5_3_6_patch1
cd CMSSW_5_3_6_patch1/src/
cmsenv
addpkg CommonTools/ParticleFlow V00-03-16
addpkg RecoParticleFlow/PFProducer V15-02-09  
addpkg CommonTools/RecoUtils V00-00-09
addpkg RecoMET/METFilters V00-00-08
addpkg RecoMET/METAnalyzers V00-00-08
addpkg CommonTools/RecoAlgos V00-03-23 
#cvs co -d EGamma/EGammaAnalysisTools -r V00-00-08 UserCode/EGamma/EGammaAnalysisTools
addpkg JetMETCorrections/Type1MET V04-06-09
addpkg PhysicsTools/PatAlgos V08-09-14-00
addpkg PhysicsTools/PatUtils V03-09-22
addpkg TopQuarkAnalysis/TopPairBSM tlbsm_53x_v1_004
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
cvs co -d MyAna/bpkHitFit UserCode/NTUHEP/bpkHitFit

#bprimeKit HEAD version 
cvs co -d MyAna/bprimeKit UserCode/NTUHEP/bprimeKit
#cvs co -r Summer12_Pro2_53x_patch3 -d MyAna/bprimeKit UserCode/NTUHEP/bprimeKit
#cvs co -d MyAna/ExcitedQuarkAnalysis UserCode/NTUHEP/ExcitedQuarkAnalysis
cvs co -d MyAna/bpkHitFitAnalysis UserCode/NTUHEP/bpkHitFitAnalysis

#do a little cleaning for this analysis
rm MyAna/bprimeKit/interface/doHitFitForExcitedQuark.h
rm MyAna/bprimeKit/src/doHitFitForExcitedQuark.cc
sed -i 's/<flags/#<flags/g' MyAna/bprimeKit/BuildFile.xml

cvs co -r HEAD UserCode/NTUHEP/jetTools.py.txt 
cp UserCode/NTUHEP/jetTools.py.txt PhysicsTools/PatAlgos/python/tools/jetTools.py 

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

###
### How to use the bprimeKit with the double b tagger:
###
### The jetTools.py in PhysicsTools/PatAlgos/python/tools/ needs to be replaced:
### Check out the modified jetTools.py from 
### https://twiki.cern.ch/twiki/pub/CMS/InclusiveVertexFinderRecipes/jetTools.py.txt (it is already saved to the /UserCode/NTUHEP/ directory)
### Copy to PhysicsTools/PatAlgos/python/tools/ 
### Compile and run. 
###
