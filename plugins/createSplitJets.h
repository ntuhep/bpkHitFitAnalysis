#ifndef __CREATESPLITJETS__
#define __CREATESPLITJETS__

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TChain;
class TFile;

class EvtInfoBranches;
class JetInfoBranches;
class LepInfoBranches;

using std::string;

class createSplitJets : public edm::EDAnalyzer {

public:
   explicit createSplitJets( const edm::ParameterSet& );
   ~createSplitJets();

   //       createSplitJets(const edm::ParameterSet&,
   //                         EvtInfoBranches &evt, JetInfoBranches &jets1, JetInfoBranches &jets2);


private:
   virtual void beginJob() ;
   virtual void analyze( const edm::Event&, const edm::EventSetup& );
   virtual void endJob() ;

   bool shouldSplit( int jet );
   void insertSplitJets( int idx_w, int idx_std, bool isMC );
   void insertJetCopy( int idx );

   TChain*          chain;
   int              maxEvents;
   string           inFile;
   string           outFile;
   TFile*            newfile;
   TTree*            newtree;

   EvtInfoBranches  _evt;
   JetInfoBranches  _inStdJets; //standard jet collection (AK5PFchs)
   JetInfoBranches  _inWJets;   //Wide jet collection to be split
   JetInfoBranches  _outJets;

   const string     _inStdJetCollection;
   const string     _inWJetCollection;
   const string     _outJetCollection;


   //W-jet tagging
   double _wMassLo;
   double _wMassHi;
   double _massDropMax;

   //jet matching
   double _matchDRMax;
   double _bTagCut;

   bool   _debug;

};

#endif
