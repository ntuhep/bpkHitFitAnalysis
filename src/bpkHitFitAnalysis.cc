// -*- C++ -*-
//
// Package:    bpkHitFitAnalysis
// Class:      bpkHitFitAnalysis
// 
/**\class bpkHitFitAnalysis bpkHitFitAnalysis.cc MyAna/bpkHitFitAnalysis/src/bpkHitFitAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ulysses Grundler,598 R-016,+41227679822,
//         Created:  Thu Aug  4 16:05:23 CEST 2011
// Second Author: Yeng-Ming Tzeng, B13 2-054, +41764872910, 
// $Id: bpkHitFitAnalysis.cc,v 1.1 2012/11/20 13:46:38 grundler Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TH1.h"
#include "TFile.h"
#include "TKey.h"
#include "TChain.h"
#include "time.h"

#include "MyAna/bprimeKit/interface/format.h"
#include "MyAna/bprimeKit/interface/bpkUtils.h"
#include "MyAna/bprimeKit/interface/objectSelector.h"
#include "MyAna/bprimeKit/interface/eventSelector.h"
#include "MyAna/bpkHitFit/interface/HitFitInfoBranches.h"
#include "MyAna/bpkHitFit/interface/doHitFit.h"

#include "MyAna/bpkHitFitAnalysis/interface/jetMetSystematics.h"
#include "MyAna/bpkHitFitAnalysis/interface/BTagSFUtil-tprime.h"

#include "MyAna/bpkHitFitAnalysis/interface/GetBTag_SF_EFF.h"

static const int maxChannels = 2;
static const int maxJetsToModify = 8;//8 is maximum HitFit will take, so this should be plenty
//
// class declaration
//

class bpkHitFitAnalysis : public edm::EDAnalyzer {
   public:
      explicit bpkHitFitAnalysis(const edm::ParameterSet&);
      ~bpkHitFitAnalysis();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  std::vector<int> prioritizeBtags(std::vector<int> jets);
  bool jetIsBTagged(int ijet);

      // ----------member data ---------------------------
  bool               debug;

  TChain*            chain;
  int                maxEvents;
  std::string        inFile;
  std::string        outFile;
  TFile              *newfile;
  TTree              *newtree;
  int                nAutoSave;

  EvtInfoBranches    EvtInfo;
  VertexInfoBranches VtxInfo;
  LepInfoBranches    LepInfo;
  JetInfoBranches    JetInfo;
  GenInfoBranches    GenInfo;
  HitFitInfoBranches HitFitInfo;

  doHitFit           *hitfit;
  eventSelector      *eSelector;
  jetMetSystematics  *jmsScaler;

  edm::ParameterSet        HitFitParameters;
  edm::ParameterSet        SelectionParameters;
  edm::ParameterSet        JetMetSystParameters;
  edm::ParameterSet        BTagSFUtilParameters;

  std::string              lepcollection;
  std::string              jetcollection;

  bool                     doCutFlow;
  bool                     skimNtuple;
  bool                     runHitFit;
  int                      priorityBtags;
  bool                     doJetMetSystematics;
  bool                     modifyBTags;
  int                      nJetsToModify;
  std::vector<std::string> stripBranches;

  const std::vector<int>   channels;

  std::vector<std::string> cutLevels;
  TH1F* h_cutflow[maxChannels];

  std::string        bTagAlgo;
  double             bTagCut;
  double             btag_sigma;

//   TH1F* BTAG_SF;
//   TH1F* BTAG_EFF;
//   TH1F* BTAG_SF_LIGHT;
//   TH1F* BTAG_EFF_LIGHT;
//   TH1F* BTAG_EFF_LIGHT_PASS;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
bpkHitFitAnalysis::bpkHitFitAnalysis(const edm::ParameterSet& iConfig)
  : channels(iConfig.getUntrackedParameter< std::vector<int> >("Channels"))
{
   //now do what ever initialization is needed
  inFile  = iConfig.getUntrackedParameter<std::string>("InputFile");
  maxEvents = iConfig.getUntrackedParameter<int>("MaxEvents",-1);
  outFile = iConfig.getUntrackedParameter<std::string>("OutputFile");
  nAutoSave = iConfig.getUntrackedParameter<int>("NAutoSave",0);
  debug = iConfig.getUntrackedParameter<bool>("Debug",false);
  HitFitParameters = iConfig.getParameter<edm::ParameterSet>("HitFitParameters");
  SelectionParameters = iConfig.getParameter<edm::ParameterSet>("SelectionParameters");
  JetMetSystParameters = iConfig.getParameter<edm::ParameterSet>("JetMetSystematicsParameters");
  BTagSFUtilParameters = iConfig.getParameter<edm::ParameterSet>("BTagSFUtilParameters");

  lepcollection = iConfig.getUntrackedParameter<std::string>("LeptonCollection");
  jetcollection = iConfig.getUntrackedParameter<std::string>("JetCollection");

  doCutFlow = iConfig.getUntrackedParameter<bool>("CutFlow",false);
  skimNtuple = iConfig.getUntrackedParameter<bool>("Skim",false);
  runHitFit = iConfig.getUntrackedParameter<bool>("RunHitFit",false);
  priorityBtags = iConfig.getUntrackedParameter<int>("PriorityBTags",0);
  doJetMetSystematics = iConfig.getUntrackedParameter<bool>("DoJetMetSystematics",false);
  modifyBTags = iConfig.getUntrackedParameter<bool>("ModifyBTags",false);
  nJetsToModify = iConfig.getUntrackedParameter<int>("NJetsToModify",8);
  stripBranches = iConfig.getUntrackedParameter< std::vector<std::string> >("StripBranches");

  bTagAlgo = iConfig.getUntrackedParameter<std::string>("BTagAlgorithm","CSV");
  bTagCut  = iConfig.getUntrackedParameter<double>("BTagCut",0.679);
  btag_sigma  = iConfig.getUntrackedParameter<double>("BTagUtilitySigma",0.0);

}


bpkHitFitAnalysis::~bpkHitFitAnalysis()
{
  if(debug) std::cout << "destruct\n";
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete chain;

  delete hitfit;
  delete eSelector;
  delete jmsScaler;

  newfile->Close();
  delete newfile;

  if(debug) std::cout << "destruction done\n";
}


//
// member functions
//
// ------------ method called once each job just before starting event loop  ------------
void 
bpkHitFitAnalysis::beginJob()
{
  if(debug) std::cout << "Starting beginJob\n";
   chain = new TChain("bprimeKit/root");
   chain->Add(inFile.c_str());

   if(maxEvents<0 || maxEvents>chain->GetEntries())
     maxEvents = chain->GetEntries();

   EvtInfo.Register(chain);
   VtxInfo.Register(chain);
   LepInfo.Register(chain,lepcollection);
   JetInfo.Register(chain,jetcollection);
   GenInfo.Register(chain);
   UInt_t found;
   if(runHitFit) chain->SetBranchStatus("HitFitInfo*",0,&found);
   for(unsigned i=0; i<stripBranches.size(); i++) {
     chain->SetBranchStatus(stripBranches[i].c_str(),0,&found);
   }

   eSelector = new eventSelector(SelectionParameters,EvtInfo,LepInfo,JetInfo,VtxInfo);
   cutLevels = eSelector->getCutLevels();
   if(doJetMetSystematics) {
     jmsScaler = new jetMetSystematics(JetMetSystParameters,EvtInfo,JetInfo, LepInfo);
     if(debug) std::cout << "Initialized jetMetSystematics\n";
   }

   newfile = new TFile(outFile.c_str(),"recreate");
   if(skimNtuple || runHitFit) {
     newtree = chain->CloneTree(0);
     if(nAutoSave!=0) newtree->SetAutoSave(nAutoSave);

     if(runHitFit) {
       hitfit = new doHitFit(HitFitParameters,EvtInfo,LepInfo,JetInfo,GenInfo);
       HitFitInfo.RegisterTree(newtree);
     }

     if(debug) std::cout << "New tree ready\n";
   }

   //prepare histograms
	for(int i=0; i<(int)channels.size(); i++) {
		char cName[4];
	   	if (ELECTRON == channels[i]) sprintf(cName,"el");
   		else if(MUON == channels[i]) sprintf(cName,"mu");
	   	else                         sprintf(cName,"NA");//don't think this should happen.

    	char hName[120];
	    sprintf(hName,"cutflow_%s",cName);

    	if(doCutFlow) {
    		const int nCutLevels = (int)cutLevels.size();
	    	h_cutflow[i] = new TH1F(hName,hName,nCutLevels,-0.5,nCutLevels-0.5);
    	}
	}
    if(debug) std::cout << "Finishing beginJob\n";
	
// 	BTAG_SF = new TH1F("BTag SF", "", 100, 0.9, 1.0);
// 	BTAG_EFF = new TH1F("BTag EFF", "", 100, 1., 2.);
// 	BTAG_SF_LIGHT = new TH1F("BTag SF LIGHT", "", 100, 0.8, 1.2);
// 	BTAG_EFF_LIGHT = new TH1F("BTag EFF LIGHT", "", 100, 0., 0.1);
// 	BTAG_EFF_LIGHT_PASS = new TH1F("BTag EFF LIGHT_PASS", "", 100, 0., 0.1);
}



// ------------ method called for each event  ------------
void
bpkHitFitAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    //using namespace reco;
    if(debug) std::cout << "Starting analyze\n";

    BTagSFUtil* btsf = new BTagSFUtil(BTagSFUtilParameters);
	//btsf->setSeed(time(NULL));
	btsf->setSeed(1);

	for(int entry=0; entry<maxEvents; entry++) {//loop over entries
		if((entry%1000)==0) printf("Loading event %i of %i.\n",entry,maxEvents);
     	chain->GetEntry(entry);
		if(runHitFit) HitFitInfo.clear();

     	if(debug) std::cout << "Entry " << entry << ": Run " << EvtInfo.RunNo << " Event " << EvtInfo.EvtNo << std::endl;

     	if(doJetMetSystematics) jmsScaler->scale();

     	std::vector< std::pair<float,float> > jetmom;
     	std::vector<bool> jetisbtag;
	    for(int ijet=0; ijet<maxJetsToModify; ijet++) {
			jetisbtag.push_back(jetIsBTagged(ijet));
		 	std::pair<float,float> thisjetmom(JetInfo.Et[ijet],JetInfo.Eta[ijet]);
		 	jetmom.push_back(thisjetmom);
	 	}
     	if(EvtInfo.McFlag && modifyBTags) {
	     	btsf->setSeed(entry*1e+12+EvtInfo.RunNo*1e+6+EvtInfo.EvtNo);
		    //get et, eta of jetsjets
    		for(int ijet=0; ijet<JetInfo.Size; ijet++) {
			std::pair<float,float> thisjetmom(JetInfo.Et[ijet],JetInfo.Eta[ijet]);
			jetmom.push_back(thisjetmom);
			if(ijet<nJetsToModify) jetisbtag[ijet] = jetIsBTagged(ijet);
	    	}
		    btsf->readDB(iSetup,jetmom);
			//BTag_Utility btag_uti(JetInfo, btag_sigma);
		    for(int ijet=0; ijet<JetInfo.Size; ijet++) {
				if(ijet>=nJetsToModify) break;
// 				double btag_sf = btag_uti.GetSF(ijet);
// 				if(abs(JetInfo.GenFlavor[ijet]) == 5 || abs(JetInfo.GenFlavor[ijet]) ==4){BTAG_SF->Fill(btag_sf);}
// 				else {BTAG_SF_LIGHT->Fill(btag_sf);}
// 				double btag_eff = btag_uti.GetEFF(ijet);
// 				if(abs(JetInfo.GenFlavor[ijet]) == 5 || abs(JetInfo.GenFlavor[ijet]) ==4) {BTAG_EFF->Fill(btag_eff);}
// 				else {BTAG_EFF_LIGHT->Fill(btag_eff);	
// 					if(JetInfo.Pt[ijet] > 30 && fabs(JetInfo.Eta[ijet]) < 2.4 && JetInfo.JetIDLOOSE[ijet] ){ BTAG_EFF_LIGHT_PASS->Fill(btag_eff);}
// 				}
				
// 				btag_sf = btag_sf + btag_uti.Btag_SF_ERR;
// 				if(btag_sf == -999 || btag_eff == -999 || btag_uti.Btag_SF_ERR == -999) std::cout << "SF = -999, Please Check" << std::endl;
				double btag_sf = btsf->getSF("BTAG" + bTagAlgo + "M",ijet);
				double btag_eff = btsf->BtagEff_[0];//getSF("BTAG" + bTagAlgo + "Meff",ijet);
				//double bmistag_sf = btsf->getSF("MISTAG" + bTagAlgo + "M",ijet);
				//double bmistag_eff = btsf->getSF("MISTAG" + bTagAlgo + "Meff",ijet);
				if(debug) {
	   				std::cout << "Jet " << ijet << ": Et,Eta,PdgId,tag= " << jetmom[ijet].first << "," << jetmom[ijet].second << "," << JetInfo.GenFlavor[ijet] << "," << jetisbtag[ijet] << std::endl;
					//std::cout << " btageff_sf=" << btag_sf << " bmistag_sf=" << bmistag_sf << " bmistag_eff=" << bmistag_eff << std::endl;
					std::cout << " btag_sf=" << btag_sf << " btag_eff=" << btag_eff <<std::endl;
		 		}
     			//if(btag_sf<0 || btag_eff<0 || bmistag_sf<0 || bmistag_eff<0) {
				if(btag_sf<0 || btag_eff<0) {
		        std::cout << "modifyBTagsWithSF: SF<0, something is wrong(maybe just out of range where SF measured).  Doing nothing\n";
    		    continue;
     			}
				//btsf->modifyBTagsWithSF(jetisbtag[ijet],JetInfo.GenFlavor[ijet],btag_sf,btag_eff,bmistag_sf,bmistag_eff);
				bool _jetisbtag = jetisbtag[ijet];
			    btsf->modifyBTagsWithSF(_jetisbtag,JetInfo.GenFlavor[ijet],btag_sf,btag_eff);
				if(_jetisbtag != jetisbtag[ijet]) std::cout << "jet btag been modified" << std::endl;
				jetisbtag[ijet] = _jetisbtag;
				if(debug) std::cout << " after modification, tag= " << jetisbtag[ijet] << std::endl;
			}
     	}

		for(int ich=0; ich<(int)channels.size(); ich++) {//loop over channels
	        eSelector->setChannel(channels[ich]);
	
    	    if(doCutFlow) {//cutflow
				eSelector->reset();
				int code = eSelector->passCode();
				if(debug) std::cout << "Passes " << cutLevels[code] << std::endl;
				for(int i=0; i<(int)cutLevels.size(); i++) {
	    			if(debug) std::cout << "Checking against " << cutLevels[i] << " level" << std::endl;
					if(code >= i) h_cutflow[ich]->Fill(i);
		   			else break;
		 		}
		 		if(debug) eSelector->printEventStatus();
			}//cutflow

       		if(skimNtuple) {//skim
	   			if(eSelector->passes()) {
					if(runHitFit) {
						std::vector<int> selectedjets = eSelector->getGoodJets();
			    	 	if(priorityBtags>0) selectedjets = prioritizeBtags(eSelector->getGoodJets());
				     	hitfit->runHitFit(eSelector->getPrimaryLepton(),eSelector->getGoodJets(), jetisbtag);
						hitfit->fillHitFitInfo(HitFitInfo);
		   			}
	  				newtree->Fill();
		 		}//event passed
    	    }//skim
        	else if(runHitFit) {
				if(eSelector->passes()) {
	    			std::vector<int> selectedjets = eSelector->getGoodJets();
			  		if(priorityBtags>0) selectedjets = prioritizeBtags(eSelector->getGoodJets());
 		   			hitfit->runHitFit(eSelector->getPrimaryLepton(),eSelector->getGoodJets(), jetisbtag);
					hitfit->fillHitFitInfo(HitFitInfo);
	 			}
		 		newtree->Fill();
			}
     	}//loop over channels
	}//loop over entries

	if(debug) std::cout << "Finishing analyze\n";
}

// ------------ Rearrange jet list to prioritize b-tags  ------------
std::vector<int>
bpkHitFitAnalysis::prioritizeBtags(std::vector<int> jets)
{
  if(debug) {
    std::cout << "Prioritize b-tags\n";
    std::cout << "  Starting with jets:\n";
    for(int j=0; j<(int)jets.size(); j++) {
      double tagvalue=0;
      if(bTagAlgo.compare("CSV")==0)  tagvalue= JetInfo.CombinedSVBJetTags[jets[j]];
      else if(bTagAlgo.compare("TCHE")==0) tagvalue= JetInfo.TrackCountHiEffBJetTags[jets[j]];
      std::cout << "\tIndex: " << jets[j] 
		<< " pt=" << JetInfo.Pt[jets[j]] 
		<< " Btag(" << bTagAlgo << ")=" << tagvalue
		<< std::endl;
    }
  }
  if(jets.size()<5) {
    if(debug) std::cout << "  Less than 5 jets, don't bother\n";
    return jets;
  }

  std::pair<size_t,double> idxBtag;
  std::vector< std::pair<size_t,double> > jetsIdxBtag;

  for(int j=0; j<5; j++) {
    double tagvalue=0;
    if(bTagAlgo.compare("CSV")==0)  tagvalue= JetInfo.CombinedSVBJetTags[jets[j]];
    else if(bTagAlgo.compare("TCHE")==0) tagvalue= JetInfo.TrackCountHiEffBJetTags[jets[j]];
    idxBtag.first = jets[j];
    idxBtag.second = tagvalue;//JetInfo.TrackCountHiEffBJetTags[jets[j]];
    jetsIdxBtag.push_back(idxBtag);
  }

  std::stable_sort(jetsIdxBtag.begin(),jetsIdxBtag.end(),::IndexedQuantityGreaterThan<double>);

  std::vector<int> finalJets;
  for(int i=0; i<priorityBtags; i++) {
    int idx = jetsIdxBtag[i].first;
    finalJets.push_back(idx);
  }

  for(int j=0; j<(int)jets.size(); j++) {
    bool alreadyIn=false;
    for(int k=0; k<priorityBtags; k++) {
      if(jets[j] == finalJets[k]) {
	alreadyIn=true;
	break;
      }
    }
    if(alreadyIn) continue;
    finalJets.push_back(jets[j]);
  }

  if(debug) {
    std::cout << "  Finishing with jets:\n";
    for(int j=0; j<(int)finalJets.size(); j++) {
      double tagvalue=0;
      if(bTagAlgo.compare("CSV")==0)  tagvalue= JetInfo.CombinedSVBJetTags[jets[j]];
      else if(bTagAlgo.compare("TCHE")==0) tagvalue= JetInfo.TrackCountHiEffBJetTags[jets[j]];
      std::cout << "\tIndex: " << finalJets[j] 
		<< " pt=" << JetInfo.Pt[finalJets[j]] 
		<< " Btag=" << tagvalue
		<< std::endl;
    }
  }

  return finalJets;

}

bool
bpkHitFitAnalysis::jetIsBTagged(int ijet) {

  if(bTagAlgo.compare("CSV")==0)  return (JetInfo.CombinedSVBJetTags[ijet] > bTagCut);
  if(bTagAlgo.compare("TCHE")==0) return (JetInfo.TrackCountHiEffBJetTags[ijet] > bTagCut);

  std::cout << "WARNING: No valid b-tagger found, jet not tagged\n";
  return false;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
bpkHitFitAnalysis::endJob() 
{
  if(debug) std::cout << "Starting endJob\n";

  for(int ich=0; ich<(int)channels.size(); ich++) {//loop over channels
    if(doCutFlow) h_cutflow[ich]->Write();
  }

  if(debug) std::cout << "Cut flows written to file\n";

  if(skimNtuple || runHitFit) {
    newfile->mkdir("bprimeKit");
    newfile->cd("bprimeKit");

    TKey* key = (TKey*)newfile->GetListOfKeys()->FindObject("root");
    if(key) {
       key->Delete();
       delete key;
    }
    newtree->Write();
  }

  if(debug) std::cout << "New tree written\n";

  if(doCutFlow) {
    std::cout << "Event selection: cutflow\n";
    for(int ich=0; ich<(int)channels.size(); ich++) {//loop channels
      if (ELECTRON == channels[ich]) std::cout << "  Electron Channel:\n";
      else if(MUON == channels[ich]) std::cout << "  Muon Channel:\n";
      else                         std::cout << "  NA:\n"; //shouldn't happen.
      for(int i=0; i<h_cutflow[ich]->GetNbinsX(); i++) {
	if(i>=(int)cutLevels.size()) continue;
	std::cout << "\t" << cutLevels[i] << ":\t" << h_cutflow[ich]->GetBinContent(i+1) << std::endl;
      }
    }
  }

  if(debug) std::cout << "Finishing endJob\n";
  //BTAG_SF->Write();BTAG_EFF->Write();BTAG_SF_LIGHT->Write();BTAG_EFF_LIGHT->Write();BTAG_EFF_LIGHT_PASS->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(bpkHitFitAnalysis);
