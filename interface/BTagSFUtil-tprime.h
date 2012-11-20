#ifndef BTagSFUtil_h
#define BTagSFUtil_h

#include <Riostream.h>
#include "TRandom3.h"
#include "TMath.h"

#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class BTagSFUtil{

 public:
    
  BTagSFUtil( const edm::ParameterSet& );//std::string algo, int seed=0 );

  void setSeed(int seed);

  void readDB(const edm::EventSetup& iSetup, std::vector< std::pair<float,float> > jets); //pair is Et,Eta
  void modifyBTagsWithSF( bool& isBTagged, int pdgIdPart,float Btag_SF = 0.98, float Btag_eff = 1.0, float Bmistag_SF = 1.0, float Bmistag_eff = 1.0);
  void modifyBTagsWithSF(bool& isBTagged, int pdgIdPart, double Btag_SF, double Btag_eff);
  bool applySF(bool& isBTagged, float Btag_SF, float Btag_eff);

  void modifyBTagsWithSF(bool& isBTagged_l, bool& isBTagged_m, int pdgIdPart, 
			 float Btag_SF_l=0.95,  float Btag_SF_m=0.94,  float Btag_eff_l=1.,    float Btag_eff_m=1., 
			 float Bmistag_SF_l=1., float Bmistag_SF_m=1., float Bmistag_eff_l=1., float Bmistag_eff_m=1.);
  void applySF(bool& isBTagged_l, bool& isBTagged_m, float Btag_SF_l, float Btag_SF_m, float Btag_eff_l, float Btag_eff_m);

  float getSF(std::string s, int ijet);

  const std::vector<double> BtagCut_;
  const std::vector<double> BtagEff_;//value not stored in DB, so store it here

 private:

  void setupMaps();

  std::map<std::string,PerformanceResult::ResultType> measureMap;

  edm::ESHandle<BtagPerformance> perfH;
  std::vector<std::string> measureName;
  std::vector<std::string> measureType;

  const size_t maxNJets;// = 2;
  float Btageff_SF_;

  std::map< std::string, std::vector<float> > ScaleFactors;
  std::map< std::string, std::vector<float> > ScaleFactorsEff;

  TRandom3* rand_;
  std::string algo_;

  bool debug_;

};

#endif
