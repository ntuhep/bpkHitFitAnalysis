#include "MyAna/bpkHitFitAnalysis/interface/BTagSFUtil-tprime.h"

BTagSFUtil::BTagSFUtil( const edm::ParameterSet& iConfig )//std::string algo, int seed ) 
  : BtagCut_(iConfig.getUntrackedParameter< std::vector<double> >("BTagCuts")),
    BtagEff_(iConfig.getUntrackedParameter< std::vector<double> >("BTagEffs")),
    maxNJets(iConfig.getUntrackedParameter<int>("MaxNJets",8)),
    Btageff_SF_(iConfig.getUntrackedParameter<double>("BTagEffSF",1.)),
    algo_(iConfig.getUntrackedParameter<std::string>("BTagAlgorithm","CSV")),
    debug_(iConfig.getUntrackedParameter<bool>("Debug",false))
{
  int seed = iConfig.getUntrackedParameter<int>("Seed",0);
  rand_ = new TRandom3(seed);

  setupMaps();

}

void BTagSFUtil::setSeed(int seed) {
  rand_->SetSeed(seed);
}

void BTagSFUtil::readDB(const edm::EventSetup& iSetup, std::vector< std::pair<float,float> > jets) {

  ScaleFactors.clear();
  ScaleFactorsEff.clear();

  for( size_t iMeasure = 0; iMeasure < measureName.size(); iMeasure++ ) {
    if(debug_) std::cout << "iMeasure = " << iMeasure << " Name = " <<  measureName[iMeasure] << " Type = " << measureType[iMeasure] << " Map = " << measureMap[measureType[iMeasure]] << std::endl;
    //Setup our measurement
    iSetup.get<BTagPerformanceRecord>().get( measureName[ iMeasure ],perfH);
    const BtagPerformance & perf = *(perfH.product());
    BinningPointByMap measurePoint;

    float scaler = 1.;
    std::string suffix = "";
    if ( measureType[ iMeasure ] == "BTAGLEFF" || measureType[ iMeasure ] == "BTAGBEFF" ) {
      if(debug_) std::cout << "efficiency\n";
      suffix = "eff";
    }
    else {
      if(debug_) std::cout << "scale factor\n";
      scaler = Btageff_SF_;
    }

    std::vector<float> jetsf;
    size_t nJets = (jets.size()>maxNJets) ? maxNJets : jets.size();
    for(size_t ijet=0; ijet<nJets; ijet++) {
      measurePoint.reset();
      measurePoint.insert(BinningVariables::JetEt, (jets[ijet]).first );                // pass in the et of the jet
      measurePoint.insert(BinningVariables::JetAbsEta, fabs(  (jets[ijet]).second ) );  // pass in the absolute eta of the jet

      jetsf.push_back(scaler*perf.getResult( measureMap[measureType[iMeasure]], measurePoint ));
      if(debug_) std::cout << "jet " << ijet << " result = " << jetsf[ijet] << std::endl;
    // Extract the mistag eff value
    }

    if ( measureType[ iMeasure ] == "BTAGLEFF" || measureType[ iMeasure ] == "BTAGBEFF" ) {
      ScaleFactorsEff[ measureName[iMeasure] + suffix ] = jetsf;
    }
    else {
      ScaleFactors[ measureName[iMeasure] ] = jetsf;
    }

  }


}

void BTagSFUtil::modifyBTagsWithSF(bool& isBTagged, int pdgIdPart, float Btag_SF, float Btag_eff, float Bmistag_SF, float Bmistag_eff) {

  bool newBTag = isBTagged;

  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 ||  abs( pdgIdPart ) == 4) { 

    double bctag_eff = Btag_eff;
    if(abs(pdgIdPart)==4) bctag_eff = Btag_eff/5; // take ctag eff as one 5th of Btag eff
    newBTag = applySF(isBTagged, Btag_SF, bctag_eff);

//     float coin = rand_->Uniform(1.);
    
//     //cout << "Sf / coin = "<< Btageff_SF << " " << coin  << endl;

//     if( isBTagged ){ 
//       if( coin > Btageff_SF ) {isBTagged=false;} //turn medium off 
//     }
    

    // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    newBTag = applySF(isBTagged, Bmistag_SF, Bmistag_eff);

//     // no need to upgrade if is medium tagged
//     if( isBTagged ) return;

//     float mistagPercent = ( Btagmistag_SF*Btagmistag_eff - Btagmistag_eff ) / ( 1. - Btagmistag_eff );    
//     float coin = rand_->Uniform(1.);

//     // for light quarks, the jet has to be upgraded:

//     if( coin < mistagPercent ) {isBTagged = true;}

  } //if light quark

  isBTagged = newBTag;

}
void BTagSFUtil::modifyBTagsWithSF(bool& isBTagged, int pdgIdPart, double Btag_SF, double Btag_eff) {
	bool newBTag = isBTagged;
	if (pdgIdPart != 0) newBTag = applySF(isBTagged, Btag_SF, Btag_eff);
	if (debug_ && pdgIdPart == 0) std::cout << "Jet flavor is 0, take it as Light Quark" << std::endl;
	if (newBTag != isBTagged) std::cout << "Modified tag SF/EFF/pdgId : " << Btag_SF << "," << Btag_eff << "," << pdgIdPart << std::endl;
	isBTagged = newBTag;
}

bool BTagSFUtil::applySF(bool& isBTagged, float Btag_SF, float Btag_eff){
  
  bool newBTag = isBTagged;

  if (Btag_SF == 1) return newBTag; //no correction needed 

  //throw die
  float coin = rand_->Uniform(1.);    
  
  if(Btag_SF > 1){  // use this if SF>1

    if( !isBTagged ) {

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );

      //upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }

  }else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {newBTag = false;}

  }
  return newBTag;
}

void BTagSFUtil::modifyBTagsWithSF(bool& isBTagged_l, bool& isBTagged_m, int pdgIdPart, 
				   float Btag_SF_l,  float Btag_SF_m,  float Btag_eff_l,    float Btag_eff_m,
				   float Bmistag_SF_l, float Bmistag_SF_m, float Bmistag_eff_l, float Bmistag_eff_m) {

  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 ||  abs( pdgIdPart ) == 4) { 

    double bctag_eff_l = Btag_eff_l;
    double bctag_eff_m = Btag_eff_m;
    if(abs(pdgIdPart)==4) {
      bctag_eff_l = Btag_eff_l/5; // take ctag eff as one 5th of Btag eff
      bctag_eff_m = Btag_eff_m/5; // take ctag eff as one 5th of Btag eff
    }
    applySF(isBTagged_l, isBTagged_m, Btag_SF_l, Btag_SF_m, bctag_eff_l, bctag_eff_m);

  // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    applySF(isBTagged_l, isBTagged_m, Bmistag_SF_l, Bmistag_SF_m, Bmistag_eff_l, Bmistag_eff_m);
    
  } //if light quark

}

void BTagSFUtil::applySF(bool& isBTagged_l, bool& isBTagged_m, float Btag_SF_l, float Btag_SF_m, float Btag_eff_l, float Btag_eff_m){
  
  bool newBTag_l = isBTagged_l;
  bool newBTag_m = isBTagged_m;

  if (Btag_SF_l == 1 && Btag_SF_m == 1) return; //no correction needed 

  //throw die
  float coin = rand_->Uniform(1.);
  
  if(Btag_SF_m > 1){  // use this if SF>1

    if( !isBTagged_m ) {

      float btag_sf = 1.;
      float btag_eff = 1.;

      if(isBTagged_l) {
	btag_sf = Btag_SF_m;
	btag_eff = Btag_eff_m;
      }
      else {
	btag_sf = Btag_SF_l;
	btag_eff = Btag_eff_l;
      }

      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - btag_sf) / (1.0 - (btag_sf/btag_eff) );

      //upgrade to tagged
      if( coin < mistagPercent ) {
	if(!isBTagged_l) newBTag_l = true;
	else if(!isBTagged_m) newBTag_m = true;
      }
    }

  } else{  // use this if SF<1
      
    //downgrade tagged to untagged
    if( isBTagged_m ) {
      if(coin > Btag_SF_m ) {newBTag_m = false; newBTag_l = false;}
    }
    else if( isBTagged_l && !isBTagged_m ) {
      if(coin > Btag_SF_l) {newBTag_l = false;}
    }

  }

  isBTagged_l = newBTag_l;
  isBTagged_m = newBTag_m;
}

void BTagSFUtil::setupMaps() {
  // This is needed for the DB

  measureMap["BTAGBEFF"]=PerformanceResult::BTAGBEFF;
  measureMap["BTAGBERR"]=PerformanceResult::BTAGBERR;
  measureMap["BTAGCEFF"]=PerformanceResult::BTAGCEFF;
  measureMap["BTAGCERR"]=PerformanceResult::BTAGCERR;
  measureMap["BTAGLEFF"]=PerformanceResult::BTAGLEFF;
  measureMap["BTAGLERR"]=PerformanceResult::BTAGLERR;
  measureMap["BTAGNBEFF"]=PerformanceResult::BTAGNBEFF;
  measureMap["BTAGNBERR"]=PerformanceResult::BTAGNBERR;
  measureMap["BTAGBEFFCORR"]=PerformanceResult::BTAGBEFFCORR;
  measureMap["BTAGBERRCORR"]=PerformanceResult::BTAGBERRCORR;
  measureMap["BTAGCEFFCORR"]=PerformanceResult::BTAGCEFFCORR;
  measureMap["BTAGCERRCORR"]=PerformanceResult::BTAGCERRCORR;
  measureMap["BTAGLEFFCORR"]=PerformanceResult::BTAGLEFFCORR;
  measureMap["BTAGLERRCORR"]=PerformanceResult::BTAGLERRCORR;
  measureMap["BTAGNBEFFCORR"]=PerformanceResult::BTAGNBEFFCORR;
  measureMap["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  measureMap["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  measureMap["MUEFF"]=PerformanceResult::MUEFF;
  measureMap["MUERR"]=PerformanceResult::MUERR;
  measureMap["MUFAKE"]=PerformanceResult::MUFAKE; 
  measureMap["MUEFAKE"]=PerformanceResult::MUEFAKE;
      
  // Define which Btag and Mistag algorithm you want to use. These are not user defined and need to be exact
  measureName.push_back("MISTAG" + algo_ + "M");
  measureName.push_back("BTAG" + algo_ + "M");
  measureName.push_back("MISTAG" + algo_ + "L");
  measureName.push_back("BTAG" + algo_ + "L");
  measureName.push_back("MISTAG" + algo_ + "M");
  measureName.push_back("MISTAG" + algo_ + "L");
//   measureName.push_back("BTAG" + algo_ + "M");
//   measureName.push_back("BTAG" + algo_ + "L");

  // Tell DB you want the SF. These are not user defined and need to be exact
  measureType.push_back("BTAGLEFFCORR");
  measureType.push_back("BTAGBEFFCORR");
  measureType.push_back("BTAGLEFFCORR");
  measureType.push_back("BTAGBEFFCORR");
  measureType.push_back("BTAGLEFF");
  measureType.push_back("BTAGLEFF");
//   measureType.push_back("BTAGBEFF");
//   measureType.push_back("BTAGBEFF");

}

float BTagSFUtil::getSF(std::string s, int ijet) {

//   if(ijet < 0 || ijet>=(int)(ScaleFactors[s]).size()) {
//     std::cout << "WARNING: BTagSFUtil::getSF - " << s << " invalid jet index " << ijet << "greater than size " << (ScaleFactors[s]).size() << std::endl;
//     return -999.;
//   }

  if(debug_) std::cout << "getSF: name " << s << ", jet" << ijet;// << std::endl;

  if(s.find("eff")!=std::string::npos) {//eff
    if(debug_) std::cout << " Eff=" << (ScaleFactorsEff[s]).at(ijet) << std::endl;
    return (ScaleFactorsEff[s]).at(ijet);
  }
  else {
    if(debug_) std::cout << " SF=" << (ScaleFactors[s]).at(ijet) << std::endl;
    return (ScaleFactors[s]).at(ijet);
  }

}
