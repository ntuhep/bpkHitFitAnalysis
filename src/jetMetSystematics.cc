#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "TVector2.h"
#include "TMath.h"

#include "MyAna/bprimeKit/interface/format.h"
#include "MyAna/bprimeKit/interface/bpkUtils.h"

#include "MyAna/bpkHitFitAnalysis/interface/jetMetSystematics.h"

const string jetMetSystematics::sType[]= {
   "None",
   "JES",
   "JER",
   "MET",
   "UNC"
};

jetMetSystematics::jetMetSystematics(const edm::ParameterSet& iConfig,
                                     EvtInfoBranches &evt, JetInfoBranches &jet, LepInfoBranches &lep)
   : _jecFile(iConfig.getUntrackedParameter<std::string>("JECfile")),
     stype(iConfig.getUntrackedParameter<std::string>("Type")),
     _evt(&evt),
     _jets(&jet),
     _leps(&lep),
     _jesUncType(iConfig.getUntrackedParameter<std::string>("JESUncType","Total")),
     _jerEta(iConfig.getUntrackedParameter< std::vector<double> >("JEREtaMax")),
     _jerNominal(iConfig.getUntrackedParameter< std::vector<double> >("JERNominal")),
     _jerSigmaSym(iConfig.getUntrackedParameter< std::vector<double> >("JERSigmaSym")),
     _jerSigmaNeg(iConfig.getUntrackedParameter< std::vector<double> >("JERSigmaNeg")),
     _jerSigmaPos(iConfig.getUntrackedParameter< std::vector<double> >("JERSigmaPos")),
     _jerMinGenJetPt(iConfig.getUntrackedParameter<double>("JERMinGenJetPt",15.)),
     _skipUNC(iConfig.getUntrackedParameter<bool>("SkipUNC",false)),
     _debug(iConfig.getUntrackedParameter<bool>("Debug",false))
{

   setType(stype);

   if(_skipUNC && UNC==_type) {
      std::cout << "jetMetSystematics: Will not allow adjustments for unclustered energy\n"; 
      return;
   }

   _scale = iConfig.getUntrackedParameter<double>("Scale",0.);
   if(_debug) std::cout << "jetMetSystematics: systematic scale = " << _scale << std::endl;

   _minDRljet = iConfig.getUntrackedParameter<double>("MinDRLepJet",0.3);

   if(JES == _type) {
      ifstream f(_jecFile.c_str());
      if(f) {
         JetCorrectorParameters* jcp = new JetCorrectorParameters(_jecFile,_jesUncType);
         _jesSigma = new JetCorrectionUncertainty(*jcp);
      }
      else {
         std::cout << "jetMetSystematics ERROR: Couldn't open JES File, can't continue!\n";
         assert(false);
      }
   }
   if(_debug) std::cout << "jetMetSystematics: initialized.\n";
}

void jetMetSystematics::scale() {
   if(NOSYS == _type) return;
   if(_skipUNC && UNC==_type) return;
   if(JER==_type && !_evt->McFlag) return; //no need to do any scaling on data
   
   if(JES == _type || JER == _type) { //JES or JER
      if(_debug) std::cout << "jetMetSystematics: Scaling jet energy\n";
      TVector2 metVec(_evt->PFMETx,_evt->PFMETy);
      
      for(int j=0; j<_jets->Size; j++) {
         metVec += TVector2(_jets->PtCorrRaw[j]*cos(_jets->Phi[j]), 
                            _jets->PtCorrRaw[j]*sin(_jets->Phi[j]) );
         
         double ptscale = 1.;
         if(JES == _type) {
            ptscale = jesScale(j,_scale);
         }
         if(JER == _type) {
            ptscale = jerScale(j,_scale);
         }
         
         _jets->Pt[j] *= ptscale;
         _jets->Et[j] *= ptscale;
         _jets->PtCorrRaw[j] *= ptscale;
         _jets->PtCorrL2[j]  *= ptscale; _jets->PtCorrL3[j]    *= ptscale;
         _jets->PtCorrL7g[j] *= ptscale; _jets->PtCorrL7uds[j] *= ptscale;
         _jets->PtCorrL7c[j] *= ptscale; _jets->PtCorrL7b[j]   *= ptscale;
         _jets->Px[j] *= ptscale; _jets->Py[j] *= ptscale; _jets->Pz[j] *= ptscale; _jets->Energy[j] *= ptscale;
         _jets->Mass[j] *= ptscale;
         
         metVec -= TVector2(_jets->PtCorrRaw[j]*cos(_jets->Phi[j]), 
                            _jets->PtCorrRaw[j]*sin(_jets->Phi[j]) );

         if(metVec.Mod() > 0.) {         
            _evt->PFMET    = metVec.Mod();
            _evt->PFMETPhi = TVector2::Phi_mpi_pi(metVec.Phi());
            _evt->PFMETx   = metVec.X();
            _evt->PFMETy   = metVec.Y();
         }
         else {
            std::cout << "\t*** jetMetSystematics: protection for metVec.Mod() = "<< metVec.Mod() 
                      << "\tmetVec.X() = "<<  metVec.X() << "\tmetVec.Y() = " << metVec.Y() << std::endl;
            std::cout << "\t*** Original MET: _evt->PFMET  = "<< _evt->PFMET   
                      << "\t_evt->PFMETx = "<<  _evt->PFMETx << "\t_evt->PFMETy = " << _evt->PFMETy << std::endl;
            _evt->PFMETx = 0.;
            _evt->PFMETy = 0.;
            _evt->PFMET = 0.;
            _evt->PFMETPhi = 0.;
         }

      }
   }//JES or JER
   
   if(MET == _type) { //MET
      if(_debug) std::cout << "jetMetSystematics: Scaling MET\n";
      double unc = sqrt(pow(0.599,2)+pow(0.563,2)*_evt->PFSumEt);
      double shift = _scale*unc;
      _evt->PFMET  += shift;
      _evt->PFMETx += shift*cos(_evt->PFMETPhi);
      _evt->PFMETy += shift*sin(_evt->PFMETPhi);
   }//MET
   
   if(UNC == _type) { //Unclustered energy
      if(_debug) std::cout << "jetMetSystematics: Scaling unclustered energy\n";
      TVector2 metVec(_evt->PFMETx,_evt->PFMETy);
      for(int j=0; j<_jets->Size; j++) {
         metVec += TVector2(_jets->PtCorrRaw[j]*cos(_jets->Phi[j]), 
                            _jets->PtCorrRaw[j]*sin(_jets->Phi[j]) );
      }
      //do i need to loop over leptons, or are those covered in pfjets?
      std::vector<int> lepIdx;
      for(int l=0; l<_leps->Size; l++) {
         bool notInJet=true;
         for(int j=0; j<_jets->Size; j++) {
            double dr = ::dR(_leps->Eta[l],_jets->Eta[j],_leps->Phi[l],_jets->Phi[j]);
            // if(_debug) {
            //    std::cout << "\tjet " << j << "(" << _jets->Eta[j] << "," << _jets->Phi[j] << ") ";
            //    std::cout << "lepton " << l << "(" << _leps->Eta[l] << "," << _leps->Phi[l] << ") ";
            //    std::cout << "dR = " << dr << std::endl;
            // }
            if(dr < _minDRljet) {
               notInJet=false;
               if(_debug) std::cout << "Lepton "<< l << " and jet " << j << " overlap with dR " << dr << ", less than minimum of " << _minDRljet << ". Already covered.\n";
               break;
            }
         }
         if(notInJet) {
            if(_debug) std::cout << "Got lep " << l << " which is not in a jet\n";
            lepIdx.push_back(l);
            metVec += TVector2(_leps->Px[l],_leps->Py[l]);
         }
      }

      metVec *= _scale;

      //do i need to loop over leptons, or are those covered in pfjets?
      for(int i=0; i<(int)lepIdx.size(); i++) {
         metVec -= TVector2(_leps->Px[lepIdx[i]],_leps->Py[lepIdx[i]]);
      }

      for(int j=0; j<_jets->Size; j++) {
         metVec -= TVector2(_jets->PtCorrRaw[j]*cos(_jets->Phi[j]), 
                            _jets->PtCorrRaw[j]*sin(_jets->Phi[j]) );
      }

      _evt->PFMETx = metVec.X();
      _evt->PFMETy = metVec.Y();
      _evt->PFMET = metVec.Mod();
      if(metVec.Mod() > 0.) 
         _evt->PFMETPhi = TVector2::Phi_mpi_pi(metVec.Phi());
   }//Unclustered energy

}

double jetMetSystematics::jesScale(const int index, const double sigma) {
   if(sigma==0.) return 1.;//why are we even here?
   if(_jets->Pt[index]<=0.) return 0.;//shouldn't happen, but just in case
   //https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources

   _jesSigma->setJetPt(_jets->Pt[index]);
   _jesSigma->setJetEta(_jets->Eta[index]);
   double uncert = _jesSigma->getUncertainty(sigma>0.);
   double ptscale = (sigma > 0.) ? (1. + uncert) : (1. - uncert);
   if(_debug) std::cout << "jetMetSystematics: JES scale for jet " << index << ", with Pt=" << _jets->Pt[index] << " is " << ptscale << std::endl; 

   return ptscale;
}

double jetMetSystematics::jerScale(const int index, const double sigma) {
   //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

   if(_jets->Pt[index]<=0.) return 1.;//not sure how this could happen, but just in case
   if(_jets->GenJetPt[index] < _jerMinGenJetPt) return 1.;
   int etaB = -1;
   for(int ieta=0; ieta<(int)_jerEta.size(); ieta++) {
      if(etaB>=0) break;
      if(fabs(_jets->GenJetEta[index]) < _jerEta[ieta]) etaB=ieta; 
   }
   if(etaB<0) etaB = (int)_jerEta.size(); //must be forward

   double sf = _jerNominal[etaB];
   if(sigma<0.)      sf -= sqrt(pow(_jerSigmaSym[etaB],2)+pow(_jerSigmaNeg[etaB],2));
   else if(sigma>0.) sf += sqrt(pow(_jerSigmaSym[etaB],2)+pow(_jerSigmaPos[etaB],2));

   double deltaPt = (_jets->Pt[index] - _jets->GenJetPt[index])*sf;
   //double ptscale = std::max(0.0, (_jets->PtCorrRaw[index] + deltaPt)/_jets->PtCorrRaw[index]);
   double ptscale = std::max(0.0, (_jets->GenJetPt[index] + deltaPt))/_jets->Pt[index];
   if(_debug) std::cout << "jetMetSystematics: JER scale for jet " << index << ", with Eta=" << _jets->Eta[index] << " is " << ptscale << std::endl; 
   if(ptscale<=0.) {
      std::cout << "jetMetSystematics WARNING: JER scale is " << ptscale << ". May have unintended consequences, so not scaling.\n";
      return 1.;
   }
   return ptscale;
}

void jetMetSystematics::setType(std::string s) {
   _type = NOSYS;

   for(Systematic i=NOSYS; i<nSystematics; i=static_cast<Systematic>(i+1)) {
      if(boost::iequals(s,sType[i])) _type = i;
   }
   
   std::cout << "jetMetSystematics: systematic type = " << sType[_type] << std::endl;
}

