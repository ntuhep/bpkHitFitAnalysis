#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"

#include "MyAna/bprimeKit/interface/format.h"
#include "MyAna/bprimeKit/interface/bpkUtils.h"

#include "MyAna/bpkHitFitAnalysis/plugins/createSplitJets.h"

createSplitJets::createSplitJets( const edm::ParameterSet& iConfig ) : 
   maxEvents( iConfig.getUntrackedParameter<int>( "MaxEvents", -1 ) ),
   inFile( iConfig.getUntrackedParameter<std::string>( "InputFile" ) ),
   outFile( iConfig.getUntrackedParameter<std::string>( "OutputFile" ) ),
   _inStdJetCollection( iConfig.getUntrackedParameter<std::string>( "InputStandardJets" ) ),
   _inWJetCollection( iConfig.getUntrackedParameter<std::string>( "InputWJets" ) ),
   _outJetCollection( iConfig.getUntrackedParameter<std::string>( "OutputJets" ) ),
   _wMassLo( iConfig.getUntrackedParameter<double>( "MinWMass", 60. ) ),
   _wMassHi( iConfig.getUntrackedParameter<double>( "MaxWMass", 100. ) ),
   _massDropMax( iConfig.getUntrackedParameter<double>( "MaxMassDrop", -1. ) ),
   _matchDRMax( iConfig.getUntrackedParameter<double>( "MaxDRMatch", 0.4 ) ),
   _bTagCut( iConfig.getUntrackedParameter<double>( "BTagCut", 0.679 ) ),
   _debug( iConfig.getUntrackedParameter<bool>( "Debug", false ) )
{

}

createSplitJets::~createSplitJets()
{

   if( _debug ) { std::cout << "createSplitJets: destruction done\n"; }
}

   void
createSplitJets::beginJob()
{
   if( _debug ) { std::cout << "createSplitJets: starting beginJob\n"; }

   chain = new TChain( "bprimeKit/root" );
   chain->Add( inFile.c_str() );

   if( maxEvents < 0 || maxEvents > chain->GetEntries() )
   { maxEvents = chain->GetEntries(); }

   _evt.Register( chain );
   _inStdJets.Register( chain, _inStdJetCollection );
   _inWJets.Register( chain, _inWJetCollection );

   newfile = new TFile( outFile.c_str(), "recreate" );
   newtree = chain->CloneTree( 0 );
   _outJets.RegisterTree( newtree, _outJetCollection );
   if( _debug ) { std::cout << "createSplitJets: new tree ready\n"; }


   if( _debug ) { std::cout << "createSplitJets: finishing beginJob\n"; }
}

   void
createSplitJets::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
   if( _debug ) { std::cout << "createSplitJets: starting analyze\n"; }

   for( int entry = 0; entry < maxEvents; entry++ ) { //loop over entries
      if( ( entry % 1000 ) == 0 ) { printf( "Loading event %i of %i.\n", entry, maxEvents ); }
      chain->GetEntry( entry );
      memset( &_outJets, 0x00, sizeof( _outJets ) );

      if( _debug ) { std::cout << "Entry " << entry << ": Run " << _evt.RunNo << " Event " << _evt.EvtNo << std::endl; }

      //loop over wide jets: w-tag, put split jets in another collection
      //within wide jet loop, loop over ak5 jets: match (allow for multiple matches?), index if matched
      //loop over ak5 jets again, insert non-indexed jets into collection
      std::vector<int> matchedStdJets;
      matchedStdJets.clear();
      for( int iwjet = 0; iwjet < _inWJets.Size; iwjet++ ) { //loop over wide jets
         if( _debug ) { std::cout << "\twide jet " << iwjet << std::endl; }

         if( !shouldSplit( iwjet ) ) { continue; }

         int match = -1;
         bool hadBmatch = false;
         for( int jstdjet = 0; jstdjet < _inStdJets.Size; jstdjet++ ) { //loop over standard jets
            double dR = ::dR( _inWJets.Eta[iwjet], _inStdJets.Eta[jstdjet], _inWJets.Phi[iwjet], _inStdJets.Phi[jstdjet] );
            if( _debug ) { std::cout << "\t\tstandard jet " << jstdjet << " dR=" << dR << std::endl; }
            if( dR < _matchDRMax ) {
               if( _inStdJets.pfCombinedSecondaryVertexV2BJetTags[jstdjet] > _bTagCut ) {
                  if( _debug ) { std::cout << "\t\tGot match to b-tagged jet. Won't use\n"; }
                  hadBmatch = true;
                  continue; //don't want to actually match to b-tagged jets
               }
               match = jstdjet;
               if( _debug ) {
                  for( unsigned i = 0; i < matchedStdJets.size(); i++ )
                     if( jstdjet == matchedStdJets[i] )
                     { std::cout << "This jet was already matched to another wide jet\n"; }
               }
               matchedStdJets.push_back( jstdjet );
               break;
            }
         }//loop over standard jets
         if( match == -1 ) {
            if( !hadBmatch ) { std::cout << "No match found? Can that be right?\n"; }
            continue;
         }

         if( _debug ) { std::cout << "\tGot match at " << match << ". Now splitting jets\n"; }
         insertSplitJets( iwjet, match, _evt.McFlag == 1 );
      }//loop over wide jets

      if( _debug ) { std::cout << "Done looping over wide jets, now to standard jets\n"; }

      for( int istdjet = 0; istdjet < _inStdJets.Size; istdjet++ ) { //loop over standard jets
         if( _debug ) { std::cout << "\tStandard jet " << istdjet; }
         bool wasMatched = false;
         for( unsigned jmatched = 0; jmatched < matchedStdJets.size(); jmatched++ ) { //loop over matched jets
            if( matchedStdJets[jmatched] == istdjet ) {
               wasMatched = true;
               break;
            }
         }
         if( _debug ) { std::cout << ( wasMatched ? " was" : " wasn't" ) << " matched\n"; }
         if( wasMatched ) { continue; } //already added split jets, don't want this too.

         if( _debug ) { std::cout << "\tInserting jet\n"; }
         insertJetCopy( istdjet );
      }//loop over standard jets

      newtree->Fill();
   }//loop over entries

   if( _debug ) { std::cout << "createSplitJets: finishing analyze\n"; }
}

   void
createSplitJets::endJob()
{
   if( _debug ) { std::cout << "createSplitJets: starting endJob\n"; }

   newfile->mkdir( "bprimeKit" );
   newfile->cd( "bprimeKit" );
   newtree->Write();
   if( _debug ) { std::cout << "createSplitJets: new tree written\n"; }

   if( _debug ) { std::cout << "createSplitJets: finishing endJob\n"; }
}

   bool
createSplitJets::shouldSplit( int jet )
{
   int sSubjets = _inWJets.SubjetsIdxStart[jet];
   int nSubjets = _inWJets.NSubjets[jet];

   if( _debug ) {
      std::cout << "\tSplit jet " << jet << "? mass,sSub,nSub,1,2,..="
         << _inWJets.Mass[jet] << ","
         << sSubjets << "," << nSubjets << ",";
      for( int i = sSubjets; i < sSubjets + nSubjets; i++ )
      { std::cout << _inWJets.SubjetMass->at( i ) << ","; }
      std::cout << std::endl;
   }

   if( nSubjets != 2 ) {
      std::cout << "WARNING: don't have 2 subjets. Getting out of here.\n";
      return false;
   }

   double mass = _inWJets.Mass[jet];
   if( mass < _wMassLo ) { return false; }
   if( mass > _wMassHi ) { return false; }

   if( _massDropMax < 0. ) { return true; } //let <0 mean don't cut on mass drop
   double massDrop = ( _inWJets.SubjetMass->at( sSubjets ) > _inWJets.SubjetMass->at( sSubjets + 1 ) ) ? _inWJets.SubjetMass->at( sSubjets ) / mass : _inWJets.SubjetMass->at( sSubjets + 1 ) / mass;
   if( massDrop > _massDropMax ) { return false; }

   return true;
}

   void
createSplitJets::insertSplitJets( int idx_w, int idx_std, bool isMC )
{
   //insert the daughter jets of wide jet at idx_w, matched to standard jet at idx_std
   if( _debug ) std::cout << "\t\tInsert split jets from wide jet " << idx_w
      << " matched to standard jet " << idx_std
         << " isMC? " << isMC
         << " to index " << _outJets.Size << std::endl;

   const int sSubjets = _inWJets.SubjetsIdxStart[idx_w];
   const int nSubjets = _inWJets.NSubjets[idx_w];

   //split wide jet
   TLorentzVector jW, jS[nSubjets];
   jW.SetPtEtaPhiM( 0, 0, 0, 0 );
   for( int iSub = 0; iSub < nSubjets; iSub++ ) {
      jS[iSub].SetPtEtaPhiM( _inWJets.SubjetPt->at( sSubjets + iSub ),
            _inWJets.SubjetEta->at( sSubjets + iSub ),
            _inWJets.SubjetPhi->at( sSubjets + iSub ),
            _inWJets.SubjetMass->at( sSubjets + iSub ) );
      jW += jS[iSub];
   }
   double rescale = _inWJets.Pt[idx_w] / jW.Pt();

   for( int iSub = 0; iSub < nSubjets; iSub++ ) {
      jS[iSub].SetPtEtaPhiM( _inWJets.SubjetPt->at( sSubjets + iSub )*rescale,
            _inWJets.SubjetEta->at( sSubjets + iSub ),
            _inWJets.SubjetPhi->at( sSubjets + iSub ),
            _inWJets.SubjetMass->at( sSubjets + iSub )*rescale );
   }

   //insert split jets into collection
   for( int iS = 0; iS < nSubjets; iS++ ) { //loop over 2 daughters
      _outJets.Index[_outJets.Size] = _outJets.Size;
      _outJets.Unc[_outJets.Size] = 1 + iS;  //use to indicate jet is split jet

      _outJets.Et  [_outJets.Size] = _inWJets.SubjetEt->at( sSubjets + iS ) * rescale;
      _outJets.Pt  [_outJets.Size] = _inWJets.SubjetPt->at( sSubjets + iS ) * rescale;
      _outJets.Phi [_outJets.Size] = _inWJets.SubjetPhi->at( sSubjets + iS );
      _outJets.Eta [_outJets.Size] = _inWJets.SubjetEta->at( sSubjets + iS );
      _outJets.Mass[_outJets.Size] = _inWJets.SubjetMass->at( sSubjets + iS );

      _outJets.Px    [_outJets.Size] = jS[iS].Px();
      _outJets.Py    [_outJets.Size] = jS[iS].Py();
      _outJets.Pz    [_outJets.Size] = jS[iS].Pz();
      _outJets.Energy[_outJets.Size] = jS[iS].E();

      //copy variables from the Std jet 'parent', nothing better available
      //only JetIDLOOSE used in jet selection
      _outJets.NTracks      [_outJets.Size] = _inStdJets.NTracks      [idx_std];
      _outJets.JetIDLOOSE   [_outJets.Size] = _inStdJets.JetIDLOOSE   [idx_std];
      _outJets.NConstituents[_outJets.Size] = _inStdJets.NConstituents[idx_std];
      _outJets.NCH          [_outJets.Size] = _inStdJets.NCH          [idx_std];
      _outJets.CEF          [_outJets.Size] = _inStdJets.CEF          [idx_std];
      _outJets.NHF          [_outJets.Size] = _inStdJets.NHF          [idx_std];
      _outJets.NEF          [_outJets.Size] = _inStdJets.NEF          [idx_std];
      _outJets.CHF          [_outJets.Size] = _inStdJets.CHF          [idx_std];
      //_outJets.CombinedSVBJetTags[_outJets.Size] = _inStdJets.CombinedSVBJetTags[idx_std];
      //_outJets.Area         [_outJets.Size] = _inStdJets.Area         [idx_std];
      _outJets.pfCombinedSecondaryVertexV2BJetTags[_outJets.Size] = _inWJets.SubjetCombinedSVBJetTags->at( sSubjets + iS );
      _outJets.Area         [_outJets.Size] = _inWJets.SubjetArea->at( sSubjets + iS );

      if( isMC ) {
         _outJets.GenFlavor[_outJets.Size] = _inStdJets.GenFlavor[idx_std];
         _outJets.GenPdgID [_outJets.Size] = _inStdJets.GenPdgID [idx_std];

         //   approx
         _outJets.GenJetEta[_outJets.Size] = _outJets.Eta[_outJets.Size];
         _outJets.GenJetPhi[_outJets.Size] = _outJets.Phi[_outJets.Size];
         _outJets.GenJetPt[_outJets.Size]  = _inWJets.GenJetPt[idx_std] * ( _inWJets.SubjetPt->at( sSubjets + iS ) / _inWJets.Pt[idx_w] );
      }

      double scf = _outJets.Pt[_outJets.Size] / _inWJets.Pt[idx_w];

      //_outJets.PtCorrRaw  [_outJets.Size] = _inWJets.PtCorrRaw  [idx_w]*scf;
      _outJets.PtCorrRaw  [_outJets.Size] = _inWJets.SubjetPtUncorr->at( sSubjets + iS );
      _outJets.PtCorrL2   [_outJets.Size] = _inWJets.PtCorrL2   [idx_w] * scf;
      _outJets.PtCorrL3   [_outJets.Size] = _inWJets.PtCorrL3   [idx_w] * scf;
      _outJets.PtCorrL7g  [_outJets.Size] = _inWJets.PtCorrL7g  [idx_w] * scf;
      _outJets.PtCorrL7uds[_outJets.Size] = _inWJets.PtCorrL7uds[idx_w] * scf;
      _outJets.PtCorrL7c  [_outJets.Size] = _inWJets.PtCorrL7c  [idx_w] * scf;
      _outJets.PtCorrL7b  [_outJets.Size] = _inWJets.PtCorrL7b  [idx_w] * scf;

      _outJets.Size++;
   }//loop over 2 daughters

}

   void
createSplitJets::insertJetCopy( int idx )
{
   //insert copy of standard jet into new collection
   if( _debug ) std::cout << "\t\tInsert copy of standard jet " << idx
      << " to index " << _outJets.Size << std::endl;

   _outJets.Index[_outJets.Size] = _outJets.Size;
   _outJets.NTracks[_outJets.Size] = _inStdJets.NTracks[idx];
   _outJets.Et[_outJets.Size] = _inStdJets.Et[idx];
   _outJets.Pt[_outJets.Size] = _inStdJets.Pt[idx];
   _outJets.Unc[_outJets.Size] = 0; //set to indicate from standard jet
   _outJets.Eta[_outJets.Size] = _inStdJets.Eta[idx];
   _outJets.Phi[_outJets.Size] = _inStdJets.Phi[idx];
   _outJets.JetIDLOOSE[_outJets.Size] = _inStdJets.JetIDLOOSE[idx];
   _outJets.JetCharge[_outJets.Size] = _inStdJets.JetCharge[idx];
   _outJets.NConstituents[_outJets.Size] = _inStdJets.NConstituents[idx];
   _outJets.NCH[_outJets.Size] = _inStdJets.NCH[idx];
   _outJets.CEF[_outJets.Size] = _inStdJets.CEF[idx];
   _outJets.NHF[_outJets.Size] = _inStdJets.NHF[idx];
   _outJets.NEF[_outJets.Size] = _inStdJets.NEF[idx];
   _outJets.CHF[_outJets.Size] = _inStdJets.CHF[idx];
   _outJets.PtCorrRaw[_outJets.Size] = _inStdJets.PtCorrRaw[idx];
   _outJets.PtCorrL2[_outJets.Size] = _inStdJets.PtCorrL2[idx];
   _outJets.PtCorrL3[_outJets.Size] = _inStdJets.PtCorrL3[idx];
   _outJets.PtCorrL7g[_outJets.Size] = _inStdJets.PtCorrL7g[idx];
   _outJets.PtCorrL7uds[_outJets.Size] = _inStdJets.PtCorrL7uds[idx];
   _outJets.PtCorrL7c[_outJets.Size] = _inStdJets.PtCorrL7c[idx];
   _outJets.PtCorrL7b[_outJets.Size] = _inStdJets.PtCorrL7b[idx];

   _outJets.combinedSecondaryVertexBJetTags              [_outJets.Size] = _inStdJets.combinedSecondaryVertexBJetTags             [idx] ;
   _outJets.pfJetBProbabilityBJetTags                    [_outJets.Size] = _inStdJets.pfJetBProbabilityBJetTags                   [idx] ;
   _outJets.pfJetProbabilityBJetTags                     [_outJets.Size] = _inStdJets.pfJetProbabilityBJetTags                    [idx] ;
   _outJets.pfTrackCountingHighPurBJetTags               [_outJets.Size] = _inStdJets.pfTrackCountingHighPurBJetTags              [idx] ;
   _outJets.pfTrackCountingHighEffBJetTags               [_outJets.Size] = _inStdJets.pfTrackCountingHighEffBJetTags              [idx] ;
   _outJets.pfSimpleSecondaryVertexHighEffBJetTags       [_outJets.Size] = _inStdJets.pfSimpleSecondaryVertexHighEffBJetTags      [idx] ;
   _outJets.pfSimpleSecondaryVertexHighPurBJetTags       [_outJets.Size] = _inStdJets.pfSimpleSecondaryVertexHighPurBJetTags      [idx] ;
   _outJets.pfCombinedSecondaryVertexV2BJetTags          [_outJets.Size] = _inStdJets.pfCombinedSecondaryVertexV2BJetTags         [idx] ;
   _outJets.pfCombinedInclusiveSecondaryVertexV2BJetTags [_outJets.Size] = _inStdJets.pfCombinedInclusiveSecondaryVertexV2BJetTags[idx] ;
   _outJets.pfCombinedSecondaryVertexSoftLeptonBJetTags  [_outJets.Size] = _inStdJets.pfCombinedSecondaryVertexSoftLeptonBJetTags [idx] ;
   _outJets.pfCombinedMVABJetTags                        [_outJets.Size] = _inStdJets.pfCombinedMVABJetTags                       [idx] ;


   _outJets.GenJetPt[_outJets.Size] = _inStdJets.GenJetPt[idx];
   _outJets.GenJetEta[_outJets.Size] = _inStdJets.GenJetEta[idx];
   _outJets.GenJetPhi[_outJets.Size] = _inStdJets.GenJetPhi[idx];
   _outJets.GenPt[_outJets.Size] = _inStdJets.GenPt[idx];
   _outJets.GenEta[_outJets.Size] = _inStdJets.GenEta[idx];
   _outJets.GenPhi[_outJets.Size] = _inStdJets.GenPhi[idx];
   _outJets.GenPdgID[_outJets.Size] = _inStdJets.GenPdgID[idx];
   _outJets.GenFlavor[_outJets.Size] = _inStdJets.GenFlavor[idx];
   _outJets.GenMCTag[_outJets.Size] = _inStdJets.GenMCTag[idx];
   _outJets.Px[_outJets.Size] = _inStdJets.Px[idx];
   _outJets.Py[_outJets.Size] = _inStdJets.Py[idx];
   _outJets.Pz[_outJets.Size] = _inStdJets.Pz[idx];
   _outJets.Energy[_outJets.Size] = _inStdJets.Energy[idx];

   _outJets.Mass[_outJets.Size] = _inStdJets.Mass[idx];
   _outJets.Area[_outJets.Size] = _inStdJets.Area[idx];

   //    _outJets.MassD1[_outJets.Size] = _inStdJets.MassD1[idx];
   //    _outJets.MassD2[_outJets.Size] = _inStdJets.MassD2[idx];
   //    _outJets.PtD1[_outJets.Size] = _inStdJets.PtD1[idx];
   //    _outJets.PtD2[_outJets.Size] = _inStdJets.PtD2[idx];
   //    _outJets.EtD1[_outJets.Size] = _inStdJets.EtD1[idx];
   //    _outJets.EtD2[_outJets.Size] = _inStdJets.EtD2[idx];
   //    _outJets.EtaD1[_outJets.Size] = _inStdJets.EtaD1[idx];
   //    _outJets.EtaD2[_outJets.Size] = _inStdJets.EtaD2[idx];
   //    _outJets.PhiD1[_outJets.Size] = _inStdJets.PhiD1[idx];
   //    _outJets.PhiD2[_outJets.Size] = _inStdJets.PhiD2[idx];

   _outJets.Size++;
}

DEFINE_FWK_MODULE( createSplitJets );
