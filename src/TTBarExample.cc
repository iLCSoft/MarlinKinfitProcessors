#include "TTBarExample.h"
#include <iostream>
#include <vector>
#include <string>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>
#endif

#include "UTIL/LCRelationNavigator.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
#include "NeutrinoFitObject.h"
#include "MomentumConstraint.h"
#include "OPALFitterGSL.h"
#include "TwoB4JPairing.h"
#include "FourJetPairing.h"
#include "MassConstraint.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


TTBarExample aTTBarExample ;


TTBarExample::TTBarExample() : Processor("TTBarExample") {
  
  // modify processor description
  _description = "TTBarExample does a 6C fit on 6jet ttbar events (Px, Py, Pz, E, M12 = M34 = M_W (for all permutations assuming the b-jets are tagged))" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "LightJetCollectionName" , 
			   "Name of the Light Jet collection"  ,
			   _lightjetcolName ,
			   std::string("DurhamLightJets") ) ;
                           
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "BJetCollectionName" , 
			   "Name of the B Jet collection"  ,
			   _bjetcolName ,
			   std::string("DurhamBJets") ) ;
                           
  registerProcessorParameter( "ECM" ,
                              "Center-of-Mass Energy",
                              _ecm,
                              (float)500.);

}


void TTBarExample::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TTBarExample::processRunHeader( LCRunHeader* ) { 
  _nRun++ ;
} 

void TTBarExample::processEvent( LCEvent * evt ) { 

    
  streamlog_out(MESSAGE) 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      << std::endl ;
  // this gets called for every event 
  // usually the working horse ...

#ifdef MARLIN_USE_AIDA
  
  // define a histogram pointer
  static AIDA::IHistogram1D* hRecTopMassBest ;    
  static AIDA::IHistogram1D* hRecTopMassAll ;    
  static AIDA::IHistogram1D* hRecTopMassNoFitBest ;    
  static AIDA::IHistogram1D* hRecTopMassNoFitAll ;    
  static AIDA::IHistogram1D* hFitProbBest ;    
  static AIDA::IHistogram1D* hFitProbAll ;    
             
  streamlog_out(MESSAGE) 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      << std::endl ;
  
  if( isFirstEvent() ) { 
    
    hRecTopMassBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassBest", "M_W", 200, 0., 200. ) ; 
    hRecTopMassAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassAll", "M_W", 200, 0., 200. ) ; 
    hRecTopMassNoFitBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassNoFitBest", "M_W", 200, 0., 200. ) ; 
    hRecTopMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassNoFitAll", "M_W", 200, 0., 200. ) ; 
    hFitProbBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ; 
    hFitProbAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProbAll", "fit probability", 100, 0., 1. ) ; 

  }

#endif
   
  streamlog_out(MESSAGE) 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      << std::endl ;
  
  
  HepLorentzVector lvec;
  
    
  // fill histogram from LCIO data :

  //////////////////   JETS ///////////////////////////
   
     LCCollection* lightjetcol = evt->getCollection( _lightjetcolName ) ;
     LCCollection* bjetcol = evt->getCollection( _bjetcolName ) ;
     if (lightjetcol != 0 && bjetcol != 0) {
  
       int nlightJETS = lightjetcol->getNumberOfElements()  ;
       streamlog_out(MESSAGE) 
                      << " found " << nlightJETS
                      << " light jets in event" << evt->getEventNumber() 
                      << "  in run "           << evt->getRunNumber() 
                      << std::endl ;
                   
       int nbJETS = bjetcol->getNumberOfElements()  ;
       streamlog_out(MESSAGE) 
                      << " found " << nbJETS
                      << " b jets in event" << evt->getEventNumber() 
                      << "  in run "           << evt->getRunNumber() 
                      << std::endl ;
                   
       
  // original fit objects - save for next permutation
       JetFitObject* j1 = 0;
       JetFitObject* j2 = 0;
       JetFitObject* j3 = 0;
       JetFitObject* j4 = 0;
       // these are assumed to be the b-jets
       JetFitObject* j5 = 0;
       JetFitObject* j6 = 0;
       
       double erre = 1.0;        //   100%/sqrt(E)
       double errtheta = 0.01;   //   10mrad
       double errphi = 0.01;     //   10mrad
       
       
       for(int i=0; i< nlightJETS ; i++){
         
          ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( lightjetcol->getElementAt( i ) ) ;
               
          if (j) {
             streamlog_out(MESSAGE) 
                       << " found jet in event " << evt->getEventNumber() 
                       << "  in run "          << evt->getRunNumber() 
                       << std::endl ;
             lvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy()); 
             erre *= std::sqrt(lvec.e());
             if (i == 0) {
               j1 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               streamlog_out(MESSAGE) 
                       << " start four-vector of first  jet: " << *j1 
                       << std::endl ;
             }
             else if (i == 1) {
               j2 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               streamlog_out(MESSAGE) 
                       << " start four-vector of second  jet: " << *j2 
                       << std::endl ;
             }
             else if (i == 2) {
               j3 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               streamlog_out(MESSAGE) 
                       << " start four-vector of third  jet: " << *j3 
                       << std::endl ;
             }
             else if (i == 3) {
               j4 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               streamlog_out(MESSAGE) 
                       << " start four-vector of forth  jet: " << *j4 
                       << std::endl ;
             }
           
          }
       }
       
       for(int i=0; i< nbJETS ; i++){
         
          ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( bjetcol->getElementAt( i ) ) ;
               
          if (j) {
             streamlog_out(MESSAGE) 
                       << " found b-jet in event " << evt->getEventNumber() 
                       << "  in run "          << evt->getRunNumber() 
                       << std::endl ;
             lvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy()); 
             erre *= std::sqrt(lvec.e());
             if (i == 0) {
               j5 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               streamlog_out(MESSAGE) 
                       << " start four-vector of first b-jet: " << *j5 
                       << std::endl ;
             }
             else if (i == 1) {
               j6 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(), erre, errtheta, errphi);
               streamlog_out(MESSAGE) 
                       << " start four-vector of second b-jet: " << *j6 
                       << std::endl ;
             }
           
          }
       }
       
       const int NJETS = 6;
       
       // these get changed by the fit -> reset after each permutation!
       JetFitObject fitjets[NJETS] = {*j1, *j2, *j3, *j4, *j5, *j6};
 
       // these point allways to the fitjets array, which gets reset.
       JetFitObject *jets[NJETS];
       for (int i = 0; i < NJETS; ++i) jets[i] = &fitjets[i];

       TwoB4JPairing pairing (jets);
       JetFitObject *permutedjets[NJETS];
       
       double bestprob = 0.;
       double bestmass1 = 0, bestmass2 = 0;
       double beststartmass1 = 0, beststartmass2 = 0;
       double startmass1 = 0, startmass2 = 0;
 
       for (int iperm = 0; iperm < pairing.getNPerm(); iperm++) {
     
         streamlog_out(MESSAGE) 
                       << " ================================================= "  
                       << std::endl ;
         streamlog_out(MESSAGE) 
                       << " iperm = " << iperm 
                       << std::endl ;

         // important: (re-)set fitjets array!
         fitjets[0] = *j1;
         fitjets[1] = *j2;
         fitjets[2] = *j3;
         fitjets[3] = *j4;
         fitjets[4] = *j5;
         fitjets[5] = *j6;

         pairing.nextPermutation (permutedjets);
         for (int i = 0; i < NJETS; ++i) {
            streamlog_out(MESSAGE) 
                       << "start four-vector of jet " << i << ": " << *(permutedjets[i])
                       << std::endl ;
         }              
        
         MomentumConstraint pxc (1, 0);
         for (int i = 0; i < NJETS; ++i)
            pxc.addToFOList (*(permutedjets[i]));
        
         MomentumConstraint pyc (0, 1);
         for (int i = 0; i < NJETS; ++i)
            pyc.addToFOList (*(permutedjets[i]));
        
         MomentumConstraint pzc (0, 0, 1);
         for (int i = 0; i < NJETS; ++i)
            pzc.addToFOList (*(permutedjets[i]));
            
         streamlog_out(MESSAGE) 
                   << "ECM = " << _ecm
                       << std::endl ;
         MomentumConstraint ec(0, 0, 0, 1, _ecm);
         for (int i = 0; i < NJETS; ++i)
            ec.addToFOList (*(permutedjets[i]));
        
         streamlog_out(MESSAGE) 
                    << "Value of pxc before fit: " << pxc.getValue()
                    << std::endl ;
         streamlog_out(MESSAGE) 
                    << "Value of pyc before fit: " << pyc.getValue()
                    << std::endl ;
         streamlog_out(MESSAGE) 
                    << "Value of pzc before fit: " << pzc.getValue()
                    << std::endl ;
         streamlog_out(MESSAGE) 
                    << "Value of ec before fit: " << ec.getValue()
                    << std::endl ;
  
         MassConstraint w1(80.4);
         MassConstraint w2(80.4);
         w1.addToFOList (*(permutedjets[0]), 1);
         w1.addToFOList (*(permutedjets[1]), 1);
         w2.addToFOList (*(permutedjets[2]), 1);
         w2.addToFOList (*(permutedjets[3]), 1);
         
         streamlog_out(MESSAGE) 
                       << "start mass of W 1: " << w1.getMass(1)
                       << std::endl ;
         streamlog_out(MESSAGE) 
                       << "start mass of W 2: " << w2.getMass(1)
                       << std::endl ;
                       
         // this is just a cheap way to monitor the resulting top mass:
         MassConstraint t1(175.);
         t1.addToFOList (*(permutedjets[0]), 1);
         t1.addToFOList (*(permutedjets[1]), 1);
         t1.addToFOList (*(permutedjets[4]), 1);
         MassConstraint t2(175.);
         t2.addToFOList (*(permutedjets[2]), 1);
         t2.addToFOList (*(permutedjets[3]), 1);
         t2.addToFOList (*(permutedjets[5]), 1);
         startmass1 = t1.getMass(1);
         startmass2 = t2.getMass(1);
         streamlog_out(MESSAGE) 
                       << "start mass of top 1: " << startmass1
                       << std::endl ;
         streamlog_out(MESSAGE) 
                       << "start mass of top 2: " << startmass2
                       << std::endl ;
#ifdef MARLIN_USE_AIDA
         hRecTopMassNoFitAll->fill( startmass1 ) ;
         hRecTopMassNoFitAll->fill( startmass2 ) ;
#endif        
       
         OPALFitterGSL fitter;
         for (int i = 0; i < NJETS; ++i)
            fitter.addFitObject (*(permutedjets[i]));
         fitter.addConstraint (pxc);
         fitter.addConstraint (pyc);
         fitter.addConstraint (pzc);
         fitter.addConstraint (ec);
         fitter.addConstraint (w1);
         fitter.addConstraint (w2);

         double prob = fitter.fit();
         streamlog_out(MESSAGE) 
                       << "fit probability = " << prob 
                       << std::endl ;
         streamlog_out(MESSAGE) 
                       << "error code: " << fitter.getError() 
                       << std::endl ;
         for (int i = 0; i < NJETS; ++i) {
            streamlog_out(MESSAGE) 
                       << "final four-vector of jet " << i << ": " << *(permutedjets[i])
                       << std::endl ;
         }              
         
         streamlog_out(MESSAGE) 
                       << "final mass of W 1: " << w1.getMass(1)
                       << std::endl ;
         streamlog_out(MESSAGE) 
                       << "final mass of W 2: " << w2.getMass(1)
                       << std::endl ;
         streamlog_out(MESSAGE) 
                       << "final mass of top 1: " << t1.getMass(1)
                       << std::endl ;
         streamlog_out(MESSAGE) 
                       << "final mass of top 2: " << t2.getMass(1)
                       << std::endl ;
         if (fitter.getError() == 0) {
#ifdef MARLIN_USE_AIDA
           hFitProbAll->fill( prob ) ;
           hRecTopMassAll->fill( t1.getMass(1)) ;
           hRecTopMassAll->fill( t2.getMass(1)) ;
#endif       
           if (prob > bestprob) {
             bestprob = prob;
             bestmass1 = t1.getMass(1);
             bestmass2 = t2.getMass(1);
             beststartmass1 = startmass1;
             beststartmass2 = startmass2;
           }
         }
         else {
         streamlog_out(MESSAGE) 
                       << "FIT ERROR = " << fitter.getError() << ", not filling histograms!"
                       << std::endl ;
         }

       }

#ifdef MARLIN_USE_AIDA
       if (bestprob > 0) {
         hFitProbBest->fill( bestprob ) ;
         hRecTopMassBest->fill( bestmass1 ) ;
         hRecTopMassBest->fill( bestmass2 ) ;
         hRecTopMassNoFitBest->fill( beststartmass1 ) ;
         hRecTopMassNoFitBest->fill( beststartmass2 ) ;

       } 
#endif       

       delete j1;
       delete j2;
       delete j3;
       delete j4;
       delete j5;
       delete j6;
     }
    


  _nEvt ++ ;
}



void TTBarExample::check( LCEvent* ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TTBarExample::end(){ 
  
}

