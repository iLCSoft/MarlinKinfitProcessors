#include "WW5CFit.h"
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
//#include <root/TLorentzVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
#include "ISRPhotonFitObject.h"
#include "MomentumConstraint.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "TextTracer.h"
#include "NewtonFitterGSL.h"
#include "FourJetPairing.h"
#include "MassConstraint.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


WW5CFit aWW5CFit ;

// function to define the jet energy resolution (in GeV)
double WW5CFit::JetEnergyResolution(double E)
{
  // examples here derived by Benjamin Hermberg from e+e- -> udsc:
  
  // 1) comparing jet-level to quark-level energies 
  //    (using MarlinReco/Analysis/RecoMCTruthLink/QuarkJetPairing.cc)
  // double result = std::sqrt(pow(0.6908,2)*(E)+(pow(0.02596,2)*pow(E,2))); 
  
  // 2) 120%/sqrt(E), gives best convergence of 5C fit on e+e- -> udsc
  double result = 1.2*std::sqrt(E);  
  return result;      
}

WW5CFit::WW5CFit() : Processor("WW5CFit") {
  
  // modify processor description
  _description = "WW5CFit does a 5C fit on 4 jet events (Px, Py, Pz, E, M12 = M34 (for all three permutations))" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "JetCollectionName" , 
			   "Name of the Jet collection"  ,
			   _jetcolName ,
			   std::string("Durham2Jets") ) ;
                           
  registerProcessorParameter( "ECM" ,
                              "Center-of-Mass Energy in GeV",
                              _ecm,
                              (float)500.);
                              
  registerProcessorParameter( "FitISR" ,
                              "0: Fit hypothesis without ISR   1: Fit hypothesis including ISR",
                              _fitISR,
                              (int) 1);
                              
  registerProcessorParameter( "ISRPzMax" ,
                              "Maximum possible energy for a single ISR photon",
                              _isrpzmax,
                              (float)225.);

  registerProcessorParameter( "fitter" ,
                              "0 = OPALFitter, 1 = NewFitter, 2 = NewtonFitter",
                              _ifitter,
                              (int)0);
                              
  registerProcessorParameter( "traceall" ,
                              "set true if every event should be traced",
                              _traceall,
                              (bool)false);
                              
  registerProcessorParameter( "ievttrace" ,
                              "number of individual event to be traced",
                              _ievttrace,
                              (int)0);
                              
}


void WW5CFit::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  b = (double) 0.00464564*( std::log(_ecm*_ecm*3814714.)-1. );
  //= 2*alpha/pi*( ln(s/m_e^2)-1 )
  ISRPzMaxB = std::pow((double)_isrpzmax,b);
  
}

void WW5CFit::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void WW5CFit::processEvent( LCEvent * evt ) { 

    
    message<DEBUG>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
  // this gets called for every event 
  // usually the working horse ...

#ifdef MARLIN_USE_AIDA
    
  // define a histogram pointer
  static AIDA::IHistogram1D* hRecWMassBest ;    
  static AIDA::IHistogram1D* hRecWMassAll ;    
  static AIDA::IHistogram1D* hRecWMassNoFitBest ;    
  static AIDA::IHistogram1D* hRecWMassNoFitAll ;    
  static AIDA::IHistogram1D* hTestWMassNoFitAll ;    
  static AIDA::IHistogram1D* hFitProbBest ;    
  static AIDA::IHistogram1D* hFitProbAll ;    
  static AIDA::IHistogram1D* hNItBest ;    
  static AIDA::IHistogram1D* hNItAll ;    
  static AIDA::IHistogram1D* hPhotonEnergy ;    
  static AIDA::IHistogram1D* hJetMass ;    
  static AIDA::IHistogram1D* hFitError;    
             
    message<DEBUG>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
  
  if( isFirstEvent() ) { 
    
    hRecWMassBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassBest", "M_W", 200, 0., 200. ) ; 
    hRecWMassAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassAll", "M_W", 200, 0., 200. ) ; 
    hRecWMassNoFitBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassNoFitBest", "M_W", 200, 0., 200. ) ; 
    hRecWMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassNoFitAll", "M_W", 200, 0., 200. ) ; 
    hTestWMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hTestWMassNoFitAll", "M_W", 200, 0., 200. ) ; 
    hFitProbBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ; 
    hFitProbAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProbAll", "fit probability", 100, 0., 1. ) ; 
    hNItBest = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hNItBest", "number of iterations", 200, 0., 200. ) ; 
    hNItAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hNItAll", "number of iterations", 200, 0., 200. ) ; 
    hPhotonEnergy = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPhotonEnergy", "ISR photon energy", 200, 0., 400. ) ; 
    hJetMass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hJetMass", "Jet Mass", 200, 0., 100. ) ; 
    if (_ifitter == 1) {
      hFitError = 
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 100, -0.5, 99.5 ) ; 
    }
    else {    
      hFitError = 
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 10, -0.5, 9.5 ) ; 
    }    

  }

#endif   
   
  message<DEBUG>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
  
  
  HepLorentzVector lvec;
  
    
  // fill histogram from LCIO data :

  //////////////////   JETS ///////////////////////////
   
     LCCollection* jetcol = evt->getCollection( _jetcolName ) ;
     if (jetcol != 0) {
  
       int nJETS = jetcol->getNumberOfElements()  ;
       message<DEBUG>( log() 
                      << " found " << nJETS
                      << " jets in event " << evt->getEventNumber() 
                      << "  in run "          << evt->getRunNumber() 
                      ) ;
                      
       if (nJETS != 4) return;               
                   
       float yminus = jetcol ->parameters().getFloatVal( "YMinus");              
       message<DEBUG>( log()  << " yminus = " << yminus ) ;
                     
       float yplus = jetcol ->parameters().getFloatVal( "YPlus");              
       message<DEBUG>( log()  << " yplus = " << yplus ) ;
                            
  // original fit objects - save for next permutation
       JetFitObject* j1 = 0;
       JetFitObject* j2 = 0;
       JetFitObject* j3 = 0;
       JetFitObject* j4 = 0;
        
       // angular resolutions - optimised for best convergence of 5C fit on e+e- -> udsc
       double errtheta = 0.1;   //   100mrad
       double errphi = 0.1;     //   100mrad
       
       
       ReconstructedParticle* jrps[4];
       
       for(int i=0; i< nJETS ; i++){
         
          ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i ) ) ;
          
               
          if (j) {
             jrps[i] = j;
             message<DEBUG>( log() 
                       << " found jet in event " << evt->getEventNumber() 
                       << "  in run "          << evt->getRunNumber() 
                       ) ;
             lvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy()); 
             if (i == 0) { 
               j1 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(),
                  JetEnergyResolution(lvec.e()), errtheta, errphi, lvec.m());
               j1->setName("Jet1");
               message<DEBUG>( log()  << " start four-vector of first  jet: " << *j1  ) ;
	     }
             else if (i == 1) { 
               j2 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(),
                  JetEnergyResolution(lvec.e()), errtheta, errphi, lvec.m());
               j2->setName("Jet2");
               message<DEBUG>( log() << " start four-vector of second  jet: " << *j2  ) ;
	     }
             else if (i == 2) {
               j3 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(),
                  JetEnergyResolution(lvec.e()), errtheta, errphi, lvec.m());
               j3->setName("Jet3");
               message<DEBUG>( log() << " start four-vector of third  jet: " << *j3  ) ;
	     }
             else if (i == 3) { 
               j4 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(),
                  JetEnergyResolution(lvec.e()), errtheta, errphi, lvec.m());
               j4->setName("Jet4");
               message<DEBUG>( log() << " start four-vector of forth  jet: " << *j4  ) ;
	     }
#ifdef MARLIN_USE_AIDA
             hJetMass->fill(j->getMass());
#endif           
          }
       }
#ifdef MARLIN_USE_AIDA
       double en, px, py, pz, mass; 
       if (jrps[0] && jrps[1] && jrps[2] && jrps[3]) {
         en = jrps[0]->getEnergy()     +jrps[1]->getEnergy();
         px = jrps[0]->getMomentum()[0]+jrps[1]->getMomentum()[0];
         py = jrps[0]->getMomentum()[1]+jrps[1]->getMomentum()[1];
         pz = jrps[0]->getMomentum()[2]+jrps[1]->getMomentum()[2];
         mass = en*en-px*px-py*py-pz*pz;
         if (mass >= 0) hTestWMassNoFitAll->fill( std::sqrt(mass) ) ;
         en = jrps[2]->getEnergy()     +jrps[3]->getEnergy();
         px = jrps[2]->getMomentum()[0]+jrps[3]->getMomentum()[0];
         py = jrps[2]->getMomentum()[1]+jrps[3]->getMomentum()[1];
         pz = jrps[2]->getMomentum()[2]+jrps[3]->getMomentum()[2];
         mass = en*en-px*px-py*py-pz*pz;
         if (mass >= 0) hTestWMassNoFitAll->fill( std::sqrt(mass) ) ;
         en = jrps[0]->getEnergy()     +jrps[2]->getEnergy();
         px = jrps[0]->getMomentum()[0]+jrps[2]->getMomentum()[0];
         py = jrps[0]->getMomentum()[1]+jrps[2]->getMomentum()[1];
         pz = jrps[0]->getMomentum()[2]+jrps[2]->getMomentum()[2];
         mass = en*en-px*px-py*py-pz*pz;
         if (mass >= 0) hTestWMassNoFitAll->fill( std::sqrt(mass) ) ;
         en = jrps[1]->getEnergy()     +jrps[3]->getEnergy();
         px = jrps[1]->getMomentum()[0]+jrps[3]->getMomentum()[0];
         py = jrps[1]->getMomentum()[1]+jrps[3]->getMomentum()[1];
         pz = jrps[1]->getMomentum()[2]+jrps[3]->getMomentum()[2];
         mass = en*en-px*px-py*py-pz*pz;
         if (mass >= 0) hTestWMassNoFitAll->fill( std::sqrt(mass) ) ;
         en = jrps[0]->getEnergy()     +jrps[3]->getEnergy();
         px = jrps[0]->getMomentum()[0]+jrps[3]->getMomentum()[0];
         py = jrps[0]->getMomentum()[1]+jrps[3]->getMomentum()[1];
         pz = jrps[0]->getMomentum()[2]+jrps[3]->getMomentum()[2];
         mass = en*en-px*px-py*py-pz*pz;
         if (mass >= 0) hTestWMassNoFitAll->fill( std::sqrt(mass) ) ;
         en = jrps[1]->getEnergy()     +jrps[2]->getEnergy();
         px = jrps[1]->getMomentum()[0]+jrps[2]->getMomentum()[0];
         py = jrps[1]->getMomentum()[1]+jrps[2]->getMomentum()[1];
         pz = jrps[1]->getMomentum()[2]+jrps[2]->getMomentum()[2];
         mass = en*en-px*px-py*py-pz*pz;
         if (mass >= 0) hTestWMassNoFitAll->fill( std::sqrt(mass) ) ;
       }  
#endif           
       
       const int NJETS = 4;
       message<DEBUG>( log()  << "*j1" << *j1  << "*j2" << *j2  << "*j3" << *j3  << "*j4" << *j4  ) ;
       
       // these get changed by the fit -> reset after each permutation!
       JetFitObject fitjets[NJETS] = {*j1, *j2, *j3, *j4};
       for (int i = 0; i < NJETS; ++i)
         message<DEBUG>( log()  << "fitjets[ " << i << "]: " << fitjets[i]  ) ;
 
       // these point allways to the fitjets array, which gets reset.
       JetFitObject *jets[NJETS];
       for (int i = 0; i < NJETS; ++i) jets[i] = &fitjets[i];
       for (int i = 0; i < NJETS; ++i)
         message<DEBUG>( log()  << "start four-vector of jets[ " << i << "]: " << *(jets[i])  ) ;

       FourJetPairing pairing (jets);
       JetFitObject *permutedjets[NJETS];

       double bestprob = 0.;
       int bestnit = 0;
       double bestmass1 = 0., bestmass2 = 0.;
       double beststartmass1 = 0., beststartmass2 = 0.;
       double startmass1 = 0., startmass2 = 0.;
       double bestphotonenergy = 0.;

       for (int iperm = 0; iperm < pairing.getNPerm(); iperm++) {

         message<DEBUG>( log() 
                       << " ================================================= "  
                       ) ;
         message<DEBUG>( log() 
                       << " iperm = " << iperm 
                       ) ;

         // important: (re-)set fitjets array!
         fitjets[0] = *j1;
         fitjets[1] = *j2;
         fitjets[2] = *j3;
         fitjets[3] = *j4;

         pairing.nextPermutation (permutedjets);
         for (int i = 0; i < NJETS; ++i)
            message<DEBUG>( log()  << "start four-vector of jet " << i << ": " << *(permutedjets[i])  ) ;
                              
         //MomentumConstraint pxc (1, 0);
         // crossing angle 14 mrad = 7/500
         MomentumConstraint pxc (0, 1, 0, 0, 7.0);
         pxc.setName("sum(p_x)");
         for (int i = 0; i < NJETS; ++i)
            pxc.addToFOList (*(permutedjets[i]));
        
         MomentumConstraint pyc (0, 0, 1);
         pyc.setName("sum(p_y)");
         for (int i = 0; i < NJETS; ++i)
            pyc.addToFOList (*(permutedjets[i]));
        
         MomentumConstraint pzc (0, 0, 0, 1);
         pzc.setName("sum(p_z)");
         for (int i = 0; i < NJETS; ++i)
            pzc.addToFOList (*(permutedjets[i]));
        
         message<DEBUG>( log() << "ECM = " << _ecm  ); 
	 MomentumConstraint ec(1, 0, 0, 0, _ecm);
         ec.setName("sum(E)");
         for (int i = 0; i < NJETS; ++i)
            ec.addToFOList (*(permutedjets[i]));
        
            message<DEBUG>( log()  << "Value of pxc before fit: " << pxc.getValue() ) ;
	    message<DEBUG>( log()  << "Value of pyc before fit: " << pyc.getValue() ) ;
	    message<DEBUG>( log()  << "Value of pzc before fit: " << pzc.getValue() ) ;
	    message<DEBUG>( log()  << "Value of ec before fit: " << ec.getValue() ) ;
                      

         // ISR Photon initialized with missing p_z
         ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);
        
	 if(_fitISR){
            message<DEBUG>( log()  << "start four-vector of ISR photon: " << *(photon) ) ;
                      
            pxc.addToFOList (*(photon));
            pyc.addToFOList (*(photon));
            pzc.addToFOList (*(photon));
            ec.addToFOList  (*(photon));
         }

         MassConstraint w(0.);
         w.addToFOList (*(permutedjets[0]), 1);
         w.addToFOList (*(permutedjets[1]), 1);
         w.addToFOList (*(permutedjets[2]), 2);
         w.addToFOList (*(permutedjets[3]), 2);

         startmass1 = w.getMass(1);
         startmass2 = w.getMass(2);
         
	 message<DEBUG>( log() << "start mass of W 1: " << startmass1 ) ;
	 message<DEBUG>( log() << "start mass of W 2: " << startmass2 ) ;
                      
#ifdef MARLIN_USE_AIDA
         hRecWMassNoFitAll->fill( 0.5*(startmass1 + startmass2) ) ;
         //hRecWMassNoFitAll->fill( startmass2 ) ;
#endif
         BaseFitter *pfitter;
         if (_ifitter == 1) {
           pfitter = new NewFitterGSL();
           if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<NewFitterGSL*>(pfitter))->setDebug (1);
         }
         else if (_ifitter == 2) {
           pfitter = new NewtonFitterGSL();
           if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<NewtonFitterGSL*>(pfitter))->setDebug (1);
         }
         else {
           // OPALFitter has no method setDebug !
           pfitter = new OPALFitterGSL();
         }
         BaseFitter &fitter = *pfitter;
  
         TextTracer tracer (std::cout);
         if (evt->getEventNumber()== _ievttrace || _traceall) fitter.setTracer (tracer);
         
         for (int i = 0; i < NJETS; ++i)
            fitter.addFitObject (*(permutedjets[i]));
         if(_fitISR){
            fitter.addFitObject (*(photon));
         }
         fitter.addConstraint (pxc);
         fitter.addConstraint (pyc);
         fitter.addConstraint (pzc);
         fitter.addConstraint (ec);
         fitter.addConstraint (w);

         double prob = fitter.fit();
         double chi2 = fitter.getChi2();
         int nit = fitter.getIterations();

         message<DEBUG>( log() << "fit probability = " << prob ) ;  
         message<DEBUG>( log() << "fit chi2 = " << chi2  ) ; 
         message<DEBUG>( log() << "error code: " << fitter.getError() ) ;
                      
         for (int i = 0; i < NJETS; ++i) {
            message<DEBUG>( log()  << "final four-vector of jet " << i << ": " << *(permutedjets[i])) ;
	 }
         if(_fitISR){
            message<DEBUG>( log()  << "final four-vector of ISR photon: " << *(photon) ) ;
	 }

         message<DEBUG>( log()  << "final mass of W 1: " << w.getMass(1) ) ;
	 message<DEBUG>( log()  << "final mass of W 2: " << w.getMass(2) ) ;
                      
         hFitError->fill( fitter.getError() ) ;
         if (fitter.getError() == 0) {
#ifdef MARLIN_USE_AIDA
           hFitProbAll->fill( prob ) ;
           hNItAll->fill( nit ) ;
           hRecWMassAll->fill( 0.5*(w.getMass(1)+w.getMass(2)) ) ;
           //hRecWMassAll->fill( w.getMass(2)) ;
#endif
//           if (prob > bestprob && w.getMass(1) > 50 && w.getMass(1) < 110) {
           if (prob > bestprob) {
             bestprob = prob;
             bestnit  = nit;
             bestmass1 = w.getMass(1);
             bestmass2 = w.getMass(2);
             beststartmass1 = startmass1;
             beststartmass2 = startmass2;
             bestphotonenergy = photon->getE();
           }
         }
         else {
         message<DEBUG>( log() << "FIT ERROR = " << fitter.getError() << ", not filling histograms!"  ) ;
	 }
         delete photon;
         message<DEBUG>( log() << "end permutation ") ;
       }

#ifdef MARLIN_USE_AIDA
       if (bestprob > 0) {
         hFitProbBest->fill( bestprob ) ;
         hNItBest->fill( bestnit ) ;
         hRecWMassBest->fill( 0.5*(bestmass1+bestmass2) ) ;
         //hRecWMassBest->fill( bestmass2 ) ;
         hRecWMassNoFitBest->fill( 0.5*(beststartmass1+beststartmass2) ) ;
         //hRecWMassNoFitBest->fill( beststartmass2 ) ;
         hPhotonEnergy->fill( _fitISR ? bestphotonenergy : 0. );
       } 
#endif

       delete j1;
       delete j2;
       delete j3;
       delete j4;
     }



  _nEvt ++ ;
}



void WW5CFit::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void WW5CFit::end(){ 

}

