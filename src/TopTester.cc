#include "TopTester.h"
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

#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
//#include "PConstraint.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "TextTracer.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"
#include "TopEventILC.h"

using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


TopTester aTopTester ;


TopTester::TopTester() : Processor("TopTester") {
  
  // modify processor description
  _description = "TopTester does a 7C fit on 6 jet events (Px, Py, Pz, E, M34 = M56 = 80.4 GeV, M134 = M256)" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "ECM" ,
                              "Center-of-Mass Energy in GeV",
                              _ecm,
                              (float)500.);
                              
  registerProcessorParameter( "ntoy" ,
                              "number of toy events",
                              _ntoy,
                              (int)200);
                              
  registerProcessorParameter( "semileptonic" ,
                              "set true if semi-leptonic ttbar events should be used",
                              _semileptonic,
                              (bool)false);
                              
  registerProcessorParameter( "leptonasjet" ,
                              "set true if lepton should be parametrised at JetFitObject",
                              _leptonasjet,
                              (bool)false);
                              
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
                              
                              
  topevent = new TopEventILC();   
}


void TopTester::init() { 

  // usually a good idea to
  printParameters() ;
  topevent->leptonic = _semileptonic;
  topevent->leptonasjet = _leptonasjet;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TopTester::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 

void TopTester::processEvent( LCEvent * evt ) { 

    
    message<ERROR>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
                      
  // this gets called for every event 
  // usually the working horse ...

#ifdef MARLIN_USE_AIDA
    
  // define a histogram pointer
  static AIDA::IHistogram1D* hRecTop1Mass;    
  static AIDA::IHistogram1D* hRecTop2Mass;    
  static AIDA::IHistogram1D* hRecTop1MassNoFitOK;    
  static AIDA::IHistogram1D* hRecTop1MassNoFitAll;    
  static AIDA::IHistogram1D* hRecTop2MassNoFitOK;    
  static AIDA::IHistogram1D* hRecTop2MassNoFitAll;    
  static AIDA::IHistogram1D* hRecTopMass;    
  static AIDA::IHistogram1D* hRecTopMassNoFitOK;    
  static AIDA::IHistogram1D* hRecTopMassNoFitAll;    
  static AIDA::IHistogram1D* hRecW1Mass;    
  static AIDA::IHistogram1D* hRecW2Mass;    
  static AIDA::IHistogram1D* hRecW1MassNoFitOK;    
  static AIDA::IHistogram1D* hRecW1MassNoFitAll;    
  static AIDA::IHistogram1D* hRecW2MassNoFitOK;    
  static AIDA::IHistogram1D* hRecW2MassNoFitAll;    
  static AIDA::IHistogram1D* hRecWMass;    
  static AIDA::IHistogram1D* hRecWMassNoFitOK;    
  static AIDA::IHistogram1D* hRecWMassNoFitAll;    
  static AIDA::IHistogram1D* hFitProb;    
  static AIDA::IHistogram1D* hNIt;    
  static AIDA::IHistogram1D* hFitError;
      
  static AIDA::IHistogram1D* hPullEJetOK;    
  static AIDA::IHistogram1D* hPullELepOK;    
  static AIDA::IHistogram1D* hPullENeuOK;    
  static AIDA::IHistogram1D* hPullThJetOK;    
  static AIDA::IHistogram1D* hPullThLepOK;    
  static AIDA::IHistogram1D* hPullThNeuOK;    
  static AIDA::IHistogram1D* hPullPhJetOK;    
  static AIDA::IHistogram1D* hPullPhLepOK;    
  static AIDA::IHistogram1D* hPullPhNeuOK;
      
  static AIDA::IHistogram1D* hPullEJetMea;    
  static AIDA::IHistogram1D* hPullELepMea;    
  static AIDA::IHistogram1D* hPullENeuMea;    
  static AIDA::IHistogram1D* hPullThJetMea;    
  static AIDA::IHistogram1D* hPullThLepMea;    
  static AIDA::IHistogram1D* hPullThNeuMea;    
  static AIDA::IHistogram1D* hPullPhJetMea;    
  static AIDA::IHistogram1D* hPullPhLepMea;    
  static AIDA::IHistogram1D* hPullPhNeuMea; 
     
  static AIDA::IHistogram1D* hPullEJetTrue;    
  static AIDA::IHistogram1D* hPullELepTrue;    
  static AIDA::IHistogram1D* hPullENeuTrue;    
  static AIDA::IHistogram1D* hPullThJetTrue;    
  static AIDA::IHistogram1D* hPullThLepTrue;    
  static AIDA::IHistogram1D* hPullThNeuTrue;    
  static AIDA::IHistogram1D* hPullPhJetTrue;    
  static AIDA::IHistogram1D* hPullPhLepTrue;    
  static AIDA::IHistogram1D* hPullPhNeuTrue;    
               
  if( isFirstEvent() ) { 
    
    hRecTopMass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMass", "M_top", 200, 125., 225. ) ; 

    hRecTop1Mass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop1Mass", "M_top1", 200, 125., 225. ) ; 

    hRecTop2Mass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop2Mass", "M_top2", 200, 125., 225. ) ; 

    hRecTopMassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassNoFitOK", "prefit average M_top, fit OK", 200, 125., 225. ) ; 

    hRecTopMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTopMassNoFitAll", "prefit average M_top, all", 200, 125., 225. ) ; 

    hRecTop1MassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop1MassNoFitOK", "prefit M_top1, fit OK", 200, 125., 225. ) ; 

    hRecTop1MassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop1MassNoFitAll", "prefit M_top1, all", 200, 125., 225. ) ; 

    hRecTop2MassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop2MassNoFitOK", "prefit M_top2, fit OK", 200, 125., 225. ) ; 

    hRecTop2MassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecTop2MassNoFitAll", "prefit M_top2, all", 200, 125., 225. ) ; 

    hRecWMass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass", "M_W", 200, 0., 200. ) ; 

    hRecW1Mass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecW1Mass", "M_W1", 200, 0., 200. ) ; 

    hRecW2Mass = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecW2Mass", "M_W2", 200, 0., 200. ) ; 

    hRecW1MassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass1NoFitOK", "prefit M_W1, fit OK", 200, 0., 200. ) ; 

    hRecW1MassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass1NoFitAll", "prefit M_W1, all", 200, 0., 200. ) ; 

    hRecW2MassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass2NoFitOK", "prefit M_W2, fit OK", 200, 0., 200. ) ; 

    hRecW2MassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMass2NoFitAll", "prefit M_W2, all", 200, 0., 200. ) ; 

    hRecWMassNoFitOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassNoFitOK", "prefit average M_W, fit OK", 200, 0., 200. ) ; 

    hRecWMassNoFitAll = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecWMassNoFitAll", "prefit average M_W, all", 200, 0., 200. ) ; 

    hFitProb = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ; 

    hNIt = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hNIt", "number of iterations", 201, -0.5, 200.5 ) ; 

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

    hPullEJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetOK", "pull of jet energy", 100, -5., 5. ) ;    
    hPullELepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepOK", "pull of lepton energy", 100, -5., 5. ) ;    
    hPullENeuOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullENeuOK", "pull of neutrino energy", 100, -5., 5. ) ;    
    hPullThJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetOK", "pull of jet theta", 100, -5., 5. ) ;    
    hPullThLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThLepOK", "pull of lepton theta", 100, -5., 5. ) ;    
    hPullThNeuOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThNeuOK", "pull of neutrino theta", 100, -5., 5. ) ;    
    hPullPhJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetOK", "pull of jet phi", 100, -5., 5. ) ;    
    hPullPhLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhLepOK", "pull of lepton phi", 100, -5., 5. ) ;    
    hPullPhNeuOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhNeuOK", "pull of neutrino phi", 100, -5., 5. ) ;    

    hPullEJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetMea", "pull of jet energy, #sigma_mea", 100, -5., 5. ) ;    
    hPullELepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepMea", "pull of lepton energy, #sigma_mea", 100, -5., 5. ) ;    
    hPullENeuMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullENeuMea", "pull of neutrino energy, #sigma_mea", 100, -5., 5. ) ;    
    hPullThJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetMea", "pull of jet theta, #sigma_mea", 100, -5., 5. ) ;    
    hPullThLepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThLepMea", "pull of lepton theta, #sigma_mea", 100, -5., 5. ) ;    
    hPullThNeuMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThNeuMea", "pull of neutrino theta, #sigma_mea", 100, -5., 5. ) ;    
    hPullPhJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetMea", "pull of jet phi, #sigma_mea", 100, -5., 5. ) ;    
    hPullPhLepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhLepMea", "pull of lepton phi, #sigma_mea", 100, -5., 5. ) ;    
    hPullPhNeuMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhNeuMea", "pull of neutrino phi, #sigma_mea", 100, -5., 5. ) ;    

    hPullEJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetTrue", "pull of jet energy vs true", 100, -5., 5. ) ;    
    hPullELepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepTrue", "pull of lepton energy vs true", 100, -5., 5. ) ;    
    hPullENeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullENeuTrue", "pull of neutrino energy vs true", 100, -5., 5. ) ;    
    hPullThJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetTrue", "pull of jet theta vs true", 100, -5., 5. ) ;    
    hPullThLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThLepTrue", "pull of lepton theta vs true", 100, -5., 5. ) ;    
    hPullThNeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThNeuTrue", "pull of neutrino theta vs true", 100, -5., 5. ) ;    
    hPullPhJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetTrue", "pull of jet phi vs true", 100, -5., 5. ) ;    
    hPullPhLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhLepTrue", "pull of lepton phi vs true", 100, -5., 5. ) ;    
    hPullPhNeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhNeuTrue", "pull of neutrino phi vs true", 100, -5., 5. ) ;    
  }

#endif   
      
   
  for (int ievt = 0; ievt < _ntoy; ievt++) { 
  
     message<MESSAGE>( log()  << "start to process toy event number " << ievt ) ;
     
     double startmassW1 = 0., startmassW2 = 0.;
     double startmasstop1 = 0., startmasstop2 = 0.;
     
     int debug = 0;
     if (ievt == _ievttrace || _traceall) debug = 4;
     topevent->setDebug (debug);
          
     topevent->genEvent();
       
     startmassW1 = topevent->getW1Mass();
     startmassW2 = topevent->getW2Mass();
     message<DEBUG>( log()  << "start mass of W 1: " << startmassW1 ) ;
     message<DEBUG>( log()  << "start mass of W 2: " << startmassW2 ) ;
     startmasstop1 = topevent->getTop1Mass();
     startmasstop2 = topevent->getTop2Mass();
     message<DEBUG>( log()  << "start mass of top 1: " << startmasstop1 ) ;
     message<DEBUG>( log()  << "start mass of top 2: " << startmasstop2 ) ;
                     
#ifdef MARLIN_USE_AIDA
     hRecTop1MassNoFitAll->fill( startmasstop1 ) ;
     hRecTop2MassNoFitAll->fill( startmasstop2 ) ;
     hRecTopMassNoFitAll->fill( 0.5*(startmasstop1+startmasstop2) ) ;
     hRecW1MassNoFitAll->fill( startmassW1 ) ;
     hRecW2MassNoFitAll->fill( startmassW2 ) ;
     hRecWMassNoFitAll->fill( 0.5*(startmassW1+startmassW2) ) ;
#endif

     BaseFitter *pfitter;
     if (_ifitter == 1) { 
       pfitter = new NewFitterGSL();
       (dynamic_cast<NewFitterGSL*>(pfitter))->setDebug (debug);
     }  
     else if (_ifitter == 2) { 
       pfitter = new NewtonFitterGSL();
       (dynamic_cast<NewtonFitterGSL*>(pfitter))->setDebug (debug);
     }  
     else {
       // OPALFitter has no method setDebug !
       pfitter = new OPALFitterGSL();
     }
     BaseFitter &fitter = *pfitter;
  
     TextTracer tracer (std::cout);
     if (ievt == _ievttrace || _traceall) fitter.setTracer (tracer);
           
     int ierr = topevent->fitEvent(fitter);
  
     double prob = fitter.getProbability();
     double chi2 = fitter.getChi2();
     int nit = fitter.getIterations();

     message<DEBUG>( log() << "fit probability = " << prob ) ;  
     message<DEBUG>( log() << "fit chi2 = " << chi2  ) ; 
     message<DEBUG>( log() << "error code: " << ierr ) ;
                                  
     message<DEBUG>( log()  << "final mass of W 1: " << topevent->getW1Mass() ) ;
     message<DEBUG>( log()  << "final mass of W 2: " << topevent->getW2Mass() ) ;
     message<DEBUG>( log()  << "final mass of top 1: " << topevent->getTop1Mass() ) ;
     message<DEBUG>( log()  << "final mass of top 2: " << topevent->getTop2Mass() ) ;
       
     bool usesigma_evt = true;
                  
#ifdef MARLIN_USE_AIDA
     hFitError->fill( ierr ) ;
     if (ierr == 0) {
       hFitProb->fill( prob ) ;
       hNIt->fill( nit ) ;
       hRecTop1MassNoFitOK->fill( startmasstop1 ) ;
       hRecTop2MassNoFitOK->fill( startmasstop2 ) ;
       hRecTopMassNoFitOK->fill( 0.5*(startmasstop1+startmasstop2) ) ;
       hRecW1MassNoFitOK->fill( startmassW1 ) ;
       hRecW2MassNoFitOK->fill( startmassW2 ) ;
       hRecWMassNoFitOK->fill( 0.5*(startmassW1+startmassW2) ) ;
       hRecTopMass->fill( 0.5*(topevent->getTopMass(1)+topevent->getTopMass(2)) ) ;
       hRecTop1Mass->fill( topevent->getTopMass(1) ) ;
       hRecTop2Mass->fill( topevent->getTopMass(2) ) ;
       hRecW1Mass->fill( topevent->getW1Mass() ) ;
       hRecW2Mass->fill( topevent->getW2Mass() ) ;
       hRecWMass->fill( 0.5*(topevent->getW1Mass()+topevent->getW2Mass()) ) ;
       
       for (int ifo = 0; ifo < 6; ifo++){
         double errfit, errmea, start, sigma; 
         double pull[3], pullmea[3], pulltrue[3];
         bool usesigma[3];
         for (int ipar = 0; ipar < 3; ipar++) {
           errfit = topevent->getFittedFitObject(ifo)->getError(ipar);  // should make difference for NewtonFitter
           errmea = topevent->getStartFitObject(ifo)->getError(ipar);  // SmearedFO are not fitted: original errors
           start = topevent->getStartFitObject(ifo)->getParam(ipar);  // SmearedFO are not fitted: original values
           sigma = errmea*errmea-errfit*errfit;
           if (sigma > 0) {
             sigma = sqrt(sigma);
             usesigma[ipar] = true;
           }
           else {
             message<WARNING>( log() << " SIGMA <= 0, taking only measured errors for pull for ifo = " << ifo 
                                     << " in evt " << ievt << ", errmea =  " << errmea << ", errfit = " << errfit ) ;
             usesigma[ipar] = false;
             usesigma_evt = false;
           }  
           pull[ipar] = (start - topevent->getFittedFitObject(ifo)->getParam(ipar))/sigma;
           pullmea[ipar] = (start - topevent->getFittedFitObject(ifo)->getParam(ipar))/errmea;
           pulltrue[ipar] = (start - topevent->getTrueFitObject(ifo)->getParam(ipar))/errmea;
         }  
         if ( !_semileptonic || ifo < 4 ) {
           if (usesigma[0]) hPullEJetOK->fill (pull[0]); hPullEJetMea->fill (pullmea[0]); hPullEJetTrue->fill (pulltrue[0]);
           if (usesigma[1]) hPullThJetOK->fill(pull[1]); hPullThJetMea->fill(pullmea[1]); hPullThJetTrue->fill(pulltrue[1]);
           if (usesigma[2]) hPullPhJetOK->fill(pull[2]); hPullPhJetMea->fill(pullmea[2]); hPullPhJetTrue->fill(pulltrue[2]);
         }
         else if ( _semileptonic && ifo ==  4 ) {
           if (usesigma[0]) hPullELepOK->fill (pull[0]); hPullELepMea->fill (pullmea[0]); hPullELepTrue->fill (pulltrue[0]);
           if (usesigma[1]) hPullThLepOK->fill(pull[1]); hPullThLepMea->fill(pullmea[1]); hPullThLepTrue->fill(pulltrue[1]);
           if (usesigma[2]) hPullPhLepOK->fill(pull[2]); hPullPhLepMea->fill(pullmea[2]); hPullPhLepTrue->fill(pulltrue[2]);
         }
         else if ( _semileptonic && ifo ==  5 ) {
           if (usesigma[0]) hPullENeuOK->fill (pull[0]); hPullENeuMea->fill (pullmea[0]); hPullENeuTrue->fill (pulltrue[0]);
           if (usesigma[1]) hPullThNeuOK->fill(pull[1]); hPullThNeuMea->fill(pullmea[1]); hPullThNeuTrue->fill(pulltrue[1]);
           if (usesigma[2]) hPullPhNeuOK->fill(pull[2]); hPullPhNeuMea->fill(pullmea[2]); hPullPhNeuTrue->fill(pulltrue[2]);
         }
       }
     }
#endif
     if (ierr > 0) message<WARNING>( log() << "FIT ERROR = " << ierr << " in toy event " << ievt ) ;
     
     if (!usesigma_evt) break;

     _nEvt ++ ;
     
   }   
}



void TopTester::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TopTester::end(){ 

 delete topevent;

}

