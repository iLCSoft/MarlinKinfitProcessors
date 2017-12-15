#include "DijetTester.h"
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
#include "DijetEventILC.h"

using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


DijetTester aDijetTester ;


DijetTester::DijetTester() : Processor("DijetTester") {
  
  // modify processor description
  _description = "DijetTester does a 4C fit on 2 jet events (Px, Py, Pz, E)" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "ECM" ,
                              "Center-of-Mass Energy in GeV",
                              _ecm,
                              (float)500.);
                              
  registerProcessorParameter( "ntoy" ,
                              "number of toy events",
                              _ntoy,
                              (int)200);
                              
  registerProcessorParameter( "leptonic" ,
                              "set true if di-lepton events should be used",
                              _leptonic,
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
                              
                              
  dijetevent = new DijetEventILC();   
}


void DijetTester::init() { 

  // usually a good idea to
  printParameters() ;
  dijetevent->leptonic = _leptonic;
  dijetevent->leptonasjet = _leptonasjet;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void DijetTester::processRunHeader( LCRunHeader* ) { 
  _nRun++ ;
} 

void DijetTester::processEvent( LCEvent * evt ) { 

    
    message<ERROR>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
                      
  // this gets called for every event 
  // usually the working horse ...

#ifdef MARLIN_USE_AIDA
    
  // define a histogram pointer
  static AIDA::IHistogram1D* hFitProb;    
  static AIDA::IHistogram1D* hNIt;    
  static AIDA::IHistogram1D* hFitError;
      
  static AIDA::IHistogram1D* hPxConstStart;
  static AIDA::IHistogram1D* hPyConstStart;
  static AIDA::IHistogram1D* hPzConstStart;
  static AIDA::IHistogram1D* hEConstStart;
  static AIDA::IHistogram1D* hMConstStart;
      
  static AIDA::IHistogram1D* hPxConstStop;
  static AIDA::IHistogram1D* hPyConstStop;
  static AIDA::IHistogram1D* hPzConstStop;
  static AIDA::IHistogram1D* hEConstStop;
  static AIDA::IHistogram1D* hMConstStop;
      
  static AIDA::IHistogram1D* hDistEJetOK;    
  static AIDA::IHistogram1D* hDistELepOK;    
  static AIDA::IHistogram1D* hDistThJetOK;    
  static AIDA::IHistogram1D* hDistThLepOK;    
  static AIDA::IHistogram1D* hDistPhJetOK;    
  static AIDA::IHistogram1D* hDistPhLepOK; 
     
  static AIDA::IHistogram1D* hPullEJetOK;    
  static AIDA::IHistogram1D* hPullELepOK;    
  static AIDA::IHistogram1D* hPullThJetOK;    
  static AIDA::IHistogram1D* hPullThLepOK;    
  static AIDA::IHistogram1D* hPullPhJetOK;    
  static AIDA::IHistogram1D* hPullPhLepOK;    
      
  static AIDA::IHistogram1D* hDistEJetMea;    
  static AIDA::IHistogram1D* hDistELepMea;    
  static AIDA::IHistogram1D* hDistThJetMea;    
  static AIDA::IHistogram1D* hDistThLepMea;    
  static AIDA::IHistogram1D* hDistPhJetMea;    
  static AIDA::IHistogram1D* hDistPhLepMea;
      
  static AIDA::IHistogram1D* hPullEJetMea;    
  static AIDA::IHistogram1D* hPullELepMea;    
  static AIDA::IHistogram1D* hPullThJetMea;    
  static AIDA::IHistogram1D* hPullThLepMea;    
  static AIDA::IHistogram1D* hPullPhJetMea;    
  static AIDA::IHistogram1D* hPullPhLepMea;    
     
  static AIDA::IHistogram1D* hDistEJetTrue;    
  static AIDA::IHistogram1D* hDistELepTrue;    
  static AIDA::IHistogram1D* hDistThJetTrue;    
  static AIDA::IHistogram1D* hDistThLepTrue;    
  static AIDA::IHistogram1D* hDistPhJetTrue;    
  static AIDA::IHistogram1D* hDistPhLepTrue;
      
  static AIDA::IHistogram1D* hPullEJetTrue;    
  static AIDA::IHistogram1D* hPullELepTrue;    
  static AIDA::IHistogram1D* hPullThJetTrue;    
  static AIDA::IHistogram1D* hPullThLepTrue;    
  static AIDA::IHistogram1D* hPullPhJetTrue;    
  static AIDA::IHistogram1D* hPullPhLepTrue;    
               
  if( isFirstEvent() ) { 
    
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

    hPxConstStart = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPxConstStart", "start value of Px constraint", 200, -50., 50. ) ;    
    hPyConstStart = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPyConstStart", "start value of Py constraint", 200, -50., 50. ) ;    
    hPzConstStart = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPzConstStart", "start value of Pz constraint", 200, -50., 50. ) ;    
    hEConstStart = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hEConstStart", "start value of E constraint", 200, -50., 50. ) ;    
    hMConstStart = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hMConstStart", "start value of M constraint", 200, -50., 50. ) ;    

    hPxConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPxConstStop", "final value of Px constraint", 200, -50., 50. ) ;    
    hPyConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPyConstStop", "final value of Py constraint", 200, -50., 50. ) ;    
    hPzConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPzConstStop", "final value of Pz constraint", 200, -50., 50. ) ;    
    hEConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hEConstStop", "final value of E constraint", 200, -50., 50. ) ;    
    hMConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hMConstStop", "final value of M constraint", 200, -50., 50. ) ;    

    hDistEJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistEJetOK", "fitted jet energy", 200, 150., 350. ) ;    
    hDistELepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistELepOK", "fitted lepton 1/pt", 200, 0.001, 1. ) ;    
    hDistThJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThJetOK", "fitted jet theta", 100, 0., 3.2 ) ;    
    hDistThLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThLepOK", "fitted lepton theta", 100, 0., 3.2 ) ;    
    hDistPhJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhJetOK", "fitted jet phi", 100, -3.2, 3.2 ) ;    
    hDistPhLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhLepOK", "fitted lepton phi", 100, -3.2, 3.2 ) ;    

    hPullEJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetOK", "pull of jet energy", 100, -5., 5. ) ;    
    hPullELepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepOK", "pull of lepton 1/pt", 100, -5., 5. ) ;    
    hPullThJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetOK", "pull of jet theta", 100, -5., 5. ) ;    
    hPullThLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThLepOK", "pull of lepton theta", 100, -5., 5. ) ;    
    hPullPhJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetOK", "pull of jet phi", 100, -5., 5. ) ;    
    hPullPhLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhLepOK", "pull of lepton phi", 100, -5., 5. ) ;    

    hDistEJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistEJetMea", "measured jet energy", 200, 150., 350. ) ;    
    hDistELepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistELepMea", "measured lepton 1/pt", 200, 0.001, 1. ) ;    
    hDistThJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThJetMea", "measured jet theta", 100, 0., 3.2 ) ;    
    hDistThLepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThLepMea", "measured lepton theta", 100, 0., 3.2 ) ;    
    hDistPhJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhJetMea", "measured jet phi", 100, -3.2, 3.2 ) ;    
    hDistPhLepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhLepMea", "measured lepton phi", 100, -3.2, 3.2 ) ;    

    hPullEJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetMea", "pull of jet energy, #sigma_mea", 100, -5., 5. ) ;    
    hPullELepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepMea", "pull of lepton 1/pt, #sigma_mea", 100, -5., 5. ) ;    
    hPullThJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetMea", "pull of jet theta, #sigma_mea", 100, -5., 5. ) ;    
    hPullThLepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThLepMea", "pull of lepton theta, #sigma_mea", 100, -5., 5. ) ;    
    hPullPhJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetMea", "pull of jet phi, #sigma_mea", 100, -5., 5. ) ;    
    hPullPhLepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhLepMea", "pull of lepton phi, #sigma_mea", 100, -5., 5. ) ;    

    hDistEJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistEJetTrue", "true jet energy", 200, 150., 350. ) ;    
    hDistELepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistELepTrue", "true lepton 1/pt", 200, 0.001, 1. ) ;    
    hDistThJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThJetTrue", "true jet theta", 100, 0., 3.2 ) ;    
    hDistThLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThLepTrue", "true lepton theta", 100, 0., 3.2 ) ;    
    hDistPhJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhJetTrue", "true jet phi", 100, -3.2, 3.2 ) ;    
    hDistPhLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhLepTrue", "true lepton phi", 100, -3.2, 3.2 ) ;    

    hPullEJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetTrue", "pull of jet energy vs true", 100, -5., 5. ) ;    
    hPullELepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepTrue", "pull of lepton 1/pt vs true", 100, -5., 5. ) ;    
    hPullThJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetTrue", "pull of jet theta vs true", 100, -5., 5. ) ;    
    hPullThLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThLepTrue", "pull of lepton theta vs true", 100, -5., 5. ) ;    
    hPullPhJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetTrue", "pull of jet phi vs true", 100, -5., 5. ) ;    
    hPullPhLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhLepTrue", "pull of lepton phi vs true", 100, -5., 5. ) ;    
  }

#endif   
      
   
  for (int ievt = 0; ievt < _ntoy; ievt++) { 
  
     message<MESSAGE>( log()  << "start to process toy event number " << ievt ) ;
     
     int debug = 0;
     if (ievt == _ievttrace || _traceall) debug = 4;
     dijetevent->setDebug (debug);
          
     dijetevent->genEvent();
                            
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
     
     // fill start values of constraints
#ifdef MARLIN_USE_AIDA
     hPxConstStart->fill (dijetevent->getPxConstraint().getValue()); 
     hPyConstStart->fill (dijetevent->getPyConstraint().getValue()); 
     hPzConstStart->fill (dijetevent->getPzConstraint().getValue()); 
     hEConstStart->fill (dijetevent->getEConstraint().getValue()); 
     hMConstStart->fill (dijetevent->getMassConstraint().getValue()); 
#endif
     
           
     int ierr = dijetevent->fitEvent(fitter);
  
     prob = fitter.getProbability();
     double chi2 = fitter.getChi2();
     nit = fitter.getIterations();

     message<DEBUG>( log() << "fit probability = " << prob ) ;  
     message<DEBUG>( log() << "fit chi2 = " << chi2  ) ; 
     message<DEBUG>( log() << "error code: " << ierr ) ;
                                  
     //bool usesigma_evt = true;
                  
#ifdef MARLIN_USE_AIDA
     hFitError->fill( ierr ) ;
     if (ierr == 0) {
       hFitProb->fill( prob ) ;
       hNIt->fill( nit ) ;
       // fill final values of constraints
       hPxConstStop->fill (dijetevent->getPxConstraint().getValue()); 
       hPyConstStop->fill (dijetevent->getPyConstraint().getValue()); 
       hPzConstStop->fill (dijetevent->getPzConstraint().getValue()); 
       hEConstStop->fill (dijetevent->getEConstraint().getValue()); 
       hMConstStop->fill (dijetevent->getMassConstraint().getValue()); 
       
       message<DEBUG>( log() << "looping over FOs " ) ;
       for (int ifo = 0; ifo < 2; ifo++){
         message<DEBUG>( log() << "ifo =  " << ifo ) ;
         double errfit, errmea, start, sigma; 
         double pull[3], pullmea[3], pulltrue[3];
         bool usesigma[3];
         for (int ipar = 0; ipar < 3; ipar++) {
           message<DEBUG>( log() << "ipar =  " << ipar ) ;
           errfit = dijetevent->getFittedFitObject(ifo)->getError(ipar);  // should make difference for NewtonFitter
           errmea = dijetevent->getStartFitObject(ifo)->getError(ipar);  // SmearedFO are not fitted: original errors
           start = dijetevent->getStartFitObject(ifo)->getParam(ipar);  // SmearedFO are not fitted: original values
           sigma = errmea*errmea-errfit*errfit;
           message<DEBUG>( log() << " sigma =  " << sigma << " for ifo " << ifo 
                                 << " in evt " << ievt << ", errmea =  " << errmea << ", errfit = " << errfit ) ;
           if (sigma > 0) {
             sigma = sqrt(sigma);
             usesigma[ipar] = true;
           }
           else {
             message<WARNING>( log() << " SIGMA <= 0, taking only measured errors for pull for ifo = " << ifo 
                                     << " in evt " << ievt << ", errmea =  " << errmea << ", errfit = " << errfit ) ;
             usesigma[ipar] = false;
             //usesigma_evt = false;
           }  
           pull[ipar] = (dijetevent->getFittedFitObject(ifo)->getParam(ipar) - start)/sigma;
           pullmea[ipar] = (dijetevent->getFittedFitObject(ifo)->getParam(ipar) - start)/errmea;
           pulltrue[ipar] = (start - dijetevent->getTrueFitObject(ifo)->getParam(ipar))/errmea;
           message<DEBUG>( log() << " pull =  " << pull[ipar] << " for ifo " << ifo 
                                 << " in evt " << ievt << ", pullmea =  " << pullmea[ipar] 
                                 << ", delta = " << start - dijetevent->getFittedFitObject(ifo)->getParam(ipar) ) ;
         }  
         if ( !_leptonic) {
           if (usesigma[0]) {hDistEJetOK->fill (dijetevent->getFittedFitObject(ifo)->getParam(0));  hPullEJetOK->fill (pull[0]); }
           if (usesigma[1]) {hDistThJetOK->fill (dijetevent->getFittedFitObject(ifo)->getParam(1)); hPullThJetOK->fill(pull[1]); }
           if (usesigma[2]) {hDistPhJetOK->fill (dijetevent->getFittedFitObject(ifo)->getParam(2)); hPullPhJetOK->fill(pull[2]); }
           hDistEJetMea->fill (dijetevent->getStartFitObject(ifo)->getParam(0));  hPullEJetMea->fill (pullmea[0]); 
           hDistThJetMea->fill (dijetevent->getStartFitObject(ifo)->getParam(1)); hPullThJetMea->fill (pullmea[1]); 
           hDistPhJetMea->fill (dijetevent->getStartFitObject(ifo)->getParam(2)); hPullPhJetMea->fill (pullmea[2]); 
           hDistEJetTrue->fill (dijetevent->getTrueFitObject(ifo)->getParam(0));  hPullEJetTrue->fill (pulltrue[0]);
           hDistThJetTrue->fill (dijetevent->getTrueFitObject(ifo)->getParam(1)); hPullThJetTrue->fill (pulltrue[1]);
           hDistPhJetTrue->fill (dijetevent->getTrueFitObject(ifo)->getParam(2)); hPullPhJetTrue->fill (pulltrue[2]);
         }
         else if ( _leptonic) {
           if (usesigma[0]) {hDistELepOK->fill (dijetevent->getFittedFitObject(ifo)->getParam(0));  hPullELepOK->fill (pull[0]); } 
           if (usesigma[1]) {hDistThLepOK->fill (dijetevent->getFittedFitObject(ifo)->getParam(1)); hPullThLepOK->fill(pull[1]); } 
           if (usesigma[2]) {hDistPhLepOK->fill (dijetevent->getFittedFitObject(ifo)->getParam(2)); hPullPhLepOK->fill(pull[2]); } 
           hDistELepMea->fill (dijetevent->getStartFitObject(ifo)->getParam(0));  hPullELepMea->fill (pullmea[0]); 
           hDistThLepMea->fill (dijetevent->getStartFitObject(ifo)->getParam(1)); hPullThLepMea->fill(pullmea[1]); 
           hDistPhLepMea->fill (dijetevent->getStartFitObject(ifo)->getParam(2)); hPullPhLepMea->fill(pullmea[2]); 
           hDistELepTrue->fill (dijetevent->getTrueFitObject(ifo)->getParam(0));  hPullELepTrue->fill (pulltrue[0]);
           hDistThLepTrue->fill (dijetevent->getTrueFitObject(ifo)->getParam(1)); hPullThLepTrue->fill(pulltrue[1]);
           hDistPhLepTrue->fill (dijetevent->getTrueFitObject(ifo)->getParam(2)); hPullPhLepTrue->fill(pulltrue[2]);
         }
       }
     }
#endif
     if (ierr > 0) message<WARNING>( log() << "FIT ERROR = " << ierr << " in toy event " << ievt ) ;
     
     //if (!usesigma_evt) break;

     _nEvt ++ ;
     
   }   
}



void DijetTester::check( LCEvent* ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void DijetTester::end(){ 

 delete dijetevent;

}

