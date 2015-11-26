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
      
  static AIDA::IHistogram1D* hPxConstStart;
  static AIDA::IHistogram1D* hPyConstStart;
  static AIDA::IHistogram1D* hPzConstStart;
  static AIDA::IHistogram1D* hEConstStart;
  static AIDA::IHistogram1D* hMW1ConstStart;
  static AIDA::IHistogram1D* hMW2ConstStart;
  static AIDA::IHistogram1D* hMConstStart;
      
  static AIDA::IHistogram1D* hPxConstStop;
  static AIDA::IHistogram1D* hPyConstStop;
  static AIDA::IHistogram1D* hPzConstStop;
  static AIDA::IHistogram1D* hEConstStop;
  static AIDA::IHistogram1D* hMW1ConstStop;
  static AIDA::IHistogram1D* hMW2ConstStop;
  static AIDA::IHistogram1D* hMConstStop;
      
  static AIDA::IHistogram1D* hDistEJetOK;    
  static AIDA::IHistogram1D* hDistELepOK;    
  static AIDA::IHistogram1D* hDistENeuOK;    
  static AIDA::IHistogram1D* hDistThJetOK;    
  static AIDA::IHistogram1D* hDistThLepOK;    
  static AIDA::IHistogram1D* hDistThNeuOK;    
  static AIDA::IHistogram1D* hDistPhJetOK;    
  static AIDA::IHistogram1D* hDistPhLepOK; 
  static AIDA::IHistogram1D* hDistPhNeuOK; 
     
  static AIDA::IHistogram1D* hPullEJetOK;    
  static AIDA::IHistogram1D* hPullELepOK;    
  static AIDA::IHistogram1D* hPullENeuOK;    
  static AIDA::IHistogram1D* hPullThJetOK;    
  static AIDA::IHistogram1D* hPullThLepOK;    
  static AIDA::IHistogram1D* hPullThNeuOK;    
  static AIDA::IHistogram1D* hPullPhJetOK;    
  static AIDA::IHistogram1D* hPullPhLepOK;    
  static AIDA::IHistogram1D* hPullPhNeuOK;
      
  static AIDA::IHistogram1D* hDistEJetMea;    
  static AIDA::IHistogram1D* hDistELepMea;    
  static AIDA::IHistogram1D* hDistENeuMea;    
  static AIDA::IHistogram1D* hDistThJetMea;    
  static AIDA::IHistogram1D* hDistThLepMea;    
  static AIDA::IHistogram1D* hDistThNeuMea;    
  static AIDA::IHistogram1D* hDistPhJetMea;    
  static AIDA::IHistogram1D* hDistPhLepMea;
  static AIDA::IHistogram1D* hDistPhNeuMea;
      
  static AIDA::IHistogram1D* hPullEJetMea;    
  static AIDA::IHistogram1D* hPullELepMea;    
  static AIDA::IHistogram1D* hPullENeuMea;    
  static AIDA::IHistogram1D* hPullThJetMea;    
  static AIDA::IHistogram1D* hPullThLepMea;    
  static AIDA::IHistogram1D* hPullThNeuMea;    
  static AIDA::IHistogram1D* hPullPhJetMea;    
  static AIDA::IHistogram1D* hPullPhLepMea;    
  static AIDA::IHistogram1D* hPullPhNeuMea; 
     
  static AIDA::IHistogram1D* hDistEJetTrue;    
  static AIDA::IHistogram1D* hDistELepTrue;    
  static AIDA::IHistogram1D* hDistENeuTrue;    
  static AIDA::IHistogram1D* hDistThJetTrue;    
  static AIDA::IHistogram1D* hDistThLepTrue;    
  static AIDA::IHistogram1D* hDistThNeuTrue;    
  static AIDA::IHistogram1D* hDistPhJetTrue;    
  static AIDA::IHistogram1D* hDistPhLepTrue;
  static AIDA::IHistogram1D* hDistPhNeuTrue;
      
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
    hMW1ConstStart = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hMW1ConstStart", "start value of MW1 constraint", 200, -50., 50. ) ;    
    hMW2ConstStart = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hMW2ConstStart", "start value of MW2 constraint", 200, -50., 50. ) ;    
    hMConstStart = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hMConstStart", "start value of top mass constraint", 200, -50., 50. ) ;    

    hPxConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPxConstStop", "final value of Px constraint", 200, -5., 5. ) ;    
    hPyConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPyConstStop", "final value of Py constraint", 200, -5., 5. ) ;    
    hPzConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPzConstStop", "final value of Pz constraint", 200, -5., 5. ) ;    
    hEConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hEConstStop", "final value of E constraint", 200, -5., 5. ) ;    
    hMW1ConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hMW1ConstStop", "final value of MW1 constraint", 200, -5., 5. ) ;    
    hMW2ConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hMW2ConstStop", "final value of MW2 constraint", 200, -5., 5. ) ;    
    hMConstStop = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hMConstStop", "final value of top mass constraint", 200, -5., 5. ) ;    

    hDistEJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistEJetOK", "fitted - true jet energy", 200, -50., 50. ) ;    
    hDistELepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistELepOK", "fitted - true lepton 1/pt", 200,  -0.001, 0.001 ) ;    
    hDistENeuOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistENeuOK", "fitted - true neutrino energy", 200,  -50., 50. ) ;    
    hDistThJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThJetOK", "fitted - true jet theta", 200, -1., 1.) ;    
    hDistThLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThLepOK", "fitted - true lepton theta", 200,  -0.01, 0.01) ;    
    hDistThNeuOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThNeuOK", "fitted - true neutrino theta", 200,  -1., 1.) ;    
    hDistPhJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhJetOK", "fitted - true jet phi", 200, -1., 1. ) ;    
    hDistPhLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhLepOK", "fitted - true lepton phi", 200, -0.01, 0.01 ) ;    
    hDistPhNeuOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhNeuOK", "fitted - true neutrino phi", 200, -1., 1. ) ;    

    hPullEJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetOK", "pull of jet energy", 100, -5., 5. ) ;    
    hPullELepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepOK", "pull of lepton 1/pt", 100, -5., 5. ) ;    
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

    hDistEJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistEJetMea", "measured - true jet energy", 200, -50., 50. ) ;    
    hDistELepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistELepMea", "measured - true lepton 1/pt", 200, -0.001, 0.001 ) ;    
    hDistENeuMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistENeuMea", "measured - true neutrino 1/pt", 200, -50., 50. ) ;    
    hDistThJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThJetMea", "measured - true jet theta", 200, -1., 1. ) ;    
    hDistThLepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThLepMea", "measured - true lepton theta", 200, -0.01, 0.01 ) ;    
    hDistThNeuMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThNeuMea", "measured - true neutrino theta", 200, -1., 1. ) ;    
    hDistPhJetMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhJetMea", "measured - true jet phi", 200, -1., 1. ) ;    
    hDistPhLepMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhLepMea", "measured - true lepton phi", 200, -0.01, 0.01 ) ;    
    hDistPhNeuMea = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhNeuMea", "measured - true neutrino phi", 200, -1., 1. ) ;    

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

    hDistEJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistEJetTrue", "true jet energy", 250, 0., 500. ) ;    
    hDistELepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistELepTrue", "true lepton 1/pt", 200, 0.001, 1. ) ;    
    hDistENeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistENeuTrue", "true neutrino energy", 250, 0., 500. ) ;    
    hDistThJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThJetTrue", "true jet theta", 100, 0., 3.2 ) ;    
    hDistThLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThLepTrue", "true lepton theta", 100, 0., 3.2 ) ;    
    hDistThNeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThNeuTrue", "true neutrino theta", 100, 0., 3.2 ) ;    
    hDistPhJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhJetTrue", "true jet phi", 100, -3.2, 3.2 ) ;    
    hDistPhLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhLepTrue", "true lepton phi", 100, -3.2, 3.2 ) ;    
    hDistPhNeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhNeuTrue", "true neutrino phi", 100, -3.2, 3.2 ) ;    

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
           
     // fill start values of constraints
#ifdef MARLIN_USE_AIDA
     hPxConstStart->fill (topevent->getPxConstraint().getValue()); 
     hPyConstStart->fill (topevent->getPyConstraint().getValue()); 
     hPzConstStart->fill (topevent->getPzConstraint().getValue()); 
     hEConstStart->fill (topevent->getEConstraint().getValue()); 
     hMW1ConstStart->fill (topevent->getW1Constraint().getValue()); 
     hMW2ConstStart->fill (topevent->getW2Constraint().getValue()); 
     hMConstStart->fill (topevent->getTopConstraint().getValue()); 
#endif

     int ierr = topevent->fitEvent(fitter);
  
     double prob = fitter.getProbability();
     if (prob < 0.002) message<WARNING>( log() << "fit probability = " << prob << " in TOY EVENT " << ievt ) ;  
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
       
       // fill final values of constraints
       hPxConstStop->fill (topevent->getPxConstraint().getValue()); 
       hPyConstStop->fill (topevent->getPyConstraint().getValue()); 
       hPzConstStop->fill (topevent->getPzConstraint().getValue()); 
       hEConstStop->fill (topevent->getEConstraint().getValue()); 
       hMW1ConstStop->fill (topevent->getW1Constraint().getValue()); 
       hMW2ConstStop->fill (topevent->getW2Constraint().getValue()); 
       hMConstStop->fill (topevent->getTopConstraint().getValue()); 
       
       message<DEBUG>( log() << "looping over FOs " ) ;
       for (int ifo = 0; ifo < 6; ifo++){
         double errfit, errmea, start, sigma; 
         double pull[3], pullmea[3], pulltrue[3];
         bool usesigma[3];
         for (int ipar = 0; ipar < 3; ipar++) {
           errfit = topevent->getFittedFitObject(ifo)->getError(ipar);  // should make difference for NewtonFitter
           errmea = topevent->getStartFitObject(ifo)->getError(ipar);  // SmearedFO are not fitted: original errors
           start = topevent->getStartFitObject(ifo)->getParam(ipar);  // SmearedFO are not fitted: original values
           sigma = errmea*errmea-errfit*errfit;
           message<DEBUG>( log() << " sigma =  " << sigma << " for ifo " << ifo 
                                 << " in evt " << ievt << ", errmea =  " << errmea << ", errfit = " << errfit ) ;
           if (sigma > 0) {
             sigma = sqrt(sigma);
             usesigma[ipar] = true;
           }
           else {
             message<WARNING>( log() << " SIGMA <= 0, taking only fitted errors for pull for ifo = " << ifo 
                                     << " in evt " << ievt << ", errmea =  " << errmea << ", errfit = " << errfit ) ;
             sigma = errfit*errfit;
             usesigma[ipar] = false;
             usesigma_evt = false;
           }  
           pull[ipar] = (start - topevent->getFittedFitObject(ifo)->getParam(ipar))/sigma;
           pullmea[ipar] = (start - topevent->getFittedFitObject(ifo)->getParam(ipar))/errmea;
           pulltrue[ipar] = (start - topevent->getTrueFitObject(ifo)->getParam(ipar))/errmea;
           message<DEBUG>( log() << " pull =  " << pull[ipar] << " for ifo " << ifo 
                                 << " in evt " << ievt << ", pullmea =  " << pullmea[ipar] 
                                 << ", delta = " << start - topevent->getFittedFitObject(ifo)->getParam(ipar) ) ;
         }  
         if ( !_semileptonic || ifo < 4 ) {
           hDistEJetOK->fill  (topevent->getFittedFitObject(ifo)->getParam(0) - topevent->getTrueFitObject(ifo)->getParam(0)); hPullEJetOK->fill (pull[0]); 
           hDistThJetOK->fill (topevent->getFittedFitObject(ifo)->getParam(1) - topevent->getTrueFitObject(ifo)->getParam(1)); hPullThJetOK->fill(pull[1]); 
           hDistPhJetOK->fill (topevent->getFittedFitObject(ifo)->getParam(2) - topevent->getTrueFitObject(ifo)->getParam(2)); hPullPhJetOK->fill(pull[2]); 
           hDistEJetMea->fill  (topevent->getStartFitObject(ifo)->getParam(0) - topevent->getTrueFitObject(ifo)->getParam(0));  hPullEJetMea->fill (pullmea[0]); 
           hDistThJetMea->fill (topevent->getStartFitObject(ifo)->getParam(1) - topevent->getTrueFitObject(ifo)->getParam(1)); hPullThJetMea->fill (pullmea[1]); 
           hDistPhJetMea->fill (topevent->getStartFitObject(ifo)->getParam(2) - topevent->getTrueFitObject(ifo)->getParam(2)); hPullPhJetMea->fill (pullmea[2]); 
           hDistEJetTrue->fill (topevent->getTrueFitObject(ifo)->getParam(0));  hPullEJetTrue->fill (pulltrue[0]);
           hDistThJetTrue->fill(topevent->getTrueFitObject(ifo)->getParam(1)); hPullThJetTrue->fill (pulltrue[1]);
           hDistPhJetTrue->fill(topevent->getTrueFitObject(ifo)->getParam(2)); hPullPhJetTrue->fill (pulltrue[2]);
         }
         else if ( _semileptonic && ifo ==  4 ) {
           hDistELepOK->fill  (topevent->getFittedFitObject(ifo)->getParam(0) - topevent->getTrueFitObject(ifo)->getParam(0)); hPullELepOK->fill (pull[0]);  
           hDistThLepOK->fill (topevent->getFittedFitObject(ifo)->getParam(1) - topevent->getTrueFitObject(ifo)->getParam(1)); hPullThLepOK->fill(pull[1]);  
           hDistPhLepOK->fill (topevent->getFittedFitObject(ifo)->getParam(2) - topevent->getTrueFitObject(ifo)->getParam(2)); hPullPhLepOK->fill(pull[2]);  
           hDistELepMea->fill  (topevent->getStartFitObject(ifo)->getParam(0) - topevent->getTrueFitObject(ifo)->getParam(0)); hPullELepMea->fill (pullmea[0]); 
           hDistThLepMea->fill (topevent->getStartFitObject(ifo)->getParam(1) - topevent->getTrueFitObject(ifo)->getParam(1)); hPullThLepMea->fill(pullmea[1]); 
           hDistPhLepMea->fill (topevent->getStartFitObject(ifo)->getParam(2) - topevent->getTrueFitObject(ifo)->getParam(2)); hPullPhLepMea->fill(pullmea[2]); 
           hDistELepTrue->fill (topevent->getTrueFitObject(ifo)->getParam(0)); hPullELepTrue->fill (pulltrue[0]);
           hDistThLepTrue->fill(topevent->getTrueFitObject(ifo)->getParam(1)); hPullThLepTrue->fill(pulltrue[1]);
           hDistPhLepTrue->fill(topevent->getTrueFitObject(ifo)->getParam(2)); hPullPhLepTrue->fill(pulltrue[2]);
         }
         else if ( _semileptonic && ifo ==  5 ) {
           hDistENeuOK->fill  (topevent->getFittedFitObject(ifo)->getParam(0) - topevent->getTrueFitObject(ifo)->getParam(0)); hPullENeuOK->fill (pull[0]);  
           hDistThNeuOK->fill (topevent->getFittedFitObject(ifo)->getParam(1) - topevent->getTrueFitObject(ifo)->getParam(1)); hPullThNeuOK->fill(pull[1]);  
           hDistPhNeuOK->fill (topevent->getFittedFitObject(ifo)->getParam(2) - topevent->getTrueFitObject(ifo)->getParam(2)); hPullPhNeuOK->fill(pull[2]);  
           hDistENeuMea->fill  (topevent->getStartFitObject(ifo)->getParam(0) - topevent->getTrueFitObject(ifo)->getParam(0)); hPullENeuMea->fill (pullmea[0]); 
           hDistThNeuMea->fill (topevent->getStartFitObject(ifo)->getParam(1) - topevent->getTrueFitObject(ifo)->getParam(1)); hPullThNeuMea->fill(pullmea[1]); 
           hDistPhNeuMea->fill (topevent->getStartFitObject(ifo)->getParam(2) - topevent->getTrueFitObject(ifo)->getParam(2)); hPullPhNeuMea->fill(pullmea[2]); 
           hDistENeuTrue->fill (topevent->getTrueFitObject(ifo)->getParam(0)); hPullENeuTrue->fill (pulltrue[0]);
           hDistThNeuTrue->fill(topevent->getTrueFitObject(ifo)->getParam(1)); hPullThNeuTrue->fill(pulltrue[1]);
           hDistPhNeuTrue->fill(topevent->getTrueFitObject(ifo)->getParam(2)); hPullPhNeuTrue->fill(pulltrue[2]);
         }
       }
     }
#endif
     if (ierr > 0) message<WARNING>( log() << "FIT ERROR = " << ierr << " in toy event " << ievt ) ;
     
     //if (!usesigma_evt) break;

     _nEvt ++ ;
     
   }   
}



void TopTester::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TopTester::end(){ 

 delete topevent;

}

