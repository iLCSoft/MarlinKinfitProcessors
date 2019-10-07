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
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "TextTracer.h"
#include "NewtonFitterGSL.h"
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
                              
  registerProcessorParameter( "softmasses" ,
                              "set true if mass constraints should be soft",
                              _softmasses,
                              (bool)false);
                              
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
  topevent->softmasses = _softmasses;
  topevent->leptonic = _semileptonic;
  topevent->leptonasjet = _leptonasjet;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TopTester::processRunHeader( LCRunHeader* ) { 
  _nRun++ ;
} 

void TopTester::processEvent( LCEvent * evt ) { 

    
    streamlog_out(ERROR) 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      << std::endl ;
                      
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
     
  static AIDA::IHistogram1D* hDistEJetSmear;    
  static AIDA::IHistogram1D* hDistELepSmear;    
  static AIDA::IHistogram1D* hDistENeuSmear;    
  static AIDA::IHistogram1D* hDistThJetSmear;    
  static AIDA::IHistogram1D* hDistThLepSmear;    
  static AIDA::IHistogram1D* hDistThNeuSmear;    
  static AIDA::IHistogram1D* hDistPhJetSmear;    
  static AIDA::IHistogram1D* hDistPhLepSmear;
  static AIDA::IHistogram1D* hDistPhNeuSmear;
      
  static AIDA::IHistogram1D* hPullEJetSmear;    
  static AIDA::IHistogram1D* hPullELepSmear;    
  static AIDA::IHistogram1D* hPullENeuSmear;    
  static AIDA::IHistogram1D* hPullThJetSmear;    
  static AIDA::IHistogram1D* hPullThLepSmear;    
  static AIDA::IHistogram1D* hPullThNeuSmear;    
  static AIDA::IHistogram1D* hPullPhJetSmear;    
  static AIDA::IHistogram1D* hPullPhLepSmear;    
  static AIDA::IHistogram1D* hPullPhNeuSmear;    
               
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
      createHistogram1D( "hDistEJetOK", "fitted - start jet energy", 200, -50., 50. ) ;    
    hDistELepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistELepOK", "fitted - start lepton 1/pt", 200,  -0.001, 0.001 ) ;    
    hDistENeuOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistENeuOK", "fitted - start neutrino energy", 200,  -50., 50. ) ;    
    hDistThJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThJetOK", "fitted - start jet theta", 200, -1., 1.) ;    
    hDistThLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThLepOK", "fitted - start lepton theta", 200,  -0.01, 0.01) ;    
    hDistThNeuOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThNeuOK", "fitted - start neutrino theta", 200,  -1., 1.) ;    
    hDistPhJetOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhJetOK", "fitted - start jet phi", 200, -1., 1. ) ;    
    hDistPhLepOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhLepOK", "fitted - start lepton phi", 200, -0.01, 0.01 ) ;    
    hDistPhNeuOK = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhNeuOK", "fitted - start neutrino phi", 200, -1., 1. ) ;    

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

    hDistEJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistEJetTrue", "fitted - true jet energy", 200, -50., 50. ) ;    
    hDistELepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistELepTrue", "fitted - true lepton 1/pt", 200, -0.001, 0.001 ) ;    
    hDistENeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistENeuTrue", "fitted - true neutrino energy", 200, -50., 50. ) ;    
    hDistThJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThJetTrue", "fitted - true jet theta", 200, -1., 1. ) ;    
    hDistThLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThLepTrue", "fitted - true lepton theta", 200, -0.01, 0.01 ) ;    
    hDistThNeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThNeuTrue", "fitted - true neutrino theta", 200, -1., 1. ) ;    
    hDistPhJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhJetTrue", "fitted - true jet phi", 200, -1., 1. ) ;    
    hDistPhLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhLepTrue", "fitted - true lepton phi", 200, -0.01, 0.01 ) ;    
    hDistPhNeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhNeuTrue", "fitted - true neutrino phi", 200, -1., 1. ) ;    

    hPullEJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetTrue",  "(fit - truth)/#sigma_fit of jet energy", 100, -5., 5. ) ;    
    hPullELepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepTrue",  "(fit - truth)/#sigma_fit of lepton 1/pt", 100, -5., 5. ) ;    
    hPullENeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullENeuTrue",  "(fit - truth)/#sigma_fit of neutrino energy", 100, -5., 5. ) ;    
    hPullThJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetTrue", "(fit - truth)/#sigma_fit of jet theta", 100, -5., 5. ) ;    
    hPullThLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThLepTrue", "(fit - truth)/#sigma_fit of lepton theta", 100, -5., 5. ) ;    
    hPullThNeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThNeuTrue", "(fit - truth)/#sigma_fit of neutrino theta", 100, -5., 5. ) ;    
    hPullPhJetTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetTrue", "(fit - truth)/#sigma_fit of jet phi", 100, -5., 5. ) ;    
    hPullPhLepTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhLepTrue", "(fit - truth)/#sigma_fit of lepton phi", 100, -5., 5. ) ;    
    hPullPhNeuTrue = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhNeuTrue", "(fit - truth)/#sigma_fit of neutrino phi", 100, -5., 5. ) ;    

    hDistEJetSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistEJetSmear", "start - true jet energy", 200, -50., 50. ) ;    
    hDistELepSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistELepSmear", "start - true lepton 1/pt", 200, -0.001, 0.001 ) ;    
    hDistENeuSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistENeuSmear", "start - true neutrino energy", 200, -50., 50. ) ;    
    hDistThJetSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThJetSmear", "start - true jet theta", 200, -1., 1. ) ;    
    hDistThLepSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThLepSmear", "start - true lepton theta", 200, -0.01, 0.01 ) ;    
    hDistThNeuSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistThNeuSmear", "start - true neutrino theta", 200, -1., 1. ) ;    
    hDistPhJetSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhJetSmear", "start - true jet phi", 200, -1., 1. ) ;    
    hDistPhLepSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhLepSmear", "start - true lepton phi", 200, -0.01, 0.01 ) ;    
    hDistPhNeuSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hDistPhNeuSmear", "start - true neutrino phi", 100, -3.2, 3.2 ) ;    

    hPullEJetSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetSmear", "(start - true)/#sigma_mea of jet energy", 100, -5., 5. ) ;    
    hPullELepSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepSmear", "(start - true)/#sigma_mea of lepton energy", 100, -5., 5. ) ;    
    hPullENeuSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullENeuSmear", "(start - true)/#sigma_mea of neutrino energy", 100, -5., 5. ) ;    
    hPullThJetSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetSmear", "(start - true)/#sigma_mea of jet theta", 100, -5., 5. ) ;    
    hPullThLepSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThLepSmear", "(start - true)/#sigma_mea of lepton theta", 100, -5., 5. ) ;    
    hPullThNeuSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThNeuSmear", "(start - true)/#sigma_mea of neutrino theta", 100, -5., 5. ) ;    
    hPullPhJetSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetSmear", "(start - true)/#sigma_mea of jet phi", 100, -5., 5. ) ;    
    hPullPhLepSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhLepSmear", "(start - true)/#sigma_mea of lepton phi", 100, -5., 5. ) ;    
    hPullPhNeuSmear = 
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhNeuSmear", "(start - true)/#sigma_mea of neutrino phi", 100, -5., 5. ) ;    
  }

#endif   
      
   
  for (int ievt = 0; ievt < _ntoy; ievt++) { 
  
     streamlog_out(MESSAGE)  << "start to process toy event number " << ievt << std::endl ;
     
     // double startmassW1 = 0., startmassW2 = 0.;
     double startmasstop1 = 0., startmasstop2 = 0.;
     
     int debug = 0;
     if (ievt == _ievttrace || _traceall) debug = 10;
     topevent->setDebug (debug);
          
     topevent->genEvent();
       
     startmassW1 = topevent->getW1Mass();
     startmassW2 = topevent->getW2Mass();
     streamlog_out(MESSAGE)  << "start mass of W 1: " << startmassW1 << std::endl ;
     streamlog_out(MESSAGE)  << "start mass of W 2: " << startmassW2 << std::endl ;
     startmasstop1 = topevent->getTop1Mass();
     startmasstop2 = topevent->getTop2Mass();
     streamlog_out(MESSAGE)  << "start mass of top 1: " << startmasstop1 << std::endl ;
     streamlog_out(MESSAGE)  << "start mass of top 2: " << startmasstop2 << std::endl ;
                     
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
       // OPALFitter now also has method setDebug !
       pfitter = new OPALFitterGSL();
       (dynamic_cast<OPALFitterGSL*>(pfitter))->setDebug (debug);
     }
     BaseFitter &fitter = *pfitter;
  
     TextTracer tracer (std::cout);
     if (ievt == _ievttrace || _traceall) fitter.setTracer (tracer);
           
     // fill start values of constraints
#ifdef MARLIN_USE_AIDA
     hPxConstStart->fill (topevent->getPxConstraint()->getValue()); 
     hPyConstStart->fill (topevent->getPyConstraint()->getValue()); 
     hPzConstStart->fill (topevent->getPzConstraint()->getValue()); 
     hEConstStart->fill (topevent->getEConstraint()->getValue()); 
     hMW1ConstStart->fill (topevent->getW1Constraint()->getValue()); 
     hMW2ConstStart->fill (topevent->getW2Constraint()->getValue()); 
     hMConstStart->fill (topevent->getTopConstraint()->getValue()); 
#endif

     int ierr = topevent->fitEvent(fitter);
  
     prob = fitter.getProbability();
     if (prob < 0.002) streamlog_out(WARNING) << "fit probability = " << prob << " in TOY EVENT " << ievt << std::endl ;  
     double chi2 = fitter.getChi2();
     nit = fitter.getIterations();

     streamlog_out(MESSAGE) << "fit probability = " << prob << std::endl ;  
     streamlog_out(MESSAGE) << "fit chi2 = " << chi2  << std::endl ; 
     streamlog_out(MESSAGE) << "error code: " << ierr << std::endl ;
                                  
     streamlog_out(MESSAGE)  << "final mass of W 1: " << topevent->getW1Mass() << std::endl ;
     streamlog_out(MESSAGE)  << "final mass of W 2: " << topevent->getW2Mass() << std::endl ;
     streamlog_out(MESSAGE)  << "final mass of top 1: " << topevent->getTop1Mass() << std::endl ;
     streamlog_out(MESSAGE)  << "final mass of top 2: " << topevent->getTop2Mass() << std::endl ;
       
     //bool usesigma_evt = true;
                  
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
       hPxConstStop->fill (topevent->getPxConstraint()->getValue()); 
       hPyConstStop->fill (topevent->getPyConstraint()->getValue()); 
       hPzConstStop->fill (topevent->getPzConstraint()->getValue()); 
       hEConstStop->fill (topevent->getEConstraint()->getValue()); 
       hMW1ConstStop->fill (topevent->getW1Constraint()->getValue()); 
       hMW2ConstStop->fill (topevent->getW2Constraint()->getValue()); 
       hMConstStop->fill (topevent->getTopConstraint()->getValue()); 
       
       streamlog_out(MESSAGE) << "looping over FOs " << std::endl ;
       for (int ifo = 0; ifo < 6; ifo++){
         double errfit, errmea, sigma; 
         double truth, start, fitted;
         double pull[3], pulltrue[3], pullsmear[3];
         double dist[3], disttrue[3], distsmear[3];
         //bool usesigma[3];
         for (int ipar = 0; ipar < 3; ipar++) {
           errfit = topevent->getFittedFitObject(ifo)->getError(ipar);  // should make difference for NewtonFitter
           errmea = topevent->getStartFitObject(ifo)->getError(ipar);  // SmearedFO are not fitted: original errors
           truth = topevent->getTrueFitObject(ifo)->getParam(ipar);  // TrueFO are not fitted: truth values           
           start = topevent->getStartFitObject(ifo)->getParam(ipar);  // SmearedFO are not fitted: original values           
           fitted = topevent->getFittedFitObject(ifo)->getParam(ipar);  // fitted values           
           if ( topevent->getFittedFitObject(ifo)->isParamMeasured(ipar) ) {
             sigma = errmea*errmea-errfit*errfit;
           }
           else {
             sigma = errfit*errfit;     
           }  
           
           streamlog_out(MESSAGE) << " sigma =  " << sigma << " for ifo " << ifo 
                                   << " in evt " << ievt << ", errmea =  " << errmea << ", errfit = " << errfit << std::endl ;
           if (sigma > 0) {
             sigma = sqrt(sigma);
             //usesigma[ipar] = true;
           }
           else {
             streamlog_out(WARNING) << " SIGMA <= 0, taking only fitted errors for pull for ifo = " << ifo 
                                     << " in evt " << ievt << ", errmea =  " << errmea << ", errfit = " << errfit
                                     << " for ierr =  " << ierr << std::endl ;
             streamlog_out(WARNING) << "fit probability = " << prob << std::endl ;  
             streamlog_out(WARNING) << "fit chi2 = " << chi2  << std::endl ; 
             streamlog_out(WARNING) << "error code: " << ierr << std::endl ;
                                          
             streamlog_out(WARNING)  << "final mass of W 1: " << topevent->getW1Mass() << std::endl ;
             streamlog_out(WARNING)  << "final mass of W 2: " << topevent->getW2Mass() << std::endl ;
             streamlog_out(WARNING)  << "final mass of top 1: " << topevent->getTop1Mass() << std::endl ;
             streamlog_out(WARNING)  << "final mass of top 2: " << topevent->getTop2Mass() << std::endl ;
             sigma = errfit*errfit;
             //usesigma[ipar] = false;
             //usesigma_evt = false;
           } 
           
           dist[ipar] = (fitted - start);
           disttrue[ipar] = (fitted - truth);
           distsmear[ipar] = (start - truth);
           pull[ipar] = dist[ipar]/sigma;
           pulltrue[ipar] = disttrue[ipar]/errfit;
           pullsmear[ipar] = (start - truth)/errmea;
           streamlog_out(MESSAGE) << " pull =  " << pull[ipar] << " for ifo " << ifo                                
                                   << " in evt " << ievt << ", pulltrue =  " << pulltrue[ipar]                         
                                   << ", delta = " << start - topevent->getFittedFitObject(ifo)->getParam(ipar) << std::endl ;  
                                   
         }  
         if ( !_semileptonic || ifo < 4 ) {
           hDistEJetOK->fill    (dist[0]);      hPullEJetOK->fill (pull[0]); 
           hDistThJetOK->fill   (dist[1]);      hPullThJetOK->fill(pull[1]); 
           hDistPhJetOK->fill   (dist[2]);      hPullPhJetOK->fill(pull[2]); 
           hDistEJetTrue->fill  (disttrue[0]);  hPullEJetTrue->fill (pulltrue[0]); 
           hDistThJetTrue->fill (disttrue[1]);  hPullThJetTrue->fill (pulltrue[1]); 
           hDistPhJetTrue->fill (disttrue[2]);  hPullPhJetTrue->fill (pulltrue[2]); 
           hDistEJetSmear->fill (distsmear[0]); hPullEJetSmear->fill (pullsmear[0]);
           hDistThJetSmear->fill(distsmear[1]); hPullThJetSmear->fill (pullsmear[1]);
           hDistPhJetSmear->fill(distsmear[2]); hPullPhJetSmear->fill (pullsmear[2]);
         }
         else if ( _semileptonic && ifo ==  4 ) {
           hDistELepOK->fill    (dist[0]);      hPullELepOK->fill (pull[0]);  
           hDistThLepOK->fill   (dist[1]);      hPullThLepOK->fill(pull[1]);  
           hDistPhLepOK->fill   (dist[2]);      hPullPhLepOK->fill(pull[2]);  
           hDistELepTrue->fill  (disttrue[0]);  hPullELepTrue->fill (pulltrue[0]); 
           hDistThLepTrue->fill (disttrue[1]);  hPullThLepTrue->fill(pulltrue[1]); 
           hDistPhLepTrue->fill (disttrue[2]);  hPullPhLepTrue->fill(pulltrue[2]); 
           hDistELepSmear->fill (distsmear[0]); hPullELepSmear->fill (pullsmear[0]);
           hDistThLepSmear->fill(distsmear[1]); hPullThLepSmear->fill(pullsmear[1]);
           hDistPhLepSmear->fill(distsmear[2]); hPullPhLepSmear->fill(pullsmear[2]);
         }
         else if ( _semileptonic && ifo ==  5 ) {
           hDistENeuOK->fill    (dist[0]);      hPullENeuOK->fill (pull[0]); 
           hDistThNeuOK->fill   (dist[1]);      hPullThNeuOK->fill(pull[1]); 
           hDistPhNeuOK->fill   (dist[2]);      hPullPhNeuOK->fill(pull[2]); 
           hDistENeuTrue->fill  (disttrue[0]);  hPullENeuTrue->fill (pulltrue[0]); 
           hDistThNeuTrue->fill (disttrue[1]);  hPullThNeuTrue->fill(pulltrue[1]); 
           hDistPhNeuTrue->fill (disttrue[2]);  hPullPhNeuTrue->fill(pulltrue[2]); 
           hDistENeuSmear->fill (distsmear[0]); hPullENeuSmear->fill (pullsmear[0]);
           hDistThNeuSmear->fill(distsmear[1]); hPullThNeuSmear->fill(pullsmear[1]);
           hDistPhNeuSmear->fill(distsmear[2]); hPullPhNeuSmear->fill(pullsmear[2]);
         }
       }
     }
#endif
     if (ierr > 0) streamlog_out(WARNING) << "FIT ERROR = " << ierr << " in toy event " << ievt << std::endl ;
     
     //if (!usesigma_evt) break;

     _nEvt ++ ;
     
   }   
}



void TopTester::check( LCEvent * ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TopTester::end(){ 

 delete topevent;

}

