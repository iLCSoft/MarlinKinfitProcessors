#include "ZHllqq5CFit.h"
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
#include <EVENT/Track.h>
#include <EVENT/LCGenericObject.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <EVENT/LCFloatVec.h>


#include <GeometryUtil.h>
#include "TLorentzVector.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
#include "LeptonFitObject.h"
#include "ISRPhotonFitObject.h"
#include "MomentumConstraint.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "TextTracer.h"
#include "NewtonFitterGSL.h"
#include "FourJetZHPairing.h"
#include "MassConstraint.h"
#include "SoftGaussParticleConstraint.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCCollectionVec.h"
#include <EVENT/LCCollection.h>
using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


ZHllqq5CFit aZHllqq5CFit ;

// function to define the jet energy resolution (in GeV)
double ZHllqq5CFit::JetEnergyResolution(double E)
{

  // examples here derived by Benjamin Hermberg from e+e- -> udsc:
  // 1) default is 120%/sqrt(E), gives best convergence of 5C fit on e+e- -> udsc
  double result = _errene*std::sqrt(E);

  // 2) comparing jet-level to quark-level energies
  //    (using MarlinReco/Analysis/RecoMCTruthLink/QuarkJetPairing.cc)
  if (_errene == 0 ) result = std::sqrt(pow(0.6908,2)*(E)+(pow(0.02596,2)*pow(E,2)));

  return result;
}

ZHllqq5CFit::ZHllqq5CFit() : Processor("ZHllqq5CFit"),
			     m_Bfield(0.f),
			     c(0.),
			     mm2m(0.),
			     eV2GeV(0.),
			     eB(0.)
{

  // modify processor description
  _description = "ZHllqq5CFit does a 5C fit on 4 jet events (Px, Py, Pz, E, M12 = MZ (for all six permutations))" ;


  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "JetCollectionName" ,
			   "Name of the Jet collection"  ,
			   _jetcolName ,
			   std::string("Durham2Jets") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "InputIsoLeptonsCollection",
                           "Name of collection with the selected isolated lepton",
                           _colLeptons,
                           std::string("ISOLeptons") );

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

 registerProcessorParameter( "SigmaEnergyScaleFactor",
                        		 "Scale Factor t be applied to jet energy uncertainty",
                        		 sigmaScaleFactor,
                        		 (float) 1.0 );

  registerProcessorParameter( "errene" ,
                              "assumed energy resolution for jets as x/sqrt(E) - if 0, then parametrisation is used",
                              _errene,
                              (double)1.2);

  registerProcessorParameter( "errtheta" ,
                              "assumed theta resolution for jet axis",
                              _errtheta,
                              (double)0.1);

  registerProcessorParameter( "errphi" ,
                              "assumed phi resolution for jet axis",
                              _errphi,
                              (double)0.1);

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
  registerInputCollection( LCIO::MCPARTICLE,
                                "MCParticleCollection" , //name
                                "Name of the MCParticle collection"  , //description
                                _colMCP , //my collection name
                                std::string("MCParticle")   );//take this string from input slcio

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                                "FitOutputColection",
                      	        " Output Fit Colection" ,
                      	         _OutCol,
                                std::string("FitReco")) ;

 registerOutputCollection( LCIO::LCFLOATVEC,
                                "PullsEOutputColection",
                                " Output Pull E  Colection" ,
                                 _OutPullECol,
                                std::string("PullEnergy")) ;

  registerOutputCollection( LCIO::LCFLOATVEC,
                                "PullsThetaOutputColection",
                                " Output Pull Theta  Colection" ,
                                 _OutPullThetaCol,
                                std::string("PullTheta")) ;

 registerOutputCollection( LCIO::LCFLOATVEC,
                                "PullsPhiOutputColection",
                                " Output Pull Phi  Colection" ,
                                 _OutPullPhiCol,
                                std::string("PullPhi")) ;



 registerProcessorParameter("outputFilename",
                               "name of output file",
                               _outfile,
                                std::string("")
                                );
}


void ZHllqq5CFit::init() {

  // usually a good idea to
  printParameters() ;
  trackcounter = 0;
  _nRun = 0 ;
  _nEvt = 0 ;

  b = (double) 0.00464564*( std::log(_ecm*_ecm*3814714.)-1. );
  //= 2*alpha/pi*( ln(s/m_e^2)-1 )
  ISRPzMaxB = std::pow((double)_isrpzmax,b);

  m_Bfield = MarlinUtil::getBzAtOrigin();
  streamlog_out(DEBUG) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
  c = 2.99792458e8;
  mm2m = 1e-3;
  eV2GeV = 1e-9;
  eB = m_Bfield * c * mm2m * eV2GeV;

  _fout = new TFile(_outfile.c_str(),"recreate");

   ZHTree = new TTree("ZHTree","ZHTree");
   ZHTree->Branch("Hmass_before_fit_best",&Hmass_before_fit,"Hmass_before_fit/F") ;
   ZHTree->Branch("Hmass_No_fit_best",&Hmass_NoFit,"Hmass_NoFit/F") ;
   ZHTree->Branch("Hmass_after_fit_best",&Hmass_after_fit,"Hmass_after_fit/F") ;
   ZHTree->Branch("Error_Code",&Error_code,"Error_code/F") ;
   ZHTree->Branch("hpull_jet1_E_best",&hpull_jet1_E,"hpull_jet1_E/F") ;
   ZHTree->Branch("hpull_jet2_E_best",&hpull_jet2_E,"hpull_jet2_E/F") ;
   ZHTree->Branch("hpull_jet1_th_best",&hpull_jet1_th,"hpull_jet1_th/F") ;
   ZHTree->Branch("hpull_jet2_th_best",&hpull_jet2_th,"hpull_jet2_th/F") ;
   ZHTree->Branch("hpull_jet1_phi_best",&hpull_jet1_phi,"hpull_jet1_phi/F") ;
   ZHTree->Branch("hpull_jet2_phi_best",&hpull_jet2_phi,"hpull_jet2_phi/F") ;
   ZHTree->Branch("hpull_lepton1_InvpT_best",&hpull_lepton1_InvpT,"hpull_lepton1_InvpT/F") ;
   ZHTree->Branch("hpull_lepton2_InvpT_best",&hpull_lepton2_InvpT,"hpull_lepton2_InvpT/F") ;
   ZHTree->Branch("hpull_lepton1_th_best",&hpull_lepton1_th,"hpull_lepton1_th/F") ;
   ZHTree->Branch("hpull_lepton2_th_best",&hpull_lepton2_th,"hpull_lepton2_th/F") ;
   ZHTree->Branch("hpull_lepton1_phi_best",&hpull_lepton1_phi,"hpull_lepton1_phi/F") ;
   ZHTree->Branch("hpull_lepton2_phi_best",&hpull_lepton2_phi,"hpull_lepton2_phi/F") ;
   ZHTree->Branch("jetmatch",&jetmatch,"jetmatch/I") ;
   ZHTree->Branch("jetmatchth",&jetmatchth,"jetmatchth/I") ;
   ZHTree->Branch("jetmatchphi",&jetmatchphi,"jetmatchphi/I") ;


}

void ZHllqq5CFit::processRunHeader( LCRunHeader* ) {
  _nRun++ ;
}
void ZHllqq5CFit::compcorrect() //finds the jet with cross checking
{
      if(std::abs(delta_theta[bestjet_phi])<std::abs(delta_phi[bestjet_th])){
      // if(std::abs(diff_besttheta)<std::abs(diff_bestphi)){
      bestjet=bestjet_phi;
      // }
      // else{
        //    bestjet=bestjet_th;
      // }
      }
      if(std::abs(delta_theta[bestjet_phi])>std::abs(delta_phi[bestjet_th])){
      bestjet=bestjet_th;
      }
}

void ZHllqq5CFit::SetZero()
{

  Px=0.;
  Px2=0.;
  Py=0.;
  Py2=0.;
  Pz=0.;
  Pz2=0.;
  pT2=0.;
  P=0.;
  P2=0.;
  SigPx2=0.;
  SigPxSigPy=0.;
  SigPxSigPz=0.;
  SigPy2=0.;
  SigPySigPz=0.;
  SigPz2=0.;
  SigE2=0.;
  dth_dpx=0.;
  dth_dpy=0.;
  dth_dpz=0.;
  dphi_dpx=0.;
  dphi_dpy=0.;
  JetResE=0.;
  JetResTheta=0.;
  JetResPhi=0.;
  LeptonResE=0.;
  LeptonResTheta=0.;
  LeptonResPhi=0.;
  invmass = 0;  

  bestprob = 0.;
  bestnit = 0;
  //double bestmass1 = 0., bestmass2 = 0.;
  beststartmassZ = 0., beststartmassH = 0.;
  startmassZ = 0., startmassH = 0.;
  bestphotonenergy = 0.;
   besterr = 999;
   bestzvalue = 10000.;
   chi2startmassZ = 0.;
   chi2startmassH = 0.;
  // for (int i=0 ; i<3; i++){
  // Zmomentum [i]=0.;
  // Hmomentum [i]=0.;
  // ISRmomentum [i]=0.;
  // }setZero
  


  memset(Zmomentum, 0, sizeof(Zmomentum));
  memset(Hmomentum, 0, sizeof(Hmomentum));
  memset(ISRmomentum, 0, sizeof(ISRmomentum));
  Z_Energy=0.;
  H_Energy=0.;
  chi2best=0.;
  errorcode=0.;
     streamlog_out(DEBUG)  << "Values set to zero" <<std::endl;

}


void ZHllqq5CFit::processEvent( LCEvent * evt ) { //event start
  
  
  streamlog_out(MESSAGE)
    << " processing event " << evt->getEventNumber()
    << "  in run "          << evt->getRunNumber()
    << std::endl ;
  // this gets called for every event
  // usually the working horse ...
  
  int debug = 0;
  if ( evt->getEventNumber() == _ievttrace || _traceall) debug = 10;
  
#ifdef MARLIN_USE_AIDA
  
  // define a histogram pointer
  static AIDA::IHistogram1D* hRecHMassBest ;
  static AIDA::IHistogram1D* hRecHMassAll ;
  static AIDA::IHistogram1D* hRecHMassNoFitBest ;
  static AIDA::IHistogram1D* hRecZMassNoFitBest ;
  static AIDA::IHistogram1D* hRecHMassNoFitAll ;
  static AIDA::IHistogram1D* hTestHMassNoFitAll ;
  static AIDA::IHistogram1D* hRecHMassNoFitFail ;
  static AIDA::IHistogram1D* hRecZMassNoFitFail ;
  static AIDA::IHistogram1D* hFitProbBest ;
  static AIDA::IHistogram1D* hFitProbAll ;
  static AIDA::IHistogram1D* hNItBest ;
  static AIDA::IHistogram1D* hNItAll ;
  static AIDA::IHistogram1D* hPhotonEnergy ;
  static AIDA::IHistogram1D* hJetMass ;
  static AIDA::IHistogram1D* hFitError;
  static AIDA::IHistogram1D* hFitErrorBest;
  static AIDA::IHistogram1D* hPullEJetOK;
  static AIDA::IHistogram1D* hPullThJetOK;
  static AIDA::IHistogram1D* hPullPhJetOK;
  static AIDA::IHistogram1D* hPullELeptonOK;
  static AIDA::IHistogram1D* hPullThLeptonOK;
  static AIDA::IHistogram1D* hPullPhLeptonOK;
  static AIDA::IHistogram1D* hPullEJetBest;
  static AIDA::IHistogram1D* hPullThJetBest;
  static AIDA::IHistogram1D* hPullPhJetBest;
  static AIDA::IHistogram1D* hPullELeptonBest;
  static AIDA::IHistogram1D* hPullThLeptonBest;
  static AIDA::IHistogram1D* hPullPhLeptonBest;


  if( isFirstEvent() ) {
    
    hRecHMassBest =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassBest", "M_H", 200, 0., 250. ) ;
    hRecHMassAll =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassAll", "M_H", 200, 0., 250. ) ;
    hRecHMassNoFitBest =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassNoFitBest", "M_H", 500, 0., 500. ) ;
    hRecZMassNoFitBest =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecZMassNoFitBest", "M_Z", 200, 0., 250. ) ;
    hRecHMassNoFitAll =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassNoFitAll", "M_H", 200, 0., 250. ) ;
    hTestHMassNoFitAll =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hTestHMassNoFitAll", "M_H", 200, 0., 250. ) ;
    hRecHMassNoFitFail =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassNoFitFail", "M_H", 200, 0., 250. ) ;
    hRecZMassNoFitFail =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecZMassNoFitFail", "M_Z", 200, 0., 250. ) ;
    hFitProbBest =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ;
    hFitProbAll =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProbAll", "fit probability", 100, 0., 1. ) ;
    hNItBest =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hNItBest", "number of iterations", 200, 0, 200 ) ;
    hNItAll =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hNItAll", "number of iterations", 200, 0, 200 ) ;
    hPhotonEnergy =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPhotonEnergy", "ISR photon energy", 250, 0., 250. ) ;
    hJetMass =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hJetMass", "Jet Mass", 200, 0., 100. ) ;
    if (_ifitter == 1) {
      hFitError =
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 100, -1.5, 99.5 ) ;
      hFitErrorBest =
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitErrorBest", "Error flag", 100, -1.5, 98.5 ) ;
    }
    else {
      hFitError =
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 11, -1.5, 9.5 ) ;
      hFitErrorBest =
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitErrorBest", "Error flag", 11, -1.5, 9.5 ) ;
    }
    hPullEJetOK =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetOK", "pull of jet energy", 100, -5., 5. ) ;
    hPullThJetOK =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetOK", "pull of jet theta", 100, -5., 5. ) ;
    hPullPhJetOK =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetOK", "pull of jet phi", 100, -5., 5. ) ;
    hPullEJetBest =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJetBest", "pull of jet energy", 100, -5., 5. ) ;
    hPullThJetBest =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJetBest", "pull of jet theta", 100, -5., 5. ) ;
    hPullPhJetBest =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJetBest", "pull of jet phi", 100, -5., 5. ) ;
   
  }

#endif





  HepLorentzVector leptonvec;
  HepLorentzVector jetvec;


  // fill histogram from LCIO data :

  //////////////////   JETS ///////////////////////////

     LCCollection* jetcol = evt->getCollection( _jetcolName ) ;
     LCCollection* lepcol = evt->getCollection( _colLeptons ) ;
     if (jetcol != 0 && lepcol != 0) {//jetcol is not null

       int nJETS = jetcol->getNumberOfElements()  ;
       int nLEPS = lepcol->getNumberOfElements()  ;
       streamlog_out(MESSAGE)
                      << " found " << nJETS

                      << "  in run "          << evt->getRunNumber()
                      << std::endl ;

       if (nJETS != 2) return;
       if (nLEPS != 2) return;

       float yminus = jetcol ->parameters().getFloatVal( "YMinus");
       streamlog_out(MESSAGE)  << " yminus = " << yminus << std::endl ;

       float yplus = jetcol ->parameters().getFloatVal( "YPlus");
       streamlog_out(MESSAGE)  << " yplus = " << yplus << std::endl ;

       //float besttheta =100.;
       //float bestphi =100.;
       //bestjet_th=0;
       //bestjet_phi=0;
       //bestjet=1000000;
       //float diff_besttheta;
       //float diff_bestphi;


       double Omega = 0.;
       double Omega_uncert = 0.;
       double TanLambda = 0.;
       double TanLambda_err = 0.;
       double theta = 0.;
       double theta_err = 0.;
       double phi = 0.;
       double phi_err = 0.;
       double invers_pT = 0.;
       double invers_pT_err = 0.;

       float pxc_before_ISR = 0.;
       float pyc_before_ISR = 0.;
       float pzc_before_ISR = 0.;
       float ec_before_ISR = 0.;

       float pxc_before_fit = 0.;
       float pyc_before_fit = 0.;
       float pzc_before_fit = 0.;
       float ec_before_fit = 0.;

       float pxc_after_fit = 0.;
       float pyc_after_fit = 0.;
       float pzc_after_fit = 0.;
       float ec_after_fit = 0.;


// original fit objects - save for next permutation
     JetFitObject* j1 = 0;
     JetFitObject* j2 = 0;
     //orginal lepton fit objects
     LeptonFitObject* l1 = 0;
     LeptonFitObject* l2 = 0;

      // if (nSLDB == 1 && nSLDC == 0){
      //    NuEne[0]=ENuplus;
      //    NuEne[1]=ENuminus;
      // }
     SetZero();

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////
 ////Set JetFitObjects
 ////
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////


 ReconstructedParticle* jrps[2];
 for (int i_jet = 0; i_jet < nJETS; i_jet++)
   {
     ReconstructedParticle *j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i_jet ) );
     ReconstructedParticle *jet2= dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( 1 ) );
     if (j) {
       jrps[i_jet] = j;
       streamlog_out(MESSAGE)
	 << " found jet in event " << evt->getEventNumber()
	 << "  in run "          << evt->getRunNumber()
	 << std::endl ;

       Px=j->getMomentum()[0];
       Px2=std::pow(Px,2);
       Py=j->getMomentum()[1];
       Py2=std::pow(Py,2);
       Pz=j->getMomentum()[2];
       Pz2=std::pow(Pz,2);
       pT2=Px2+Py2;
       P=std::sqrt(std::pow(Px,2)+std::pow(Py,2)+std::pow(Pz,2));
       
       SigPx2     = j->getCovMatrix()[ 0 ];
       SigPxSigPy = j->getCovMatrix()[ 1 ];
       SigPxSigPz = j->getCovMatrix()[ 3 ];
       SigPy2     = j->getCovMatrix()[ 2 ];
       SigPySigPz = j->getCovMatrix()[ 4 ];
       SigPz2     = j->getCovMatrix()[ 5 ];
       SigE2      = j->getCovMatrix()[ 9 ];
       
       dth_dpx    = ( Px * Pz ) / ( std::pow( P , 2 ) * std::sqrt( pT2 ) );
       dth_dpy    = ( Py * Pz ) / ( std::pow( P , 2 ) * std::sqrt( pT2 ) );
       dth_dpz    = -std::sqrt( pT2 ) / std::pow( P , 2 );
       dphi_dpx   = -Py / pT2;
       dphi_dpy   = Px / pT2 ;
       
       JetResE    = std::sqrt( SigE2 ) * sigmaScaleFactor;
       JetResTheta= std::sqrt( std::fabs( SigPx2 * std::pow( dth_dpx , 2 ) + SigPy2 * std::pow( dth_dpy , 2 ) + SigPz2 * std::pow( dth_dpz , 2 ) + 2 * ( SigPxSigPy * dth_dpx * dth_dpy ) + 2 * ( SigPySigPz * dth_dpy * dth_dpz ) + 2 * ( SigPxSigPz * dth_dpx * dth_dpz ) ) );
       JetResPhi  = std::sqrt( std::fabs( SigPx2 * std::pow( dphi_dpx , 2 ) + SigPy2 * std::pow( dphi_dpy , 2 ) + 2 * ( SigPxSigPy * dphi_dpx * dphi_dpy ) ) );
       
       streamlog_out(DEBUG4)  << "SigPx2= " << SigPx2 <<std::endl;
       streamlog_out(DEBUG4)  << "SigPxSigPy= " << SigPxSigPy <<std::endl;
       streamlog_out(DEBUG4)  << "SigPxSigPz= " << SigPxSigPz <<std::endl;
       streamlog_out(DEBUG4)  << "SigPy2= " <<  SigPy2<<std::endl;
       streamlog_out(DEBUG4)  << "SigPySigPz= " << SigPySigPz <<std::endl;
       streamlog_out(DEBUG4)  << "SigPz2= " << SigPz2 <<std::endl;
       streamlog_out(DEBUG4)  << "SigE2= " << SigE2 <<std::endl;
       
       streamlog_out(DEBUG)  << "Px= " << Px <<std::endl;
       streamlog_out(DEBUG)  << "Py= " << Py <<std::endl;
       streamlog_out(DEBUG)  << "Pz= " << Pz <<std::endl;
       
       streamlog_out(DEBUG)  << "JetResE= " << JetResE <<std::endl;
       streamlog_out(DEBUG)  << "JetResTheta= " << JetResTheta <<std::endl;
       streamlog_out(DEBUG)  << "JetResPhi= " << JetResPhi <<std::endl;
       
       if(SigPx2==0. || SigPxSigPy==0. || SigPxSigPz==0. || SigPy2==0. || SigPySigPz==0. || SigPz2==0. || SigE2==0.){
	 streamlog_out(WARNING) << "Covariance Matrix is singular"<<std::endl;
	 streamlog_out(WARNING) << "Setting theta and phi Resolution back to previous values "<<std::endl;
	 JetResTheta=0.1;
	 JetResPhi=0.1;
	 nCo++;
       }

       jetvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy());
       if (i_jet == 0 ) {
	 j1 = new JetFitObject (jetvec.e(), jetvec.theta(), jetvec.phi(),
				JetResE*sigmaScaleFactor, JetResTheta, JetResPhi, jetvec.m());
	 j1->setName("Jet1");
	 streamlog_out(DEBUG)  << " start four-vector of first  jet: " << *j1  << std::endl ;
       }
       else if (i_jet == 1) {
	 j2 = new JetFitObject (jetvec.e(), jetvec.theta(), jetvec.phi(),
				JetResE*sigmaScaleFactor, JetResTheta, JetResPhi, jetvec.m());
	 j2->setName("Jet2");
	 streamlog_out(DEBUG) << " start four-vector of second  jet: " << *j2  << std::endl ;
      }
#ifdef MARLIN_USE_AIDA
       hJetMass->fill(j->getMass());
#endif
       
  
   // invariante Masse der b quarks.

       double px = j->getMomentum()[0]+jet2->getMomentum()[0];
       double py = j->getMomentum()[1]+jet2->getMomentum()[1];
       double pz = j->getMomentum()[2]+jet2->getMomentum()[2];

 
   //    pT2=Px2+Py2;
    //   P=std::sqrt(std::pow(Px,2)+std::pow(Py,2)+std::pow(Pz,2));
       float E1 = j -> getEnergy();
       float E2 = jet2 -> getEnergy();
      // float P2 = std::sqrt(std::pow(Px_2,2)+std::pow(Py_2,2)+std::pow(Pz_2,2));
       invmass = std::sqrt((E1+E2)*(E1+E2)-(px*px)-(py*py)-(pz*pz));      
      
       //streamlog_out(WARNING)  << "ivariant mass= " << invmass <<std::endl;
       if (i_jet == 0){

      // streamlog_out(WARNING)  << "setting ivariant mass (step1)= " << invmass <<std::endl;
       ml[0] = invmass;}

     }//end of   if (reco particle j!=0)
 }//end loop over nJets

#ifdef MARLIN_USE_AIDA
 double en, px, py, pz, mass;
 if (jrps[0] && jrps[1]) {
   en = jrps[0]->getEnergy()     +jrps[1]->getEnergy();
   px = jrps[0]->getMomentum()[0]+jrps[1]->getMomentum()[0];
   py = jrps[0]->getMomentum()[1]+jrps[1]->getMomentum()[1];
   pz = jrps[0]->getMomentum()[2]+jrps[1]->getMomentum()[2];
   mass = en*en-px*px-py*py-pz*pz;
   if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
 }
#endif


 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////
 ////Set LeptonFitObjects
 ////
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////

 for(int i_lep = 0; i_lep < nLEPS ; i_lep++)
   {
     ReconstructedParticle* l = dynamic_cast<ReconstructedParticle*>( lepcol->getElementAt( i_lep ) ) ;
     leptonvec = HepLorentzVector ((l->getMomentum())[0],(l->getMomentum())[1],(l->getMomentum())[2],l->getEnergy());
     TrackVec tckvec = l->getTracks();
 
       Px=l->getMomentum()[0];
       Px2=std::pow(Px,2);
       Py=l->getMomentum()[1];
       Py2=std::pow(Py,2);
       Pz=l->getMomentum()[2];
       Pz2=std::pow(Pz,2);
 
       ReconstructedParticle* lep2 = dynamic_cast<ReconstructedParticle*>( lepcol->getElementAt( 1 ) ) ;
      // float Px_2=lep2->getMomentum()[0];
      // float Py_2=lep2->getMomentum()[1];
      // float Pz_2=lep2->getMomentum()[2];
       
       double px = l->getMomentum()[0]+lep2->getMomentum()[0];
       double py = l->getMomentum()[1]+lep2->getMomentum()[1];
       double pz = l->getMomentum()[2]+lep2->getMomentum()[2];

 
       pT2=Px2+Py2;
       P=std::sqrt(std::pow(Px,2)+std::pow(Py,2)+std::pow(Pz,2));
       float E1 = l -> getEnergy();
       float E2 = lep2 -> getEnergy();
      // float P2 = std::sqrt(std::pow(Px_2,2)+std::pow(Py_2,2)+std::pow(Pz_2,2));
       invmass = std::sqrt((E1+E2)*(E1+E2)-(px*px)-(py*py)-(pz*pz));      
      
       //streamlog_out(WARNING)  << "ivariant mass= " << invmass <<std::endl;
       if (i_lep == 0){

      // streamlog_out(WARNING)  << "setting ivariant mass (step1)= " << invmass <<std::endl;
       ml[1] = invmass;}
/*      
       SigPx2     = l->getCovMatrix()[ 0 ];
       SigPxSigPy = l->getCovMatrix()[ 1 ];
       SigPxSigPz = l->getCovMatrix()[ 3 ];
       SigPy2     = l->getCovMatrix()[ 2 ];
       SigPySigPz = l->getCovMatrix()[ 4 ];
       SigPz2     = l->getCovMatrix()[ 5 ];
       SigE2      = l->getCovMatrix()[ 9 ];
       
       dth_dpx    = ( Px * Pz ) / ( std::pow( P , 2 ) * std::sqrt( pT2 ) );
       dth_dpy    = ( Py * Pz ) / ( std::pow( P , 2 ) * std::sqrt( pT2 ) );
       dth_dpz    = -std::sqrt( pT2 ) / std::pow( P , 2 );
       dphi_dpx   = -Py / pT2;
       dphi_dpy   = Px / pT2 ;
       
       LeptonResE    = std::sqrt( SigE2 ) * sigmaScaleFactor;
       LeptonResTheta= std::sqrt( std::fabs( SigPx2 * std::pow( dth_dpx , 2 ) + SigPy2 * std::pow( dth_dpy , 2 ) + SigPz2 * std::pow( dth_dpz , 2 ) + 2 * ( SigPxSigPy * dth_dpx * dth_dpy ) + 2 * ( SigPySigPz * dth_dpy * dth_dpz ) + 2 * ( SigPxSigPz * dth_dpx * dth_dpz ) ) );
       LeptonResPhi  = std::sqrt( std::fabs( SigPx2 * std::pow( dphi_dpx , 2 ) + SigPy2 * std::pow( dphi_dpy , 2 ) + 2 * ( SigPxSigPy * dphi_dpx * dphi_dpy ) ) );
       
       streamlog_out(DEBUG4)  << "SigPx2= " << SigPx2 <<std::endl;
       streamlog_out(DEBUG4)  << "SigPxSigPy= " << SigPxSigPy <<std::endl;
       streamlog_out(DEBUG4)  << "SigPxSigPz= " << SigPxSigPz <<std::endl;
       streamlog_out(DEBUG4)  << "SigPy2= " <<  SigPy2<<std::endl;
       streamlog_out(DEBUG4)  << "SigPySigPz= " << SigPySigPz <<std::endl;
       streamlog_out(DEBUG4)  << "SigPz2= " << SigPz2 <<std::endl;
       streamlog_out(DEBUG4)  << "SigE2= " << SigE2 <<std::endl;
       
       streamlog_out(DEBUG)  << "Px= " << Px <<std::endl;
       streamlog_out(DEBUG)  << "Py= " << Py <<std::endl;
       streamlog_out(DEBUG)  << "Pz= " << Pz <<std::endl;
       
       streamlog_out(DEBUG)  << "LeptonResE= " << LeptonResE <<std::endl;
       streamlog_out(DEBUG)  << "LeptonResTheta= " << LeptonResTheta <<std::endl;
       streamlog_out(DEBUG)  << "LeptonResPhi= " << LeptonResPhi <<std::endl;
*/  
      theta_err = LeptonResTheta;
      phi_err = LeptonResPhi;
      invers_pT_err = LeptonResE;
     if ( tckvec.size() != 1 ){
	 streamlog_out(DEBUG)  << "Number of tracks for lepton[" << i_lep <<"] is not exactly ONE!!! (nTracks = " << tckvec.size() << " ) " << std::endl ;
	 invers_pT = 1/(std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2)));
	 invers_pT_err = 2*std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2))*0.00001;
	 theta = leptonvec.theta();
	 theta_err = 0.1;
	 phi = leptonvec.phi();
	 phi_err = 0.1;
	 trackcounter ++;
       }
     else
       {
	 streamlog_out(DEBUG)  << "Number of tracks for lepton[" << i_lep <<"] is exactly ONE!!!" << std::endl;
	 Omega = tckvec[0]->getOmega();
	 Omega_uncert = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[5]) );
	 streamlog_out(DEBUG)  << "Omega = " << Omega << std::endl;
	 streamlog_out(DEBUG)  << "Omega_uncert = " << Omega_uncert << std::endl;
	 TanLambda = tckvec[0]->getTanLambda();
	 TanLambda_err = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[14]) );
	 streamlog_out(DEBUG)  << "TanLambda_err = " << TanLambda_err << std::endl;
	 theta = 2 * atan( 1. ) - atan( TanLambda );
	 theta_err = TanLambda_err / ( 1 + pow( TanLambda , 2 ) );
	 streamlog_out(DEBUG)  << "theta = " << theta << std::endl;
	 streamlog_out(DEBUG)  << "theta_err = " << theta_err << std::endl;
	 phi = tckvec[0]->getPhi();
	 phi_err = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[2]) );
	 streamlog_out(DEBUG)  << "phi = " << phi << std::endl;
	 streamlog_out(DEBUG)  << "phi_err = " << phi_err << std::endl;
	 invers_pT = Omega / eB;
	 invers_pT_err = std::fabs( 1. / eB ) * Omega_uncert;
	 streamlog_out(DEBUG)  << "invers_pT = " << invers_pT << std::endl;
	 streamlog_out(DEBUG)  << "invers_pT_err = " << invers_pT_err << std::endl;
       }
     streamlog_out(DEBUG)  << "Lepton fit object from leptonvec: "
			   << 1/(std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2))) <<" +- " << 2*std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2))*0.00001 << " , "
			   << leptonvec.theta() <<" +- " << 0.1 << " , "
			   << leptonvec.phi() <<" +- " << 0.1 << std::endl ;

     streamlog_out(DEBUG)  << "Lepton fit object from track:     "
			   << std::fabs( tckvec[0]->getOmega() / eB ) <<" +- " << std::fabs( 1. / eB ) * std::sqrt( tckvec[0]->getCovMatrix()[5] ) << " , "
			   << 2 * atan( 1. ) - atan( tckvec[0]->getTanLambda() ) <<" +- " << std::abs( std::sqrt( tckvec[0]->getCovMatrix()[14]) ) / ( 1 + pow( tckvec[0]->getTanLambda() , 2 ) ) << " , "
			   << tckvec[0]->getPhi() <<" +- " << std::abs( std::sqrt( tckvec[0]->getCovMatrix()[2] ) ) << std::endl ;

    // lepton[i_lep] = new LeptonFitObject (invers_pT , theta , phi , invers_pT_err , theta_err , phi_err, leptonvec.m());
     if (i_lep == 0 )
       {
	 l1 = new LeptonFitObject (invers_pT , theta , phi , invers_pT_err , theta_err , phi_err, leptonvec.m());
	 l1->setName("Lepton1");
	 streamlog_out(DEBUG)  << " start four-vector of first  lepton: " << *l1  << std::endl ;
       }
     else if (i_lep == 1)
       {
	 l2 = new LeptonFitObject (invers_pT , theta , phi , invers_pT_err , theta_err , phi_err, leptonvec.m());
         l2->setName("Lepton2");
         streamlog_out(DEBUG)  << " start four-vector of second  lepton: " << *l2  << std::endl ;
       }
   }



 const int NJETS = 2;
 streamlog_out(MESSAGE)  << "*j1" << *j1  << "*j2" << *j2 << std::endl ;
 const int NLEPTONS = 2;
 streamlog_out(MESSAGE)  << "*l1" << *l1  << "*l2" << *l2 << std::endl ;

 ////these don't get changed by the fit -> to obtain start values later!
 JetFitObject startjets[NJETS] = {*j1, *j2};
 for (int i = 0; i < NJETS; ++i)streamlog_out(DEBUG)  << "startjets[" << i << "]: " << startjets[i]  << std::endl;

 LeptonFitObject startleptons[NLEPTONS] = {*l1, *l2};
 for (int i = 0; i < NLEPTONS; ++i) streamlog_out(DEBUG)  << "startleptons[" << i << "]: " << startleptons[i]  << std::endl;

 // these get changed by the fit
 JetFitObject fitjets[NJETS] = {*j1, *j2};
 for (int i = 0; i < NJETS; ++i)
   streamlog_out(DEBUG)  << "fitjets[" << i << "]: " << fitjets[i]  << std::endl ;

 LeptonFitObject fitleptons[NLEPTONS] = {*l1, *l2};
 for (int i = 0; i < NLEPTONS; ++i)
   streamlog_out(DEBUG)  << "fitleptons[" << i << "]: " << fitleptons[i]  << std::endl ;

 // these point allways to the fitjets array, which gets reset.
 JetFitObject *jets[NJETS];
 for (int i = 0; i < NJETS; ++i)
   {
     jets[i] = &fitjets[i];
     streamlog_out(MESSAGE)  << "start four-vector of jet " << i << ": " << *(jets[i])  << std::endl ;
   }

 LeptonFitObject *leptons[NLEPTONS];
 for (int i = 0; i < NLEPTONS; ++i)
   {

     leptons[i] = &fitleptons[i];
     streamlog_out(MESSAGE)  << "start four-vector of leptons " << i << ": " << *(leptons[i])  << std::endl ;
   }
 

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////
 ////set constraints befor fit
 ////
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////

 float target_p_due_crossing_angle = _ecm * 0.007; // crossing angle = 14 mrad
 MomentumConstraint pxc ( 0 , 1 , 0 , 0 , target_p_due_crossing_angle );//Factor for: (energy sum, px sum, py sum,pz sum,target value of sum)
 
 pxc.setName("sum(p_x)");
 for (int i = 0; i < NJETS; ++i) pxc.addToFOList(*(jets[i]));
 for (int i = 0; i < NLEPTONS; ++i) pxc.addToFOList(*(leptons[i]));

 MomentumConstraint pyc (0, 0, 1, 0, 0);
 pyc.setName("sum(p_y)");
 for (int i = 0; i < NJETS; ++i) pyc.addToFOList(*(jets[i]));
 for (int i = 0; i < NLEPTONS; ++i) pyc.addToFOList(*(leptons[i]));

 MomentumConstraint pzc (0, 0, 0, 1, 0);
 pzc.setName("sum(p_z)");
 for (int i = 0; i < NJETS; ++i) pzc.addToFOList(*(jets[i]));
 for (int i = 0; i < NLEPTONS; ++i) pzc.addToFOList(*(leptons[i]));

 E_lab= 2 * sqrt( std::pow( 0.548579909e-3 , 2 ) + std::pow( _ecm / 2 , 2 ) + std::pow( target_p_due_crossing_angle , 2 ) + 0. + 0.);
 MomentumConstraint ec(1, 0, 0, 0, E_lab);
 ec.setName("sum(E)");
 for (int i = 0; i < NJETS; ++i) ec.addToFOList(*(jets[i]));
 for (int i = 0; i < NLEPTONS; ++i) ec.addToFOList(*(leptons[i]));


 streamlog_out(DEBUG)  << "Value of E_lab before adding ISR: " << E_lab << std::endl ;
 streamlog_out(DEBUG)  << "Value of target_p_due_crossing_angle before adding ISR: " << target_p_due_crossing_angle << std::endl ;
 streamlog_out(DEBUG)  << "Value of pxc before adding ISR: " << pxc.getValue() << std::endl ;
 streamlog_out(DEBUG)  << "Value of pyc before adding ISR: " << pyc.getValue() << std::endl ;
 streamlog_out(DEBUG)  << "Value of pzc before adding ISR: " << pzc.getValue() << std::endl ;
 streamlog_out(DEBUG)  << "Value of ec before adding ISR: " << ec.getValue() << std::endl ;
 pxc_before_ISR = pxc.getValue();
 pyc_before_ISR = pyc.getValue();
 pzc_before_ISR = pzc.getValue();
 ec_before_ISR = ec.getValue();


 // ISR Photon initialized with missing p_z
 ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);
 float ISRstartPx = photon->getPx();
 float ISRstartPy = photon->getPy();
 float ISRstartPz = photon->getPz();
 //ISRPhotonFitObject(double px, double py, double pz,
 //double b_, double PzMaxB_, double PzMinB_ = 0.);

 if(_fitISR){
   streamlog_out(MESSAGE)  << "start four-vector of ISR photon: " << *(photon) << std::endl ;
   
   pxc.addToFOList (*(photon));
   pyc.addToFOList (*(photon));
   pzc.addToFOList (*(photon));
   ec.addToFOList  (*(photon));
 }

 streamlog_out(DEBUG)  << "Value of E_lab before fit: " << E_lab << std::endl ;
 streamlog_out(DEBUG)  << "Value of target_p_due_crossing_angle before fit: " << target_p_due_crossing_angle << std::endl ;
 streamlog_out(DEBUG)  << "Value of pxc after adding ISR before fit: " << pxc.getValue() << std::endl ;
 streamlog_out(DEBUG)  << "Value of pyc after adding ISR before fit: " << pyc.getValue() << std::endl ;
 streamlog_out(DEBUG)  << "Value of pzc after adding ISR before fit: " << pzc.getValue() << std::endl ;
 streamlog_out(DEBUG)  << "Value of ec after adding ISR before fit: " << ec.getValue() << std::endl ;
 pxc_before_fit = pxc.getValue();
 pyc_before_fit = pyc.getValue();
 pzc_before_fit = pzc.getValue();
 ec_before_fit = ec.getValue();
 
  MassConstraint z(91.2);
// SoftGaussMassConstraint z(91.2,2.4952/2);
 z.addToFOList (*(leptons[0]), 1);
 z.addToFOList (*(leptons[1]), 1);
 
 MassConstraint h(125.);
 h.addToFOList (*(jets[0]), 1);
 h.addToFOList (*(jets[1]), 1);
 streamlog_out(MESSAGE)  << "final mass of Z: " << z.getMass(1) << std::endl ;
 
 startmassZ = z.getMass(1);
 startmassH = h.getMass(1);
 
 streamlog_out(MESSAGE) << "start mass of Z: " << startmassZ << std::endl ;
 streamlog_out(MESSAGE) << "start mass of H: " << startmassH << std::endl ;
 
 Hmass_NoFit=startmassH;
 
 streamlog_out(DEBUG4) << "before AIDA"  << std::endl ;
#ifdef MARLIN_USE_AIDA
 hRecHMassNoFitAll->fill( startmassH ) ;
 //hRecHMassNoFitAll->fill( startmassH ) ;
#endif

 streamlog_out(DEBUG4) << "after AIDA"  << std::endl ;
 BaseFitter *pfitter;

 if (_ifitter == 1) {
   pfitter = new NewFitterGSL();
   if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<NewFitterGSL*>(pfitter))->setDebug (debug);
   
   streamlog_out(DEBUG4) << "ifitter is 1"  << std::endl ;
 }
 else if (_ifitter == 2) {
   pfitter = new NewtonFitterGSL();
   if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<NewtonFitterGSL*>(pfitter))->setDebug (debug);
   
   streamlog_out(DEBUG4) << "ifitter is 2"  << std::endl ;
 }
 else {
   // OPALFitter has no method setDebug !
   pfitter = new OPALFitterGSL();
   
   streamlog_out(DEBUG4) << "ifitter is not 1 or 2"  << std::endl ;
   if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<OPALFitterGSL*>(pfitter))->setDebug (debug);
 }
 BaseFitter &fitter = *pfitter;
 
 TextTracer tracer (std::cout);
 if (evt->getEventNumber()== _ievttrace || _traceall) fitter.setTracer (tracer);
 

 for (int i = 0; i < NJETS; ++i) fitter.addFitObject(*(jets[i]));
 for (int i = 0; i < NLEPTONS; ++i) fitter.addFitObject(*(leptons[i]));

 if(_fitISR){
   
   streamlog_out(DEBUG4) << "isr added to fit"  << std::endl ;
   fitter.addFitObject (*(photon));
 }
 fitter.addConstraint (pxc);
 fitter.addConstraint (pyc);
 fitter.addConstraint (pzc);
 fitter.addConstraint (ec);
 fitter.addConstraint (z);
 
 streamlog_out(DEBUG4) << "constraints added"  << std::endl ;
 // don't constrain Higgs mass, just use constraints for convenient mass calculation
 //fitter.addConstraint (h);
 // initial value of Z mass constraint
 if (fabs(startmassZ-91.2) + fabs(startmassH-125.) < bestzvalue) {
   chi2startmassZ = startmassZ;
   chi2startmassH = startmassH;
   bestzvalue = fabs(startmassZ-91.2) + fabs(startmassH-125.);
   
   streamlog_out(DEBUG4) << "best z value is this..." <<  bestzvalue << std::endl ;
 }
 
 
 prob = fitter.fit();
 double chi2 = fitter.getChi2();
 nit = fitter.getIterations();
 
 streamlog_out(DEBUG4) << "fit probability = " << prob << std::endl ;
 streamlog_out(DEBUG4) << "fit chi2 = " << chi2  << std::endl ;
 streamlog_out(DEBUG4) << "error code: " << fitter.getError() << std::endl ;
 
 for (int i = 0; i < NJETS; ++i) {
   streamlog_out(MESSAGE)  << "final four-vector of jet " << i << ": " << *(jets[i]) << std::endl ;
   streamlog_out(MESSAGE)  << "final px of jet " << i << ": " << (jets[i]) << std::endl ;
 }
 for (int i = 0; i < NLEPTONS; ++i) {
   streamlog_out(MESSAGE)  << "final four-vector of lepton " << i << ": " << *(leptons[i]) << std::endl ;
   streamlog_out(DEBUG)  << "final px of lepton " << i << ": " << (leptons[i]) << std::endl ;
 }

 if(_fitISR){
   streamlog_out(MESSAGE)  << "final four-vector of ISR photon: " << *(photon) << std::endl ;
 }
 
 int ierr = fitter.getError();
 streamlog_out(MESSAGE)  << "fitter error: " << ierr << std::endl;
 hFitError->fill( ierr ) ;
 if ((besterr > 0 && ierr < besterr) || ( besterr < 0 && ierr == 0)) besterr = ierr;
 
 streamlog_out(DEBUG)  << "Value of pxc after fit: " << pxc.getValue() << std::endl ;
 streamlog_out(DEBUG)  << "Value of pyc after fit: " << pyc.getValue() << std::endl ;
 streamlog_out(DEBUG)  << "Value of pzc after fit: " << pzc.getValue() << std::endl ;
 streamlog_out(DEBUG)  << "Value of ec after fit: " << ec.getValue() << std::endl ;
 pxc_after_fit = pxc.getValue();
 pyc_after_fit = pyc.getValue();
 pzc_after_fit = pzc.getValue();
 ec_after_fit = ec.getValue();

 bestphotonenergy = photon->getE();

 ISRmomentum[0]= photon->getPx();
 ISRmomentum[1]= photon->getPx();
 ISRmomentum[2]= photon->getPx();

 streamlog_out(DEBUG)  << "Value of ISRmomentum x after fit: " << ISRmomentum[0] << std::endl ;
 streamlog_out(DEBUG)  << "Value of ISRmomentum y after fit: " << ISRmomentum[1] << std::endl ;
 streamlog_out(DEBUG)  << "Value of ISRmomentum z after fit: " << ISRmomentum[2] << std::endl ;



 Zmomentum[0] = leptons[0]->getPx() + leptons[1]->getPx();
 Zmomentum[1] = leptons[0]->getPy() + leptons[1]->getPy();
 Zmomentum[2] = leptons[0]->getPz() + leptons[1]->getPz();



 streamlog_out(DEBUG)  << "Value of Zmomx after fit: " << Zmomentum[0] << std::endl ;
 streamlog_out(DEBUG)  << "Value of Zmomy after fit: " << Zmomentum[1] << std::endl ;
 streamlog_out(DEBUG)  << "Value of Zmomz after fit: " << Zmomentum[2] << std::endl ;




 Hmomentum[0] = jets[0]->getPx() + jets[1]->getPx();
 Hmomentum[1] = jets[0]->getPy() + jets[1]->getPy();
 Hmomentum[2] = jets[0]->getPz() + jets[1]->getPz();
 Z_Energy = leptons[0]->getE() + leptons[1]->getE();
 H_Energy = jets[0]->getE() + jets[1]->getE();
 chi2best= fitter.getChi2();
 errorcode= fitter.getError();


 float J1momentum[3];
 float J0momentum[3];
 float L1momentum[3];
 float L0momentum[3];

// float Jpull0[3];
// float Jpull1[3];


EVENT::LCFloatVec*  PullE = new LCFloatVec();
EVENT::LCFloatVec*  PullTheta = new LCFloatVec();
EVENT::LCFloatVec*  PullPhi = new LCFloatVec();
/*
EVENT::LCFloatVec*  PullE = new LCFloatVec();
EVENT::LCFloatVec*  PullTheta = new LCFloatVec();
EVENT::LCFloatVec*  PullPhi = new LCFloatVec();
*/

//streamlog_out(DEBUG)  << "Creating LCGenericObjectImpl(0,3,0) " << std::endl ;
//LCGenericObjectImpl*  Jpull0 = new LCGenericObjectImpl(0,3,0);
//streamlog_out(DEBUG)  << "Creating LCGenericObjectImpl(0,3,0) " << std::endl ;
//LCGenericObjectImpl*  Jpull1 = new LCGenericObjectImpl(0,3,0);
//streamlog_out(DEBUG)  << "Done with both "<< std::endl ;

 ReconstructedParticle* lepton0 = dynamic_cast<ReconstructedParticle*>( lepcol->getElementAt( 0 ) ) ;
 ReconstructedParticle* lepton1 = dynamic_cast<ReconstructedParticle*>( lepcol->getElementAt( 1 ) ) ;



 L1momentum[0] = leptons[1]->getPx();
 L1momentum[1] = leptons[1]->getPy();
 L1momentum[2] = leptons[1]->getPz();

 L0momentum[0] = leptons[0]->getPx();
 L0momentum[1] = leptons[0]->getPy();
 L0momentum[2] = leptons[0]->getPz();

 J1momentum[0] = jets[1]->getPx();
 J1momentum[1] = jets[1]->getPy();
 J1momentum[2] = jets[1]->getPz();

 J0momentum[0] = jets[0]->getPx();
 J0momentum[1] = jets[0]->getPy();
 J0momentum[2] = jets[0]->getPz();


 streamlog_out(DEBUG)  << "Value of L0 momentum after fit: " << L0momentum[0] << "," << L0momentum[1] << ","  << L0momentum[2]  << std::endl ;
 streamlog_out(DEBUG)  << "Value of L1 momentum after fit: " << L1momentum[0] << "," << L1momentum[1] << ","  << L1momentum[2]  << std::endl ;
 streamlog_out(DEBUG)  << "Value of J0 momentum after fit: " << J0momentum[0] << "," << J0momentum[1] << ","  << J0momentum[2]  << std::endl ;
 streamlog_out(DEBUG)  << "Value of J1 momentum after fit: " << J1momentum[0] << "," << J1momentum[1] << ","  << J1momentum[2]  << std::endl ;



 // ierr == -1 only means that error calculation for fitted parameters failed!
 if (ierr <= 0) { //if fitter.getError() <=0
#ifdef MARLIN_USE_AIDA
   hFitProbAll->fill( prob ) ;
   hNItAll->fill( nit ) ;
   hRecHMassAll->fill( h.getMass(1) ) ;
#endif
   double pullJet[3][2];
   double pullLepton[3][2];
   // require successfull error calculation for pulls!
   if (ierr == 0) {
     for (int ifo = 0; ifo < 2; ifo++){
       double start, fitted;
       double errfit, errmea, sigma;
       for (int ipar = 0; ipar < 3; ipar++) {
	 fitted = jets[ifo]->getParam(ipar);
	 start = startjets[ifo].getParam(ipar);
	 errfit = jets[ifo]->getError(ipar);
	 errmea = startjets[ifo].getError(ipar);
	 sigma = errmea*errmea-errfit*errfit;
	 if (sigma > 0) {
	   sigma = sqrt(sigma);
	   pullJet[ipar][ifo] = (fitted - start)/sigma;
	 }
	 else {
	   pullJet[ipar][ifo] = -4.5;
	 }
       }
     }
     for (int ifo = 0; ifo < 2; ifo++) {
       double start, fitted;
       double errfit, errmea, sigma;
       for (int ipar = 0; ipar < 3; ipar++) {
	 fitted = leptons[ifo]->getParam(ipar);
	 start = startleptons[ifo].getParam(ipar);
	// startPT =  sqrt(px**2+py**2)
         errfit = leptons[ifo]->getError(ipar);
	 errmea = startleptons[ifo].getError(ipar);
	 sigma = errmea*errmea-errfit*errfit;
	 if (sigma > 0) {
	   sigma = sqrt(sigma);
	   pullLepton[ipar][ifo] = (fitted - start)/sigma;
	 }
	 else {
	   pullLepton[ipar][ifo] = -4.5;
	 }
       }
     }
   }


// Jpulls [Energy,theta,phi] Jpull0[0]=pullJet[0][0];




/*

streamlog_out(DEBUG)  << "Setting Values in GenObj1 " << std::endl ;

  Jpull0[0] =pullJet[0][0];
  Jpull0[1] =pullJet[1][0];
  Jpull0[2] =pullJet[2][0];

  Jpull1[0] =pullJet[0][1];
  Jpull1[1] =pullJet[1][1];
  Jpull1[2] =pullJet[2][1];
*/


/*
  PullE -> setFloatVal(0,pullJet[0][0]);
  PullTheta -> setFloatVal(0,pullJet[1][0]);
  PullPhi -> setFloatVal(0,pullJet[2][0]);

streamlog_out(DEBUG)  << "Setting Values in GenObj2 " << std::endl ;

  PullE -> setFloatVal(1,pullJet[0][1]);
  PullTheta -> setFloatVal(1,pullJet[1][1]);
  PullPhi -> setFloatVal(1,pullJet[2][1]);
*/
 PullE -> push_back(pullJet[0][0]);
 PullTheta -> push_back(pullJet[1][0]);
 PullPhi -> push_back(pullJet[2][0]);

 PullE -> push_back(pullJet[0][1]);
 PullTheta -> push_back(pullJet[1][1]);
 PullPhi -> push_back(pullJet[2][1]);

 PullE -> push_back(pullLepton[0][0]);
 PullTheta -> push_back(pullLepton[1][0]);
 PullPhi -> push_back(pullLepton[2][0]);

 PullE -> push_back(pullLepton[0][1]);
 PullTheta -> push_back(pullLepton[1][1]);
 PullPhi -> push_back(pullLepton[2][1]);



   //if (prob > bestprob && h.getMass(1) > 70 && w.getMass(1) < 150) {
   if (prob > bestprob) {
     bestprob = prob;
     streamlog_out(DEBUG)  << "BESTPROB: " << bestprob << std::endl ;
     bestnit  = nit;
     bestmassZ = z.getMass(1);
     bestmassH = h.getMass(1);
     beststartmassZ = startmassZ;
     beststartmassH = startmassH;
     bestphotonenergy = photon->getE();
     ISRmomentum [0]= photon->getPx();
     ISRmomentum [1]= photon->getPx();
     ISRmomentum [2]= photon->getPx();
     Zmomentum[0] = leptons[0]->getPx() + leptons[1]->getPx();
     Zmomentum[1] = leptons[0]->getPy() + leptons[1]->getPy();
     Zmomentum[2] = leptons[0]->getPz() + leptons[1]->getPz();
     Hmomentum[0] = jets[0]->getPx() + jets[1]->getPx();
     Hmomentum[1] = jets[0]->getPy() + jets[1]->getPy();
     Hmomentum[2] = jets[0]->getPz() + jets[1]->getPz();
     Z_Energy = leptons[0]->getE() + leptons[1]->getE();
     H_Energy = jets[0]->getE() + jets[1]->getE();
     chi2best= fitter.getChi2();
     errorcode= fitter.getError();
     if (ierr == 0) { //if  fitter.getError() is = 0
       hpull_jet1_E=pullJet[0][0];
       hpull_jet2_E=pullJet[0][1];
       hpull_jet1_th=pullJet[1][0];
       hpull_jet2_th=pullJet[1][1];
       hpull_jet1_phi=pullJet[2][0];
       hpull_jet2_phi=pullJet[2][1];
       hpull_lepton1_InvpT=pullLepton[0][0];
       hpull_lepton2_InvpT=pullLepton[0][1];
       hpull_lepton1_th=pullLepton[1][0];
       hpull_lepton2_th=pullLepton[1][1];
       hpull_lepton1_phi=pullLepton[2][0];
       hpull_lepton2_phi=pullLepton[2][1];
     }//end if  fitter.getError() is = 0
     else {//if  fitter.getError() is not = 0
       streamlog_out(DEBUG) << " ERROR CALCULATION FAILED for best permutation "
			    << " in event " << evt->getEventNumber()
			    << std::endl ;
     }
   }
 }//end-if fitter.getError() <=0
 else {
   streamlog_out(DEBUG4) << "FIT ERROR = " << fitter.getError()
			 << " in event " << evt->getEventNumber()
			 << ", not filling histograms!"  << std::endl ;
   streamlog_out(MESSAGE)  << "start mass of Z: " << startmassZ << std::endl ;
   streamlog_out(DEBUG4)  << "start mass of H: " << startmassH << std::endl ;
   streamlog_out(DEBUG4)  << "final mass of Z: " << z.getMass(1) << std::endl ;
   streamlog_out(DEBUG4)  << "final mass of H: " << h.getMass(1) << std::endl ;
   
 }
 delete photon;

//*******write something

LCCollectionVec *PullEOutCol = new LCCollectionVec(LCIO::LCFLOATVEC);
LCCollectionVec *PullThetaOutCol = new LCCollectionVec(LCIO::LCFLOATVEC);
LCCollectionVec *PullPhiOutCol = new LCCollectionVec(LCIO::LCFLOATVEC);





 LCCollectionVec *OutputCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
 ReconstructedParticleImpl* ISRfitrec = new ReconstructedParticleImpl;
 ReconstructedParticleImpl* Zfitrec = new ReconstructedParticleImpl;
 ReconstructedParticleImpl* Hfitrec = new ReconstructedParticleImpl;
 ReconstructedParticleImpl* L1fitrec = new ReconstructedParticleImpl;
 ReconstructedParticleImpl* L0fitrec = new ReconstructedParticleImpl; 
 ReconstructedParticleImpl* J1fitrec = new ReconstructedParticleImpl;
 ReconstructedParticleImpl* J0fitrec = new ReconstructedParticleImpl;
// ReconstructedParticleImpl* Jpull1fitrec = new ReconstructedParticleImpl;
//ReconstructedParticleImpl* Jpull0fitrec = new ReconstructedParticleImpl;



 ISRfitrec->setMomentum(ISRmomentum);
 streamlog_out(DEBUG) << "  ISRmomentum :   "
                      << ISRfitrec->getMomentum()[0] << "," << ISRfitrec->getMomentum()[1]<<","<< ISRfitrec->getMomentum()[2] << std::endl ;
 

 ISRfitrec->setEnergy(bestphotonenergy);
 streamlog_out(DEBUG) << " Energy ISR:   "
		      << ISRfitrec->getEnergy() << std::endl ;
 ISRfitrec->setType (22);
 streamlog_out(DEBUG) << " IS ISR:   "
		      << ISRfitrec->getType() << std::endl ;
 OutputCol->addElement(ISRfitrec);


 ISRfitrec->setCharge(bestprob);

 
 Zfitrec->setMomentum(Zmomentum);
 streamlog_out(DEBUG) << "  Zmomentum :   "
		      << Zfitrec->getMomentum()[0] << "," << Zfitrec->getMomentum()[1]<<","<< Zfitrec->getMomentum()[2] << std::endl ;
 
 Zfitrec->setEnergy(Z_Energy);
 streamlog_out(DEBUG) << " Energy Z:   "
		      << Zfitrec->getEnergy() << std::endl ;
 Zfitrec->setMass(bestmassZ);
 streamlog_out(DEBUG) << " Mass Z:   "
		      << Zfitrec->getMass() << std::endl ;


 Zfitrec -> setCharge(ml[1]);

 Zfitrec->setType (23);
 streamlog_out(DEBUG) << " IS Z :   "
		      << Zfitrec->getType() << std::endl ;
 OutputCol->addElement(Zfitrec);
 
 
 Hfitrec->setMomentum(Hmomentum);
 Hfitrec->setEnergy(H_Energy);
 streamlog_out(DEBUG) << " Energy H:   "
		      << Hfitrec->getEnergy() << std::endl ;
 Hfitrec->setMass(bestmassH);
 streamlog_out(DEBUG) << " Mass H:   "
		      << Hfitrec->getMass() << std::endl ;

 Hfitrec->setCharge(ml[0]);
 Hfitrec->setGoodnessOfPID(bestprob);

 Hfitrec->setType (25);
 streamlog_out(DEBUG) << " IS H:   "
		      << Hfitrec->getType() << std::endl ;

 OutputCol->addElement(Hfitrec); 
 
//L0fitrec->setCovMatrix(leptons[0] -> getCovMatrix())
L0fitrec->setMomentum(L0momentum);
 L0fitrec->setType (13);
 L0fitrec->setEnergy(leptons[0]->getE());
 //L0fitrec->setMomentum(L0momentum);

 L0fitrec -> setCovMatrix(lepton0 -> getCovMatrix());
 //L0fitrec -> setCharge(ml[0]);
// L0fitrec -> setParticleIDUsed(lepton0 -> id());
// L0fitrec -> setCharge(lepton0 -> getCharge());
// L0fitrec -> setCharge(lepton0 -> getCharge());

 OutputCol->addElement(L0fitrec);

 L1fitrec->setMomentum(L1momentum);
 L1fitrec->setType (13);
 L1fitrec->setEnergy(leptons[1]->getE());
// L1fitrec->setMomentum(L1momentum);
 L1fitrec -> setCovMatrix(lepton1 -> getCovMatrix());
 OutputCol->addElement(L1fitrec);

 J0fitrec->setMomentum(J0momentum);
 J0fitrec->setType (5);
 J0fitrec->setEnergy(jets[0]->getE());
 //J0fitrec->setMomentum(J0momentum);
 OutputCol->addElement(J0fitrec);
 
 J1fitrec->setMomentum(J1momentum);
 J1fitrec->setType (5);
 J1fitrec->setEnergy(jets[1]->getE());
 //J1fitrec->setMomentum(J1momentum);
 OutputCol->addElement(J1fitrec);


//Jpull0fitrec -> setMomentum(Jpull0);
//Jpull0fitrec -> setType(8888);
//OutputCol->addElement(Jpull0fitrec);
//streamlog_out(DEBUG)  << "Adding Obj1 to Col " << std::endl ;

PullEOutCol->addElement(PullE);
PullThetaOutCol->addElement(PullTheta);
PullPhiOutCol->addElement(PullPhi);



//Jpull1fitrec -> setMomentum(Jpull1);
//Jpull1fitrec -> setType(8888);
//OutputCol->addElement(Jpull1fitrec);
//
//streamlog_out(DEBUG)  << "Adding Obj1 to Col " << std::endl ;
//PullsOutputCol->addElement(Jpull1);

 OutputCol->parameters().setValue("fitprob", (float)prob);

// OutputCol->parameters().setValue("pull0energy",(float)Jpull0[0]); 
// OutputCol->parameters().setValue("pull0theta", (float)Jpull0[1]);
//OutputCol->parameters().setValue("pull0phi", (float)Jpull0[2]);
// OutputCol->parameters().setValue("pull1energy", (float)Jpull1[0]);
// OutputCol->parameters().setValue("pull1theta", (float)Jpull1[1]);
// OutputCol->parameters().setValue("pull1phi",(float)Jpull1[2]);


 OutputCol->parameters().setValue("bestchisq", (float)chi2best);
 streamlog_out(DEBUG) << " chi2:   " << chi2best << std::endl ;
 OutputCol->parameters().setValue("best_prob", (float)bestprob);
 streamlog_out(DEBUG) << " prob:   " << bestprob << std::endl ;
 OutputCol->parameters().setValue("error_code", (float)errorcode);
 streamlog_out(DEBUG) << "Error Code:   " << errorcode << std::endl ;
 
 Hmass_before_fit=beststartmassH;
 Hmass_after_fit=bestmassH;
 Error_code=errorcode;
 
 streamlog_out(DEBUG) << "==============  end of permutations for event " << evt->getEventNumber() <<  " ==============" << std::endl ;
 streamlog_out(DEBUG)  << "min chi2 start mass of Z: " << chi2startmassZ << std::endl ;
 streamlog_out(DEBUG)  << "min chi2 start mass of H: " << chi2startmassH << std::endl ;
 streamlog_out(DEBUG)  << "best start mass of Z: " << beststartmassZ << std::endl ;
 streamlog_out(DEBUG)  << "best start mass of H: " << beststartmassH << std::endl ;
 streamlog_out(DEBUG)  << "best mass of Z: " << bestmassZ << std::endl ;
 streamlog_out(DEBUG)  << "best mass of H: " << bestmassH << std::endl ;
 streamlog_out(DEBUG)  << "Error Code: " << errorcode << std::endl ;
 
 
#ifdef MARLIN_USE_AIDA
 hFitErrorBest->fill( besterr ) ;
 if (bestprob > 0) {
   hFitProbBest->fill( bestprob ) ;
   hNItBest->fill( bestnit ) ;
   hRecHMassBest->fill( bestmassH ) ;
   hRecHMassNoFitBest->fill( beststartmassH ) ;
   hRecZMassNoFitBest->fill( beststartmassZ ) ;
   hPhotonEnergy->fill( _fitISR ? bestphotonenergy : 0. );
 }
 // if none of the permutations converged, fill here startmass for permutation with smallest initial chi2!
 else {
   hRecHMassNoFitFail->fill( chi2startmassH ) ;
   hRecZMassNoFitFail->fill( chi2startmassZ ) ;
 }
#endif
 //streamlog_out(MESSAGE)  << "final mass of Z: " << z.getMass(1) << std::endl ;
 
 delete j1;
 delete j2;
 delete l1;
 delete l2;
 
 ZHTree->Fill();
 evt->addCollection( OutputCol, _OutCol.c_str() );
 
 evt->addCollection( PullEOutCol, _OutPullECol.c_str() );
 evt->addCollection( PullThetaOutCol, _OutPullThetaCol.c_str() );
 evt->addCollection( PullPhiOutCol, _OutPullPhiCol.c_str() );
 //}// end if # SLD-B =1
 // else{
    //   streamlog_out(ERROR) << "BREAK*********" << std::endl;
    // }
     }//end if jetcol is not null


	streamlog_out(DEBUG4) << "next event******************" << std::endl;

  _nEvt ++ ;
}//event end



void ZHllqq5CFit::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZHllqq5CFit::end(){

    	streamlog_out(ERROR) << "# of events: " << _nEvt << std::endl;
    	streamlog_out(ERROR) << "# of nucorrection: " << correction<< std::endl;
    	streamlog_out(ERROR) << "# of Covariance failed: " << nCo<< std::endl;
	streamlog_out(ERROR) << "# of tack != 1: " << trackcounter<< std::endl;

  _fout->Write(0);
  _fout->Close();


}

