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


ZHllqq5CFit::ZHllqq5CFit() : Processor("ZHllqq5CFit"),
			     m_Bfield(0.f),
			     c(0.),
			     mm2m(0.),
			     eV2GeV(0.),
			     eB(0.)
{

  // modify processor description
  _description = "ZHllqq5CFit does a 5C fit on ZH->mumubb events" ;


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
                                _sigmaScaleFactor,
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

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                                "FitPostRecoColection",
                      	        " Output Fit Colection" ,
                      	         _PostFitRecoCol,
                                std::string("FitReco")) ;

 registerOutputCollection( LCIO::LCFLOATVEC,
                                "PullsEPostFitRecoColection",
                                " Output Pull E  Colection" ,
                                 _OutPullECol,
                                std::string("PullEnergy")) ;

  registerOutputCollection( LCIO::LCFLOATVEC,
                                "PullsThetaPostFitRecoColection",
                                " Output Pull Theta  Colection" ,
                                 _OutPullThetaCol,
                                std::string("PullTheta")) ;

 registerOutputCollection( LCIO::LCFLOATVEC,
                                "PullsPhiPostFitRecoColection",
                                " Output Pull Phi  Colection" ,
                                 _OutPullPhiCol,
                                std::string("PullPhi")) ;

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

}

void ZHllqq5CFit::processRunHeader( LCRunHeader* ) {
  _nRun++ ;
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
  static AIDA::IHistogram1D* hRecHMassFitOK ;
  static AIDA::IHistogram1D* hRecZMassFitOK ;
  static AIDA::IHistogram1D* hRecHMassNoFitAll ;
  static AIDA::IHistogram1D* hRecZMassNoFitAll ;
  static AIDA::IHistogram1D* hRecHMassNoFitFail ;
  static AIDA::IHistogram1D* hRecZMassNoFitFail ;
  static AIDA::IHistogram1D* hFitProb ;
  static AIDA::IHistogram1D* hChi2 ;
  static AIDA::IHistogram1D* hNIt ;
  static AIDA::IHistogram1D* hPhotonEnergy ;
  static AIDA::IHistogram1D* hJetMass ;
  static AIDA::IHistogram1D* hFitError;
  static AIDA::IHistogram1D* hPullEJet;
  static AIDA::IHistogram1D* hPullThJet;
  static AIDA::IHistogram1D* hPullPhJet;
  static AIDA::IHistogram1D* hPullELepton;
  static AIDA::IHistogram1D* hPullThLepton;
  static AIDA::IHistogram1D* hPullPhLepton;


  if( isFirstEvent() ) {
    
    hRecHMassFitOK =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassFitOK", "M_H", 200, 0., 250. ) ;
    hRecZMassFitOK =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecZMassFitOK", "M_Z", 200, 0., 250. ) ;
    hRecHMassNoFitAll =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassNoFitAll", "M_H", 200, 0., 250. ) ;
    hRecZMassNoFitAll =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecZMassNoFitAll", "M_Z", 200, 0., 250. ) ;
    hRecHMassNoFitFail =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecHMassNoFitFail", "M_H", 200, 0., 250. ) ;
    hRecZMassNoFitFail =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hRecZMassNoFitFail", "M_Z", 200, 0., 250. ) ;
    hFitProb =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ;
    hChi2 =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hChi2", "chi2", 200, 0, 200. ) ;
    hNIt =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hNIt", "number of iterations", 200, 0, 200. ) ;
    hPhotonEnergy =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPhotonEnergy", "ISR photon energy", 100, 0., 50. ) ;
    hJetMass =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hJetMass", "Jet Mass", 200, 0., 100. ) ;
    if (_ifitter == 1) {
      hFitError =
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 100, -1.5, 99.5 ) ;
    }
    else {
      hFitError =
        AIDAProcessor::histogramFactory(this)->
        createHistogram1D( "hFitError", "Error flag", 11, -1.5, 9.5 ) ;
    }
    hPullEJet =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullEJet", "pull of jet energy", 100, -5., 5. ) ;
    hPullThJet =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThJet", "pull of jet theta", 100, -5., 5. ) ;
    hPullPhJet =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhJet", "pull of jet phi", 100, -5., 5. ) ;
    hPullELepton =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullELepton", "pull of lepton energy", 100, -5., 5. ) ;
    hPullThLepton =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullThLepton", "pull of lepton theta", 100, -5., 5. ) ;
    hPullPhLepton =
      AIDAProcessor::histogramFactory(this)->
      createHistogram1D( "hPullPhLepton", "pull of lepton phi", 100, -5., 5. ) ;
   
  }

#endif




//////////////////  get lepton and jet four-momenta ///////////////////////////

  LCCollection* jetcol = evt->getCollection( _jetcolName ) ;
//   float yminus = jetcol ->parameters().getFloatVal( "YMinus");
//   streamlog_out(MESSAGE)  << " yminus = " << yminus << std::endl ;
//   float yplus = jetcol ->parameters().getFloatVal( "YPlus");
//   streamlog_out(MESSAGE)  << " yplus = " << yplus << std::endl ;
  
  LCCollection* lepcol = evt->getCollection( _colLeptons ) ;
  
  if (!jetcol  || !lepcol) return;
  
  int nJETS = jetcol->getNumberOfElements()  ;
  int nLEPS = lepcol->getNumberOfElements()  ;
  streamlog_out(MESSAGE)
                 << " found " << nJETS << "jets and " << nLEPS

                 << " leptons in event "          << evt->getEventNumber()
                 << std::endl ;

  if (nJETS != 2) return;
  if (nLEPS != 2) return;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////Set JetFitObjects
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// original fit objects
  JetFitObject* j1 = 0;
  JetFitObject* j2 = 0;
  HepLorentzVector jetvec;
  
  double JetResE, JetResTheta, JetResPhi;

  for (int i_jet = 0; i_jet < nJETS; i_jet++) {
    ReconstructedParticle *j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i_jet ) );
    if (j) {
      streamlog_out(MESSAGE)
	 << " found jet in event " << evt->getEventNumber()
	 << "  in run "          << evt->getRunNumber()
	 << std::endl ;
         
      jetvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy());      
      
      // use simple parametrisation or jet covariances?
      if (_errene > 0) {
        JetResE    = _errene;
        JetResTheta=_errtheta;
        JetResPhi=_errphi;
      }
      else {
        // vect gets 3 vector, mag2 length squared, perp the transverse component of spacial vector
        double dth_dpx    = ( jetvec.px() * jetvec.pz() ) / ( jetvec.vect().mag2() * jetvec.perp() );
        double dth_dpy    = ( jetvec.py() * jetvec.pz() ) / ( jetvec.vect().mag2() * jetvec.perp() );
        double dth_dpz    = -jetvec.perp() / jetvec.vect().mag2();
        double dphi_dpx   = -jetvec.py() / jetvec.perp2();
        double dphi_dpy   =  jetvec.px()/ jetvec.perp2();
        
        double SigPx2     = j->getCovMatrix()[ 0 ];
        double SigPxSigPy = j->getCovMatrix()[ 1 ];
        double SigPxSigPz = j->getCovMatrix()[ 3 ];
        double SigPy2     = j->getCovMatrix()[ 2 ];
        double SigPySigPz = j->getCovMatrix()[ 4 ];
        double SigPz2     = j->getCovMatrix()[ 5 ];
        double SigE2      = j->getCovMatrix()[ 9 ];
        streamlog_out(DEBUG4)  << "SigPx2= " << SigPx2 <<std::endl;
        streamlog_out(DEBUG4)  << "SigPxSigPy= " << SigPxSigPy <<std::endl;
        streamlog_out(DEBUG4)  << "SigPxSigPz= " << SigPxSigPz <<std::endl;
        streamlog_out(DEBUG4)  << "SigPy2= " <<  SigPy2<<std::endl;
        streamlog_out(DEBUG4)  << "SigPySigPz= " << SigPySigPz <<std::endl;
        streamlog_out(DEBUG4)  << "SigPz2= " << SigPz2 <<std::endl;
        streamlog_out(DEBUG4)  << "SigE2= " << SigE2 <<std::endl;
      
        JetResE    = std::sqrt( SigE2 );
        if(SigPx2==0. || SigPxSigPy==0. || SigPxSigPz==0. || SigPy2==0. || SigPySigPz==0. || SigPz2==0. || SigE2==0.){
          streamlog_out(WARNING) << "Covariance Matrix is singular"<<std::endl;
          streamlog_out(WARNING) << "Setting theta and phi Resolution to default values "<<std::endl;
          JetResTheta=_errtheta;
          JetResPhi=_errphi;
          nCo++;
        }
        else {
          JetResTheta= std::sqrt( std::fabs( SigPx2 * std::pow( dth_dpx , 2 ) + SigPy2 * std::pow( dth_dpy , 2 ) + SigPz2 * std::pow( dth_dpz , 2 ) + 2 * ( SigPxSigPy * dth_dpx * dth_dpy ) + 2 * ( SigPySigPz * dth_dpy * dth_dpz ) + 2 * ( SigPxSigPz * dth_dpx * dth_dpz ) ) );
          JetResPhi  = std::sqrt( std::fabs( SigPx2 * std::pow( dphi_dpx , 2 ) + SigPy2 * std::pow( dphi_dpy , 2 ) + 2 * ( SigPxSigPy * dphi_dpx * dphi_dpy ) ) );
        }
      
      } 
      streamlog_out(DEBUG)  << "_errene = " << _errene <<std::endl;
      streamlog_out(DEBUG)  << "JetResE= " << JetResE <<std::endl;
      streamlog_out(DEBUG)  << "JetResTheta= " << JetResTheta <<std::endl;
      streamlog_out(DEBUG)  << "JetResPhi= " << JetResPhi <<std::endl;
      
      // ...........................................................................................

      if (i_jet == 0 ) {
        j1 = new JetFitObject (jetvec.e(), jetvec.theta(), jetvec.phi(),
                               JetResE*_sigmaScaleFactor, JetResTheta, JetResPhi, jetvec.m());
        j1->setName("Jet1");
        streamlog_out(DEBUG)  << " start four-vector of first  jet: " << *j1  << std::endl ;
      }
      else if (i_jet == 1) {
        j2 = new JetFitObject (jetvec.e(), jetvec.theta(), jetvec.phi(),
                               JetResE*_sigmaScaleFactor, JetResTheta, JetResPhi, jetvec.m());
        j2->setName("Jet2");
        streamlog_out(DEBUG) << " start four-vector of second  jet: " << *j2  << std::endl ;
      }
#ifdef MARLIN_USE_AIDA
      hJetMass->fill(j->getMass());
#endif
    }//end of   if (reco particle j!=0)
     
  }//end loop over nJets


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////Set LeptonFitObjects
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//orginal lepton fit objects
  LeptonFitObject* l1 = 0;
  LeptonFitObject* l2 = 0;
  HepLorentzVector leptonvec;
  
  for(int i_lep = 0; i_lep < nLEPS ; i_lep++) {
    ReconstructedParticle* l = dynamic_cast<ReconstructedParticle*>( lepcol->getElementAt( i_lep ) ) ;
    if (l) {
      streamlog_out(MESSAGE)
	 << " found lepton in event " << evt->getEventNumber()
	 << "  in run "          << evt->getRunNumber()
	 << std::endl ;
         
      leptonvec = HepLorentzVector ((l->getMomentum())[0],(l->getMomentum())[1],(l->getMomentum())[2],l->getEnergy());
    
      // ...........................................................................................
    
      TrackVec tckvec = l->getTracks();
      double Omega, Omega_uncert = 0.;
      double TanLambda, TanLambda_err = 0.;
      double theta, theta_err = 0.;
      double phi, phi_err = 0.;
      double invers_pT, invers_pT_err = 0.;
      if ( tckvec.size() != 1 ){
        streamlog_out(WARNING)  << "Number of tracks for lepton[" << i_lep <<"] is not exactly ONE!!! (nTracks = " << tckvec.size() << " ) " << std::endl ;
        // initialize from PFO 4-momentum and fixed resolutions
        invers_pT = 1/leptonvec.perp();
        // sigma (1/pt) = 2 *10^5 GeV^-1 \oplus 10^-3 /(pt*sintheta)
        invers_pT_err = std::sqrt( std::pow(0.00002, 2) + std::pow ( 0.001 /(leptonvec.perp2()/leptonvec.vect().mag()), 2));
        theta = leptonvec.theta();
        theta_err = 0.1;
        phi = leptonvec.phi();
        phi_err = 0.1;
      }
      else {
        // initialize from track & covariance matrix
        streamlog_out(DEBUG)  << "Number of tracks for lepton[" << i_lep <<"] is exactly ONE!!!" << std::endl;
        Omega = tckvec[0]->getOmega();
        Omega_uncert = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[5]) );
        streamlog_out(DEBUG)  << "Omega = " << Omega << std::endl;
        streamlog_out(DEBUG)  << "Omega_uncert = " << Omega_uncert << std::endl;
        invers_pT = std::fabs( Omega / eB );
        invers_pT_err = std::fabs( 1. / eB ) * Omega_uncert;
        streamlog_out(DEBUG)  << "invers_pT = " << invers_pT << std::endl;
        streamlog_out(DEBUG)  << "invers_pT_err = " << invers_pT_err << std::endl;
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
      }
 
 
      streamlog_out(MESSAGE)  << "Lepton fit object from leptonvec: "
                            << 1/leptonvec.perp() <<" +- " << std::sqrt( std::pow(0.00002, 2) + std::pow ( 0.001 /(leptonvec.perp2()/leptonvec.vect().mag()), 2)) << " , "
                            << leptonvec.theta() <<" +- " << 0.1 << " , "
                            << leptonvec.phi() <<" +- " << 0.1 << std::endl ;

      streamlog_out(MESSAGE)  << "Lepton fit object from track:     "
                            << invers_pT <<" +- " << invers_pT_err << " , "
                            << theta <<" +- " << theta_err << " , "
                            << phi <<" +- " << phi_err << std::endl ;

      // ...........................................................................................
   
      if (i_lep == 0 ) {
        l1 = new LeptonFitObject (invers_pT , theta , phi , invers_pT_err , theta_err , phi_err, leptonvec.m());
        l1->setName("Lepton1");
        streamlog_out(DEBUG)  << " start four-vector of first  lepton: " << *l1  << std::endl ;
      }
      else if (i_lep == 1) {
        l2 = new LeptonFitObject (invers_pT , theta , phi , invers_pT_err , theta_err , phi_err, leptonvec.m());
        l2->setName("Lepton2");
        streamlog_out(DEBUG)  << " start four-vector of second  lepton: " << *l2  << std::endl ;
      }
    } //end of   if (reco particle l!=0)
  } //end loop over nLeptons


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////  Assure book-keeping of fit objects: keep copy of start values and additional copy for case of permutations 
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


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

    // these point always to the fitjets array, which gets reset in case we want to check different permutations
  JetFitObject *jets[NJETS];
  for (int i = 0; i < NJETS; ++i) {
    jets[i] = &fitjets[i];
    streamlog_out(MESSAGE)  << "start four-vector of jet " << i << ": " << *(jets[i])  << std::endl ;
  }

  LeptonFitObject *leptons[NLEPTONS];
  for (int i = 0; i < NLEPTONS; ++i) {
    leptons[i] = &fitleptons[i];
    streamlog_out(MESSAGE)  << "start four-vector of leptons " << i << ": " << *(leptons[i])  << std::endl ;
  }
  // ......................................................................................
 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////set up constraints 
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

  float target_p_due_crossing_angle = _ecm * 0.007; // crossing angle = 14 mrad
  MomentumConstraint pxc ( 0 , 1 , 0 , 0 , target_p_due_crossing_angle ); //Factor for: (energy sum, px sum, py sum,pz sum,target value of sum)
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

  double E_lab= 2 * sqrt( std::pow( 0.548579909e-3 , 2 ) + std::pow( _ecm / 2 , 2 ) + std::pow( target_p_due_crossing_angle , 2 ) + 0. + 0.);
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

  // ISR Photon initialized with missing p_z
  ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);

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
 
  // important: use soft constraints for Z->leptons since natural width is NON-negligible wrt momentum resolution!
  //MassConstraint z(91.2);
  SoftGaussMassConstraint z(91.2,2.4952/2);
  z.addToFOList (*(leptons[0]), 1);
  z.addToFOList (*(leptons[1]), 1);
 
  // only used for convenient invariant mass computation, will NOT be added to the fit
  MassConstraint h(125.);
  h.addToFOList (*(jets[0]), 1);
  h.addToFOList (*(jets[1]), 1);
 
  startmassZ = z.getMass(1);
  startmassH = h.getMass(1);
  streamlog_out(MESSAGE) << "start mass of Z: " << startmassZ << std::endl ;
  streamlog_out(MESSAGE) << "start mass of H: " << startmassH << std::endl ;
  streamlog_out(MESSAGE) << "sum of massconstraints  before fit = " << z.getValue() + h.getValue() << std::endl ; 
 
#ifdef MARLIN_USE_AIDA
  hRecHMassNoFitAll->fill( startmassH ) ;
  hRecZMassNoFitAll->fill( startmassZ ) ;
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
//// initialize fitter
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
//// add fit objects to fitter
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  // don't constrain Higgs mass, just use constraints for convenient mass calculation
  //fitter.addConstraint (h);
 
  streamlog_out(DEBUG4) << "constraints added"  << std::endl ;
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
//// do it !
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

  prob = fitter.fit();
  double chi2 = fitter.getChi2();
  nit = fitter.getIterations();
  int ierr = fitter.getError();
  
  streamlog_out(MESSAGE)  << "fitter error: " << ierr << std::endl;
  hFitError->fill( ierr ) ;
 
 
  streamlog_out(MESSAGE) << "fit probability = " << prob << std::endl ;
  streamlog_out(MESSAGE) << "fit chi2 = " << chi2  << std::endl ;
  streamlog_out(MESSAGE)  << "final mass of Z: " << z.getMass(1) << std::endl ;
  streamlog_out(MESSAGE)  << "final mass of H: " << h.getMass(1) << std::endl ;
 
  for (int i = 0; i < NJETS; ++i) {
    streamlog_out(MESSAGE)  << "final four-vector of jet " << i << ": " << *(jets[i]) << std::endl ;
    streamlog_out(MESSAGE)  << "final px of jet " << i << ": " << (jets[i])[0] << std::endl ;
  }
  for (int i = 0; i < NLEPTONS; ++i) {
    streamlog_out(MESSAGE)  << "final four-vector of lepton " << i << ": " << *(leptons[i]) << std::endl ;
    streamlog_out(DEBUG)  << "final px of lepton " << i << ": " << (leptons[i])[0] << std::endl ;
  }

  if(_fitISR){
    streamlog_out(MESSAGE)  << "final four-vector of ISR photon: " << *(photon) << std::endl ;
  }
 
  streamlog_out(DEBUG)  << "Value of pxc after fit: " << pxc.getValue() << std::endl ;
  streamlog_out(DEBUG)  << "Value of pyc after fit: " << pyc.getValue() << std::endl ;
  streamlog_out(DEBUG)  << "Value of pzc after fit: " << pzc.getValue() << std::endl ;
  streamlog_out(DEBUG)  << "Value of ec after fit: " << ec.getValue() << std::endl ;
  

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
//// Calculate pulls
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double pullJet[3][2];
  double pullLepton[3][2];
  
  // ierr == -1 only means that error calculation for fitted parameters failed!
  if (ierr <= 0) { //if fitter.getError() <=0
    if ( prob > 0 ) {
#ifdef MARLIN_USE_AIDA
    hFitProb->fill( prob ) ;
    hChi2->fill( chi2 ) ;
    hNIt->fill( nit ) ;
    hRecHMassFitOK->fill( h.getMass(1) ) ;
    hRecZMassFitOK->fill( z.getMass(1) ) ;
    hPhotonEnergy->fill( _fitISR ? photon->getE() : 0. );
#endif
    }
    else {
      streamlog_out(WARNING)  << " ierr <= 0 but prob == 0, filling fail histograms; chi2 = " << chi2 << std::endl ;
      hRecHMassNoFitFail->fill( startmassH ) ;
      hRecZMassNoFitFail->fill( startmassZ ) ;
    }
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
        hPullEJet->fill(pullJet[0][ifo]);
        hPullThJet->fill(pullJet[1][ifo]);
        hPullPhJet->fill(pullJet[2][ifo]);
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
        hPullELepton->fill(pullLepton[0][ifo]);
        hPullThLepton->fill(pullLepton[1][ifo]);
        hPullPhLepton->fill(pullLepton[2][ifo]);
      }
    } // if ierr==0
    else {
      streamlog_out(WARNING) << " FIT CONVERGED, BUT ERROR CALCULATION FAILED "
			     << " in event " << evt->getEventNumber()
			     << ", not filling pull histograms!"  << std::endl ;
    }
  } // if ierr <=0
  else {
    streamlog_out(ERROR) << "FIT FAILED WITH ERROR = " << fitter.getError()
			 << " in event " << evt->getEventNumber()<< std::endl ;
    streamlog_out(MESSAGE)  << "start mass of Z: " << startmassZ << std::endl ;
    streamlog_out(MESSAGE)  << "start mass of H: " << startmassH << std::endl ;
    streamlog_out(MESSAGE)  << "final mass of Z: " << z.getMass(1) << std::endl ;
    streamlog_out(MESSAGE)  << "final mass of H: " << h.getMass(1) << std::endl ;
    
    hRecHMassNoFitFail->fill( startmassH ) ;
    hRecZMassNoFitFail->fill( startmassZ ) ;
  }

 

  // ******* write collection with fitted particles 
  // TODO: add covariance matrices, component particles etc!

  LCCollectionVec *PostFitRecoCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  ReconstructedParticleImpl* ISRfitRecoPart = new ReconstructedParticleImpl;
  ReconstructedParticleImpl* ZfitRecoPart = new ReconstructedParticleImpl;
  ReconstructedParticleImpl* HfitRecoPart = new ReconstructedParticleImpl;
  ReconstructedParticleImpl* LfitRecoPart[NLEPTONS];
  LfitRecoPart[0] = new ReconstructedParticleImpl; 
  LfitRecoPart[1] = new ReconstructedParticleImpl;
  ReconstructedParticleImpl* JfitRecoPart[NJETS];
  JfitRecoPart[0] = new ReconstructedParticleImpl;
  JfitRecoPart[1] = new ReconstructedParticleImpl;
  
  float momentum[3];
  
  // ** ISR Photon
  momentum[0]= photon->getPx();
  momentum[1]= photon->getPy();
  momentum[2]= photon->getPz();
  ISRfitRecoPart->setMomentum(momentum);
  ISRfitRecoPart->setEnergy(photon->getE());
  ISRfitRecoPart->setType (22);
  ISRfitRecoPart->setCharge (0);
  ISRfitRecoPart->setMass (0);
  PostFitRecoCol->addElement(ISRfitRecoPart);

  // ** Z boson
  // leptons[i] are the fitted FitObjects
  momentum[0]= leptons[0]->getPx()+leptons[1]->getPx();
  momentum[1]= leptons[0]->getPy()+leptons[1]->getPy();
  momentum[2]= leptons[0]->getPz()+leptons[1]->getPz();
  ZfitRecoPart->setMomentum(momentum);
  ZfitRecoPart->setEnergy(leptons[0]->getE()+leptons[1]->getE());
  ZfitRecoPart->setType (23);
  ZfitRecoPart->setCharge(0);
  ZfitRecoPart->setMass(z.getMass(1));
  PostFitRecoCol->addElement(ZfitRecoPart);
 
  // ** Higgs boson
  // jets[i] are the fitted FitObjects
  momentum[0]= jets[0]->getPx()+jets[1]->getPx();
  momentum[1]= jets[0]->getPy()+jets[1]->getPy();
  momentum[2]= jets[0]->getPz()+jets[1]->getPz();
  HfitRecoPart->setMomentum(momentum);
  HfitRecoPart->setEnergy(jets[0]->getE()+jets[1]->getE());
  HfitRecoPart->setType (25);
  HfitRecoPart->setCharge(0);
  HfitRecoPart->setMass(h.getMass(1));
  PostFitRecoCol->addElement(HfitRecoPart);
 
   // ** leptons
  // leptons[i] are the fitted FitObjects, l is the original PFO
  for(int i_lep = 0; i_lep < nLEPS ; i_lep++) {
    ReconstructedParticle* l = dynamic_cast<ReconstructedParticle*>( lepcol->getElementAt( i_lep ) ) ;
    momentum[0]= leptons[i_lep]->getPx();
    momentum[1]= leptons[i_lep]->getPy();
    momentum[2]= leptons[i_lep]->getPz();
    LfitRecoPart[i_lep]->setMomentum(momentum);
    LfitRecoPart[i_lep]->setEnergy(leptons[i_lep]->getE());
    LfitRecoPart[i_lep]->setType(l->getType());
    LfitRecoPart[i_lep]->setCharge(l->getCharge());
    LfitRecoPart[i_lep]->setMass(l->getMass());
    PostFitRecoCol->addElement(LfitRecoPart[i_lep]);
  }

   // ** jets
  // jets[i] are the fitted FitObjects, j is the original PFO
  for (int i_jet = 0; i_jet < nJETS; i_jet++) {
    ReconstructedParticle *j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i_jet ) );
    momentum[0]= jets[i_jet]->getPx();
    momentum[1]= jets[i_jet]->getPy();
    momentum[2]= jets[i_jet]->getPz();
    JfitRecoPart[i_jet]->setMomentum(momentum);
    JfitRecoPart[i_jet]->setEnergy(jets[i_jet]->getE());
    JfitRecoPart[i_jet]->setType (5);
    JfitRecoPart[i_jet]->setCharge(j->getCharge());
    JfitRecoPart[i_jet]->setMass(j->getMass());
    PostFitRecoCol->addElement(JfitRecoPart[i_jet]);
  }

  
  PostFitRecoCol->parameters().setValue("fitprob", (float)prob);
  streamlog_out(DEBUG) << " fitprob:   " << chi2 << std::endl ;
  PostFitRecoCol->parameters().setValue("chisq", (float)chi2);
  streamlog_out(DEBUG) << " chi2:   " << chi2 << std::endl ;
  PostFitRecoCol->parameters().setValue("error_code", (float)ierr);
  streamlog_out(DEBUG) << " error_code:   " << ierr << std::endl ;
 
 
  evt->addCollection( PostFitRecoCol, _PostFitRecoCol.c_str() );
  
  
  // ******* write collections with pulls  --- this needs re-thinking!!!

  LCCollectionVec *PullEOutCol = new LCCollectionVec(LCIO::LCFLOATVEC);
  LCCollectionVec *PullThetaOutCol = new LCCollectionVec(LCIO::LCFLOATVEC);
  LCCollectionVec *PullPhiOutCol = new LCCollectionVec(LCIO::LCFLOATVEC);
 
  EVENT::LCFloatVec*  PullE = new LCFloatVec();
  EVENT::LCFloatVec*  PullTheta = new LCFloatVec();
  EVENT::LCFloatVec*  PullPhi = new LCFloatVec();
  
  PullE->push_back(pullJet[0][0]);
  PullTheta->push_back(pullJet[1][0]);
  PullPhi->push_back(pullJet[2][0]);

  PullE->push_back(pullJet[0][1]);
  PullTheta->push_back(pullJet[1][1]);
  PullPhi->push_back(pullJet[2][1]);

  PullE->push_back(pullLepton[0][0]);
  PullTheta->push_back(pullLepton[1][0]);
  PullPhi->push_back(pullLepton[2][0]);

  PullE->push_back(pullLepton[0][1]);
  PullTheta->push_back(pullLepton[1][1]);
  PullPhi->push_back(pullLepton[2][1]);  
  
  PullEOutCol->addElement(PullE);
  PullThetaOutCol->addElement(PullTheta);
  PullPhiOutCol->addElement(PullPhi);
  evt->addCollection( PullEOutCol, _OutPullECol.c_str() );
  evt->addCollection( PullThetaOutCol, _OutPullThetaCol.c_str() );
  evt->addCollection( PullPhiOutCol, _OutPullPhiCol.c_str() );


  streamlog_out(DEBUG4) << "next event******************" << std::endl;

  delete j1;
  delete j2;
  delete l1;
  delete l2;
  delete photon;
   
  _nEvt ++ ;
}  



void ZHllqq5CFit::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZHllqq5CFit::end(){

  streamlog_out(MESSAGE) << "# of events: " << _nEvt << std::endl;
  streamlog_out(MESSAGE) << "# of nucorrection: " << correction<< std::endl;
  streamlog_out(MESSAGE) << "# of Covariance failed: " << nCo<< std::endl;
  streamlog_out(MESSAGE) << "# of track != 1: " << trackcounter<< std::endl;

}

