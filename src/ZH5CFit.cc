#include "ZH5CFit.h"
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

#include "TLorentzVector.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
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


ZH5CFit aZH5CFit ;

// function to define the jet energy resolution (in GeV)
double ZH5CFit::JetEnergyResolution(double E)
{

  // examples here derived by Benjamin Hermberg from e+e- -> udsc:
  // 1) default is 120%/sqrt(E), gives best convergence of 5C fit on e+e- -> udsc
  double result = _errene*std::sqrt(E);

  // 2) comparing jet-level to quark-level energies
  //    (using MarlinReco/Analysis/RecoMCTruthLink/QuarkJetPairing.cc)
  if (_errene == 0 ) result = std::sqrt(pow(0.6908,2)*(E)+(pow(0.02596,2)*pow(E,2)));

  return result;
}

ZH5CFit::ZH5CFit() : Processor("ZH5CFit") {

  // modify processor description
  _description = "ZH5CFit does a 5C fit on 4 jet events (Px, Py, Pz, E, M12 = MZ (for all six permutations))" ;


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

  registerProcessorParameter( "NuEnergyCorrection" ,
                              "0: Fit hypothesis without ISR   1: Fit hypothesis including ISR",
                              _NuE,
                              (int) 1);
  registerProcessorParameter( "useErrorFlow",
                              "If true, use covariance matrix for energy uncertainty. Otherwise 1.2/sqrt(E)",
                              _useErrorFlow,
                              (int) 1 );
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

  // registerInputCollection( 	LCIO::RECONSTRUCTEDPARTICLE,
  //         "ErrorFlowCollection",
  //         "Collection of jet with error flow information",
  //         _errorflowcollection,
  //         std::string("OutputErrorFlowJets")
  //         );
  registerInputCollection( LCIO::MCPARTICLE,
                                  "SemiLeptonicDecays" , //name
                                  "Semi-Leptonic Decays Collection"  , //description
                                  _SLDCol , //my collection name
                                  std::string("") //take this string from input slcio
                          );
  registerInputCollection( 	LCIO::MCPARTICLE,
                                "NeutrinoCorrection",
                                "Collection of Corrected/Estimated neutrino energies from SemiLeptonic decays",
                                _NuCorrector,
                                std::string("NuCorrect")
                                );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                                "FitOutputColection",
                      	        " Output Fit Colection" ,
                      	         _OutCol,
                                std::string("FitReco")) ;
  registerProcessorParameter("outputFilename",
                               "name of output file",
                               _outfile,
                                std::string("")
                                );
}


void ZH5CFit::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  b = (double) 0.00464564*( std::log(_ecm*_ecm*3814714.)-1. );
  //= 2*alpha/pi*( ln(s/m_e^2)-1 )
  ISRPzMaxB = std::pow((double)_isrpzmax,b);

  _fout = new TFile(_outfile.c_str(),"recreate");

   ZHTree = new TTree("ZHTree","ZHTree");
   ZHTree->Branch("Hmass_before_fit_best",&Hmass_before_fit,"Hmass_before_fit/F") ;
   ZHTree->Branch("Hmass_No_fit_best",&Hmass_NoFit,"Hmass_NoFit/F") ;
   ZHTree->Branch("Hmass_after_fit_best",&Hmass_after_fit,"Hmass_after_fit/F") ;
   ZHTree->Branch("Error_Code",&Error_code,"Error_code/F") ;
   ZHTree->Branch("hpull_jet1_E_best",&hpull_jet1_E,"hpull_jet1_E/F") ;
   ZHTree->Branch("hpull_jet2_E_best",&hpull_jet2_E,"hpull_jet2_E/F") ;
   ZHTree->Branch("hpull_jet3_E_best",&hpull_jet3_E,"hpull_jet3_E/F") ;
   ZHTree->Branch("hpull_jet4_E_best",&hpull_jet4_E,"hpull_jet4_E/F") ;
   ZHTree->Branch("hpull_jet1_th_best",&hpull_jet1_th,"hpull_jet1_th/F") ;
   ZHTree->Branch("hpull_jet2_th_best",&hpull_jet2_th,"hpull_jet2_th/F") ;
   ZHTree->Branch("hpull_jet3_th_best",&hpull_jet3_th,"hpull_jet3_th/F") ;
   ZHTree->Branch("hpull_jet4_th_best",&hpull_jet4_th,"hpull_jet4_th/F") ;
   ZHTree->Branch("hpull_jet1_phi_best",&hpull_jet1_phi,"hpull_jet1_phi/F") ;
   ZHTree->Branch("hpull_jet2_phi_best",&hpull_jet2_phi,"hpull_jet2_phi/F") ;
   ZHTree->Branch("hpull_jet3_phi_best",&hpull_jet3_phi,"hpull_jet3_phi/F") ;
   ZHTree->Branch("hpull_jet4_phi_best",&hpull_jet4_phi,"hpull_jet4_phi/F") ;
   ZHTree->Branch("jetmatch",&jetmatch,"jetmatch/I") ;
   ZHTree->Branch("jetmatchth",&jetmatchth,"jetmatchth/I") ;
   ZHTree->Branch("jetmatchphi",&jetmatchphi,"jetmatchphi/I") ;
}

void ZH5CFit::processRunHeader( LCRunHeader* ) {
  _nRun++ ;
}
void ZH5CFit::compcorrect() //finds the jet with cross checking
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

void ZH5CFit::SetZero()
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


void ZH5CFit::processEvent( LCEvent * evt ) { //event start


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
  static AIDA::IHistogram1D* hPullEJetBest;
  static AIDA::IHistogram1D* hPullThJetBest;
  static AIDA::IHistogram1D* hPullPhJetBest;


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





  HepLorentzVector lvec;
  HepLorentzVector jetvec;


  // fill histogram from LCIO data :

  //////////////////   JETS ///////////////////////////

     LCCollection* jetcol = evt->getCollection( _jetcolName ) ;
     if (jetcol != 0) {//jetcol is not null

       int nJETS = jetcol->getNumberOfElements()  ;
       streamlog_out(MESSAGE)
                      << " found " << nJETS
                      << " jets in event " << evt->getEventNumber()
                      << "  in run "          << evt->getRunNumber()
                      << std::endl ;

       if (nJETS != 4) return;

       float yminus = jetcol ->parameters().getFloatVal( "YMinus");
       streamlog_out(MESSAGE)  << " yminus = " << yminus << std::endl ;

       float yplus = jetcol ->parameters().getFloatVal( "YPlus");
       streamlog_out(MESSAGE)  << " yplus = " << yplus << std::endl ;



   //---------------------------neutrino collection------------------------
       LCCollection* sldcol = NULL;
       try{
       sldcol = evt->getCollection( _SLDCol );
       streamlog_out(MESSAGE) << _SLDCol << " collection available*********" << std::endl;
       }
       catch( lcio::DataNotAvailableException e )
       {
       streamlog_out(WARNING) << _SLDCol << " collection not available****" << std::endl;
       sldcol = NULL;
       }
       nSLDB=sldcol->getParameters().getIntVal("nBSLD");
       nSLDC=sldcol->getParameters().getIntVal("nCSLD");
       nSLDBC=sldcol->getParameters().getIntVal("nSLD");
       B_index={};
       B_index=sldcol->getParameters().getIntVals("BHadronIndex", B_index);
       C_index={};
       C_index=sldcol->getParameters().getIntVals("CHadronIndex", C_index);

       streamlog_out(MESSAGE) << "Number of Semi-Leptonic decay of B-Hadron: " << nSLDB << std::endl;
       streamlog_out(MESSAGE) << "Number of Semi-Leptonic decay of C-Hadron: " << nSLDC << std::endl;
       streamlog_out(MESSAGE) << "Total Number of Semi-Leptonic decays: " << nSLDBC << std::endl;


       LCCollection* nucol = NULL;
       try{
       nucol = evt->getCollection( _NuCorrector );
       streamlog_out(MESSAGE) << _NuCorrector << " collection available*********" << std::endl;
       }
       catch( lcio::DataNotAvailableException e )
       {
       streamlog_out(WARNING) << _NuCorrector << " collection not available****" << std::endl;
       nucol = NULL;
       }
       ENuplus=0.;
       ENuminus=0.;
       ENuplus=nucol->getParameters().getFloatVal("recENuPlus");
       ENuminus=nucol->getParameters().getFloatVal("recENuMinus");


       for ( int i=0; i< B_index.size(); i++){
         streamlog_out(DEBUG) << " Index of B[" << i <<"]: " << B_index[i] << std::endl;
       }
       for ( int i=0; i< C_index.size(); i++){
         streamlog_out(DEBUG) << " Index of C[" << i <<"]: " << C_index[i] << std::endl;
       }

       // TLorentzVector l_3p(0,0,0,0);
       HepLorentzVector l_3p;
       LCCollection* mccol = evt->getCollection( _colMCP ) ;
              if(nSLDB == 1 && nSLDC == 0)  {
                streamlog_out(DEBUG)  << "SLD-B is 1      ------******************" << std::endl ;
                 if(B_index.size()!=0){
                   for (int nB=0; nB<B_index.size(); nB++ ){
                     MCParticle* mcpB = dynamic_cast<MCParticle*>( mccol->getElementAt(B_index[nB])) ;
                     MCParticleVec mcpBD = mcpB->getDaughters() ;
                     for(int nBD=0; nBD<mcpBD.size(); nBD++){
                       // streamlog_out(ERROR)  << "PDG OF DAUGHTER:[" << nBD<< "], is " << mcpBD[nBD]->getPDG() <<std::endl ;
                       if((std::abs(mcpBD[nBD]->getPDG()) == 11) || (std::abs(mcpBD[nBD]->getPDG()) == 13) || (std::abs(mcpBD[nBD]->getPDG()) == 15)){
                         l_theta=0. ; l_phi=0.;
                         // l_3p = TLorentzVector(mcpBD[nBD]->getMomentum(), mcpBD[nBD]->getEnergy());
                         // l_p=std::sqrt((pow(l_3p[0],2)+pow(l_3p[1],2)+pow(l_3p[2],2)));
                         // l_theta=acos(l_3p[2]/l_p);
                         // l_phi=atan(l_3p[1]/l_3p[0]);
                         l_3p = HepLorentzVector (mcpBD[nBD]->getMomentum()[0],mcpBD[nBD]->getMomentum()[1],mcpBD[nBD]->getMomentum()[2], mcpBD[nBD]->getEnergy());
                         l_theta=l_3p.theta();
                         l_phi=l_3p.phi();
                         streamlog_out(DEBUG)  << "lepton px= " <<l_3p[0] << ", py= " <<l_3p[1]  << ", pz= " <<l_3p[2] <<", E= " << l_3p[3] <<std::endl ;
                         streamlog_out(DEBUG)  << "lepton theta= " <<l_theta<< ", phi= " <<l_phi <<std::endl ;
                         streamlog_out(DEBUG)  << "=========================================== "<<std::endl ;
                       }
                       if((std::abs(mcpBD[nBD]->getPDG()) == 12) || (std::abs(mcpBD[nBD]->getPDG()) == 14) || (std::abs(mcpBD[nBD]->getPDG()) == 16)){
                         streamlog_out(DEBUG) << " ENuplus:  "<< ENuplus << ", ENuminus: "<< ENuminus << std::endl;
                       }
                     }
                   }
                 }
             else{
               streamlog_out(DEBUG)  << "there is no B meson" << std::endl ;
             }

           }
           else{
             streamlog_out(MESSAGE)  << "# of B-SLD  != 1" << std::endl ;
           }

       float besttheta =100.;
       float bestphi =100.;
       bestjet_th=0;
       bestjet_phi=0;
       bestjet=1000000;
       float diff_besttheta;
       float diff_bestphi;
//-----------------------------------find the jet which should have the neutrino-------------------

      if(nSLDB == 1 && nSLDC == 0){// if # SLD-B =1
        for(int i=0; i< nJETS ; i++){//loop over nJets
              ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i ) ) ;
              if (j) {
                jetvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy());
                 streamlog_out(DEBUG)  << "jet["<<i<<"] px= "<<j->getMomentum()[0] << ", py= " <<j->getMomentum()[1]  << ", pz= " <<j->getMomentum()[2] <<", E= " << j->getEnergy() <<std::endl ;
                 streamlog_out(DEBUG)  << "jet["<<i<<"] Energy= "<<jetvec.e() << ", theta= " <<jetvec.theta()<< ", phi= " <<jetvec.phi()<<std::endl ;
                 delta_theta[i]=jetvec.theta()-l_3p.theta();
                 delta_phi[i]=jetvec.phi()-l_3p.phi();
                 streamlog_out(DEBUG)  << "jet["<<i<<"] difference in theta= "<<delta_theta[i] << ", difference in phi= " <<delta_phi[i]<<std::endl ;
                 if (besttheta > std::abs(delta_theta[i])){
                   besttheta= std::abs(delta_theta[i]);
                   bestjet_th=i;
                   streamlog_out(DEBUG)  << "best theta is "<< besttheta<< " at jet: " << bestjet_th<<std::endl ;
                 }
                 if (bestphi > std::abs(delta_phi[i])){
                   bestphi= std::abs(delta_phi[i]);
                   bestjet_phi=i;
                   streamlog_out(DEBUG)  << "best phi is "<< bestphi<< " at jet: " << bestjet_phi<<std::endl ;
                 }

                 streamlog_out(DEBUG)  << "=========================================== "<<std::endl ;

              }
        }
        if (bestjet_th != bestjet_phi){
           streamlog_out(DEBUG)  << "BESTJET FOR THETA: " <<bestjet_th << " and BESTJET FOR PHI:  " <<bestjet_phi<< " are different " <<std::endl ;
           diff_besttheta= std::abs(delta_theta[bestjet_phi])-std::abs(delta_theta[bestjet_th]);
           diff_bestphi=std::abs(delta_phi[bestjet_phi])-std::abs(delta_phi[bestjet_th]);

               if (abs(delta_theta[bestjet_phi]) < 0.5 && abs(delta_phi[bestjet_phi]) <0.5){  // 1
                     if ( abs(delta_theta[bestjet_th]) > 0.5 && abs(delta_phi[bestjet_th]) >0.5 ){ //2
                       bestjet=bestjet_phi;
                     }
                     else{ //3
                         streamlog_out(DEBUG)  << "all 4 combi  < 0.5" << std::endl ;
                         compcorrect();
                     }
             }
             else if (abs(delta_theta[bestjet_th]) < 0.5 && abs(delta_phi[bestjet_th]) <0.5 ){ //4
                 bestjet=bestjet_th;
             }
             else{//5
               streamlog_out(DEBUG)  << "large uncertainties"<<std::endl ;
               compcorrect();
             }
               streamlog_out(DEBUG)  << "BESTJET FINALLY IS   " <<bestjet <<std::endl;
         }
         else{
             bestjet=bestjet_th;
             streamlog_out(DEBUG)  << "BESTJET FINALLY IS   " <<bestjet <<std::endl;
         }
         jetmatch=bestjet;
         jetmatchth=bestjet_th;
         jetmatchphi=bestjet_phi;

       }//end if  # SLD-B =1




// original fit objects - save for next permutation
     JetFitObject* j1 = 0;
     JetFitObject* j2 = 0;
     JetFitObject* j3 = 0;
     JetFitObject* j4 = 0;
      // if (nSLDB == 1 && nSLDC == 0){
      //    NuEne[0]=ENuplus;
      //    NuEne[1]=ENuminus;
      // }
     SetZero();
 if (_useErrorFlow){
   streamlog_out(DEBUG)  << "Using ErrorFlow...." <<std::endl;
 }
 if (_NuE){
   streamlog_out(DEBUG)  << "Using NeutrinoCorrection...." <<std::endl;
 }
    ReconstructedParticle* jrps[4];
    for(int nuenesign=0; nuenesign<2; nuenesign++){//neutrino correction + and - loop

      streamlog_out(MESSAGE)  << "nuenesign: " << nuenesign << std::endl ;
       for(int i=0; i< nJETS ; i++){//loop over nJets

          ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i ) ) ;
          if (j) {
             jrps[i] = j;
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

             SigPx2=j->getCovMatrix()[ 0 ];
             SigPxSigPy=j->getCovMatrix()[ 1 ];
             SigPxSigPz=j->getCovMatrix()[ 2 ];
             SigPy2=j->getCovMatrix()[ 4 ];
             SigPySigPz=j->getCovMatrix()[ 5 ];
             SigPz2=j->getCovMatrix()[ 7 ];
             SigE2=j->getCovMatrix()[ 9 ];


             dth_dpx =	( Px * Pz ) / ( std::pow( P , 2 ) * std::sqrt( pT2 ) );
             dth_dpy =	( Py * Pz ) / ( std::pow( P , 2 ) * std::sqrt( pT2 ) );
             dth_dpz =	-std::sqrt( pT2 ) / std::pow( P , 2 );
             dphi_dpx=	-Py / pT2;
             dphi_dpy=	Px / pT2 ;

             JetResE	=	std::sqrt( SigE2 ) * sigmaScaleFactor;
             JetResTheta=	std::sqrt( std::fabs( SigPx2 * std::pow( dth_dpx , 2 ) + SigPy2 * std::pow( dth_dpy , 2 ) + SigPz2 * std::pow( dth_dpz , 2 ) + 2 * ( SigPxSigPy * dth_dpx * dth_dpy ) + 2 * ( SigPySigPz * dth_dpy * dth_dpz ) + 2 * ( SigPxSigPz * dth_dpx * dth_dpz ) ) );
             JetResPhi	=	std::sqrt( std::fabs( SigPx2 * std::pow( dphi_dpx , 2 ) + SigPy2 * std::pow( dphi_dpy , 2 ) + 2 * ( SigPxSigPy * dphi_dpx * dphi_dpy ) ) );


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
             lvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy());
             float NuEnergy[2]={0.,0.};
             if (_NuE){
               if (i == bestjet){
                 correction++;
                 NuEnergy[0]=ENuplus;
                 NuEnergy[1]=ENuminus;
                 streamlog_out(DEBUG)  << " added neutrino energy to the best-jet "<< NuEnergy[nuenesign] << std::endl ;
               }
             }

             if (_useErrorFlow){
               if (i == 0 ) {
                 j1 = new JetFitObject (lvec.e()+NuEnergy[nuenesign], lvec.theta(), lvec.phi(),
                    JetResE*sigmaScaleFactor, JetResTheta, JetResPhi, lvec.m());
                 j1->setName("Jet1");
                 streamlog_out(DEBUG)  << " start four-vector of first  jet: " << *j1  << std::endl ;
                }
               else if (i == 1) {
                 j2 = new JetFitObject (lvec.e()+NuEnergy[nuenesign], lvec.theta(), lvec.phi(),
                    JetResE*sigmaScaleFactor, JetResTheta, JetResPhi, lvec.m());
                 j2->setName("Jet2");
                 streamlog_out(DEBUG) << " start four-vector of second  jet: " << *j2  << std::endl ;
                }
               else if (i == 2) {
                 j3 = new JetFitObject (lvec.e()+NuEnergy[nuenesign], lvec.theta(), lvec.phi(),
                    JetResE*sigmaScaleFactor, JetResTheta, JetResPhi, lvec.m());
                 j3->setName("Jet3");
                 streamlog_out(DEBUG) << " start four-vector of third  jet: " << *j3  << std::endl ;
                }
               else if (i == 3) {
                 j4 = new JetFitObject (lvec.e()+NuEnergy[nuenesign], lvec.theta(), lvec.phi(),
                    JetResE*sigmaScaleFactor, JetResTheta, JetResPhi, lvec.m());
                 j4->setName("Jet4");
                   streamlog_out(DEBUG) << " start four-vector of forth  jet: " << *j4  << std::endl ;
                }
             }
             else{
               if (i == 0 ) {
                 j1 = new JetFitObject (lvec.e()+NuEnergy[nuenesign], lvec.theta(), lvec.phi(),
                    JetEnergyResolution(lvec.e()), _errtheta, _errphi, lvec.m());
                 j1->setName("Jet1");
                 streamlog_out(DEBUG)  << " start four-vector of first  jet: " << *j1  << std::endl ;
                }
               else if (i == 1) {
                 j2 = new JetFitObject (lvec.e()+NuEnergy[nuenesign], lvec.theta(), lvec.phi(),
                    JetEnergyResolution(lvec.e()), _errtheta, _errphi, lvec.m());
                 j2->setName("Jet2");
                 streamlog_out(DEBUG) << " start four-vector of second  jet: " << *j2  << std::endl ;
                }
               else if (i == 2) {
                 j3 = new JetFitObject (lvec.e()+NuEnergy[nuenesign], lvec.theta(), lvec.phi(),
                    JetEnergyResolution(lvec.e()), _errtheta, _errphi, lvec.m());
                 j3->setName("Jet3");
                 streamlog_out(DEBUG) << " start four-vector of third  jet: " << *j3  << std::endl ;
                }
               else if (i == 3) {
                 j4 = new JetFitObject (lvec.e()+NuEnergy[nuenesign], lvec.theta(), lvec.phi(),
                    JetEnergyResolution(lvec.e()), _errtheta, _errphi, lvec.m());
                 j4->setName("Jet4");
                   streamlog_out(DEBUG) << " start four-vector of forth  jet: " << *j4  << std::endl ;
                }
             }


#ifdef MARLIN_USE_AIDA
             hJetMass->fill(j->getMass());
#endif
          } //end of   if (reco particle j!=0)
        }//end loop over nJets

    #ifdef MARLIN_USE_AIDA
           double en, px, py, pz, mass;
           if (jrps[0] && jrps[1] && jrps[2] && jrps[3]) {
             en = jrps[0]->getEnergy()     +jrps[1]->getEnergy();
             px = jrps[0]->getMomentum()[0]+jrps[1]->getMomentum()[0];
             py = jrps[0]->getMomentum()[1]+jrps[1]->getMomentum()[1];
             pz = jrps[0]->getMomentum()[2]+jrps[1]->getMomentum()[2];
             mass = en*en-px*px-py*py-pz*pz;
             if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
             en = jrps[2]->getEnergy()     +jrps[3]->getEnergy();
             px = jrps[2]->getMomentum()[0]+jrps[3]->getMomentum()[0];
             py = jrps[2]->getMomentum()[1]+jrps[3]->getMomentum()[1];
             pz = jrps[2]->getMomentum()[2]+jrps[3]->getMomentum()[2];
             mass = en*en-px*px-py*py-pz*pz;
             if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
             en = jrps[0]->getEnergy()     +jrps[2]->getEnergy();
             px = jrps[0]->getMomentum()[0]+jrps[2]->getMomentum()[0];
             py = jrps[0]->getMomentum()[1]+jrps[2]->getMomentum()[1];
             pz = jrps[0]->getMomentum()[2]+jrps[2]->getMomentum()[2];
             mass = en*en-px*px-py*py-pz*pz;
             if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
             en = jrps[1]->getEnergy()     +jrps[3]->getEnergy();
             px = jrps[1]->getMomentum()[0]+jrps[3]->getMomentum()[0];
             py = jrps[1]->getMomentum()[1]+jrps[3]->getMomentum()[1];
             pz = jrps[1]->getMomentum()[2]+jrps[3]->getMomentum()[2];
             mass = en*en-px*px-py*py-pz*pz;
             if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
             en = jrps[0]->getEnergy()     +jrps[3]->getEnergy();
             px = jrps[0]->getMomentum()[0]+jrps[3]->getMomentum()[0];
             py = jrps[0]->getMomentum()[1]+jrps[3]->getMomentum()[1];
             pz = jrps[0]->getMomentum()[2]+jrps[3]->getMomentum()[2];
             mass = en*en-px*px-py*py-pz*pz;
             if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
             en = jrps[1]->getEnergy()     +jrps[2]->getEnergy();
             px = jrps[1]->getMomentum()[0]+jrps[2]->getMomentum()[0];
             py = jrps[1]->getMomentum()[1]+jrps[2]->getMomentum()[1];
             pz = jrps[1]->getMomentum()[2]+jrps[2]->getMomentum()[2];
             mass = en*en-px*px-py*py-pz*pz;
             if (mass >= 0) hTestHMassNoFitAll->fill( std::sqrt(mass) ) ;
           }
    #endif

           const int NJETS = 4;
           streamlog_out(MESSAGE)  << "*j1" << *j1  << "*j2" << *j2  << "*j3" << *j3  << "*j4" << *j4  << std::endl ;

           // these don't get changed by the fit -> to obtain start values later!
           JetFitObject startjets[NJETS] = {*j1, *j2, *j3, *j4};
           for (int i = 0; i < NJETS; ++i)
             streamlog_out(MESSAGE)  << "startjets[ " << i << "]: " << startjets[i]  << std::endl ;

           // these get changed by the fit -> reset after each permutation!
           JetFitObject fitjets[NJETS] = {*j1, *j2, *j3, *j4};
           for (int i = 0; i < NJETS; ++i)
             streamlog_out(MESSAGE)  << "fitjets[ " << i << "]: " << fitjets[i]  << std::endl ;

           // these point allways to the fitjets array, which gets reset.
           JetFitObject *jets[NJETS];
           for (int i = 0; i < NJETS; ++i) jets[i] = &fitjets[i];
           for (int i = 0; i < NJETS; ++i)
             streamlog_out(MESSAGE)  << "start four-vector of jets[ " << i << "]: " << *(jets[i])  << std::endl ;

           FourJetZHPairing pairing (jets);
           JetFitObject *permutedjets[NJETS];

           for (int iperm = 0; iperm < pairing.getNPerm(); iperm++) { //permutation begins

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

             pairing.nextPermutation (permutedjets);
             for (int i = 0; i < NJETS; ++i)
                streamlog_out(MESSAGE)  << "start four-vector of jet " << i << ": " << *(permutedjets[i])  << std::endl ;

             //MomentumConstraint pxc (1, 0);
             // crossing angle 14 mrad = 7/500
             MomentumConstraint pxc (0, 1, 0, 0, 3.5);//Factor for: (energy sum, px sum, py sum,pz sum,target value of sum)
             //3.5 due to crossing angle
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

                E_lab= 2*sqrt(std::pow(0.548579909e-3,2) + std::pow(250.,2) + std::pow(3.5,2) + 0. + 0.);

            MomentumConstraint ec(1, 0, 0, 0, E_lab);
            ec.setName("sum(E)");
            for (int i = 0; i < NJETS; ++i)
            ec.addToFOList (*(permutedjets[i]));

                streamlog_out(MESSAGE)  << "Value of pxc before fit: " << pxc.getValue() << std::endl ;
    	    streamlog_out(MESSAGE)  << "Value of pyc before fit: " << pyc.getValue() << std::endl ;
    	    streamlog_out(MESSAGE)  << "Value of pzc before fit: " << pzc.getValue() << std::endl ;
    	    streamlog_out(MESSAGE)  << "Value of ec before fit: " << ec.getValue() << std::endl ;


             // ISR Photon initialized with missing p_z
             ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);
                                                                                    //ISRPhotonFitObject(double px, double py, double pz,
                                                                                    //double b_, double PzMaxB_, double PzMinB_ = 0.);

    	 if(_fitISR){
                streamlog_out(MESSAGE)  << "start four-vector of ISR photon: " << *(photon) << std::endl ;

                pxc.addToFOList (*(photon));
                pyc.addToFOList (*(photon));
                pzc.addToFOList (*(photon));
                ec.addToFOList  (*(photon));
             }

             MassConstraint z(91.2);
             // SoftGaussMassConstraint z(91.2,2.4952/2);
             z.addToFOList (*(permutedjets[0]), 1);
             z.addToFOList (*(permutedjets[1]), 1);

             MassConstraint h(125.);
             h.addToFOList (*(permutedjets[2]), 1);
             h.addToFOList (*(permutedjets[3]), 1);
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

             for (int i = 0; i < NJETS; ++i)
                fitter.addFitObject (*(permutedjets[i]));
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
             //fitter.addConstraint (h);\

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
                streamlog_out(MESSAGE)  << "final four-vector of jet " << i << ": " << *(permutedjets[i]) << std::endl ;
                streamlog_out(MESSAGE)  << "final px of jet " << i << ": " << (permutedjets[i]) << std::endl ;
    	 }
             if(_fitISR){
                streamlog_out(MESSAGE)  << "final four-vector of ISR photon: " << *(photon) << std::endl ;
    	 }


             int ierr = fitter.getError();
             hFitError->fill( ierr ) ;
             if ((besterr > 0 && ierr < besterr) || ( besterr < 0 && ierr == 0)) besterr = ierr;

             // ierr == -1 only means that error calculation for fitted parameters failed!
             if (ierr <= 0) { //if fitter.getError() <=0
    #ifdef MARLIN_USE_AIDA
               hFitProbAll->fill( prob ) ;
               hNItAll->fill( nit ) ;
               hRecHMassAll->fill( h.getMass(1) ) ;
    #endif
               double pull[3][4];
               // require successfull error calculation for pulls!
               if (ierr == 0) {
                 for (int ifo = 0; ifo < 4; ifo++){
                   double start, fitted;
                   double errfit, errmea, sigma;
                   for (int ipar = 0; ipar < 3; ipar++) {
                     fitted = fitjets[ifo].getParam(ipar);
                     start = startjets[ifo].getParam(ipar);
                     errfit = fitjets[ifo].getError(ipar);
                     errmea = startjets[ifo].getError(ipar);
                     sigma = errmea*errmea-errfit*errfit;
                     if (sigma > 0) {
                      sigma = sqrt(sigma);
                      pull[ipar][ifo] = (fitted - start)/sigma;
                     }
                     else {
                      pull[ipar][ifo] = -4.5;
                     }
                   }
                   hPullEJetOK->fill (pull[0][ifo]);
                   hPullThJetOK->fill(pull[1][ifo]);
                   hPullPhJetOK->fill(pull[2][ifo]);
                 }
               }
    //           if (prob > bestprob && h.getMass(1) > 70 && w.getMass(1) < 150) {
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
                 Zmomentum [0]= permutedjets[0]->getPx()+permutedjets[1]->getPx();
                 Zmomentum [1]= permutedjets[0]->getPy()+permutedjets[1]->getPy();
                 Zmomentum [2]= permutedjets[0]->getPz()+permutedjets[1]->getPz();
                 Hmomentum [0]= permutedjets[2]->getPx()+permutedjets[3]->getPx();
                 Hmomentum [1]= permutedjets[2]->getPy()+permutedjets[3]->getPy();
                 Hmomentum [2]= permutedjets[2]->getPz()+permutedjets[3]->getPz();
                 Z_Energy= permutedjets[0]->getE()+permutedjets[1]->getE();
                 H_Energy=permutedjets[2]->getE()+permutedjets[3]->getE();
                 chi2best= fitter.getChi2();
                 errorcode= fitter.getError();
                if (ierr == 0) { //if  fitter.getError() is = 0
                   for (int ifo = 0; ifo < 4; ifo++){
                     hPullEJetBest->fill (pull[0][ifo]);
                     hPullThJetBest->fill(pull[1][ifo]);
                     hPullPhJetBest->fill(pull[2][ifo]);
                   }
                   hpull_jet1_E=pull[0][0];
                   hpull_jet2_E=pull[0][1];
                   hpull_jet3_E=pull[0][2];
                   hpull_jet4_E=pull[0][3];
                   hpull_jet1_th=pull[1][0];
                   hpull_jet2_th=pull[1][1];
                   hpull_jet3_th=pull[1][2];
                   hpull_jet4_th=pull[1][3];
                   hpull_jet1_phi=pull[2][0];
                   hpull_jet2_phi=pull[2][1];
                   hpull_jet3_phi=pull[2][2];
                   hpull_jet4_phi=pull[2][3];
                 }//end if  fitter.getError() is = 0
                 else {//if  fitter.getError() is not = 0
                   streamlog_out(DEBUG) << " ERROR CALCULATION FAILED for best permutation "
                                       << " in event " << evt->getEventNumber()
                                       << " for permutation " << iperm << std::endl ;
                   for (int ifo = 0; ifo < 4; ifo++){
                     hPullEJetBest->fill (-6.);
                     hPullThJetBest->fill(-6.);
                     hPullPhJetBest->fill(-6.);
                   }
                 }
               }
             }//end-if fitter.getError() <=0
             else {
                streamlog_out(DEBUG4) << "FIT ERROR = " << fitter.getError()
                                       << " in event " << evt->getEventNumber()
                                       << " for permutation " << iperm
                                       << ", not filling histograms!"  << std::endl ;
                                       streamlog_out(MESSAGE)  << "start mass of Z: " << startmassZ << std::endl ;
    	          streamlog_out(DEBUG4)  << "start mass of H: " << startmassH << std::endl ;
                streamlog_out(DEBUG4)  << "final mass of Z: " << z.getMass(1) << std::endl ;
    	          streamlog_out(DEBUG4)  << "final mass of H: " << h.getMass(1) << std::endl ;

    	         }
             delete photon;
             streamlog_out(DEBUG4) << "end permutation" << std::endl ;
           }//permutation ends
      }//neutrino correction + and - loop


//*******write something
      LCCollectionVec *OutputCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
      ReconstructedParticleImpl* ISRfitrec = new ReconstructedParticleImpl;
      ReconstructedParticleImpl* Zfitrec = new ReconstructedParticleImpl;
      ReconstructedParticleImpl* Hfitrec = new ReconstructedParticleImpl;
      ISRfitrec->setMomentum(ISRmomentum);
      ISRfitrec->setEnergy(bestphotonenergy);
      streamlog_out(DEBUG) << " Energy ISR:   "
                         << ISRfitrec->getEnergy() << std::endl ;
      ISRfitrec->setType (22);
      streamlog_out(DEBUG) << " IS ISR:   "
                         << ISRfitrec->getType() << std::endl ;
      OutputCol->addElement(ISRfitrec);


      Zfitrec->setMomentum(Zmomentum);
      streamlog_out(DEBUG) << "  Zmomentum :   "
			   << Zfitrec->getMomentum()[0] << "," << Zfitrec->getMomentum()[1]<<","<< Zfitrec->getMomentum()[2] << std::endl ;

      Zfitrec->setEnergy(Z_Energy);
      streamlog_out(DEBUG) << " Energy Z:   "
                         << Zfitrec->getEnergy() << std::endl ;
      Zfitrec->setMass(bestmassZ);
      streamlog_out(DEBUG) << " Mass Z:   "
                         << Zfitrec->getMass() << std::endl ;
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
      Hfitrec->setType (25);
      streamlog_out(DEBUG) << " IS H:   "
                         << Hfitrec->getType() << std::endl ;
      OutputCol->addElement(Hfitrec);


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
       delete j3;
       delete j4;

       ZHTree->Fill();
       evt->addCollection( OutputCol, _OutCol.c_str() );

    //}// end if # SLD-B =1
    // else{
    //   streamlog_out(ERROR) << "BREAK*********" << std::endl;
    // }
}//end if jetcol is not null


	streamlog_out(DEBUG4) << "next event******************" << std::endl;

  _nEvt ++ ;
}//event end



void ZH5CFit::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZH5CFit::end(){

    	streamlog_out(ERROR) << "# of events: " << _nEvt << std::endl;
    	streamlog_out(ERROR) << "# of nucorrection: " << correction<< std::endl;
    	streamlog_out(ERROR) << "# of Covariance failed: " << nCo<< std::endl;

  _fout->Write(0);
  _fout->Close();


}
