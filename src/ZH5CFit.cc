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
//#include <root/TLorentzVector.h>
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

  registerProcessorParameter( "ISRPzMax" ,
                              "Maximum possible energy for a single ISR photon",
                              _isrpzmax,
                              (float)225.);

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
                            "FitOutputColection",
                		        " Output Fit Colection" ,
                		         _OutCol,
                            std::string("FitReco")) ;
}


void ZH5CFit::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  b = (double) 0.00464564*( std::log(_ecm*_ecm*3814714.)-1. );
  //= 2*alpha/pi*( ln(s/m_e^2)-1 )
  ISRPzMaxB = std::pow((double)_isrpzmax,b);

}

void ZH5CFit::processRunHeader( LCRunHeader* ) {
  _nRun++ ;
}

void ZH5CFit::processEvent( LCEvent * evt ) {


    message<MESSAGE>( log()
		      << " processing event " << evt->getEventNumber()
		      << "  in run "          << evt->getRunNumber()
		      ) ;
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
      createHistogram1D( "hRecHMassNoFitBest", "M_H", 200, 0., 250. ) ;
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


  // fill histogram from LCIO data :

  //////////////////   JETS ///////////////////////////

     LCCollection* jetcol = evt->getCollection( _jetcolName ) ;
     if (jetcol != 0) {//jetcol is not null

       int nJETS = jetcol->getNumberOfElements()  ;
       message<MESSAGE>( log()
                      << " found " << nJETS
                      << " jets in event " << evt->getEventNumber()
                      << "  in run "          << evt->getRunNumber()
                      ) ;

       if (nJETS != 4) return;

       float yminus = jetcol ->parameters().getFloatVal( "YMinus");
       message<MESSAGE>( log()  << " yminus = " << yminus ) ;

       float yplus = jetcol ->parameters().getFloatVal( "YPlus");
       message<MESSAGE>( log()  << " yplus = " << yplus ) ;

  // original fit objects - save for next permutation
       JetFitObject* j1 = 0;
       JetFitObject* j2 = 0;
       JetFitObject* j3 = 0;
       JetFitObject* j4 = 0;

       ReconstructedParticle* jrps[4];

       for(int i=0; i< nJETS ; i++){

          ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>( jetcol->getElementAt( i ) ) ;


          if (j) {
             jrps[i] = j;
             message<MESSAGE>( log()
                       << " found jet in event " << evt->getEventNumber()
                       << "  in run "          << evt->getRunNumber()
                       ) ;
             lvec = HepLorentzVector ((j->getMomentum())[0],(j->getMomentum())[1],(j->getMomentum())[2],j->getEnergy());
             if (i == 0) {
               j1 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(),
                  JetEnergyResolution(lvec.e()), _errtheta, _errphi, lvec.m());
               j1->setName("Jet1");
               message<MESSAGE>( log()  << " start four-vector of first  jet: " << *j1  ) ;
	     }
             else if (i == 1) {
               j2 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(),
                  JetEnergyResolution(lvec.e()), _errtheta, _errphi, lvec.m());
               j2->setName("Jet2");
               message<MESSAGE>( log() << " start four-vector of second  jet: " << *j2  ) ;
	     }
             else if (i == 2) {
               j3 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(),
                  JetEnergyResolution(lvec.e()), _errtheta, _errphi, lvec.m());
               j3->setName("Jet3");
               message<MESSAGE>( log() << " start four-vector of third  jet: " << *j3  ) ;
	     }
             else if (i == 3) {
               j4 = new JetFitObject (lvec.e(), lvec.theta(), lvec.phi(),
                  JetEnergyResolution(lvec.e()), _errtheta, _errphi, lvec.m());
               j4->setName("Jet4");
               message<MESSAGE>( log() << " start four-vector of forth  jet: " << *j4  ) ;
	     }
#ifdef MARLIN_USE_AIDA
             hJetMass->fill(j->getMass());
#endif
} //end of   if (reco particle j!=0)
       }
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
       message<MESSAGE>( log()  << "*j1" << *j1  << "*j2" << *j2  << "*j3" << *j3  << "*j4" << *j4  ) ;

       // these don't get changed by the fit -> to obtain start values later!
       JetFitObject startjets[NJETS] = {*j1, *j2, *j3, *j4};
       for (int i = 0; i < NJETS; ++i)
         message<MESSAGE>( log()  << "startjets[ " << i << "]: " << startjets[i]  ) ;

       // these get changed by the fit -> reset after each permutation!
       JetFitObject fitjets[NJETS] = {*j1, *j2, *j3, *j4};
       for (int i = 0; i < NJETS; ++i)
         message<MESSAGE>( log()  << "fitjets[ " << i << "]: " << fitjets[i]  ) ;

       // these point allways to the fitjets array, which gets reset.
       JetFitObject *jets[NJETS];
       for (int i = 0; i < NJETS; ++i) jets[i] = &fitjets[i];
       for (int i = 0; i < NJETS; ++i)
         message<MESSAGE>( log()  << "start four-vector of jets[ " << i << "]: " << *(jets[i])  ) ;

       FourJetZHPairing pairing (jets);
       JetFitObject *permutedjets[NJETS];

       bestprob = 0.;
       bestnit = 0;
       //double bestmass1 = 0., bestmass2 = 0.;
       beststartmassZ = 0., beststartmassH = 0.;
       startmassZ = 0., startmassH = 0.;
       bestphotonenergy = 0.;
       int besterr = 999;
       double bestzvalue = 10000.;
       double chi2startmassZ = 0., chi2startmassH = 0.;
       LCCollectionVec *OutputCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
       ReconstructedParticleImpl* ISRfitrec = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* Zfitrec = new ReconstructedParticleImpl;
       ReconstructedParticleImpl* Hfitrec = new ReconstructedParticleImpl;
       // for (int i=0 ; i<3; i++){
       // Zmomentum [i]=0.;
       // Hmomentum [i]=0.;
       // ISRmomentum [i]=0.;
       // }
       memset(Zmomentum, 0, sizeof(Zmomentum));
       memset(Hmomentum, 0, sizeof(Hmomentum));
       memset(ISRmomentum, 0, sizeof(ISRmomentum));
       chi2best=0.;

       for (int iperm = 0; iperm < pairing.getNPerm(); iperm++) { //permutation begins

         message<MESSAGE>( log()
                       << " ================================================= "
                       ) ;
         message<MESSAGE>( log()
                       << " iperm = " << iperm
                       ) ;

         // important: (re-)set fitjets array!
         fitjets[0] = *j1;
         fitjets[1] = *j2;
         fitjets[2] = *j3;
         fitjets[3] = *j4;

         pairing.nextPermutation (permutedjets);
         for (int i = 0; i < NJETS; ++i)
            message<MESSAGE>( log()  << "start four-vector of jet " << i << ": " << *(permutedjets[i])  ) ;

         //MomentumConstraint pxc (1, 0);
         // crossing angle 14 mrad = 7/500
         MomentumConstraint pxc (0, 1, 0, 0, 3.5);//Factor for: (energy sum, px sum, py sum,pz sum,t arget value of sum)

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

         message<MESSAGE>( log() << "ECM = " << _ecm  );
	       MomentumConstraint ec(1, 0, 0, 0, _ecm);
         ec.setName("sum(E)");
         for (int i = 0; i < NJETS; ++i)
            ec.addToFOList (*(permutedjets[i]));

            message<MESSAGE>( log()  << "Value of pxc before fit: " << pxc.getValue() ) ;
	    message<MESSAGE>( log()  << "Value of pyc before fit: " << pyc.getValue() ) ;
	    message<MESSAGE>( log()  << "Value of pzc before fit: " << pzc.getValue() ) ;
	    message<MESSAGE>( log()  << "Value of ec before fit: " << ec.getValue() ) ;


         // ISR Photon initialized with missing p_z
         ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);
                                                                                //ISRPhotonFitObject(double px, double py, double pz,
                                                                                //double b_, double PzMaxB_, double PzMinB_ = 0.);

	 if(_fitISR){
            message<MESSAGE>( log()  << "start four-vector of ISR photon: " << *(photon) ) ;

            pxc.addToFOList (*(photon));
            pyc.addToFOList (*(photon));
            pzc.addToFOList (*(photon));
            ec.addToFOList  (*(photon));
         }

         MassConstraint z(91.2);
         z.addToFOList (*(permutedjets[0]), 1);
         z.addToFOList (*(permutedjets[1]), 1);

         MassConstraint h(125.);
         h.addToFOList (*(permutedjets[2]), 1);
         h.addToFOList (*(permutedjets[3]), 1);
         message<MESSAGE>( log()  << "final mass of Z: " << z.getMass(1) ) ;

         startmassZ = z.getMass(1);
         startmassH = h.getMass(1);

	 message<MESSAGE>( log() << "start mass of Z: " << startmassZ ) ;
	 message<MESSAGE>( log() << "start mass of H: " << startmassH ) ;

#ifdef MARLIN_USE_AIDA
         hRecHMassNoFitAll->fill( startmassH ) ;
         //hRecHMassNoFitAll->fill( startmassH ) ;
#endif
         BaseFitter *pfitter;
         if (_ifitter == 1) {
           pfitter = new NewFitterGSL();
           if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<NewFitterGSL*>(pfitter))->setDebug (debug);
         }
         else if (_ifitter == 2) {
           pfitter = new NewtonFitterGSL();
           if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<NewtonFitterGSL*>(pfitter))->setDebug (debug);
         }
         else {
           // OPALFitter has no method setDebug !
           pfitter = new OPALFitterGSL();
           if (evt->getEventNumber()== _ievttrace || _traceall) (dynamic_cast<OPALFitterGSL*>(pfitter))->setDebug (debug);
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
         fitter.addConstraint (z);
         // don't constrain Higgs mass, just use constraints for convenient mass calculation
         //fitter.addConstraint (h);\

         // initial value of Z mass constraint
         if (fabs(startmassZ-91.2) + fabs(startmassH-125.) < bestzvalue) {
           chi2startmassZ = startmassZ;
           chi2startmassH = startmassH;
           bestzvalue = fabs(startmassZ-91.2) + fabs(startmassH-125.);
         }


         prob = fitter.fit();
         double chi2 = fitter.getChi2();
         nit = fitter.getIterations();

         message<MESSAGE>( log() << "fit probability = " << prob ) ;
         message<MESSAGE>( log() << "fit chi2 = " << chi2  ) ;
         message<MESSAGE>( log() << "error code: " << fitter.getError() ) ;

         for (int i = 0; i < NJETS; ++i) {
            message<MESSAGE>( log()  << "final four-vector of jet " << i << ": " << *(permutedjets[i])) ;
            message<MESSAGE>( log()  << "final px of jet " << i << ": " << (permutedjets[i])) ;
	 }
         if(_fitISR){
            message<MESSAGE>( log()  << "final four-vector of ISR photon: " << *(photon) ) ;
	 }


         int ierr = fitter.getError();
         hFitError->fill( ierr ) ;
         if ((besterr > 0 && ierr < besterr) || ( besterr < 0 && ierr == 0)) besterr = ierr;

         // ierr == -1 only means that error calculation for fitted parameters failed!
         if (ierr <= 0) {
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
             chi2best= fitter.getChi2();
            if (ierr == 0) {
               for (int ifo = 0; ifo < 4; ifo++){
                 hPullEJetBest->fill (pull[0][ifo]);
                 hPullThJetBest->fill(pull[1][ifo]);
                 hPullPhJetBest->fill(pull[2][ifo]);
               }
             }
             else {
               message<WARNING>( log() << " ERROR CALCULATION FAILED for best permutation "
                                   << " in event " << evt->getEventNumber()
                                   << " for permutation " << iperm ) ;
               for (int ifo = 0; ifo < 4; ifo++){
                 hPullEJetBest->fill (-6.);
                 hPullThJetBest->fill(-6.);
                 hPullPhJetBest->fill(-6.);
               }
             }
           }
         }
         else {
            message<WARNING>( log() << "FIT ERROR = " << fitter.getError()
                                   << " in event " << evt->getEventNumber()
                                   << " for permutation " << iperm
                                   << ", not filling histograms!"  ) ;
                                   message<WARNING>( log()  << "start mass of Z: " << startmassZ) ;
	          message<WARNING>( log()  << "start mass of H: " << startmassH ) ;
            message<WARNING>( log()  << "final mass of Z: " << z.getMass(1) ) ;
	          message<WARNING>( log()  << "final mass of H: " << h.getMass(1) ) ;

	 }
         delete photon;
         message<MESSAGE>( log() << "end permutation") ;
       }//permutation ends


      ISRfitrec->setMomentum(ISRmomentum);
      ISRfitrec->setEnergy(bestphotonenergy);
      message<DEBUG>( log() << " Energy ISR:   "
                         << ISRfitrec->getEnergy()) ;
      ISRfitrec->setType (22);
      message<DEBUG>( log() << " IS ISR:   "
                         << ISRfitrec->getType()) ;
      OutputCol->addElement(ISRfitrec);


      Zfitrec->setMomentum(Zmomentum);
      message<DEBUG>( log() << "  Zmomentum :   "
                         << Zfitrec->getMomentum()[0] << "," << Zfitrec->getMomentum()[1]<<","<< Zfitrec->getMomentum()[2]);

      Zfitrec->setEnergy(permutedjets[0]->getE()+permutedjets[0]->getE());
      message<DEBUG>( log() << " Energy Z:   "
                         << Zfitrec->getEnergy()) ;
      Zfitrec->setMass(bestmassZ);
      message<DEBUG>( log() << " Mass Z:   "
                         << Zfitrec->getMass()) ;
      Zfitrec->setType (23);
      message<DEBUG>( log() << " IS Z :   "
                         << Zfitrec->getType()) ;
      OutputCol->addElement(Zfitrec);


      Hfitrec->setMomentum(Hmomentum);
      Hfitrec->setEnergy(permutedjets[2]->getE()+permutedjets[3]->getE());
      message<DEBUG>( log() << " Energy H:   "
                         << Hfitrec->getEnergy()) ;
      Hfitrec->setMass(bestmassH);
      message<DEBUG>( log() << " Mass H:   "
                         << Hfitrec->getMass()) ;
      Hfitrec->setType (25);
      message<DEBUG>( log() << " IS H:   "
                         << Hfitrec->getType()) ;
      OutputCol->addElement(Hfitrec);


	    OutputCol->parameters().setValue("bestchisq", (float)chi2best);
      message<DEBUG>( log() << " chi2:   " << chi2best );
      OutputCol->parameters().setValue("best_prob", (float)bestprob);
      message<DEBUG>( log() << " prob:   " << bestprob );

       message<WARNING>( log() << "==============  end of permutations for event " << evt->getEventNumber() <<  " ==============" ) ;
       message<WARNING>( log()  << "min chi2 start mass of Z: " << chi2startmassZ ) ;
       message<WARNING>( log()  << "min chi2 start mass of H: " << chi2startmassH ) ;
       message<WARNING>( log()  << "best start mass of Z: " << beststartmassZ ) ;
       message<WARNING>( log()  << "best start mass of H: " << beststartmassH ) ;
       message<WARNING>( log()  << "best mass of Z: " << bestmassZ ) ;
       message<WARNING>( log()  << "best mass of H: " << bestmassH ) ;


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
       //message<MESSAGE>( log()  << "final mass of Z: " << z.getMass(1) ) ;

       delete j1;
       delete j2;
       delete j3;
       delete j4;


       evt->addCollection( OutputCol, _OutCol.c_str() );
}//end if jetcol is not null



  _nEvt ++ ;
}



void ZH5CFit::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZH5CFit::end(){

}
