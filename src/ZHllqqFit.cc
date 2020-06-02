#include "ZHllqqFit.h"
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
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace CLHEP ;


ZHllqqFit aZHllqqFit ;

// function to define the jet energy resolution (in GeV)
double ZHllqqFit::JetEnergyResolution(double E)
{
	// examples here derived by Benjamin Hermberg from e+e- -> udsc:
	// 1) default is 120%/sqrt(E), gives best convergence of 5C fit on e+e- -> udsc
	double result = m_jetEnergyError*std::sqrt(E);

	// 2) comparing jet-level to quark-level energies
	//    (using MarlinReco/Analysis/RecoMCTruthLink/QuarkJetPairing.cc)
	if (m_jetEnergyError == 0 ) result = std::sqrt(pow(0.6908,2)*(E)+(pow(0.02596,2)*pow(E,2)));

	return result;
}

ZHllqqFit::ZHllqqFit() :
Processor("ZHllqqFit"),
m_includeHbb(false),
m_includeHcc(false),
m_includeHother(false),
m_nRun(0),
m_nEvt(0),
m_nRunSum(0),
m_nEvtSum(0),
m_Bfield(0.f),
c(0.),
mm2m(0.),
eV2GeV(0.),
eB(0.),
m_ECM(250.f),
m_nSLDecayBHadron(0),
m_nSLDecayCHadron(0),
m_nSLDecayTotal(0),
m_nJets(0),
m_nLeptons(0),
m_iError_wNu_bestfit(-5),
m_probability_wNu_bestfit(0.),
m_chi2_wNu_bestfit(0.),
m_n_itter_wNu_bestfit(0),
m_startmassZ_wNu_bestfit(0.),
m_startmassH_wNu_bestfit(0.),
m_beststartmassZ_wNu_bestfit(0.),
m_beststartmassH_wNu_bestfit(0.),
m_Zmass_after_fit_wNu_bestfit(0.),
m_Hmass_after_fit_wNu_bestfit(0.),
m_bestphotonenergy_wNu_bestfit(0.),
m_chi2startmassZ_wNu_bestfit(0.),
m_chi2startmassH_wNu_bestfit(0.),
m_iError_woNu(-5),
m_probability_woNu(0.),
m_chi2best_woNu(0.),
m_n_itter_woNu(0),
m_startmassZ_woNu(0.),
m_startmassH_woNu(0.),
m_beststartmassZ_woNu(0.),
m_beststartmassH_woNu(0.),
m_Zmass_after_fit_woNu(0.),
m_Hmass_after_fit_woNu(0.),
m_bestphotonenergy_woNu(0.),
m_chi2startmassZ_woNu(0.),
m_chi2startmassH_woNu(0.),
m_iError_best(-5),
m_probability_best(0.),
m_chi2_best(0.),
m_n_itter_best(0),
m_startmassZ_best(0.),
m_startmassH_best(0.),
m_beststartmassZ_best(0.),
m_beststartmassH_best(0.),
m_Zmass_after_fit_best(0.),
m_Hmass_after_fit_best(0.),
m_bestphotonenergy_best(0.),
m_chi2startmassZ_best(0.),
m_chi2startmassH_best(0.),
chi2best(0.),
errorcode(0),
E_lab(0.),
m_pxc_before_ISR_wNu(0.),
m_pyc_before_ISR_wNu(0.),
m_pzc_before_ISR_wNu(0.),
m_ec_before_ISR_wNu(0.),
m_pxc_before_fit_wNu(0.),
m_pyc_before_fit_wNu(0.),
m_pzc_before_fit_wNu(0.),
m_ec_before_fit_wNu(0.),
m_pxc_after_fit_wNu(0.),
m_pyc_after_fit_wNu(0.),
m_pzc_after_fit_wNu(0.),
m_ec_after_fit_wNu(0.),
m_pxc_before_ISR_woNu(0.),
m_pyc_before_ISR_woNu(0.),
m_pzc_before_ISR_woNu(0.),
m_ec_before_ISR_woNu(0.),
m_pxc_before_fit_woNu(0.),
m_pyc_before_fit_woNu(0.),
m_pzc_before_fit_woNu(0.),
m_ec_before_fit_woNu(0.),
m_pxc_after_fit_woNu(0.),
m_pyc_after_fit_woNu(0.),
m_pzc_after_fit_woNu(0.),
m_ec_after_fit_woNu(0.),
m_pxc_before_ISR_best(0.),
m_pyc_before_ISR_best(0.),
m_pzc_before_ISR_best(0.),
m_ec_before_ISR_best(0.),
m_pxc_before_fit_best(0.),
m_pyc_before_fit_best(0.),
m_pzc_before_fit_best(0.),
m_ec_before_fit_best(0.),
m_pxc_after_fit_best(0.),
m_pyc_after_fit_best(0.),
m_pzc_after_fit_best(0.),
m_ec_after_fit_best(0.),
m_pTFile(NULL),
m_pTTree(NULL),
m_pTTree_0(NULL),
m_pTTree_1(NULL),
m_pTTree_2(NULL),
m_pTTree_3(NULL),
h_Zmass_beforefit_woNu(NULL),
h_Hmass_beforefit_woNu(NULL),
h_Zmass_beforefit_wNu(NULL),
h_Hmass_beforefit_wNu(NULL),
h_Zmass_afterfit_woNu(NULL),
h_Hmass_afterfit_woNu(NULL),
h_Zmass_afterfit_wNu(NULL),
h_Hmass_afterfit_wNu(NULL),
h_fitError_wNu(NULL),
h_fitError_woNu(NULL),
h_ErrorCode_wNu_woNu(NULL),
h_fitProbability_wNu(NULL),
h_fitProbability_woNu(NULL),
h_fitProbability_best(NULL),
h_nJets(NULL),
h_nLeptons(NULL),
h_nLeptons_nJets(NULL),
h_ISR_mcp_fit(NULL),
h_ISR_pzc_wNu(NULL),
h_ISR_pzc_woNu(NULL),
h_ISR_pzc_best(NULL),
h_pull_jet_E(NULL),
h_pull_jet_theta(NULL),
h_pull_jet_phi(NULL),
h_pull_lepton_InvPt(NULL),
h_pull_lepton_theta(NULL),
h_pull_lepton_phi(NULL)

{

//	modify processor description
	_description = "ZHllqqFit does a fit on 2 jet events (Px, Py, Pz, E, M12 = MZ)" ;


//	register steering parameters: name, description, class-variable, default value

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"JetCollectionName" ,
				"Name of the Jet collection"  ,
				jetcollection ,
				std::string("ISOLeptons")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"HiggsDecayMode" ,
				"Higgs decay mode"  ,
				hDecayMode ,
				std::string("HdecayMode")
				);

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"ISOLeptonCollectionName" ,
				"Name of the Isolated Lepton collection"  ,
				leptoncollection ,
				std::string("Durham2Jets")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"SemiLeptonicDecays",
				"Semi-Leptonic Decays Collection",
				SLDecayCollection,
				std::string("")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"NeutrinoCorrectionB",
				"Collection of Corrected/Estimated neutrino energies from SemiLeptonic decays of B-Hadron",
				NuEnergyCollectionB,
				std::string("NuCorrectB")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"NeutrinoCorrectionC",
				"Collection of Corrected/Estimated neutrino energies from SemiLeptonic decays of B-Hadron",
				NuEnergyCollectionC,
				std::string("NuCorrectC")
				);

	registerInputCollection( LCIO::MCPARTICLE,
				"MCParticleCollection" ,
				"Name of the MCParticle collection"  ,
				MCPCcollection ,
				std::string("MCParticle")
				);

	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"ErrorFlowCollection" ,
				"Name of the ErrorFlow collection"  ,
				errorflowcollection ,
				std::string("OutputErrorFlowJets")
				);

	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"FitOutputColection",
				" Output Fit Colection" ,
				outputFitcollection,
				std::string("FitReco")
				);

	registerProcessorParameter("includeHbb",
				"Include H->bb decays",
				m_includeHbb,
				bool(false)
				);

	registerProcessorParameter("includeHcc",
				"Include H->cc decays",
				m_includeHcc,
				bool(false)
				);

	registerProcessorParameter("includeHother",
				"Include other Higgs decay modes",
				m_includeHother,
				bool(false)
				);

	registerProcessorParameter("includeISR",
				"Include ISR in fit hypothesis; false: without ISR , true: with ISR",
				m_fitISR,
				bool(true)
				);

	registerProcessorParameter("includeNuCorrection",
				"Include neutrino correction in fit hypothesis; false: without Nu Correction , true: with Nu Correction",
				m_fitNuE,
				bool(true)
				);

	registerProcessorParameter( "ECM" ,
				"Center-of-Mass Energy in GeV",
				m_ECM,
				float(250.f)
				);

	registerProcessorParameter( "useErrorFlow",
				"If true, use covariance matrix for energy uncertainty. Otherwise 1.2/sqrt(E)",
				m_useErrorFlow,
				bool(true)
				);

	registerProcessorParameter( "ISRPzMax" ,
				"Maximum possible energy for a single ISR photon",
				m_isrpzmax,
				(float)35.
				);

	registerProcessorParameter( "SigmaEnergyScaleFactor",
				"Scale Factor t be applied to jet energy uncertainty",
				sigmaScaleFactor,
				(float) 1.0
				);

	registerProcessorParameter( "errene" ,
				"assumed energy resolution for jets as x/sqrt(E) - if 0, then parametrisation is used",
				m_jetEnergyError,
				(double)1.2
				);

	registerProcessorParameter( "errtheta" ,
				"assumed theta resolution for jet axis",
				m_jetThetaError,
				(double)0.1
				);

	registerProcessorParameter( "errphi" ,
				"assumed phi resolution for jet axis",
				m_jetPhiError,
				(double)0.1
				);

	registerProcessorParameter( "fitter" ,
				"0 = OPALFitter, 1 = NewFitter, 2 = NewtonFitter",
				m_fitter,
				(int)0
				);

	registerProcessorParameter( "traceall" ,
				"set true if every event should be traced",
				m_traceall,
				(bool)false
				);

	registerProcessorParameter( "ievttrace" ,
				"number of individual event to be traced",
				m_ievttrace,
				(int)0
				);

	registerProcessorParameter("outputFilename",
				"name of output file",
				m_outputFile,
				std::string("")
				);

}

void ZHllqqFit::init()
{
//	usually a good idea to
	streamlog_out(DEBUG) << "   init called  " << std::endl;
	this->Clear();
	m_Bfield = MarlinUtil::getBzAtOrigin();
	printParameters();
	streamlog_out(DEBUG) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;

	m_nRun = 0;
	m_nEvt = 0;
	m_nRunSum = 0;
	m_nEvtSum = 0;

	b = (double) 0.00464564 * ( std::log( m_ECM * m_ECM * 3814714. ) - 1. );
//	  = 2*alpha/pi*( ln(s/m_e^2)-1 )
	ISRPzMaxB = std::pow((double)m_isrpzmax,b);

	m_pTFile = new TFile(m_outputFile.c_str(),"recreate");

	m_pTTree = new TTree("eventTree","eventTree");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_pTTree->Branch("nJets",&m_nJets,"nJets/I") ;
	m_pTTree->Branch("nLeptons",&m_nLeptons,"nLeptons/I") ;

	m_pTTree_0 = new TTree("ZHlljjTree","ZHlljjTree");
	m_pTTree_0->SetDirectory(m_pTFile);
	m_pTTree_0->Branch("nSLDecayBHadron",&m_nSLDecayBHadron,"nSLDecayBHadron/I") ;
	m_pTTree_0->Branch("nSLDecayCHadron",&m_nSLDecayCHadron,"nSLDecayCHadron/I") ;
	m_pTTree_0->Branch("nSLDecayTotal",&m_nSLDecayTotal,"nSLDecayTotal/I") ;
	m_pTTree_0->Branch("bestNuCombination",&m_bestNuCombination) ;
	m_pTTree_0->Branch("iError_wNu",&m_iError_wNu) ;
	m_pTTree_0->Branch("probability_wNu",&m_probability_wNu) ;
	m_pTTree_0->Branch("chi2_wNu",&m_chi2_wNu) ;
	m_pTTree_0->Branch("n_itter_wNu",&m_n_itter_wNu) ;
	m_pTTree_0->Branch("startmassZ_wNu",&m_startmassZ_wNu) ;
	m_pTTree_0->Branch("startmassH_wNu",&m_startmassH_wNu) ;
	m_pTTree_0->Branch("beststartmassZ_wNu",&m_beststartmassZ_wNu) ;
	m_pTTree_0->Branch("beststartmassH_wNu",&m_beststartmassH_wNu) ;
	m_pTTree_0->Branch("Zmass_after_fit_wNu",&m_Zmass_after_fit_wNu) ;
	m_pTTree_0->Branch("Hmass_after_fit_wNu",&m_Hmass_after_fit_wNu) ;
	m_pTTree_0->Branch("chi2startmassZ_wNu",&m_chi2startmassZ_wNu) ;
	m_pTTree_0->Branch("chi2startmassH_wNu",&m_chi2startmassH_wNu) ;
	m_pTTree_0->Branch("iError_wNu_bestfit",&m_iError_wNu_bestfit,"iError_wNu_bestfit/I") ;
	m_pTTree_0->Branch("probability_wNu_bestfit",&m_probability_wNu_bestfit,"probability_wNu_bestfit/F") ;
	m_pTTree_0->Branch("chi2_wNu_bestfit",&m_chi2_wNu_bestfit,"chi2_wNu_bestfit/F") ;
	m_pTTree_0->Branch("n_itter_wNu_bestfit",&m_n_itter_wNu_bestfit,"n_itter_wNu_bestfit/I") ;
	m_pTTree_0->Branch("startmassZ_wNu_bestfit",&m_startmassZ_wNu_bestfit,"startmassZ_wNu_bestfit/F") ;
	m_pTTree_0->Branch("startmassH_wNu_bestfit",&m_startmassH_wNu_bestfit,"startmassH_wNu_bestfit/F") ;
	m_pTTree_0->Branch("beststartmassZ_wNu_bestfit",&m_beststartmassZ_wNu_bestfit,"beststartmassZ_wNu_bestfit/F") ;
	m_pTTree_0->Branch("beststartmassH_wNu_bestfit",&m_beststartmassH_wNu_bestfit,"beststartmassH_wNu_bestfit/F") ;
	m_pTTree_0->Branch("Zmass_after_fit_wNu_bestfit",&m_Zmass_after_fit_wNu_bestfit,"Zmass_after_fit_wNu_bestfit/F") ;
	m_pTTree_0->Branch("Hmass_after_fit_wNu_bestfit",&m_Hmass_after_fit_wNu_bestfit,"Hmass_after_fit_wNu_bestfit/F") ;
	m_pTTree_0->Branch("chi2startmassZ_wNu_bestfit",&m_chi2startmassZ_wNu_bestfit,"chi2startmassZ_wNu_bestfit/F") ;
	m_pTTree_0->Branch("chi2startmassH_wNu_bestfit",&m_chi2startmassH_wNu_bestfit,"chi2startmassH_wNu_bestfit/F") ;
	m_pTTree_0->Branch("bestphotonenergy_wNu_bestfit",&m_bestphotonenergy_wNu_bestfit,"bestphotonenergy_wNu_bestfit/F") ;
	m_pTTree_0->Branch("iError_woNu",&m_iError_woNu,"iError_woNu/I") ;
	m_pTTree_0->Branch("probability_woNu",&m_probability_woNu,"probability_woNu/F") ;
	m_pTTree_0->Branch("chi2best_woNu",&m_chi2best_woNu,"chi2best_woNu/F") ;
	m_pTTree_0->Branch("n_itter_woNu",&m_n_itter_woNu,"n_itter_woNu/I") ;
	m_pTTree_0->Branch("startmassZ_woNu",&m_startmassZ_woNu,"startmassZ_woNu/F") ;
	m_pTTree_0->Branch("startmassH_woNu",&m_startmassH_woNu,"startmassH_woNu/F") ;
	m_pTTree_0->Branch("beststartmassZ_woNu",&m_beststartmassZ_woNu,"beststartmassZ_woNu/F") ;
	m_pTTree_0->Branch("beststartmassH_woNu",&m_beststartmassH_woNu,"beststartmassH_woNu/F") ;
	m_pTTree_0->Branch("Zmass_after_fit_woNu",&m_Zmass_after_fit_woNu,"Zmass_after_fit_woNu/F") ;
	m_pTTree_0->Branch("Hmass_after_fit_woNu",&m_Hmass_after_fit_woNu,"Hmass_after_fit_woNu/F") ;
	m_pTTree_0->Branch("chi2startmassZ_woNu",&m_chi2startmassZ_woNu,"chi2startmassZ_woNu/F") ;
	m_pTTree_0->Branch("chi2startmassH_woNu",&m_chi2startmassH_woNu,"chi2startmassH_woNu/F") ;
	m_pTTree_0->Branch("bestphotonenergy_woNu",&m_bestphotonenergy_woNu,"bestphotonenergy_woNu/F") ;
	m_pTTree_0->Branch("iError_best",&m_iError_best,"iError_best/I") ;
	m_pTTree_0->Branch("probability_best",&m_probability_best,"probability_best/F") ;
	m_pTTree_0->Branch("chi2_best",&m_chi2_best,"chi2_best/F") ;
	m_pTTree_0->Branch("n_itter_best",&m_n_itter_best,"n_itter_best/I") ;
	m_pTTree_0->Branch("startmassZ_best",&m_startmassZ_best,"startmassZ_best/F") ;
	m_pTTree_0->Branch("startmassH_best",&m_startmassH_best,"startmassH_best/F") ;
	m_pTTree_0->Branch("beststartmassZ_best",&m_beststartmassZ_best,"beststartmassZ_best/F") ;
	m_pTTree_0->Branch("beststartmassH_best",&m_beststartmassH_best,"beststartmassH_best/F") ;
	m_pTTree_0->Branch("Zmass_after_fit_best",&m_Zmass_after_fit_best,"Zmass_after_fit_best/F") ;
	m_pTTree_0->Branch("Hmass_after_fit_best",&m_Hmass_after_fit_best,"Hmass_after_fit_best/F") ;
	m_pTTree_0->Branch("chi2startmassZ_best",&m_chi2startmassZ_best,"chi2startmassZ_best/F") ;
	m_pTTree_0->Branch("chi2startmassH_best",&m_chi2startmassH_best,"chi2startmassH_best/F") ;
	m_pTTree_0->Branch("bestphotonenergy_best",&m_bestphotonenergy_best,"bestphotonenergy_best/F") ;
	m_pTTree_1 = new TTree("pulls","pulls");
	m_pTTree_1->SetDirectory(m_pTFile);
	m_pTTree_1->Branch("pull_jet_E_wNu",&m_pull_jet_E_wNu) ;
	m_pTTree_1->Branch("pull_jet_th_wNu",&m_pull_jet_th_wNu) ;
	m_pTTree_1->Branch("pull_jet_phi_wNu",&m_pull_jet_phi_wNu) ;
	m_pTTree_1->Branch("pull_lepton_InvpT_wNu",&m_pull_lepton_InvpT_wNu) ;
	m_pTTree_1->Branch("pull_lepton_th_wNu",&m_pull_lepton_th_wNu) ;
	m_pTTree_1->Branch("pull_lepton_phi_wNu",&m_pull_lepton_phi_wNu) ;
	m_pTTree_1->Branch("pull_jet_E_wNu_bestfit",&m_pull_jet_E_wNu_bestfit);
	m_pTTree_1->Branch("pull_jet_th_wNu_bestfit",&m_pull_jet_th_wNu_bestfit);
	m_pTTree_1->Branch("pull_jet_phi_wNu_bestfit",&m_pull_jet_phi_wNu_bestfit);
	m_pTTree_1->Branch("pull_lepton_InvpT_wNu_bestfit",&m_pull_lepton_InvpT_wNu_bestfit);
	m_pTTree_1->Branch("pull_lepton_th_wNu_bestfit",&m_pull_lepton_th_wNu_bestfit);
	m_pTTree_1->Branch("pull_lepton_phi_wNu_bestfit",&m_pull_lepton_phi_wNu_bestfit);
	m_pTTree_1->Branch("pull_jet_E_woNu",&m_pull_jet_E_woNu) ;
	m_pTTree_1->Branch("pull_jet_th_woNu",&m_pull_jet_th_woNu) ;
	m_pTTree_1->Branch("pull_jet_phi_woNu",&m_pull_jet_phi_woNu) ;
	m_pTTree_1->Branch("pull_lepton_InvpT_woNu",&m_pull_lepton_InvpT_woNu) ;
	m_pTTree_1->Branch("pull_lepton_th_woNu",&m_pull_lepton_th_woNu) ;
	m_pTTree_1->Branch("pull_lepton_phi_woNu",&m_pull_lepton_phi_woNu) ;
	m_pTTree_1->Branch("pull_jet_E_best",&m_pull_jet_E_best);
	m_pTTree_1->Branch("pull_jet_th_best",&m_pull_jet_th_best);
	m_pTTree_1->Branch("pull_jet_phi_best",&m_pull_jet_phi_best);
	m_pTTree_1->Branch("pull_lepton_InvpT_best",&m_pull_lepton_InvpT_best);
	m_pTTree_1->Branch("pull_lepton_th_best",&m_pull_lepton_th_best);
	m_pTTree_1->Branch("pull_lepton_phi_best",&m_pull_lepton_phi_best);
	m_pTTree_2 = new TTree("constraints","constraints");
	m_pTTree_2->SetDirectory(m_pTFile);
	m_pTTree_2->Branch("pxc_before_ISR_wNu",&m_pxc_before_ISR_wNu,"pxc_before_ISR_wNu/F") ;
	m_pTTree_2->Branch("pyc_before_ISR_wNu",&m_pyc_before_ISR_wNu,"pyc_before_ISR_wNu/F") ;
	m_pTTree_2->Branch("pzc_before_ISR_wNu",&m_pzc_before_ISR_wNu,"pzc_before_ISR_wNu/F") ;
	m_pTTree_2->Branch("ec_before_ISR_wNu",&m_ec_before_ISR_wNu,"ec_before_ISR_wNu/F") ;
	m_pTTree_2->Branch("pxc_before_fit_wNu",&m_pxc_before_fit_wNu,"pxc_before_fit_wNu/F") ;
	m_pTTree_2->Branch("pyc_before_fit_wNu",&m_pyc_before_fit_wNu,"pyc_before_fit_wNu/F") ;
	m_pTTree_2->Branch("pzc_before_fit_wNu",&m_pzc_before_fit_wNu,"pzc_before_fit_wNu/F") ;
	m_pTTree_2->Branch("ec_before_fit_wNu",&m_ec_before_fit_wNu,"ec_before_fit_wNu/F") ;
	m_pTTree_2->Branch("pxc_after_fit_wNu",&m_pxc_after_fit_wNu,"pxc_after_fit_wNu/F") ;
	m_pTTree_2->Branch("pyc_after_fit_wNu",&m_pyc_after_fit_wNu,"pyc_after_fit_wNu/F") ;
	m_pTTree_2->Branch("pzc_after_fit_wNu",&m_pzc_after_fit_wNu,"pzc_after_fit_wNu/F") ;
	m_pTTree_2->Branch("ec_after_fit_wNu",&m_ec_after_fit_wNu,"ec_after_fit_wNu/F") ;
	m_pTTree_2->Branch("pxc_before_ISR_woNu",&m_pxc_before_ISR_woNu,"pxc_before_ISR_woNu/F") ;
	m_pTTree_2->Branch("pyc_before_ISR_woNu",&m_pyc_before_ISR_woNu,"pyc_before_ISR_woNu/F") ;
	m_pTTree_2->Branch("pzc_before_ISR_woNu",&m_pzc_before_ISR_woNu,"pzc_before_ISR_woNu/F") ;
	m_pTTree_2->Branch("ec_before_ISR_woNu",&m_ec_before_ISR_woNu,"ec_before_ISR_woNu/F") ;
	m_pTTree_2->Branch("pxc_before_fit_woNu",&m_pxc_before_fit_woNu,"pxc_before_fit_woNu/F") ;
	m_pTTree_2->Branch("pyc_before_fit_woNu",&m_pyc_before_fit_woNu,"pyc_before_fit_woNu/F") ;
	m_pTTree_2->Branch("pzc_before_fit_woNu",&m_pzc_before_fit_woNu,"pzc_before_fit_woNu/F") ;
	m_pTTree_2->Branch("ec_before_fit_woNu",&m_ec_before_fit_woNu,"ec_before_fit_woNu/F") ;
	m_pTTree_2->Branch("pxc_after_fit_woNu",&m_pxc_after_fit_woNu,"pxc_after_fit_woNu/F") ;
	m_pTTree_2->Branch("pyc_after_fit_woNu",&m_pyc_after_fit_woNu,"pyc_after_fit_woNu/F") ;
	m_pTTree_2->Branch("pzc_after_fit_woNu",&m_pzc_after_fit_woNu,"pzc_after_fit_woNu/F") ;
	m_pTTree_2->Branch("ec_after_fit_woNu",&m_ec_after_fit_woNu,"ec_after_fit_woNu/F") ;
	m_pTTree_2->Branch("pxc_before_ISR_best",&m_pxc_before_ISR_best,"pxc_before_ISR_best/F") ;
	m_pTTree_2->Branch("pyc_before_ISR_best",&m_pyc_before_ISR_best,"pyc_before_ISR_best/F") ;
	m_pTTree_2->Branch("pzc_before_ISR_best",&m_pzc_before_ISR_best,"pzc_before_ISR_best/F") ;
	m_pTTree_2->Branch("ec_before_ISR_best",&m_ec_before_ISR_best,"ec_before_ISR_best/F") ;
	m_pTTree_2->Branch("pxc_before_fit_best",&m_pxc_before_fit_best,"pxc_before_fit_best/F") ;
	m_pTTree_2->Branch("pyc_before_fit_best",&m_pyc_before_fit_best,"pyc_before_fit_best/F") ;
	m_pTTree_2->Branch("pzc_before_fit_best",&m_pzc_before_fit_best,"pzc_before_fit_best/F") ;
	m_pTTree_2->Branch("ec_before_fit_best",&m_ec_before_fit_best,"ec_before_fit_best/F") ;
	m_pTTree_2->Branch("pxc_after_fit_best",&m_pxc_after_fit_best,"pxc_after_fit_best/F") ;
	m_pTTree_2->Branch("pyc_after_fit_best",&m_pyc_after_fit_best,"pyc_after_fit_best/F") ;
	m_pTTree_2->Branch("pzc_after_fit_best",&m_pzc_after_fit_best,"pzc_after_fit_best/F") ;
	m_pTTree_2->Branch("ec_after_fit_best",&m_ec_after_fit_best,"ec_after_fit_best/F") ;
	m_pTTree_3 = new TTree("StartFitObjects","StartFitObjects");
	m_pTTree_3->SetDirectory(m_pTFile);
	m_pTTree_3->Branch("fitError",&m_iError_wNu) ;
	m_pTTree_3->Branch("jet_startPx",&m_jet_startPx);
	m_pTTree_3->Branch("jet_startPy",&m_jet_startPy);
	m_pTTree_3->Branch("jet_startPz",&m_jet_startPz);
	m_pTTree_3->Branch("jet_startE",&m_jet_startE);
	m_pTTree_3->Branch("lepton_startPx",&m_lepton_startPx);
	m_pTTree_3->Branch("lepton_startPy",&m_lepton_startPy);
	m_pTTree_3->Branch("lepton_startPz",&m_lepton_startPz);
	m_pTTree_3->Branch("lepton_startE",&m_lepton_startE);
	m_pTTree_3->Branch("jet_SigmaTheta",&m_jet_SigmaTheta);
	m_pTTree_3->Branch("jet_SigmaPhi",&m_jet_SigmaPhi);
	m_pTTree_3->Branch("jet_SigmaE",&m_jet_SigmaE);
	m_pTTree_3->Branch("lepton_SigmaTheta",&m_lepton_SigmaTheta);
	m_pTTree_3->Branch("lepton_SigmaPhi",&m_lepton_SigmaPhi);
	m_pTTree_3->Branch("lepton_SigmaInvpT",&m_lepton_SigmaInvpT);

	h_Zmass_beforefit_woNu = new TH1F("h_Zmass_beforefit_woNu", "Z mass before fit without #nu correction", 400, 0., 200.);
	h_Zmass_beforefit_woNu->SetDirectory(m_pTFile);
	h_Hmass_beforefit_woNu = new TH1F("h_Hmass_beforefit_woNu", "H mass before fit without #nu correction", 400, 0., 200.);
	h_Hmass_beforefit_woNu->SetDirectory(m_pTFile);
	h_Zmass_beforefit_wNu = new TH1F("h_Zmass_beforefit_wNu", "Z mass before fit with #nu correction", 400, 0., 200.);
	h_Zmass_beforefit_wNu->SetDirectory(m_pTFile);
	h_Hmass_beforefit_wNu = new TH1F("h_Hmass_beforefit_wNu", "H mass before fit with #nu correction", 400, 0., 200.);
	h_Hmass_beforefit_wNu->SetDirectory(m_pTFile);
	h_Zmass_afterfit_woNu = new TH1F("h_Zmass_afterfit_woNu", "Z mass after fit without #nu correction", 400, 0., 200.);
	h_Zmass_afterfit_woNu->SetDirectory(m_pTFile);
	h_Hmass_afterfit_woNu = new TH1F("h_Hmass_afterfit_woNu", "H mass after fit without #nu correction", 400, 0., 200.);
	h_Hmass_afterfit_woNu->SetDirectory(m_pTFile);
	h_Zmass_afterfit_wNu = new TH1F("h_Zmass_afterfit_wNu", "Z mass after fit with #nu correction", 400, 0., 200.);
	h_Zmass_afterfit_wNu->SetDirectory(m_pTFile);
	h_Hmass_afterfit_wNu = new TH1F("h_Hmass_afterfit_wNu", "H mass after fit with #nu correction", 400, 0., 200.);
	h_Hmass_afterfit_wNu->SetDirectory(m_pTFile);
	h_fitError_wNu = new TH1I("h_fitError_wNu", "fit error with #nu correction", 5, -1.5, 3.5);
	h_fitError_wNu->SetDirectory(m_pTFile);
	h_fitError_woNu = new TH1I("h_fitError_woNu", "fit error without #nu correction", 5, -1.5, 3.5);
	h_fitError_woNu->SetDirectory(m_pTFile);
	h_ErrorCode_wNu_woNu = new TH2I("h_ErrorCode_wNu_woNu", "Error Code; Error Code_{with #nu}; Error Code_{without #nu}", 5, -1.5, 3.5, 5, -1.5, 3.5);
	h_ErrorCode_wNu_woNu->SetDirectory(m_pTFile);
	h_fitProbability_wNu = new TH1F("h_fitProbability_wNu", "fit probability with #nu correction", 100, 0., 1.);
	h_fitProbability_wNu->SetDirectory(m_pTFile);
	h_fitProbability_woNu = new TH1F("h_fitProbability_woNu", "fit probability without #nu correction", 100, 0., 1.);
	h_fitProbability_woNu->SetDirectory(m_pTFile);
	h_fitProbability_best = new TH1F("h_fitProbability_best", "best fit probability", 100, 0., 1.);
	h_fitProbability_best->SetDirectory(m_pTFile);	
	h_nJets = new TH1I("h_nJets", "number of found jets per event; n_{Jets}; n_{events}", 8, -0.5, 7.5);
	h_nJets->SetDirectory(m_pTFile);
	h_nLeptons = new TH1I("h_nLeptons", "number of found isolated leptons per event; n_{ISOLeptons}; n_{events}", 8, -0.5, 7.5);
	h_nLeptons->SetDirectory(m_pTFile);
	h_nLeptons_nJets = new TH2I("h_nLeptons_nJets", "nISOleptons vs nJets per event; n_{ISOLeptons}; n_{Jets}", 8, -0.5, 7.5, 8, -0.5, 7.5);
	h_nLeptons_nJets->SetDirectory(m_pTFile);
	h_ISR_mcp_fit = new TH2F("h_ISR_mcp_fit", "ISR energy MCP vs FIT; ISR_{MCP}; ISR_{FIT}", 100, 0., m_isrpzmax * 1.5, 100, 0., m_isrpzmax * 1.5);
	h_ISR_mcp_fit->SetDirectory(m_pTFile);
	h_ISR_pzc_wNu = new TH2F("h_ISR_pzc_wNu", "ISR momentum vs p_{z} constraint with #nu correction; #Sigma p_{z} [GeV]; p_{z}^{ISR} [GeV]", 80, -40., 40., 80, -40., 40.);
	h_ISR_pzc_wNu->SetDirectory(m_pTFile);
	h_ISR_pzc_woNu = new TH2F("h_ISR_pzc_woNu", "ISR momentum vs p_{z} constraint without #nu correction; #Sigma p_{z} [GeV]; p_{z}^{ISR} [GeV]", 80, -40., 40., 80, -40., 40.);
	h_ISR_pzc_woNu->SetDirectory(m_pTFile);
	h_ISR_pzc_best = new TH2F("h_ISR_pzc_best", "ISR momentum vs p_{z} constraint (best fit); #Sigma p_{z} [GeV]; p_{z}^{ISR} [GeV]", 80, -40., 40., 80, -40., 40.);
	h_ISR_pzc_best->SetDirectory(m_pTFile);
	h_pull_jet_E = new TH1F("h_pull_jet_E", "pull of E for jets after successful fit; pull E [GeV]; n_{jets}", 100, -5., 5.);
	h_pull_jet_E->SetDirectory(m_pTFile);
	h_pull_jet_theta = new TH1F("h_pull_jet_theta", "pull of #theta for jets after successful fit; pull #theta; n_{jets}", 100, -5., 5.);
	h_pull_jet_theta->SetDirectory(m_pTFile);
	h_pull_jet_phi = new TH1F("h_pull_jet_phi", "pull of #phi for jets after successful fit; pull #phi; n_{jets}", 100, -5., 5.);
	h_pull_jet_phi->SetDirectory(m_pTFile);
	h_pull_lepton_InvPt = new TH1F("h_pull_lepton_InvPt", "pull of #frac{1}{p_{T}} for leptons after successful fit; pull #frac{1}{p_{T}}; n_{leptons}", 100, -5., 5.);
	h_pull_lepton_InvPt->SetDirectory(m_pTFile);
	h_pull_lepton_theta = new TH1F("h_pull_lepton_theta", "pull of #theta for leptons after successful fit; pull #theta; n_{leptons}", 100, -5., 5.);
	h_pull_lepton_theta->SetDirectory(m_pTFile);
	h_pull_lepton_phi = new TH1F("h_pull_lepton_phi", "pull of #phi for leptons after successful fit; pull #phi; n_{leptons}", 100, -5., 5.);
	h_pull_lepton_phi->SetDirectory(m_pTFile);
}

void ZHllqqFit::Clear()
{
	m_nSLDecayBHadron = 0;
	m_nSLDecayCHadron = 0;
	m_nSLDecayTotal = 0;
	m_nJets = 0;
	m_nLeptons = 0;
	bestprob = 0.;
	m_bestNuCombination.clear();
	m_iError_wNu.clear();
	m_probability_wNu.clear();
	m_chi2_wNu.clear();
	m_n_itter_wNu.clear();
	m_startmassZ_wNu.clear();
	m_startmassH_wNu.clear();
	m_beststartmassZ_wNu.clear();
	m_beststartmassH_wNu.clear();
	m_Zmass_after_fit_wNu.clear();
	m_Hmass_after_fit_wNu.clear();
	m_bestphotonenergy_wNu.clear();
	m_chi2startmassZ_wNu.clear();
	m_chi2startmassH_wNu.clear();
	m_jet_startPx.clear();
	m_jet_startPy.clear();
	m_jet_startPz.clear();
	m_jet_startE.clear();
	m_jet_SigmaTheta.clear();
	m_jet_SigmaPhi.clear();
	m_jet_SigmaE.clear();
	m_lepton_startPx.clear();
	m_lepton_startPy.clear();
	m_lepton_startPz.clear();
	m_lepton_startE.clear();
	m_lepton_SigmaTheta.clear();
	m_lepton_SigmaPhi.clear();
	m_lepton_SigmaInvpT.clear();
	m_pull_jet_E_wNu.clear();
	m_pull_jet_th_wNu.clear();
	m_pull_jet_phi_wNu.clear();
	m_pull_lepton_InvpT_wNu.clear();
	m_pull_lepton_th_wNu.clear();
	m_pull_lepton_phi_wNu.clear();
	m_iError_wNu_bestfit = -5;
	m_probability_wNu_bestfit = 0.;
	m_chi2_wNu_bestfit = 0.;
	m_n_itter_wNu_bestfit = 0;
	m_startmassZ_wNu_bestfit = 0.;
	m_startmassH_wNu_bestfit = 0.;
	m_beststartmassZ_wNu_bestfit = 0.;
	m_beststartmassH_wNu_bestfit = 0.;
	m_Zmass_after_fit_wNu_bestfit = 0.;
	m_Hmass_after_fit_wNu_bestfit = 0.;
	m_bestphotonenergy_wNu_bestfit = 0.;
	m_chi2startmassZ_wNu_bestfit = 0.;
	m_chi2startmassH_wNu_bestfit = 0.;
	m_pull_jet_E_wNu_bestfit.clear();
	m_pull_jet_th_wNu_bestfit.clear();
	m_pull_jet_phi_wNu_bestfit.clear();
	m_pull_lepton_InvpT_wNu_bestfit.clear();
	m_pull_lepton_th_wNu_bestfit.clear();
	m_pull_lepton_phi_wNu_bestfit.clear();
	m_iError_woNu = -5;
	m_probability_woNu = 0.;
	m_chi2best_woNu = 0.;
	m_n_itter_woNu = 0;
	m_startmassZ_woNu = 0.;
	m_startmassH_woNu = 0.;
	m_beststartmassZ_woNu = 0.;
	m_beststartmassH_woNu = 0.;
	m_Zmass_after_fit_woNu = 0.;
	m_Hmass_after_fit_woNu = 0.;
	m_bestphotonenergy_woNu = 0.;
	m_chi2startmassZ_woNu = 0.;
	m_chi2startmassH_woNu = 0.;
	m_pull_jet_E_woNu.clear();
	m_pull_jet_th_woNu.clear();
	m_pull_jet_phi_woNu.clear();
	m_pull_lepton_InvpT_woNu.clear();
	m_pull_lepton_th_woNu.clear();
	m_pull_lepton_phi_woNu.clear();
	m_iError_best = -5;
	m_probability_best = 0.;
	m_chi2_best = 0.;
	m_n_itter_best = 0;
	m_startmassZ_best = 0.;
	m_startmassH_best = 0.;
	m_beststartmassZ_best = 0.;
	m_beststartmassH_best = 0.;
	m_Zmass_after_fit_best = 0.;
	m_Hmass_after_fit_best = 0.;
	m_bestphotonenergy_best = 0.;
	m_chi2startmassZ_best = 0.;
	m_chi2startmassH_best = 0.;
	m_pull_jet_E_best.clear();
	m_pull_jet_th_best.clear();
	m_pull_jet_phi_best.clear();
	m_pull_lepton_InvpT_best.clear();
	m_pull_lepton_th_best.clear();
	m_pull_lepton_phi_best.clear();
	m_pxc_before_ISR_wNu = 0.;
	m_pyc_before_ISR_wNu = 0.;
	m_pzc_before_ISR_wNu = 0.;
	m_ec_before_ISR_wNu = 0.;
	m_pxc_before_fit_wNu = 0.;
	m_pyc_before_fit_wNu = 0.;
	m_pzc_before_fit_wNu = 0.;
	m_ec_before_fit_wNu = 0.;
	m_pxc_after_fit_wNu = 0.;
	m_pyc_after_fit_wNu = 0.;
	m_pzc_after_fit_wNu = 0.;
	m_ec_after_fit_wNu = 0.;
	m_pxc_before_ISR_woNu = 0.;
	m_pyc_before_ISR_woNu = 0.;
	m_pzc_before_ISR_woNu = 0.;
	m_ec_before_ISR_woNu = 0.;
	m_pxc_before_fit_woNu = 0.;
	m_pyc_before_fit_woNu = 0.;
	m_pzc_before_fit_woNu = 0.;
	m_ec_before_fit_woNu = 0.;
	m_pxc_after_fit_woNu = 0.;
	m_pyc_after_fit_woNu = 0.;
	m_pzc_after_fit_woNu = 0.;
	m_ec_after_fit_woNu = 0.;
	m_pxc_before_ISR_best = 0.;
	m_pyc_before_ISR_best = 0.;
	m_pzc_before_ISR_best = 0.;
	m_ec_before_ISR_best = 0.;
	m_pxc_before_fit_best = 0.;
	m_pyc_before_fit_best = 0.;
	m_pzc_before_fit_best = 0.;
	m_ec_before_fit_best = 0.;
	m_pxc_after_fit_best = 0.;
	m_pyc_after_fit_best = 0.;
	m_pzc_after_fit_best = 0.;
	m_ec_after_fit_best = 0.;
	
}

void ZHllqqFit::processRunHeader()
{
	_nRun++ ;
}

void ZHllqqFit::compcorrect() //finds the jet with cross checking
{
	if( std::abs(delta_theta[bestjet_phi]) < std::abs(delta_phi[bestjet_th]) )
	{
		bestjet=bestjet_phi;
	}
	if( std::abs(delta_theta[bestjet_phi]) > std::abs(delta_phi[bestjet_th]) )
	{
		bestjet=bestjet_th;
	}
}

void ZHllqqFit::Setvalues()
{

	Px=0.;
	Px2=0.;
	Py=0.;
	Py2=0.;
	Pz=0.;
	Pz2=0.;
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
	beststartmassZ = 0., beststartmassH = 0.;
	startmassZ = 0., startmassH = 0.;
	bestphotonenergy = 0.;
 	besterr = 999;
 	bestzvalue = 10000.;
 	chi2startmassZ = 0.;
 	chi2startmassH = 0.;
	memset(Zmomentum, 0, sizeof(Zmomentum));
	memset(Hmomentum, 0, sizeof(Hmomentum));
	memset(ISRmomentum, 0, sizeof(ISRmomentum));
	Z_Energy=0.;
	H_Energy=0.;
	chi2best=0.;
	errorcode=0.;
	streamlog_out(DEBUG)  << "Values set to defaults" <<std::endl;

}

void ZHllqqFit::processEvent( EVENT::LCEvent *pLCEvent )
{


	this->Clear();
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	streamlog_out(DEBUG) << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(DEBUG) << "///////////////////////////////////////////////// processing event " << m_nEvt << " in run " << m_nRun << " /////////////////////////////////////////////////////////////////////" << std::endl ;
	streamlog_out(DEBUG) << "/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;

	LCCollection *inputJetCollection{};
	LCCollection *HiggsDecayMode{};
	LCCollection *inputLeptonCollection{};
	LCCollection *inputSLDecayCollection{};
	LCCollection *inputNuEnergyB{};
	LCCollection *inputNuEnergyC{};
	std::vector<int> B_index{};
	std::vector<int> C_index{};
	std::vector<float> BHadENuPlus{};
	std::vector<float> BHadENuMinus{};
	std::vector<float> CHadENuPlus{};
	std::vector<float> CHadENuMinus{};
	int m_isDecayedTob = 0;
	int m_isDecayedToc = 0;
	int m_isDecayedToother = 0;
	float ISR_mcp = 0.;
	try
	{
		inputJetCollection = pLCEvent->getCollection( jetcollection );
		HiggsDecayMode = pLCEvent->getCollection( hDecayMode );
		m_isDecayedTob = HiggsDecayMode->getParameters().getIntVal("isDecayedTob");
		m_isDecayedToc = HiggsDecayMode->getParameters().getIntVal("isDecayedToc");
		m_isDecayedToother = HiggsDecayMode->getParameters().getIntVal("isDecayedToother");
		ISR_mcp = HiggsDecayMode->getParameters().getFloatVal("ISREnergy");
		if ( m_isDecayedTob == 1 )
		{
			if ( !m_includeHbb ) return;
		}
		if ( m_isDecayedToc == 1 )
		{
			if ( !m_includeHcc ) return;
		}
		if ( m_isDecayedToother == 1 )
		{
			if ( !m_includeHother ) return;
		}
		inputLeptonCollection = pLCEvent->getCollection( leptoncollection );
		if (m_fitNuE)
		{
			inputSLDecayCollection = pLCEvent->getCollection( SLDecayCollection );
			inputNuEnergyB = pLCEvent->getCollection( NuEnergyCollectionB );
			inputNuEnergyC = pLCEvent->getCollection( NuEnergyCollectionC );
			m_nSLDecayBHadron = inputSLDecayCollection->getParameters().getIntVal("nBSLD");
			m_nSLDecayCHadron = inputSLDecayCollection->getParameters().getIntVal("nCSLD");
			m_nSLDecayTotal = inputSLDecayCollection->getParameters().getIntVal("nSLD");
			B_index = inputSLDecayCollection->getParameters().getIntVals("BHadronIndex", B_index);
			C_index = inputSLDecayCollection->getParameters().getIntVals("CHadronIndex", C_index);
			BHadENuPlus = inputNuEnergyB->getParameters().getFloatVals("recEnergyENuPlusSLDB", BHadENuPlus);
			BHadENuMinus = inputNuEnergyB->getParameters().getFloatVals("recEnergyENuMinusSLDB", BHadENuMinus);
			CHadENuPlus = inputNuEnergyC->getParameters().getFloatVals("recEnergyENuPlusSLDC", CHadENuPlus);
			CHadENuMinus = inputNuEnergyC->getParameters().getFloatVals("recEnergyENuMinusSLDC", CHadENuMinus);
			streamlog_out(DEBUG)  << "Found " << m_nSLDecayBHadron << " semileptonic decay of B-hadron in event " << m_nEvt << std::endl;
			streamlog_out(DEBUG)  << "Found " << m_nSLDecayCHadron << " semileptonic decay of C-hadron in event " << m_nEvt << std::endl;
		}
		m_nJets = inputJetCollection->getNumberOfElements();
		m_nLeptons = inputLeptonCollection->getNumberOfElements();
		streamlog_out(DEBUG) << " found " << m_nJets << " jets in event: " << m_nEvt << " in run: " << m_nRun << std::endl;
		streamlog_out(DEBUG) << " found " << m_nLeptons << " isolated leptons in event: " << m_nEvt << " in run: " << m_nRun << std::endl;
		h_nJets->Fill(m_nJets);
		h_nLeptons->Fill(m_nLeptons);
		h_nLeptons_nJets->Fill( m_nLeptons , m_nJets );
		m_pTTree->Fill();
		if (m_nJets != 2 || m_nLeptons != 2) return;
		float yminus = inputJetCollection->parameters().getFloatVal( "YMinus" );
		streamlog_out(DEBUG)  << " Yminus = " << yminus << std::endl;
		float yplus = inputJetCollection->parameters().getFloatVal( "YPlus" );
		streamlog_out(DEBUG)  << " Yplus = " << yplus << std::endl;
		
//		const int m_nSLDB = m_nSLDecayBHadron;
//		const int m_nSLDC = m_nSLDecayCHadron;
		std::vector<int> JetmatchSLDB{};//const int m_nSLDB = m_nSLDecayBHadron];
		std::vector<int> JetmatchSLDC{};//const int m_nSLDC = m_nSLDecayCHadron];
		for (unsigned int i_SLDB = 0; i_SLDB < B_index.size(); i_SLDB++) JetmatchSLDB.push_back(this->FindMatchingJettoSLD(pLCEvent, B_index[i_SLDB]));
		for (unsigned int i_SLDC = 0; i_SLDC < C_index.size(); i_SLDC++) JetmatchSLDC.push_back(this->FindMatchingJettoSLD(pLCEvent, C_index[i_SLDC]));

		const int n_B_perm = pow( 3 , m_nSLDecayBHadron );
//		std::vector<int> E_nu_B_permut[n_B_perm]{};
		std::vector<std::vector<int>> E_nu_B_permut;
		std::vector<int> rowB;
		for (int E_nu_B_permutation = 0; E_nu_B_permutation < n_B_perm; E_nu_B_permutation++)
		{
			int s = E_nu_B_permutation;
			int m = 0;
			int r = 0;
			for (int i = 0; i < m_nSLDecayBHadron; i++)
			{
				m = s / 3;
				r = s - 3 * m;
				s = m;
				rowB.push_back(r);
//				E_nu_B_permut[E_nu_B_permutation].push_back(r);
			}
			E_nu_B_permut.push_back(rowB);
			rowB.clear();
		}

		int n_C_perm = pow( 3 , m_nSLDecayCHadron );
//		std::vector<int> E_nu_C_permut[n_C_perm]{};
		std::vector<std::vector<int>> E_nu_C_permut;
		std::vector<int> rowC;
		for (int E_nu_C_permutation = 0; E_nu_C_permutation < n_C_perm; E_nu_C_permutation++)
		{
			int s = E_nu_C_permutation;
			int m = 0;
			int r = 0;
			for (int i = 0; i < m_nSLDecayCHadron; i++)
			{
				m = s / 3;
				r = s - 3 * m;
				s = m;
				rowC.push_back(r);
//				E_nu_C_permut[E_nu_C_permutation].push_back(r);
			}
			E_nu_C_permut.push_back(rowC);
			rowC.clear();
		}
		float ENeutrinoJet0_BSLD = 0.;
		float ENeutrinoJet1_BSLD = 0.;
		float ENeutrinoJet0_CSLD = 0.;
		float ENeutrinoJet1_CSLD = 0.;
//		float ENeutrinoJet0 = 0.;
//		float ENeutrinoJet1 = 0.;
		float permuted_ENu_B = 0;
		float permuted_ENu_C = 0;
		TLorentzVector BNutlv(0,0,0,0);
		TLorentzVector Jet0_BNutlv(0,0,0,0);
		TLorentzVector Jet1_BNutlv(0,0,0,0);
		TLorentzVector CNutlv(0,0,0,0);
		TLorentzVector Jet0_CNutlv(0,0,0,0);
		TLorentzVector Jet1_CNutlv(0,0,0,0);
		TLorentzVector Jet0_Nutlv(0,0,0,0);
		TLorentzVector Jet1_Nutlv(0,0,0,0);
		std::vector<float> FitResultwNu{};
		std::vector<float> FitResultwoNu{};
		int ierr = 0.;
		float bestfitprob_wNu = 0.;
		float bestfitprob_woNu = 0.;

		std::vector<int> B_permutation{};
		std::vector<int> C_permutation{};
		int best_B_Nu = 0;
		int best_C_Nu = 0;
		float m_pull_jet1_E_wNu_bestfit = 0.;
		float m_pull_jet2_E_wNu_bestfit = 0.;
		float m_pull_jet1_th_wNu_bestfit = 0.;
		float m_pull_jet2_th_wNu_bestfit = 0.;
		float m_pull_jet1_phi_wNu_bestfit = 0.;
		float m_pull_jet2_phi_wNu_bestfit = 0.;
		float m_pull_lepton1_InvpT_wNu_bestfit = 0.;
		float m_pull_lepton2_InvpT_wNu_bestfit = 0.;
		float m_pull_lepton1_th_wNu_bestfit = 0.;
		float m_pull_lepton2_th_wNu_bestfit = 0.;
		float m_pull_lepton1_phi_wNu_bestfit = 0.;
		float m_pull_lepton2_phi_wNu_bestfit = 0.;
		float ISRpx_wNu_bestfit = 0.;
		float ISRpy_wNu_bestfit = 0.;
		float ISRpz_wNu_bestfit = 0.;
		float Zpx_wNu_bestfit = 0.;
		float Zpy_wNu_bestfit = 0.;
		float Zpz_wNu_bestfit = 0.;
		float ZE_wNu_bestfit = 0.;
		float Hpx_wNu_bestfit = 0.;
		float Hpy_wNu_bestfit = 0.;
		float Hpz_wNu_bestfit = 0.;
		float HE_wNu_bestfit = 0.;
		float ISRpx_woNu = 0.;
		float ISRpy_woNu = 0.;
		float ISRpz_woNu = 0.;
		float Zpx_woNu = 0.;
		float Zpy_woNu = 0.;
		float Zpz_woNu = 0.;
		float ZE_woNu = 0.;
		float Hpx_woNu = 0.;
		float Hpy_woNu = 0.;
		float Hpz_woNu = 0.;
		float HE_woNu = 0.;
		float ISRpx_best = 0.;
		float ISRpy_best = 0.;
		float ISRpz_best = 0.;
		float Zpx_best = 0.;
		float Zpy_best = 0.;
		float Zpz_best = 0.;
		float ZE_best = 0.;
		float Hpx_best = 0.;
		float Hpy_best = 0.;
		float Hpz_best = 0.;
		float HE_best = 0.;

		float jet0_startPx = 0.;
		float jet1_startPx = 0.;
		float jet0_startPy = 0.;
		float jet1_startPy = 0.;
		float jet0_startPz = 0.;
		float jet1_startPz = 0.;
		float jet0_startE = 0.;
		float jet1_startE = 0.;
		float lepton0_startPx = 0.;
		float lepton1_startPx = 0.;
		float lepton0_startPy = 0.;
		float lepton1_startPy = 0.;
		float lepton0_startPz = 0.;
		float lepton1_startPz = 0.;
		float lepton0_startE = 0.;
		float lepton1_startE = 0.;
		float jet0_SigmaTheta = 0.;
		float jet1_SigmaTheta = 0.;
		float jet0_SigmaPhi = 0.;
		float jet1_SigmaPhi = 0.;
		float jet0_SigmaE = 0.;
		float jet1_SigmaE = 0.;
		float lepton0_SigmaTheta = 0.;
		float lepton1_SigmaTheta = 0.;
		float lepton0_SigmaPhi = 0.;
		float lepton1_SigmaPhi = 0.;
		float lepton0_SigmaInvpT = 0.;
		float lepton1_SigmaInvpT = 0.;

		for (int i_perm_B = 0; i_perm_B < n_B_perm; i_perm_B++)
		{
			ENeutrinoJet0_BSLD = 0.;
			ENeutrinoJet1_BSLD = 0.;
			permuted_ENu_B = 0;
			BNutlv = TLorentzVector(0.,0.,0.,0.);
			for (unsigned int i_sldB = 0; i_sldB < E_nu_B_permut[i_perm_B].size(); i_sldB++)
			{
				EVENT::MCParticle *BNuPlus = dynamic_cast<EVENT::MCParticle*>(inputNuEnergyB->getElementAt( 2 * i_sldB ));
				EVENT::MCParticle *BNuMinus = dynamic_cast<EVENT::MCParticle*>(inputNuEnergyB->getElementAt( 2 * i_sldB + 1 ));
				if ( E_nu_B_permut[i_perm_B][i_sldB] == 0 )
				{
					permuted_ENu_B = BHadENuMinus[i_sldB];
					BNutlv = TLorentzVector(BNuMinus->getMomentum(),BNuMinus->getEnergy());
				}
				else if ( E_nu_B_permut[i_perm_B][i_sldB] == 1 )
				{
					permuted_ENu_B = BHadENuPlus[i_sldB];
					BNutlv = TLorentzVector(BNuPlus->getMomentum(),BNuPlus->getEnergy());
				}
				else
				{
					permuted_ENu_B = 0.;
					BNutlv = TLorentzVector(0.,0.,0.,0.);
				}
				if ( JetmatchSLDB[i_sldB] == 0 )
				{
					ENeutrinoJet0_BSLD += permuted_ENu_B;
					Jet0_BNutlv += BNutlv;
				}
				else if ( JetmatchSLDB[i_sldB] == 1 )
				{
					ENeutrinoJet1_BSLD += permuted_ENu_B;
					Jet1_BNutlv += BNutlv;
				}
			}
			for (int i_perm_C = 0; i_perm_C < n_C_perm; i_perm_C++)
			{
				ENeutrinoJet0_CSLD = 0;
				ENeutrinoJet1_CSLD = 0;
				permuted_ENu_C = 0;
				CNutlv = TLorentzVector(0.,0.,0.,0.);
				for (unsigned int i_sldC = 0; i_sldC < E_nu_C_permut[i_perm_C].size(); i_sldC++)
				{
					EVENT::MCParticle *CNuPlus = dynamic_cast<EVENT::MCParticle*>(inputNuEnergyC->getElementAt( 2 * i_sldC ));
					EVENT::MCParticle *CNuMinus = dynamic_cast<EVENT::MCParticle*>(inputNuEnergyC->getElementAt( 2 * i_sldC + 1 ));
					if ( E_nu_C_permut[i_perm_C][i_sldC] == 0 )
					{
						permuted_ENu_C = CHadENuMinus[i_sldC];
						CNutlv = TLorentzVector(CNuMinus->getMomentum(),CNuMinus->getEnergy());
					}
					else if ( E_nu_C_permut[i_perm_C][i_sldC] == 1 )
					{
						permuted_ENu_C = CHadENuPlus[i_sldC];
						CNutlv = TLorentzVector(CNuPlus->getMomentum(),CNuPlus->getEnergy());
					}
					else
					{
						permuted_ENu_C = 0.;
					}
					if ( JetmatchSLDC[i_sldC] == 0 )
					{
						ENeutrinoJet0_CSLD += permuted_ENu_C;
						Jet0_CNutlv += CNutlv;
					}
					else if ( JetmatchSLDC[i_sldC] == 1 )
					{
						ENeutrinoJet1_CSLD += permuted_ENu_C;
						Jet1_CNutlv += CNutlv;
					}
				}
//				ENeutrinoJet0 = ENeutrinoJet0_BSLD + ENeutrinoJet0_CSLD;
//				ENeutrinoJet1 = ENeutrinoJet1_BSLD + ENeutrinoJet1_CSLD;
				Jet0_Nutlv = Jet0_BNutlv + Jet0_CNutlv;
				Jet1_Nutlv = Jet1_BNutlv + Jet1_CNutlv;
//				if ( ENeutrinoJet0 == 0. && ENeutrinoJet1 == 0. ) continue;
				streamlog_out(DEBUG) << "**************************************************************************************************************************************************************************" << std::endl;
				streamlog_out(DEBUG) << "***************** Perform fit with Jet0_Nutlv = (" << Jet0_Nutlv.Px() << " , " << Jet0_Nutlv.Py() << " , " << Jet0_Nutlv.Pz() << " , " << Jet0_Nutlv.E() << ") and Jet1_Nutlv = (" << Jet1_Nutlv.Px() << " , " << Jet1_Nutlv.Py() << " , " << Jet1_Nutlv.Pz() << " , " << Jet1_Nutlv.E() << ") *****************" << std::endl;
				streamlog_out(DEBUG) << "**************************************************************************************************************************************************************************" << std::endl;
//				FitResultwNu = this->performFIT( pLCEvent , ENeutrinoJet0 , ENeutrinoJet1 );
				FitResultwNu = this->performFIT( pLCEvent , Jet0_Nutlv , Jet1_Nutlv );
				Jet0_BNutlv = TLorentzVector(0.,0.,0.,0.);
				Jet0_CNutlv = TLorentzVector(0.,0.,0.,0.);
				Jet1_BNutlv = TLorentzVector(0.,0.,0.,0.);
				Jet1_CNutlv = TLorentzVector(0.,0.,0.,0.);
				Jet0_Nutlv = TLorentzVector(0.,0.,0.,0.);
				Jet1_Nutlv = TLorentzVector(0.,0.,0.,0.);
				m_iError_wNu.push_back(FitResultwNu[0]);
				ierr = FitResultwNu[0];
				h_fitError_wNu->Fill(ierr);
				if ( ierr == 0 )//&& ( ENeutrinoJet0 != 0 || ENeutrinoJet1 !=0 ) )
				{
					m_probability_wNu.push_back(FitResultwNu[1]);
					m_n_itter_wNu.push_back(FitResultwNu[2]);
					m_startmassZ_wNu.push_back(FitResultwNu[3]);
					m_startmassH_wNu.push_back(FitResultwNu[4]);
					m_beststartmassZ_wNu.push_back(FitResultwNu[5]);
					m_beststartmassH_wNu.push_back(FitResultwNu[6]);
					m_Zmass_after_fit_wNu.push_back(FitResultwNu[7]);
					m_Hmass_after_fit_wNu.push_back(FitResultwNu[8]);
					m_chi2startmassZ_wNu.push_back(FitResultwNu[9]);
					m_chi2startmassH_wNu.push_back(FitResultwNu[10]);
					m_chi2_wNu.push_back(FitResultwNu[11]);
					m_bestphotonenergy_wNu.push_back(FitResultwNu[12]);
					m_pull_jet_E_wNu.push_back(FitResultwNu[13]);
					m_pull_jet_E_wNu.push_back(FitResultwNu[14]);
					m_pull_jet_th_wNu.push_back(FitResultwNu[15]);
					m_pull_jet_th_wNu.push_back(FitResultwNu[16]);
					m_pull_jet_phi_wNu.push_back(FitResultwNu[17]);
					m_pull_jet_phi_wNu.push_back(FitResultwNu[18]);
					m_pull_lepton_InvpT_wNu.push_back(FitResultwNu[19]);
					m_pull_lepton_InvpT_wNu.push_back(FitResultwNu[20]);
					m_pull_lepton_th_wNu.push_back(FitResultwNu[21]);
					m_pull_lepton_th_wNu.push_back(FitResultwNu[22]);
					m_pull_lepton_phi_wNu.push_back(FitResultwNu[23]);
					m_pull_lepton_phi_wNu.push_back(FitResultwNu[24]);
				}
				m_iError_wNu_bestfit = FitResultwNu[0];
				if ( m_iError_wNu_bestfit == 0 && bestfitprob_wNu <= FitResultwNu[1] )
				{
					best_B_Nu = i_perm_B;
					best_C_Nu = i_perm_C;
					m_probability_wNu_bestfit = FitResultwNu[1];
					bestfitprob_wNu = m_probability_wNu_bestfit;
					m_n_itter_wNu_bestfit = FitResultwNu[2];
					m_startmassZ_wNu_bestfit = FitResultwNu[3];
					m_startmassH_wNu_bestfit = FitResultwNu[4];
					m_beststartmassZ_wNu_bestfit = FitResultwNu[5];
					m_beststartmassH_wNu_bestfit = FitResultwNu[6];
					m_Zmass_after_fit_wNu_bestfit = FitResultwNu[7];
					m_Hmass_after_fit_wNu_bestfit = FitResultwNu[8];
					m_chi2startmassZ_wNu_bestfit = FitResultwNu[9];
					m_chi2startmassH_wNu_bestfit = FitResultwNu[10];
					m_chi2_wNu_bestfit = FitResultwNu[11];
					m_bestphotonenergy_wNu_bestfit = FitResultwNu[12];
					m_pull_jet1_E_wNu_bestfit = FitResultwNu[13];
					m_pull_jet2_E_wNu_bestfit = FitResultwNu[14];
					m_pull_jet1_th_wNu_bestfit = FitResultwNu[15];
					m_pull_jet2_th_wNu_bestfit = FitResultwNu[16];
					m_pull_jet1_phi_wNu_bestfit = FitResultwNu[17];
					m_pull_jet2_phi_wNu_bestfit = FitResultwNu[18];
					m_pull_lepton1_InvpT_wNu_bestfit = FitResultwNu[19];
					m_pull_lepton2_InvpT_wNu_bestfit = FitResultwNu[20];
					m_pull_lepton1_th_wNu_bestfit = FitResultwNu[21];
					m_pull_lepton2_th_wNu_bestfit = FitResultwNu[22];
					m_pull_lepton1_phi_wNu_bestfit = FitResultwNu[23];
					m_pull_lepton2_phi_wNu_bestfit = FitResultwNu[24];
					ISRpx_wNu_bestfit = FitResultwNu[25];
					ISRpy_wNu_bestfit = FitResultwNu[26];
					ISRpz_wNu_bestfit = FitResultwNu[27];
					Zpx_wNu_bestfit = FitResultwNu[28];
					Zpy_wNu_bestfit = FitResultwNu[29];
					Zpz_wNu_bestfit = FitResultwNu[30];
					ZE_wNu_bestfit = FitResultwNu[31];
					Hpx_wNu_bestfit = FitResultwNu[32];
					Hpy_wNu_bestfit = FitResultwNu[33];
					Hpz_wNu_bestfit = FitResultwNu[34];
					HE_wNu_bestfit = FitResultwNu[35];
					m_pxc_before_ISR_wNu = FitResultwNu[36];
					m_pyc_before_ISR_wNu = FitResultwNu[37];
					m_pzc_before_ISR_wNu = FitResultwNu[38];
					m_ec_before_ISR_wNu = FitResultwNu[39];
					m_pxc_before_fit_wNu = FitResultwNu[40];
					m_pyc_before_fit_wNu = FitResultwNu[41];
					m_pzc_before_fit_wNu = FitResultwNu[42];
					m_ec_before_fit_wNu = FitResultwNu[43];
					m_pxc_after_fit_wNu = FitResultwNu[44];
					m_pyc_after_fit_wNu = FitResultwNu[45];
					m_pzc_after_fit_wNu = FitResultwNu[46];
					m_ec_after_fit_wNu = FitResultwNu[47];
				}
				jet0_startPx = FitResultwNu[48];
				jet1_startPx = FitResultwNu[49];
				jet0_startPy = FitResultwNu[50];
				jet1_startPy = FitResultwNu[51];
				jet0_startPz = FitResultwNu[52];
				jet1_startPz = FitResultwNu[53];
				jet0_startE = FitResultwNu[54];
				jet1_startE = FitResultwNu[55];
				lepton0_startPx = FitResultwNu[56];
				lepton1_startPx = FitResultwNu[57];
				lepton0_startPy = FitResultwNu[58];
				lepton1_startPy = FitResultwNu[59];
				lepton0_startPz = FitResultwNu[60];
				lepton1_startPz = FitResultwNu[61];
				lepton0_startE = FitResultwNu[62];
				lepton1_startE = FitResultwNu[63];
				jet0_SigmaTheta = FitResultwNu[64];
				jet1_SigmaTheta = FitResultwNu[65];
				jet0_SigmaPhi = FitResultwNu[66];
				jet1_SigmaPhi = FitResultwNu[67];
				jet0_SigmaE = FitResultwNu[68];
				jet1_SigmaE = FitResultwNu[69];
				lepton0_SigmaTheta = FitResultwNu[70];
				lepton1_SigmaTheta = FitResultwNu[71];
				lepton0_SigmaPhi = FitResultwNu[72];
				lepton1_SigmaPhi = FitResultwNu[73];
				lepton0_SigmaInvpT = FitResultwNu[74];
				lepton1_SigmaInvpT = FitResultwNu[75];
				m_jet_startPx.push_back(jet0_startPx);
				m_jet_startPx.push_back(jet1_startPx);
				m_jet_startPy.push_back(jet0_startPy);
				m_jet_startPy.push_back(jet1_startPy);
				m_jet_startPz.push_back(jet0_startPz);
				m_jet_startPz.push_back(jet1_startPz);
				m_jet_startE.push_back(jet0_startE);
				m_jet_startE.push_back(jet1_startE);
				m_lepton_startPx.push_back(lepton0_startPx);
				m_lepton_startPx.push_back(lepton1_startPx);
				m_lepton_startPy.push_back(lepton0_startPy);
				m_lepton_startPy.push_back(lepton1_startPy);
				m_lepton_startPz.push_back(lepton0_startPz);
				m_lepton_startPz.push_back(lepton1_startPz);
				m_lepton_startE.push_back(lepton0_startE);
				m_lepton_startE.push_back(lepton1_startE);
				m_jet_SigmaTheta.push_back(jet0_SigmaTheta);
				m_jet_SigmaTheta.push_back(jet1_SigmaTheta);
				m_jet_SigmaPhi.push_back(jet0_SigmaPhi);
				m_jet_SigmaPhi.push_back(jet1_SigmaPhi);
				m_jet_SigmaE.push_back(jet0_SigmaE);
				m_jet_SigmaE.push_back(jet1_SigmaE);
				m_lepton_SigmaTheta.push_back(lepton0_SigmaTheta);
				m_lepton_SigmaTheta.push_back(lepton1_SigmaTheta);
				m_lepton_SigmaPhi.push_back(lepton0_SigmaPhi);
				m_lepton_SigmaPhi.push_back(lepton1_SigmaPhi);
				m_lepton_SigmaInvpT.push_back(lepton0_SigmaInvpT);
				m_lepton_SigmaInvpT.push_back(lepton1_SigmaInvpT);
				streamlog_out(DEBUG) << "size of FitResult with neutrino = " << FitResultwNu.size() << endl;
//				streamlog_out(DEBUG) << "Fit best probability with neutrino correction = " << FitResultwNu[1] << endl;
				FitResultwNu.clear();
			}
		}
		if ( m_iError_wNu_bestfit == 0 )
		{
			m_pull_jet_E_wNu_bestfit.push_back(m_pull_jet1_E_wNu_bestfit);
			m_pull_jet_E_wNu_bestfit.push_back(m_pull_jet2_E_wNu_bestfit);
			m_pull_jet_th_wNu_bestfit.push_back(m_pull_jet1_th_wNu_bestfit);
			m_pull_jet_th_wNu_bestfit.push_back(m_pull_jet2_th_wNu_bestfit);
			m_pull_jet_phi_wNu_bestfit.push_back(m_pull_jet1_phi_wNu_bestfit);
			m_pull_jet_phi_wNu_bestfit.push_back(m_pull_jet2_phi_wNu_bestfit);
			m_pull_lepton_InvpT_wNu_bestfit.push_back(m_pull_lepton1_InvpT_wNu_bestfit);
			m_pull_lepton_InvpT_wNu_bestfit.push_back(m_pull_lepton2_InvpT_wNu_bestfit);
			m_pull_lepton_th_wNu_bestfit.push_back(m_pull_lepton1_th_wNu_bestfit);
			m_pull_lepton_th_wNu_bestfit.push_back(m_pull_lepton2_th_wNu_bestfit);
			m_pull_lepton_phi_wNu_bestfit.push_back(m_pull_lepton1_phi_wNu_bestfit);
			m_pull_lepton_phi_wNu_bestfit.push_back(m_pull_lepton2_phi_wNu_bestfit);
		}
		std::vector<int> best_E_nu_permut{};
		int ss = best_C_Nu;
		int mm = 0;
		int rr = 0;
		for (int i = 0; i < m_nSLDecayCHadron; i++)
		{
			mm = ss / 3;
			rr = ss - 3 * mm;
			ss = mm;
			if (rr == 1)
			{
				best_E_nu_permut.push_back( i + 1 );
			}
			else if ( rr == 0 )
			{
				best_E_nu_permut.push_back( -1 - i );
			}
			else
			{
				best_E_nu_permut.push_back( 0 );
			}
		}
		ss = best_B_Nu;
		mm = 0;
		rr = 0;
		for (int i = 0; i < m_nSLDecayBHadron; i++)
		{
			mm = ss / 3;
			rr = ss - 3 * mm;
			ss = mm;
			if (rr == 1)
			{
				best_E_nu_permut.push_back( m_nSLDecayCHadron + i + 1 );
			}
			else if ( rr == 0 )
			{
				best_E_nu_permut.push_back( -1 - i - m_nSLDecayCHadron);
			}
			else
			{
				best_E_nu_permut.push_back( 0 );
			}
		}
		for (int i = 0; i < m_nSLDecayBHadron + m_nSLDecayCHadron ; i++)
		{
			m_bestNuCombination.push_back(best_E_nu_permut[m_nSLDecayBHadron + m_nSLDecayCHadron - i - 1]);
		}
		Jet0_Nutlv = TLorentzVector(0.,0.,0.,0.);
		Jet1_Nutlv = TLorentzVector(0.,0.,0.,0.);

		streamlog_out(DEBUG) << "***********************************************************************************************************************************************" << std::endl;
		streamlog_out(DEBUG) << "********************************************* Perform fit without neutrino correction *********************************************************" << std::endl;
		streamlog_out(DEBUG) << "***********************************************************************************************************************************************" << std::endl;
//		FitResultwoNu = this->performFIT( pLCEvent , 0. , 0. );
		FitResultwoNu = this->performFIT( pLCEvent , Jet0_Nutlv , Jet1_Nutlv );
		m_iError_woNu = FitResultwoNu[0];
		ierr = FitResultwoNu[0];
		h_fitError_woNu->Fill(m_iError_woNu);
		if ( ierr == 0 )
		{
			m_probability_woNu = FitResultwoNu[1];
			bestfitprob_woNu = m_probability_woNu;
			m_n_itter_woNu = FitResultwoNu[2];
			m_startmassZ_woNu = FitResultwoNu[3];
			m_startmassH_woNu = FitResultwoNu[4];
			m_beststartmassZ_woNu = FitResultwoNu[5];
			m_beststartmassH_woNu = FitResultwoNu[6];
			m_Zmass_after_fit_woNu = FitResultwoNu[7];
			m_Hmass_after_fit_woNu = FitResultwoNu[8];
			m_chi2startmassZ_woNu = FitResultwoNu[9];
			m_chi2startmassH_woNu = FitResultwoNu[10];
			m_chi2best_woNu = FitResultwoNu[11];
			m_bestphotonenergy_woNu = FitResultwoNu[12];
			m_pull_jet_E_woNu.push_back(FitResultwoNu[13]);
			m_pull_jet_E_woNu.push_back(FitResultwoNu[14]);
			m_pull_jet_th_woNu.push_back(FitResultwoNu[15]);
			m_pull_jet_th_woNu.push_back(FitResultwoNu[16]);
			m_pull_jet_phi_woNu.push_back(FitResultwoNu[17]);
			m_pull_jet_phi_woNu.push_back(FitResultwoNu[18]);
			m_pull_lepton_InvpT_woNu.push_back(FitResultwoNu[19]);
			m_pull_lepton_InvpT_woNu.push_back(FitResultwoNu[20]);
			m_pull_lepton_th_woNu.push_back(FitResultwoNu[21]);
			m_pull_lepton_th_woNu.push_back(FitResultwoNu[22]);
			m_pull_lepton_phi_woNu.push_back(FitResultwoNu[23]);
			m_pull_lepton_phi_woNu.push_back(FitResultwoNu[24]);
			ISRpx_woNu = FitResultwoNu[25];
			ISRpy_woNu = FitResultwoNu[26];
			ISRpz_woNu = FitResultwoNu[27];
			Zpx_woNu = FitResultwoNu[28];
			Zpy_woNu = FitResultwoNu[29];
			Zpz_woNu = FitResultwoNu[30];
			ZE_woNu = FitResultwoNu[31];
			Hpx_woNu = FitResultwoNu[32];
			Hpy_woNu = FitResultwoNu[33];
			Hpz_woNu = FitResultwoNu[34];
			HE_woNu = FitResultwoNu[35];
			m_pxc_before_ISR_woNu = FitResultwNu[36];
			m_pyc_before_ISR_woNu = FitResultwNu[37];
			m_pzc_before_ISR_woNu = FitResultwNu[38];
			m_ec_before_ISR_woNu = FitResultwNu[39];
			m_pxc_before_fit_woNu = FitResultwNu[40];
			m_pyc_before_fit_woNu = FitResultwNu[41];
			m_pzc_before_fit_woNu = FitResultwNu[42];
			m_ec_before_fit_woNu = FitResultwNu[43];
			m_pxc_after_fit_woNu = FitResultwNu[44];
			m_pyc_after_fit_woNu = FitResultwNu[45];
			m_pzc_after_fit_woNu = FitResultwNu[46];
			m_ec_after_fit_woNu = FitResultwNu[47];
		}
		streamlog_out(DEBUG) << "size of FitResult without neutrino = " << FitResultwoNu.size() << endl;
		streamlog_out(DEBUG) << "Fit probability without neutrino correction = " << FitResultwoNu[3] << endl;
		h_ErrorCode_wNu_woNu->Fill( m_iError_wNu_bestfit , m_iError_woNu );
//		m_iError_best = m_iError_wNu_bestfit;
		
		if ( m_iError_wNu_bestfit == 0 && bestfitprob_woNu < bestfitprob_wNu )
		{
			m_iError_best = m_iError_wNu_bestfit;
			m_probability_best = m_probability_wNu_bestfit;
			m_n_itter_best = m_n_itter_wNu_bestfit;
			m_startmassZ_best = m_startmassZ_wNu_bestfit;
			m_startmassH_best = m_startmassH_wNu_bestfit;
			m_beststartmassZ_best = m_beststartmassZ_wNu_bestfit;
			m_beststartmassH_best = m_beststartmassH_wNu_bestfit;
			m_Zmass_after_fit_best = m_Zmass_after_fit_wNu_bestfit;
			m_Hmass_after_fit_best = m_Hmass_after_fit_wNu_bestfit;
			m_chi2startmassZ_best = m_chi2startmassZ_wNu_bestfit;
			m_chi2startmassH_best = m_chi2startmassH_wNu_bestfit;
			m_chi2_best = m_chi2_wNu_bestfit;
			m_bestphotonenergy_best = m_bestphotonenergy_wNu_bestfit;
			m_pull_jet_E_best.push_back(m_pull_jet_E_wNu_bestfit[0]);
			m_pull_jet_E_best.push_back(m_pull_jet_E_wNu_bestfit[1]);
			h_pull_jet_E->Fill(m_pull_jet_E_wNu_bestfit[0]);
			h_pull_jet_E->Fill(m_pull_jet_E_wNu_bestfit[1]);

			m_pull_jet_th_best.push_back(m_pull_jet_th_wNu_bestfit[0]);
			m_pull_jet_th_best.push_back(m_pull_jet_th_wNu_bestfit[1]);
			h_pull_jet_theta->Fill(m_pull_jet_th_wNu_bestfit[0]);
			h_pull_jet_theta->Fill(m_pull_jet_th_wNu_bestfit[1]);

			m_pull_jet_phi_best.push_back(m_pull_jet_phi_wNu_bestfit[0]);
			m_pull_jet_phi_best.push_back(m_pull_jet_phi_wNu_bestfit[1]);
			h_pull_jet_phi->Fill(m_pull_jet_phi_wNu_bestfit[0]);
			h_pull_jet_phi->Fill(m_pull_jet_phi_wNu_bestfit[1]);

			m_pull_lepton_InvpT_best.push_back(m_pull_lepton_InvpT_wNu_bestfit[0]);
			m_pull_lepton_InvpT_best.push_back(m_pull_lepton_InvpT_wNu_bestfit[1]);
			h_pull_lepton_InvPt->Fill(m_pull_lepton_InvpT_wNu_bestfit[0]);
			h_pull_lepton_InvPt->Fill(m_pull_lepton_InvpT_wNu_bestfit[1]);
			
			m_pull_lepton_th_best.push_back(m_pull_lepton_th_wNu_bestfit[0]);
			m_pull_lepton_th_best.push_back(m_pull_lepton_th_wNu_bestfit[1]);
			h_pull_lepton_theta->Fill(m_pull_lepton_th_wNu_bestfit[0]);
			h_pull_lepton_theta->Fill(m_pull_lepton_th_wNu_bestfit[1]);
			
			m_pull_lepton_phi_best.push_back(m_pull_lepton_phi_wNu_bestfit[0]);
			m_pull_lepton_phi_best.push_back(m_pull_lepton_phi_wNu_bestfit[1]);
			h_pull_lepton_phi->Fill(m_pull_lepton_phi_wNu_bestfit[0]);
			h_pull_lepton_phi->Fill(m_pull_lepton_phi_wNu_bestfit[1]);

			ISRpx_best = ISRpx_wNu_bestfit;
			ISRpy_best = ISRpy_wNu_bestfit;
			ISRpz_best = ISRpz_wNu_bestfit;
			Zpx_best = Zpx_wNu_bestfit;
			Zpy_best = Zpy_wNu_bestfit;
			Zpz_best = Zpz_wNu_bestfit;
			ZE_best = ZE_wNu_bestfit;
			Hpx_best = Hpx_wNu_bestfit;
			Hpy_best = Hpy_wNu_bestfit;
			Hpz_best = Hpz_wNu_bestfit;
			HE_best = HE_wNu_bestfit;
			m_pxc_before_ISR_best = m_pxc_before_ISR_wNu;
			m_pyc_before_ISR_best = m_pyc_before_ISR_wNu;
			m_pzc_before_ISR_best = m_pzc_before_ISR_wNu;
			m_ec_before_ISR_best = m_ec_before_ISR_wNu;
			m_pxc_before_fit_best = m_pxc_before_fit_wNu;
			m_pyc_before_fit_best = m_pyc_before_fit_wNu;
			m_pzc_before_fit_best = m_pzc_before_fit_wNu;
			m_ec_before_fit_best = m_ec_before_fit_wNu;
			m_pxc_after_fit_best = m_pxc_after_fit_wNu;
			m_pyc_after_fit_best = m_pyc_after_fit_wNu;
			m_pzc_after_fit_best = m_pzc_after_fit_wNu;
			m_ec_after_fit_best = m_ec_after_fit_wNu;
		}
		else if ( m_iError_woNu == 0 )
		{
			m_iError_best = m_iError_woNu;
			m_probability_best = m_probability_woNu;
			m_n_itter_best = m_n_itter_woNu;
			m_startmassZ_best = m_startmassZ_woNu;
			m_startmassH_best = m_startmassH_woNu;
			m_beststartmassZ_best = m_beststartmassZ_woNu;
			m_beststartmassH_best = m_beststartmassH_woNu;
			m_Zmass_after_fit_best = m_Zmass_after_fit_woNu;
			m_Hmass_after_fit_best = m_Hmass_after_fit_woNu;
			m_chi2startmassZ_best = m_chi2startmassZ_woNu;
			m_chi2startmassH_best = m_chi2startmassH_woNu;
			m_chi2_best = m_chi2best_woNu;
			m_bestphotonenergy_best = m_bestphotonenergy_woNu;
			m_pull_jet_E_best.push_back(m_pull_jet_E_woNu[0]);
			m_pull_jet_E_best.push_back(m_pull_jet_E_woNu[1]);
			h_pull_jet_E->Fill(m_pull_jet_E_woNu[0]);
			h_pull_jet_E->Fill(m_pull_jet_E_woNu[1]);

			m_pull_jet_th_best.push_back(m_pull_jet_th_woNu[0]);
			m_pull_jet_th_best.push_back(m_pull_jet_th_woNu[1]);
			h_pull_jet_theta->Fill(m_pull_jet_th_woNu[0]);
			h_pull_jet_theta->Fill(m_pull_jet_th_woNu[1]);

			m_pull_jet_phi_best.push_back(m_pull_jet_phi_woNu[0]);
			m_pull_jet_phi_best.push_back(m_pull_jet_phi_woNu[1]);
			h_pull_jet_phi->Fill(m_pull_jet_phi_woNu[0]);
			h_pull_jet_phi->Fill(m_pull_jet_phi_woNu[1]);

			m_pull_lepton_InvpT_best.push_back(m_pull_lepton_InvpT_woNu[0]);
			m_pull_lepton_InvpT_best.push_back(m_pull_lepton_InvpT_woNu[1]);
			h_pull_lepton_InvPt->Fill(m_pull_lepton_InvpT_woNu[0]);
			h_pull_lepton_InvPt->Fill(m_pull_lepton_InvpT_woNu[1]);

			m_pull_lepton_th_best.push_back(m_pull_lepton_th_woNu[0]);
			m_pull_lepton_th_best.push_back(m_pull_lepton_th_woNu[1]);
			h_pull_lepton_theta->Fill(m_pull_lepton_th_woNu[0]);
			h_pull_lepton_theta->Fill(m_pull_lepton_th_woNu[1]);

			m_pull_lepton_phi_best.push_back(m_pull_lepton_phi_woNu[0]);
			m_pull_lepton_phi_best.push_back(m_pull_lepton_phi_woNu[1]);
			h_pull_lepton_phi->Fill(m_pull_lepton_phi_woNu[0]);
			h_pull_lepton_phi->Fill(m_pull_lepton_phi_woNu[1]);

			m_bestNuCombination.push_back(0);
			ISRpx_best = ISRpx_woNu;
			ISRpy_best = ISRpy_woNu;
			ISRpz_best = ISRpz_woNu;
			Zpx_best = Zpx_woNu;
			Zpy_best = Zpy_woNu;
			Zpz_best = Zpz_woNu;
			ZE_best = ZE_woNu;
			Hpx_best = Hpx_woNu;
			Hpy_best = Hpy_woNu;
			Hpz_best = Hpz_woNu;
			HE_best = HE_woNu;
			m_pxc_before_ISR_best = m_pxc_before_ISR_woNu;
			m_pyc_before_ISR_best = m_pyc_before_ISR_woNu;
			m_pzc_before_ISR_best = m_pzc_before_ISR_woNu;
			m_ec_before_ISR_best = m_ec_before_ISR_woNu;
			m_pxc_before_fit_best = m_pxc_before_fit_woNu;
			m_pyc_before_fit_best = m_pyc_before_fit_woNu;
			m_pzc_before_fit_best = m_pzc_before_fit_woNu;
			m_ec_before_fit_best = m_ec_before_fit_woNu;
			m_pxc_after_fit_best = m_pxc_after_fit_woNu;
			m_pyc_after_fit_best = m_pyc_after_fit_woNu;
			m_pzc_after_fit_best = m_pzc_after_fit_woNu;
			m_ec_after_fit_best = m_ec_after_fit_woNu;
		}
//		if ( m_iError_woNu == 0 && m_iError_wNu_bestfit == 0 )
		if ( m_iError_woNu == 0 )
		{
			h_Zmass_beforefit_woNu->Fill(m_startmassZ_woNu);
			h_Hmass_beforefit_woNu->Fill(m_startmassH_woNu);
			h_Zmass_afterfit_woNu->Fill(m_Zmass_after_fit_woNu);
			h_Hmass_afterfit_woNu->Fill(m_Hmass_after_fit_woNu);
			h_fitProbability_woNu->Fill(m_probability_woNu);
			h_ISR_pzc_woNu->Fill( -1 * m_pzc_before_ISR_woNu , ISRpz_woNu );
		}
		if ( m_iError_wNu_bestfit == 0 )
		{
			h_Zmass_beforefit_wNu->Fill(m_startmassZ_wNu_bestfit);
			h_Hmass_beforefit_wNu->Fill(m_startmassH_wNu_bestfit);
			h_fitProbability_wNu->Fill(m_probability_wNu_bestfit);
			h_ISR_pzc_wNu->Fill( -1 * m_pzc_before_ISR_wNu , ISRpz_wNu_bestfit );
		}
		if ( m_iError_best == 0 )
		{
			streamlog_out(DEBUG) << "Fit best probability = " << FitResultwoNu[3] << endl;
			h_Zmass_afterfit_wNu->Fill(m_Zmass_after_fit_best);
			h_Hmass_afterfit_wNu->Fill(m_Hmass_after_fit_best);
			h_fitProbability_best->Fill(m_probability_best);
			h_ISR_pzc_best->Fill( -1 * m_pzc_before_ISR_best , ISRpz_best );
			h_ISR_mcp_fit->Fill(ISR_mcp,m_bestphotonenergy_best);

			LCCollectionVec *OutputCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
			ReconstructedParticleImpl* ISRfitrec = new ReconstructedParticleImpl;
			ReconstructedParticleImpl* Zfitrec = new ReconstructedParticleImpl;
			ReconstructedParticleImpl* Hfitrec = new ReconstructedParticleImpl;

			ISRmomentum[0] = ISRpx_best;
			ISRmomentum[1] = ISRpy_best;
			ISRmomentum[2] = ISRpz_best;
			ISRfitrec->setMomentum(ISRmomentum);
			ISRfitrec->setEnergy(m_bestphotonenergy_best);
			ISRfitrec->setType (22);
			streamlog_out(DEBUG) << " Energy ISR:   " << ISRfitrec->getEnergy() << std::endl ;
			streamlog_out(DEBUG) << " IS ISR:   " << ISRfitrec->getType() << std::endl ;
			OutputCol->addElement(ISRfitrec);

			Zmomentum[0] = Zpx_best;
			Zmomentum[1] = Zpy_best;
			Zmomentum[2] = Zpz_best;
			Zfitrec->setMomentum(Zmomentum);
			Zfitrec->setEnergy(ZE_best);
			Zfitrec->setMass(m_Zmass_after_fit_best);
			Zfitrec->setType (23);
			streamlog_out(DEBUG) << "  Zmomentum :   " << Zfitrec->getMomentum()[0] << "," << Zfitrec->getMomentum()[1]<<","<< Zfitrec->getMomentum()[2] << std::endl ;
			streamlog_out(DEBUG) << " Energy Z:   "	<< Zfitrec->getEnergy() << std::endl ;
			streamlog_out(DEBUG) << " Mass Z:   " << Zfitrec->getMass() << std::endl ;
			streamlog_out(DEBUG) << " IS Z :   " << Zfitrec->getType() << std::endl ;
			OutputCol->addElement(Zfitrec);

			Hmomentum[0] = Hpx_best;
			Hmomentum[1] = Hpy_best;
			Hmomentum[2] = Hpz_best;
			Hfitrec->setMomentum(Hmomentum);
			Hfitrec->setEnergy(HE_best);
			Hfitrec->setMass(m_Hmass_after_fit_best);
			Hfitrec->setType (25);
			streamlog_out(DEBUG) << " Energy H:   "	<< Hfitrec->getEnergy() << std::endl ;
			streamlog_out(DEBUG) << " Mass H:   " << Hfitrec->getMass() << std::endl ;
			streamlog_out(DEBUG) << " IS H:   " << Hfitrec->getType() << std::endl ;
			OutputCol->addElement(Hfitrec);
			
			pLCEvent->addCollection( OutputCol, outputFitcollection.c_str() );


			OutputCol->parameters().setValue("bestchisq", (float)m_chi2_best);
			streamlog_out(DEBUG) << " chi2:   " << m_chi2_best << std::endl ;
			OutputCol->parameters().setValue("best_prob", (float)m_probability_best);
			streamlog_out(DEBUG) << " prob:   " << m_probability_best << std::endl ;
			OutputCol->parameters().setValue("error_code", (int)m_iError_best);
			streamlog_out(DEBUG) << "Error Code:   " << m_iError_best << std::endl ;
	
			m_pTTree_0->Fill();
			m_pTTree_1->Fill();
			m_pTTree_2->Fill();
		}
		m_pTTree_3->Fill();
		FitResultwoNu.clear();
		m_nEvtSum++;
		m_nEvt++ ;
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "Input collections not found in event " << m_nEvt << std::endl;
	}

}

int ZHllqqFit::FindMatchingJettoSLD(EVENT::LCEvent *pLCEvent, int had_index)
{
	LCCollection *inputJetCollection = pLCEvent->getCollection( jetcollection );
	LCCollection *inputMCPCollection = pLCEvent->getCollection( MCPCcollection );
	int nJets = inputJetCollection->getNumberOfElements();
	MCParticle* mcpHadron = dynamic_cast<MCParticle*>( inputMCPCollection->getElementAt(had_index));
	MCParticleVec mcpDaughters = mcpHadron->getDaughters();
	HepLorentzVector lepton_3momentum;
	HepLorentzVector jetvec;
	float lepton_theta = 0.;
	float lepton_phi = 0.;
	float match_theta = 100.;
	float match_phi = 100.;
	
	int jetIndex_thetaMatch = -1;
	int jetIndex_phiMatch = -1;
	int jetIndex_bestMatch = -1;

	for(unsigned int i_MCPD = 0; i_MCPD < mcpDaughters.size(); i_MCPD++)
	{
		if((std::abs(mcpDaughters[i_MCPD]->getPDG()) == 11) || (std::abs(mcpDaughters[i_MCPD]->getPDG()) == 13) || (std::abs(mcpDaughters[i_MCPD]->getPDG()) == 15))
		{
			lepton_3momentum = HepLorentzVector( mcpDaughters[i_MCPD]->getMomentum()[0] , mcpDaughters[i_MCPD]->getMomentum()[1] , mcpDaughters[i_MCPD]->getMomentum()[2] , mcpDaughters[i_MCPD]->getEnergy() );
			lepton_theta = lepton_3momentum.theta();
			lepton_phi = lepton_3momentum.phi();
		}
	}
	for (int i_jet = 0; i_jet < nJets; i_jet++)
	{
		ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( inputJetCollection->getElementAt( i_jet ) ) ;
		if (jet)
		{
			jetvec = HepLorentzVector ( (jet->getMomentum())[0] , (jet->getMomentum())[1] , (jet->getMomentum())[2] , jet->getEnergy() );
			delta_theta[i_jet] = jetvec.theta() - lepton_theta;
			delta_phi[i_jet] = jetvec.phi() - lepton_phi;
			if ( std::abs(delta_theta[i_jet]) < match_theta )
			{
				match_theta = std::abs(delta_theta[i_jet]);
				jetIndex_thetaMatch = i_jet;
			}
			if ( std::abs(delta_phi[i_jet]) < match_phi )
			{
				match_phi = std::abs(delta_phi[i_jet]);
				jetIndex_phiMatch = i_jet;
			}
		}
		
	}
	if ( jetIndex_thetaMatch != jetIndex_phiMatch )
	{
		if (abs(delta_theta[jetIndex_phiMatch]) < 0.5 && abs(delta_phi[jetIndex_phiMatch]) < 0.5)   // 1
		{
			if ( abs(delta_theta[jetIndex_thetaMatch]) > 0.5 && abs(delta_phi[jetIndex_thetaMatch]) > 0.5 ) //2
			{
				jetIndex_bestMatch = jetIndex_phiMatch;
			}
			else  //3
			{
				if( std::abs(delta_theta[jetIndex_phiMatch]) < std::abs(delta_phi[jetIndex_thetaMatch]) )
				{
					jetIndex_bestMatch = jetIndex_phiMatch;
				}
				if( std::abs(delta_theta[jetIndex_phiMatch]) > std::abs(delta_phi[jetIndex_thetaMatch]) )
				{
					jetIndex_bestMatch = jetIndex_thetaMatch;
				}
			}
		}
		else if (abs(delta_theta[jetIndex_thetaMatch]) < 0.5 && abs(delta_phi[jetIndex_thetaMatch]) <0.5 ) //4
		{
			jetIndex_bestMatch = jetIndex_thetaMatch;
		}
		else //5
		{
			if( std::abs(delta_theta[jetIndex_phiMatch]) < std::abs(delta_phi[jetIndex_thetaMatch]) )
			{
				jetIndex_bestMatch = jetIndex_phiMatch;
			}
			if( std::abs(delta_theta[jetIndex_phiMatch]) > std::abs(delta_phi[jetIndex_thetaMatch]) )
			{
				jetIndex_bestMatch = jetIndex_thetaMatch;
			}
		}
		streamlog_out(DEBUG)  << "BESTJET FINALLY IS   " << jetIndex_bestMatch <<std::endl;
	}
	else
	{
		jetIndex_bestMatch = jetIndex_thetaMatch;
		streamlog_out(DEBUG)  << "BESTJET FINALLY IS   " << jetIndex_bestMatch <<std::endl;
	}

	return jetIndex_bestMatch;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////					perform FIT
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> ZHllqqFit::performFIT(EVENT::LCEvent *pLCEvent, TLorentzVector Jet0_Nutlv, TLorentzVector Jet1_Nutlv)
{
	Setvalues();
	std::vector<float> FitResult{};
	LCCollection *inputJetCollection = pLCEvent->getCollection( jetcollection );
	LCCollection *inputLeptonCollection = pLCEvent->getCollection( leptoncollection );

	const int nJets = 2;
	const int nLeptons = 2;
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
	float jet0_Px,jet0_Py,jet0_Pz,jet0_E;
	float jet1_Px,jet1_Py,jet1_Pz,jet1_E;
	float lepton0_Px,lepton0_Py,lepton0_Pz,lepton0_E;
	float lepton1_Px,lepton1_Py,lepton1_Pz,lepton1_E;
	float jet0_SigmaTheta,jet0_SigmaPhi,jet0_SigmaE;
	float jet1_SigmaTheta,jet1_SigmaPhi,jet1_SigmaE;
	float lepton0_SigmaTheta,lepton0_SigmaPhi,lepton0_SigmaInvpT;
	float lepton1_SigmaTheta,lepton1_SigmaPhi,lepton1_SigmaInvpT;

	JetFitObject *jet[nJets];
	LeptonFitObject *lepton[nLeptons];
	int nSingularCovMatrix = 0;
	HepLorentzVector jetvec;
	HepLorentzVector Nuvec;
	HepLorentzVector leptonvec;
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////					Set JetFitObjects
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int i_jet = 0; i_jet < nJets; i_jet++)
	{
		ReconstructedParticle *j = dynamic_cast<ReconstructedParticle*>( inputJetCollection->getElementAt( i_jet ) );
		if ( i_jet == 0 ) Nuvec = HepLorentzVector( Jet0_Nutlv.Px() , Jet0_Nutlv.Py() , Jet0_Nutlv.Pz() , Jet0_Nutlv.E() );
		if ( i_jet == 1 ) Nuvec = HepLorentzVector( Jet1_Nutlv.Px() , Jet1_Nutlv.Py() , Jet1_Nutlv.Pz() , Jet1_Nutlv.E() );

		Px =	j->getMomentum()[0];
		Px2 =	std::pow( Px , 2 );
		Py =	j->getMomentum()[1];
		Py2 =	std::pow( Py , 2 );
		Pz =	j->getMomentum()[2];
		Pz2 =	std::pow( Pz , 2 );
		P =	std::sqrt( std::pow( Px , 2 ) + std::pow( Py , 2 ) + std::pow( Pz , 2 ) );
		P2 =	std::pow( P , 2 );
		
		if ( i_jet == 0 )
		{
			jet0_Px = Px;
			jet0_Py = Py;
			jet0_Pz = Pz;
			jet0_E = j->getEnergy();
		}
		else if ( i_jet == 1 )
		{
			jet1_Px = Px;
			jet1_Py = Py;
			jet1_Pz = Pz;
			jet1_E = j->getEnergy();
		}

		SigPx2 =	j->getCovMatrix()[0];
		SigPxSigPy =	j->getCovMatrix()[1];
		SigPxSigPz =	j->getCovMatrix()[2];
		SigPy2 =	j->getCovMatrix()[4];
		SigPySigPz =	j->getCovMatrix()[5];
		SigPz2 =	j->getCovMatrix()[7];
		SigE2 =		j->getCovMatrix()[9];

		dth_dpx =	( Px * Pz ) / ( std::pow( P , 3 ) * std::sqrt( 1 - ( Pz2 / P2 ) ) );
		dth_dpy =	( Py * Pz ) / ( std::pow( P , 3 ) * std::sqrt( 1 - ( Pz2 / P2 ) ) );
		dth_dpz =	std::sqrt( 1 - ( Pz2 / P2 ) ) / P;

		dphi_dpx =	( -Py ) / ( Px2 + Py2 );
		dphi_dpy =	( Px ) / ( Px2 + Py2 );

		JetResE =	std::sqrt( SigE2 ) * sigmaScaleFactor;
		JetResTheta =	std::sqrt( std::fabs( SigPx2 * std::pow( dth_dpx , 2 ) + SigPy2 * std::pow( dth_dpy , 2 ) + SigPz2 * std::pow( dth_dpz , 2 ) + 2 * ( SigPxSigPy * dth_dpx * dth_dpy ) + 2 * ( SigPySigPz * dth_dpy * dth_dpz ) + 2 * ( SigPxSigPz * dth_dpx * dth_dpz ) ) );
		JetResPhi =	std::sqrt( std::fabs( SigPx2 * std::pow( dphi_dpx , 2 ) + SigPy2 * std::pow( dphi_dpy , 2 ) + 2 * ( SigPxSigPy * dphi_dpx * dphi_dpy ) ) );
		streamlog_out(DEBUG) << "JET[" << i_jet << "]: JetResEnergy = " << JetResE << " , JetResTheta = " << JetResTheta << " , JetResPhi = " << JetResPhi << std::endl;

		if ( SigPx2 == 0. || SigPxSigPy == 0. || SigPxSigPz == 0. || SigPy2 == 0. || SigPySigPz == 0. || SigPz2 == 0. || SigE2 == 0. ) // Check CovMatrix singularity
		{
			streamlog_out(WARNING) << "Covariance Matrix is singular"<<std::endl;
			streamlog_out(WARNING) << "Setting theta and phi Resolution back to default values "<<std::endl;
			JetResTheta = m_jetThetaError;
			JetResPhi = m_jetPhiError;
			nSingularCovMatrix++;
		}
		jetvec = HepLorentzVector ( ( j->getMomentum() )[0] , ( j->getMomentum() )[1] , ( j->getMomentum() )[2] , j->getEnergy() );
		jetvec += Nuvec;
		if ( i_jet == 0 )
		{
			jet0_SigmaTheta = JetResTheta;
			jet0_SigmaPhi = JetResPhi;
			jet0_SigmaE = JetResE;
		}
		else if ( i_jet == 1 )
		{
			jet1_SigmaTheta = JetResTheta;
			jet1_SigmaPhi = JetResPhi;
			jet1_SigmaE = JetResE;
		}

		if (m_useErrorFlow)
		{
			jet[i_jet] = new JetFitObject ( jetvec.e(), jetvec.theta() , jetvec.phi(),
					JetResE * sigmaScaleFactor , JetResTheta , JetResPhi , jetvec.m() );
			streamlog_out(DEBUG)  << " start four-vector of jet[" << i_jet << "]: " << *jet[i_jet]  << std::endl ;
			if ( i_jet == 0 )
			{
				jet[i_jet]->setName("Jet0");
			}
			else if (i_jet == 1)
			{
				jet[i_jet]->setName("Jet1");
			}
		}
		else
		{
			jet[i_jet] = new JetFitObject ( jetvec.e(), jetvec.theta() , jetvec.phi(),
					JetEnergyResolution( jetvec.e() ) , m_jetThetaError , m_jetPhiError , jetvec.m() );
			if ( i_jet == 0 )
			{
				jet[i_jet]->setName("Jet0");
			}
			else
			{
				jet[i_jet]->setName("Jet1");
			}
			streamlog_out(DEBUG)  << " start four-vector of jet[" << i_jet << "]: " << *jet[i_jet]  << std::endl ;
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////					Set LeptonFitObjects
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(int i_lep = 0; i_lep < nLeptons ; i_lep++)
	{
		ReconstructedParticle* l = dynamic_cast<ReconstructedParticle*>( inputLeptonCollection->getElementAt( i_lep ) ) ;
		leptonvec = HepLorentzVector ((l->getMomentum())[0],(l->getMomentum())[1],(l->getMomentum())[2],l->getEnergy());
		TrackVec tckvec = l->getTracks();
		if ( tckvec.size() != 1 )
		{
			streamlog_out(DEBUG)  << "Number of tracks for lepton[" << i_lep <<"] is not exactly ONE!!! (nTracks = " << tckvec.size() << " ) " << std::endl ;
			invers_pT = 1/(std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2)));
			invers_pT_err = 2*std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2))*0.00001;
			theta = leptonvec.theta();
			theta_err = m_jetThetaError;
			phi = leptonvec.phi();
			phi_err = m_jetPhiError;
		}
		else
		{
			streamlog_out(DEBUG)  << "Number of tracks for lepton[" << i_lep <<"] is exactly ONE!!!" << std::endl;
			Omega = tckvec[0]->getOmega();
			Omega_uncert = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[5]) );
			streamlog_out(DEBUG)  << "Omega_uncert = " << Omega_uncert << std::endl;
			TanLambda = tckvec[0]->getTanLambda();
			TanLambda_err = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[14]) );
			streamlog_out(DEBUG)  << "TanLambda_err = " << TanLambda_err << std::endl;
			theta = 2 * atan( 1. ) - atan( TanLambda );
			theta_err = TanLambda_err / ( 1 + pow( TanLambda , 2 ) );
			streamlog_out(DEBUG)  << "theta_err = " << theta_err << std::endl;
			phi = tckvec[0]->getPhi();
			phi_err = std::sqrt( std::abs(tckvec[0]->getCovMatrix()[2]) );
			streamlog_out(DEBUG)  << "phi_err = " << phi_err << std::endl;
			invers_pT = Omega / eB;
			invers_pT_err = std::fabs( 1. / eB ) * Omega_uncert;
			streamlog_out(DEBUG)  << "invers_pT_err = " << invers_pT_err << std::endl;
			if ( i_lep == 0 )
			{
				lepton0_Px = cos( phi ) / invers_pT;
				lepton0_Py = sin( phi ) / invers_pT;
				lepton0_Pz = TanLambda / invers_pT;
				lepton0_E = l->getEnergy();
			}
			else if ( i_lep == 1 )
			{
				lepton1_Px = cos( phi ) / invers_pT;
				lepton1_Py = sin( phi ) / invers_pT;
				lepton1_Pz = TanLambda / invers_pT;
				lepton1_E = l->getEnergy();
			}
			
		}
		if ( i_lep == 0 )
		{
			lepton0_SigmaTheta = theta_err;
			lepton0_SigmaPhi = phi_err;
			lepton0_SigmaInvpT = invers_pT_err;
		}
		else if ( i_lep == 1 )
		{
			lepton1_SigmaTheta = theta_err;
			lepton1_SigmaPhi = phi_err;
			lepton1_SigmaInvpT = invers_pT_err;
		}		
		streamlog_out(DEBUG)  << "Lepton fit object from leptonvec: " 
			<< 1/(std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2))) <<" +- " << 2*std::sqrt(pow(leptonvec.px(),2)+pow(leptonvec.py(),2))*0.00001 << " , " 
			<< leptonvec.theta() <<" +- " << m_jetThetaError << " , " 
			<< leptonvec.phi() <<" +- " << m_jetPhiError << std::endl ;

		streamlog_out(DEBUG)  << "Lepton fit object from track:     " 
			<< std::fabs( tckvec[0]->getOmega() / eB ) <<" +- " << std::fabs( 1. / eB ) * std::sqrt( tckvec[0]->getCovMatrix()[5] ) << " , " 
			<< 2 * atan( 1. ) - atan( tckvec[0]->getTanLambda() ) <<" +- " << std::abs( std::sqrt( tckvec[0]->getCovMatrix()[14]) ) / ( 1 + pow( tckvec[0]->getTanLambda() , 2 ) ) << " , " 
			<< tckvec[0]->getPhi() <<" +- " << std::abs( std::sqrt( tckvec[0]->getCovMatrix()[2] ) ) << std::endl ;

		lepton[i_lep] = new LeptonFitObject (invers_pT , theta , phi , invers_pT_err , theta_err , phi_err, leptonvec.m());
		if (i_lep == 0 )
		{
			lepton[i_lep]->setName("Lepton0");
		}
		else if (i_lep == 1)
		{
			lepton[i_lep]->setName("Lepton1");
		}
		streamlog_out(DEBUG)  << " start four-vector of lepton[" << i_lep <<"]: " << *lepton[i_lep]  << std::endl ;
	}

	const int NJETS = 2;
	streamlog_out(MESSAGE) << "*jet[0]: " << *jet[0] << std::endl ;
	streamlog_out(MESSAGE) << "*jet[1]: " << *jet[1] << std::endl ;

	const int NLEPTONS = 2;
	streamlog_out(MESSAGE) << "*lepton[0]: " << *lepton[0] << std::endl ; 
	streamlog_out(MESSAGE) << "*lepton[1]: " << *lepton[1] << std::endl ;

////	these don't get changed by the fit -> to obtain start values later!
	JetFitObject startjets[NJETS] = {*jet[0], *jet[1]};
	for (int i = 0; i < NJETS; ++i)	streamlog_out(MESSAGE)  << "startjets[" << i << "]: " << startjets[i]  << std::endl;

	LeptonFitObject startleptons[NLEPTONS] = {*lepton[0], *lepton[1]};
	for (int i = 0; i < NLEPTONS; ++i) streamlog_out(MESSAGE)  << "startleptons[" << i << "]: " << startleptons[i]  << std::endl;

	JetFitObject fitjets[NJETS] = {*jet[0], *jet[1]};
	for (int i = 0; i < NJETS; ++i)
		streamlog_out(MESSAGE)  << "fitjets[" << i << "]: " << fitjets[i]  << std::endl ;

	LeptonFitObject fitleptons[NLEPTONS] = {*lepton[0], *lepton[1]};
	for (int i = 0; i < NLEPTONS; ++i)
		streamlog_out(MESSAGE)  << "fitleptons[" << i << "]: " << fitleptons[i]  << std::endl ;

////	these get changed by the fit
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
////				set constraints befor fit
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	float target_p_due_crossing_angle = m_ECM * 0.007; // crossing angle = 14 mrad
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

	E_lab= 2 * sqrt( std::pow( 0.548579909e-3 , 2 ) + std::pow( m_ECM / 2 , 2 ) + std::pow( target_p_due_crossing_angle , 2 ) + 0. + 0.);
	MomentumConstraint ec(1, 0, 0, 0, E_lab);
	ec.setName("sum(E)");
	for (int i = 0; i < NJETS; ++i) ec.addToFOList(*(jets[i]));
	for (int i = 0; i < NLEPTONS; ++i) ec.addToFOList(*(leptons[i]));

	streamlog_out(MESSAGE)  << "Value of E_lab before adding ISR: " << E_lab << std::endl ;
	streamlog_out(MESSAGE)  << "Value of target_p_due_crossing_angle before adding ISR: " << target_p_due_crossing_angle << std::endl ;
	streamlog_out(MESSAGE)  << "Value of pxc before adding ISR: " << pxc.getValue() << std::endl ;
	streamlog_out(MESSAGE)  << "Value of pyc before adding ISR: " << pyc.getValue() << std::endl ;
	streamlog_out(MESSAGE)  << "Value of pzc before adding ISR: " << pzc.getValue() << std::endl ;
	streamlog_out(MESSAGE)  << "Value of ec before adding ISR: " << ec.getValue() << std::endl ;
	pxc_before_ISR = pxc.getValue();
	pyc_before_ISR = pyc.getValue();
	pzc_before_ISR = pzc.getValue();
	ec_before_ISR = ec.getValue();

////	ISR Photon initialized with missing p_z
	ISRPhotonFitObject *photon = new ISRPhotonFitObject (0., 0., -pzc.getValue(), b, ISRPzMaxB);
////	ISRPhotonFitObject(double px, double py, double pz,
////	double b_, double PzMaxB_, double PzMinB_ = 0.);
	if(m_fitISR)
	{
		streamlog_out(MESSAGE)  << "start four-vector of ISR photon: " << *(photon) << std::endl ;

		pxc.addToFOList(*(photon));
		pyc.addToFOList(*(photon));
		pzc.addToFOList(*(photon));
		ec.addToFOList(*(photon));
	}

	streamlog_out(MESSAGE)  << "Value of E_lab before fit: " << E_lab << std::endl ;
	streamlog_out(MESSAGE)  << "Value of target_p_due_crossing_angle before fit: " << target_p_due_crossing_angle << std::endl ;
	streamlog_out(MESSAGE)  << "Value of pxc after adding ISR before fit: " << pxc.getValue() << std::endl ;
	streamlog_out(MESSAGE)  << "Value of pyc after adding ISR before fit: " << pyc.getValue() << std::endl ;
	streamlog_out(MESSAGE)  << "Value of pzc after adding ISR before fit: " << pzc.getValue() << std::endl ;
	streamlog_out(MESSAGE)  << "Value of ec after adding ISR before fit: " << ec.getValue() << std::endl ;
	pxc_before_fit = pxc.getValue();
	pyc_before_fit = pyc.getValue();
	pzc_before_fit = pzc.getValue();
	ec_before_fit = ec.getValue();

	SoftGaussMassConstraint z(91.2,2.4952/2);
	z.addToFOList(*(leptons[0]), 1);
	z.addToFOList(*(leptons[1]), 1);

	MassConstraint h(125.);
	h.addToFOList(*(jets[0]), 1);
	h.addToFOList(*(jets[1]), 1);

	startmassZ = z.getMass(1);
	startmassH = h.getMass(1);

	streamlog_out(MESSAGE) << "start mass of Z: " << startmassZ << std::endl ;
	streamlog_out(MESSAGE) << "start mass of H: " << startmassH << std::endl ;

	Hmass_NoFit = startmassH;

	BaseFitter *pfitter;

	int debug = 0;
	if ( pLCEvent->getEventNumber() == m_ievttrace || m_traceall) debug = 10;

	if (m_fitter == 1)
	{
		pfitter = new NewFitterGSL();
		if (pLCEvent->getEventNumber()== m_ievttrace || m_traceall) (dynamic_cast<NewFitterGSL*>(pfitter))->setDebug(debug);

		streamlog_out(DEBUG4) << "ifitter is 1"  << std::endl ;
	}
	else if (m_fitter == 2)
	{
		pfitter = new NewtonFitterGSL();
		if (pLCEvent->getEventNumber()== m_ievttrace || m_traceall) (dynamic_cast<NewtonFitterGSL*>(pfitter))->setDebug(debug);

		streamlog_out(DEBUG4) << "ifitter is 2"  << std::endl ;
	}
	else
	{
////		OPALFitter has no method setDebug !
		pfitter = new OPALFitterGSL();

		streamlog_out(DEBUG4) << "ifitter is not 1 or 2"  << std::endl ;
		if (pLCEvent->getEventNumber()== m_ievttrace || m_traceall) (dynamic_cast<OPALFitterGSL*>(pfitter))->setDebug(debug);
	}
	BaseFitter &fitter = *pfitter;

	TextTracer tracer (std::cout);
	if (pLCEvent->getEventNumber() == m_ievttrace || m_traceall) fitter.setTracer(tracer);
	for (int i = 0; i < NJETS; ++i) fitter.addFitObject(*(jets[i]));
	for (int i = 0; i < NLEPTONS; ++i) fitter.addFitObject(*(leptons[i]));

	if(m_fitISR)
	{
		fitter.addFitObject (*(photon));
		streamlog_out(DEBUG4) << "ISR added to fit"  << std::endl ;
	}

	fitter.addConstraint(pxc);
	fitter.addConstraint(pyc);
	fitter.addConstraint(pzc);
	fitter.addConstraint(ec);
	fitter.addConstraint(z);

	streamlog_out(DEBUG4) << "constraints added"  << std::endl ;

	if (fabs(startmassZ-91.2) + fabs(startmassH-125.) < bestzvalue)
	{
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

	for (int i = 0; i < NJETS; ++i)
	{
		streamlog_out(DEBUG)  << "final four-vector of jet " << i << ": " << *(jets[i]) << std::endl ;
		streamlog_out(DEBUG)  << "final px of jet " << i << ": " << (jets[i]) << std::endl ;
	}
	for (int i = 0; i < NLEPTONS; ++i)
	{
		streamlog_out(DEBUG)  << "final four-vector of lepton " << i << ": " << *(leptons[i]) << std::endl ;
		streamlog_out(DEBUG)  << "final px of lepton " << i << ": " << (leptons[i]) << std::endl ;
	}
	if(m_fitISR) streamlog_out(DEBUG)  << "final four-vector of ISR photon: " << *(photon) << std::endl;
	int ierr = fitter.getError();
	streamlog_out(MESSAGE)  << "fitter error: " << ierr << std::endl;
	if ((besterr > 0 && ierr < besterr) || ( besterr < 0 && ierr == 0)) besterr = ierr;

	streamlog_out(MESSAGE)  << "Value of pxc after fit: " << pxc.getValue() << std::endl ;
	streamlog_out(MESSAGE)  << "Value of pyc after fit: " << pyc.getValue() << std::endl ;
	streamlog_out(MESSAGE)  << "Value of pzc after fit: " << pzc.getValue() << std::endl ;
	streamlog_out(MESSAGE)  << "Value of ec after fit: " << ec.getValue() << std::endl ;
	pxc_after_fit = pxc.getValue();
	pyc_after_fit = pyc.getValue();
	pzc_after_fit = pzc.getValue();
	ec_after_fit = ec.getValue();
	

	if (ierr <= 0)
	{
		double pullJet[3][2];
		double pullLepton[3][2];
////		require successfull error calculation for pulls!
		if (ierr == 0)
		{
			for (int ifo = 0; ifo < 2; ifo++)
			{
				double start, fitted;
				double errfit, errmea, sigma;
				for (int ipar = 0; ipar < 3; ipar++)
				{
					fitted = jets[ifo]->getParam(ipar);
					start = startjets[ifo].getParam(ipar);
					errfit = jets[ifo]->getError(ipar);
					errmea = startjets[ifo].getError(ipar);
					sigma = errmea*errmea-errfit*errfit;
					if (sigma > 0)
					{
						sigma = sqrt(sigma);
						pullJet[ipar][ifo] = (fitted - start)/sigma;
					}
					else
					{
						pullJet[ipar][ifo] = -4.5;
					}
				}
			}
			for (int ifo = 0; ifo < 2; ifo++)
			{
				double start, fitted;
				double errfit, errmea, sigma;
				for (int ipar = 0; ipar < 3; ipar++)
				{
					fitted = leptons[ifo]->getParam(ipar);
					start = startleptons[ifo].getParam(ipar);
					errfit = leptons[ifo]->getError(ipar);
					errmea = startleptons[ifo].getError(ipar);
					sigma = errmea*errmea-errfit*errfit;
					if (sigma > 0)
					{
						sigma = sqrt(sigma);
						pullLepton[ipar][ifo] = (fitted - start)/sigma;
					}
					else
					{
						pullLepton[ipar][ifo] = -4.5;
					}
				}
			}
		}
		if (prob >= bestprob)
		{
			bestprob = prob;
			streamlog_out(DEBUG)  << "BESTPROB: " << bestprob << std::endl ;
			bestnit  = nit;
			bestmassZ = z.getMass(1);
			bestmassH = h.getMass(1);
			beststartmassZ = startmassZ;
			beststartmassH = startmassH;
			bestphotonenergy = photon->getE();
			ISRmomentum[0] = photon->getPx();
			ISRmomentum[1] = photon->getPy();
			ISRmomentum[2] = photon->getPz();
			Zmomentum[0] = leptons[0]->getPx() + leptons[1]->getPx();
			Zmomentum[1] = leptons[0]->getPy() + leptons[1]->getPy();
			Zmomentum[2] = leptons[0]->getPz() + leptons[1]->getPz();
			Hmomentum[0] = jets[0]->getPx() + jets[1]->getPx();
			Hmomentum[1] = jets[0]->getPy() + jets[1]->getPy();
			Hmomentum[2] = jets[0]->getPz() + jets[1]->getPz();
			Z_Energy = leptons[0]->getE() + leptons[1]->getE();
			H_Energy = jets[0]->getE() + jets[1]->getE();
			chi2best = fitter.getChi2();
			errorcode = fitter.getError();
			if (ierr == 0) //if  fitter.getError() is = 0
			{
				hpull_jet_E=pullJet[0][0];
				hpull_jet2_E=pullJet[0][1];
				hpull_jet_th=pullJet[1][0];
				hpull_jet2_th=pullJet[1][1];
				hpull_jet_phi=pullJet[2][0];
				hpull_jet2_phi=pullJet[2][1];
				hpull_lepton_InvpT=pullLepton[0][0];
				hpull_lepton2_InvpT=pullLepton[0][1];
				hpull_lepton_th=pullLepton[1][0];
				hpull_lepton2_th=pullLepton[1][1];
				hpull_lepton_phi=pullLepton[2][0];
				hpull_lepton2_phi=pullLepton[2][1];
			}//end if  fitter.getError() is = 0
			else //if  fitter.getError() is not = 0
			{
				streamlog_out(DEBUG) << " ERROR CALCULATION FAILED for best permutation in event " << pLCEvent->getEventNumber() << std::endl ;
			}
		}
	}//end-if fitter.getError() <=0
	else
	{
		streamlog_out(DEBUG4) << "FIT ERROR = " << fitter.getError()
					<< " in event " << pLCEvent->getEventNumber()
					<< ", not filling histograms!"  << std::endl ;
		streamlog_out(DEBUG4)  << "start mass of Z: " << startmassZ << std::endl ;
		streamlog_out(DEBUG4)  << "start mass of H: " << startmassH << std::endl ;
		streamlog_out(DEBUG4)  << "final mass of Z: " << z.getMass(1) << std::endl ;
		streamlog_out(DEBUG4)  << "final mass of H: " << h.getMass(1) << std::endl ;
	}

	Zmass_before_fit=beststartmassZ;
	Zmass_after_fit=bestmassZ;
	Hmass_before_fit=beststartmassH;
	Hmass_after_fit=bestmassH;
	Error_code=errorcode;

	streamlog_out(DEBUG)  << "min chi2 start mass of Z: " << chi2startmassZ << std::endl ;
	streamlog_out(DEBUG)  << "min chi2 start mass of H: " << chi2startmassH << std::endl ;
	streamlog_out(DEBUG)  << "best start mass of Z: " << beststartmassZ << std::endl ;
	streamlog_out(DEBUG)  << "best start mass of H: " << beststartmassH << std::endl ;
	streamlog_out(DEBUG)  << "best mass of Z: " << bestmassZ << std::endl ;
	streamlog_out(DEBUG)  << "best mass of H: " << bestmassH << std::endl ;
	streamlog_out(DEBUG)  << "Error Code: " << errorcode << std::endl ;

	FitResult.push_back(ierr);
	FitResult.push_back(prob);
	FitResult.push_back(nit);
	FitResult.push_back(startmassZ);
	FitResult.push_back(startmassH);
	FitResult.push_back(beststartmassZ);
	FitResult.push_back(beststartmassH);
	FitResult.push_back(Zmass_after_fit);
	FitResult.push_back(Hmass_after_fit);
	FitResult.push_back(chi2startmassZ);
	FitResult.push_back(chi2startmassH);
	FitResult.push_back(chi2best);
	FitResult.push_back(bestphotonenergy);
	FitResult.push_back(hpull_jet_E);
	FitResult.push_back(hpull_jet2_E);
	FitResult.push_back(hpull_jet_th);
	FitResult.push_back(hpull_jet2_th);
	FitResult.push_back(hpull_jet_phi);
	FitResult.push_back(hpull_jet2_phi);
	FitResult.push_back(hpull_lepton_InvpT);
	FitResult.push_back(hpull_lepton2_InvpT);
	FitResult.push_back(hpull_lepton_th);
	FitResult.push_back(hpull_lepton2_th);
	FitResult.push_back(hpull_lepton_phi);
	FitResult.push_back(hpull_lepton2_phi);
	FitResult.push_back(ISRmomentum[0]);
	FitResult.push_back(ISRmomentum[1]);
	FitResult.push_back(ISRmomentum[2]);
	FitResult.push_back(Zmomentum[0]);
	FitResult.push_back(Zmomentum[1]);
	FitResult.push_back(Zmomentum[2]);
	FitResult.push_back(Z_Energy);
	FitResult.push_back(Hmomentum[0]);
	FitResult.push_back(Hmomentum[1]);
	FitResult.push_back(Hmomentum[2]);
	FitResult.push_back(H_Energy);
	FitResult.push_back(pxc_before_ISR);
	FitResult.push_back(pyc_before_ISR);
	FitResult.push_back(pzc_before_ISR);
	FitResult.push_back(ec_before_ISR);
	FitResult.push_back(pxc_before_fit);
	FitResult.push_back(pyc_before_fit);
	FitResult.push_back(pzc_before_fit);
	FitResult.push_back(ec_before_fit);
	FitResult.push_back(pxc_after_fit);
	FitResult.push_back(pyc_after_fit);
	FitResult.push_back(pzc_after_fit);
	FitResult.push_back(ec_after_fit);
	FitResult.push_back(jet0_Px);
	FitResult.push_back(jet1_Px);
	FitResult.push_back(jet0_Py);
	FitResult.push_back(jet1_Py);
	FitResult.push_back(jet0_Pz);
	FitResult.push_back(jet1_Pz);
	FitResult.push_back(jet0_E);
	FitResult.push_back(jet1_E);
	FitResult.push_back(lepton0_Px);
	FitResult.push_back(lepton1_Px);
	FitResult.push_back(lepton0_Py);
	FitResult.push_back(lepton1_Py);
	FitResult.push_back(lepton0_Pz);
	FitResult.push_back(lepton1_Pz);
	FitResult.push_back(lepton0_E);
	FitResult.push_back(lepton1_E);
	FitResult.push_back(jet0_SigmaTheta);
	FitResult.push_back(jet1_SigmaTheta);
	FitResult.push_back(jet0_SigmaPhi);
	FitResult.push_back(jet1_SigmaPhi);
	FitResult.push_back(jet0_SigmaE);
	FitResult.push_back(jet1_SigmaE);
	FitResult.push_back(lepton0_SigmaTheta);
	FitResult.push_back(lepton1_SigmaTheta);
	FitResult.push_back(lepton0_SigmaPhi);
	FitResult.push_back(lepton1_SigmaPhi);
	FitResult.push_back(lepton0_SigmaInvpT);
	FitResult.push_back(lepton1_SigmaInvpT);
	
	delete photon;
	return FitResult;
}

void ZHllqqFit::check( LCEvent* )
{
//	nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZHllqqFit::end()
{
	streamlog_out(ERROR) << "# of events: " << m_nEvt << std::endl;
//	streamlog_out(ERROR) << "# of nucorrection: " << correction<< std::endl;
//	streamlog_out(ERROR) << "# of Covariance failed: " << nCo<< std::endl;

	m_pTFile->cd();
	m_pTTree->Write();
	m_pTTree_0->Write();
	m_pTTree_1->Write();
	m_pTTree_2->Write();
	m_pTTree_3->Write();
	h_Zmass_beforefit_woNu->Write();
	h_Hmass_beforefit_woNu->Write();
	h_Zmass_beforefit_wNu->Write();
	h_Hmass_beforefit_wNu->Write();
	h_Zmass_afterfit_woNu->Write();
	h_Hmass_afterfit_woNu->Write();
	h_Zmass_afterfit_wNu->Write();
	h_Hmass_afterfit_wNu->Write();
	h_fitError_wNu->Write();
	h_fitError_woNu->Write();
	h_ErrorCode_wNu_woNu->Write();
	h_fitProbability_wNu->Write();
	h_fitProbability_woNu->Write();
	h_fitProbability_best->Write();
	h_nJets->Write();
	h_nLeptons->Write();
	h_nLeptons_nJets->Write();
	h_ISR_mcp_fit->Write();
	h_ISR_pzc_woNu->Write();
	h_ISR_pzc_wNu->Write();
	h_ISR_pzc_best->Write();
	h_pull_jet_E->Write();
	h_pull_jet_theta->Write();
	h_pull_jet_phi->Write();
	h_pull_lepton_InvPt->Write();
	h_pull_lepton_theta->Write();
	h_pull_lepton_phi->Write();
	m_pTFile->Close();
	delete m_pTFile;

}
