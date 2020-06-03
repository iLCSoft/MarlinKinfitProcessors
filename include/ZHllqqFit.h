#ifndef ZHllqqFit_h
#define ZHllqqFit_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "TLorentzVector.h"
#include <TFile.h>
#include <TTree.h>
#include "TopEventILC.h"
#include <vector>
class TFile;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;

using namespace lcio ;
using namespace marlin ;

class ZHllqqFit : public Processor
{

	public:

		virtual Processor*  newProcessor()
		{
			return new ZHllqqFit;
		}
		ZHllqqFit() ;
		virtual ~ZHllqqFit() = default;
		ZHllqqFit(const ZHllqqFit&) = delete;
		ZHllqqFit& operator=(const ZHllqqFit&) = delete;
		virtual void				init();
		virtual void				processRunHeader() ;
		virtual void				processEvent( EVENT::LCEvent *pLCEvent ) ;
		virtual void				check( LCEvent * evt ) ;
		virtual void				end() ;
		int					FindMatchingJettoSLD(EVENT::LCEvent *pLCEvent, int had_index);
		std::vector<float>			performFIT(EVENT::LCEvent *pLCEvent, TLorentzVector Jet0_Nutlv, TLorentzVector Jet1_Nutlv);

		double					JetEnergyResolution(double E);
		void					compcorrect();
		void					Setvalues();
		void					Clear();

	private:

		std::string				jetcollection{};
		std::string				hDecayMode{};
		std::string				leptoncollection{};
		std::string				SLDecayCollection{};
		std::string				NuEnergyCollectionB{};
		std::string				NuEnergyCollectionC{};
		std::string				MCPCcollection{} ;
		std::string				errorflowcollection{};
		std::string				outputFitcollection{};
		bool					m_includeHbb;
		bool					m_includeHcc;
		bool					m_includeHother;
		typedef	std::vector<int>		IntVector;
		typedef	std::vector<float>		floatVector;
		bool					m_fitISR = true;
		bool					m_fitNuE = true;
		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		float					m_Bfield;
		double					c;
		double					mm2m;
		double					eV2GeV;
		double					eB;
		float					m_ECM;
		bool					m_useErrorFlow{};
		float					m_isrpzmax{};
		float					sigmaScaleFactor{};
		double					m_jetEnergyError{};
		double					m_jetThetaError{};
		double					m_jetPhiError{};
		int					m_fitter{};
		bool					m_traceall{};
		int					m_ievttrace{};
		std::string				m_outputFile{};


		int					m_nSLDecayBHadron;
		int					m_nSLDecayCHadron;
		int					m_nSLDecayTotal;
		int					m_nJets;
		int					m_nLeptons;

		IntVector				m_bestNuCombination{};
		IntVector				m_iError_wNu{};
		floatVector				m_probability_wNu{};
		floatVector				m_chi2_wNu{};
		IntVector				m_n_itter_wNu{};
		floatVector				m_startmassZ_wNu{};
		floatVector				m_startmassH_wNu{};
		floatVector				m_beststartmassZ_wNu{};
		floatVector				m_beststartmassH_wNu{};
		floatVector				m_Zmass_after_fit_wNu{};
		floatVector				m_Hmass_after_fit_wNu{};
		floatVector				m_bestphotonenergy_wNu{};
		floatVector				m_chi2startmassZ_wNu{};
		floatVector				m_chi2startmassH_wNu{};
		floatVector				m_jet_startPx_wNu{};
		floatVector				m_jet_startPy_wNu{};
		floatVector				m_jet_startPz_wNu{};
		floatVector				m_jet_startE_wNu{};
		floatVector				m_jet_startPx_wNu_bestfit{};
		floatVector				m_jet_startPy_wNu_bestfit{};
		floatVector				m_jet_startPz_wNu_bestfit{};
		floatVector				m_jet_startE_wNu_bestfit{};
		floatVector				m_jet_startPx_woNu{};
		floatVector				m_jet_startPy_woNu{};
		floatVector				m_jet_startPz_woNu{};
		floatVector				m_jet_startE_woNu{};
		floatVector				m_jet_startPx_best{};
		floatVector				m_jet_startPy_best{};
		floatVector				m_jet_startPz_best{};
		floatVector				m_jet_startE_best{};
		floatVector				m_lepton_startPx_wNu{};
		floatVector				m_lepton_startPy_wNu{};
		floatVector				m_lepton_startPz_wNu{};
		floatVector				m_lepton_startE_wNu{};
		floatVector				m_lepton_startPx_wNu_bestfit{};
		floatVector				m_lepton_startPy_wNu_bestfit{};
		floatVector				m_lepton_startPz_wNu_bestfit{};
		floatVector				m_lepton_startE_wNu_bestfit{};
		floatVector				m_lepton_startPx_woNu{};
		floatVector				m_lepton_startPy_woNu{};
		floatVector				m_lepton_startPz_woNu{};
		floatVector				m_lepton_startE_woNu{};
		floatVector				m_lepton_startPx_best{};
		floatVector				m_lepton_startPy_best{};
		floatVector				m_lepton_startPz_best{};
		floatVector				m_lepton_startE_best{};
		floatVector				m_jet_SigmaTheta{};
		floatVector				m_jet_SigmaPhi{};
		floatVector				m_jet_SigmaE{};
		floatVector				m_lepton_startPx{};
		floatVector				m_lepton_startPy{};
		floatVector				m_lepton_startPz{};
		floatVector				m_lepton_startE{};
		floatVector				m_lepton_SigmaTheta{};
		floatVector				m_lepton_SigmaPhi{};
		floatVector				m_lepton_SigmaInvpT{};
		floatVector				m_pull_jet_E_wNu{};
		floatVector				m_pull_jet_th_wNu{};
		floatVector				m_pull_jet_phi_wNu{};
		floatVector				m_pull_lepton_InvpT_wNu{};
		floatVector				m_pull_lepton_th_wNu{};
		floatVector				m_pull_lepton_phi_wNu{};
		int					m_iError_wNu_bestfit;
		float					m_probability_wNu_bestfit;
		float					m_chi2_wNu_bestfit;
		int					m_n_itter_wNu_bestfit;
		float					m_startmassZ_wNu_bestfit;
		float					m_startmassH_wNu_bestfit;
		float					m_beststartmassZ_wNu_bestfit;
		float					m_beststartmassH_wNu_bestfit;
		float					m_Zmass_after_fit_wNu_bestfit;
		float					m_Hmass_after_fit_wNu_bestfit;
		float					m_bestphotonenergy_wNu_bestfit;
		float					m_chi2startmassZ_wNu_bestfit;
		float					m_chi2startmassH_wNu_bestfit;
		floatVector				m_pull_jet_E_wNu_bestfit{};
		floatVector				m_pull_jet_th_wNu_bestfit{};
		floatVector				m_pull_jet_phi_wNu_bestfit{};
		floatVector				m_pull_lepton_InvpT_wNu_bestfit{};
		floatVector				m_pull_lepton_th_wNu_bestfit{};
		floatVector				m_pull_lepton_phi_wNu_bestfit{};

		int					m_iError_woNu;
		float					m_probability_woNu;
		float					m_chi2best_woNu;
		int					m_n_itter_woNu;
		float					m_startmassZ_woNu;
		float					m_startmassH_woNu;
		float					m_beststartmassZ_woNu;
		float					m_beststartmassH_woNu;
		float					m_Zmass_after_fit_woNu;
		float					m_Hmass_after_fit_woNu;
		float					m_bestphotonenergy_woNu;
		float					m_chi2startmassZ_woNu;
		float					m_chi2startmassH_woNu;
		floatVector				m_pull_jet_E_woNu{};
		floatVector				m_pull_jet_th_woNu{};
		floatVector				m_pull_jet_phi_woNu{};
		floatVector				m_pull_lepton_InvpT_woNu{};
		floatVector				m_pull_lepton_th_woNu{};
		floatVector				m_pull_lepton_phi_woNu{};

		int					m_iError_best;
		float					m_probability_best;
		float					m_chi2_best;
		int					m_n_itter_best;
		float					m_startmassZ_best;
		float					m_startmassH_best;
		float					m_beststartmassZ_best;
		float					m_beststartmassH_best;
		float					m_Zmass_after_fit_best;
		float					m_Hmass_after_fit_best;
		float					m_bestphotonenergy_best;
		float					m_chi2startmassZ_best;
		float					m_chi2startmassH_best;
		floatVector				m_pull_jet_E_best{};
		floatVector				m_pull_jet_th_best{};
		floatVector				m_pull_jet_phi_best{};
		floatVector				m_pull_lepton_InvpT_best{};
		floatVector				m_pull_lepton_th_best{};
		floatVector				m_pull_lepton_phi_best{};
		float					chi2best;
		float					errorcode;
		double					E_lab;
		float					m_pxc_before_ISR_wNu;
		float					m_pyc_before_ISR_wNu;
		float					m_pzc_before_ISR_wNu;
		float					m_ec_before_ISR_wNu;
		float					m_pxc_before_fit_wNu;
		float					m_pyc_before_fit_wNu;
		float					m_pzc_before_fit_wNu;
		float					m_ec_before_fit_wNu;
		float					m_pxc_after_fit_wNu;
		float					m_pyc_after_fit_wNu;
		float					m_pzc_after_fit_wNu;
		float					m_ec_after_fit_wNu;
		float					m_pxc_before_ISR_woNu;
		float					m_pyc_before_ISR_woNu;
		float					m_pzc_before_ISR_woNu;
		float					m_ec_before_ISR_woNu;
		float					m_pxc_before_fit_woNu;
		float					m_pyc_before_fit_woNu;
		float					m_pzc_before_fit_woNu;
		float					m_ec_before_fit_woNu;
		float					m_pxc_after_fit_woNu;
		float					m_pyc_after_fit_woNu;
		float					m_pzc_after_fit_woNu;
		float					m_ec_after_fit_woNu;
		float					m_pxc_before_ISR_best;
		float					m_pyc_before_ISR_best;
		float					m_pzc_before_ISR_best;
		float					m_ec_before_ISR_best;
		float					m_pxc_before_fit_best;
		float					m_pyc_before_fit_best;
		float					m_pzc_before_fit_best;
		float					m_ec_before_fit_best;
		float					m_pxc_after_fit_best;
		float					m_pyc_after_fit_best;
		float					m_pzc_after_fit_best;
		float					m_ec_after_fit_best;
		float					m_ISR_startPx_wNu;
		float					m_ISR_startPy_wNu;
		float					m_ISR_startPz_wNu;
		float					m_ISR_startPx_woNu;
		float					m_ISR_startPy_woNu;
		float					m_ISR_startPz_woNu;


		TFile					*m_pTFile;
	        TTree					*m_pTTree;
	        TTree					*m_pTTree_0;
	        TTree					*m_pTTree_1;
	        TTree					*m_pTTree_2;
	        TTree					*m_pTTree_3;
		TH1F					*h_Zmass_beforefit_woNu{};
		TH1F					*h_Hmass_beforefit_woNu{};
		TH1F					*h_Zmass_beforefit_wNu{};
		TH1F					*h_Hmass_beforefit_wNu{};
		TH1F					*h_Zmass_afterfit_woNu{};
		TH1F					*h_Hmass_afterfit_woNu{};
		TH1F					*h_Zmass_afterfit_wNu{};
		TH1F					*h_Hmass_afterfit_wNu{};
		TH1I					*h_fitError_wNu{};
		TH1I					*h_fitError_woNu{};
		TH2I					*h_ErrorCode_wNu_woNu{};
		TH1F					*h_fitProbability_wNu{};
		TH1F					*h_fitProbability_woNu{};
		TH1F					*h_fitProbability_best{};
		TH1I					*h_nJets{};
		TH1I					*h_nLeptons{};
		TH2I					*h_nLeptons_nJets{};
		TH2F					*h_ISRE_1mcp_fit{};
		TH2F					*h_ISRE_2mcp_fit{};
		TH2F					*h_ISRpz_1mcp_fit{};
		TH2F					*h_ISRpz_2mcp_fit{};
		TH2F					*h_ISR_pzc_wNu{};
		TH2F					*h_ISR_pzc_woNu{};
		TH2F					*h_ISR_pzc_best{};
		TH1F					*h_pull_jet_E{};
		TH1F					*h_pull_jet_theta{};
		TH1F					*h_pull_jet_phi{};
		TH1F					*h_pull_lepton_InvPt{};
		TH1F					*h_pull_lepton_theta{};
		TH1F					*h_pull_lepton_phi{};


		std::string				_name{};
		std::string				_OutputCollection{};
		double					b{};
		double					ISRPzMaxB{};

		float					prob{};
		float					bestprob{};
		float					bestnit{};
		float					bestmassZ{};
		float					bestmassH{};
		float					beststartmassZ{};
		float					beststartmassH{};
		float					bestphotonenergy{};
		float					startmassZ{};
		float					startmassH{};
		float					variable{};
		float					Zmomentum[3]{};
		float					Hmomentum[3]{};
		float					ISRmomentum[3]{};
		float					Z_Energy{};
		float					H_Energy{};
		float					momentum[3]{};
		float					energy{};
		int					_nRun{};
		int					_nEvt{};
		int					nit{};
		int					nCo{};

		int					bestperm{};
		int					errorflag{};

		float					Zmass_before_fit{};
		float					Zmass_after_fit{};
		float					Hmass_before_fit{};
		float					Hmass_after_fit{};
		float					Hmass_NoFit{};
		float					Error_code{};
		float					hpull_jet_E{};
		float					hpull_jet2_E{};
		float					hpull_jet_th{};
		float					hpull_jet2_th{};
		float					hpull_jet_phi{};
		float					hpull_jet2_phi{};
		float					hpull_lepton_InvpT{};
		float					hpull_lepton2_InvpT{};
		float					hpull_lepton_th{};
		float					hpull_lepton2_th{};
		float					hpull_lepton_phi{};
		float					hpull_lepton2_phi{};
		int					jetmatch{};
		int					jetmatchth{};
		int					jetmatchphi{};

		int					besterr{};
		double					bestzvalue{} ;
		double					chi2startmassZ{};
		double					chi2startmassH {};

		double					Px{};
		double					Px2{};
		double					Py{};
		double					Py2{};
		double					Pz{};
		double					Pz2{};
		double					P{};
		double					P2{};
		double					SigPx2{};
		double					SigPxSigPy{};
		double					SigPxSigPz{};
		double					SigPy2{};
		double					SigPySigPz{};
		double					SigPz2{};
		double					SigE2{};
		double					dth_dpx{};
		double					dth_dpy{};
		double					dth_dpz{};
		double					dphi_dpx{};
		double					dphi_dpy{};
		double					JetResE{};
		double					JetResTheta{};
		double					JetResPhi{};
		float					delta_theta[4]{};
		float					delta_phi[4]{};
		int					bestjet_th{};
		int					bestjet_phi{};
		int					bestjet{};

};

#endif
