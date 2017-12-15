#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "gear/BField.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "MassConstraintFitter.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
typedef CLHEP::HepLorentzVector LorentzVector ;
typedef CLHEP::Hep3Vector Vector3D ;

// MarlinKinfit stuff
#include "LeptonFitObject.h"
#include "JetFitObject.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"

// Marlin stuff
#include <marlin/Global.h>
// the event display

// ROOT stuff
#include "TMath.h"
#include "TMatrixD.h"

#include <cstdlib>
#include <cmath>

using namespace lcio;

MassConstraintFitter aMassConstraintFitter;

//////////////////Adapted from GammaGammaCandidateFinder.cc//////////////////////////////
MassConstraintFitter::MassConstraintFitter() : marlin::Processor("MassConstraintFitter") {

  registerProcessorParameter( "Printing" , 
			      "Print certain messages"  ,
			      _printing,
			       (int)5 );

  registerProcessorParameter("RootFile" ,
			     "Name of the output root file" ,
			     m_rootFile,
			     std::string("MassConstraintFitAnalysis.root"));

  registerProcessorParameter("parentPdg",
			     "PDG code of parent particle",
			     _parentPdg,
			     (int) 221);

  registerProcessorParameter("PhotonAngularError",
			     "Photon Angular Error (rad)",
			     _photonAngularError,
			     (double) 0.001);

  registerProcessorParameter("PhotonAngularErrorModel",
			     "Photon Angular Error Model (0 - stochastic) (1 - constant)",
			     _photonAngularErrorModel,
			     (int) 1);
  
  registerProcessorParameter("parentMass",
			     "Parent particle mass [GeV]",
			     _parentMass,
			     (double) 0.547862);
  
  registerProcessorParameter("parentCharge",
			     "Parent particle charge",
			     _parentCharge,
			     (double) 0.0);

  registerProcessorParameter("nDaughters",
			    "Number of daughter decay particles",
			     _nDaughters,
			     (int) 3);

  registerProcessorParameter("nCharged",
			    "Number of daughter decay particles with charge !=0",
			     _nCharged,
			     (int) 2);

  registerProcessorParameter("nNeutral",
			     "Number of daughter decay particles with charge ==0",
			     _nNeutral,
			     (int) 1);
 
  registerProcessorParameter("nNeutralParams",
			     "Number of parameters used in neutral particle paramterization (in FitObjects)",
			     _nNeutralParams,
			     (int) 3);
			     
  registerProcessorParameter("nChargedParams",
			     "Number of parameters used in charged particle parameterization (in FitObjects)",
			     _nChargedParams,
			     (int) 3);

  std::vector<int> daughterChargedPdgs;
  daughterChargedPdgs.push_back(211);
  registerProcessorParameter("daughterChargedPdgs",
			     "PDG codes of daughter decay particles with charge != 0",
			     _daughterChargedPdgs,
			      daughterChargedPdgs);

  std::vector<int> daughterNeutralPdgs;
  daughterNeutralPdgs.push_back(22);
  registerProcessorParameter("daughterNeutralPdgs",
			     "PDG codes of daughter decay particles with charge == 0",
			     _daughterNeutralPdgs,
			      daughterNeutralPdgs);
			     
  std::vector<float> daughterChargedMass;
  daughterChargedMass.push_back(0.13957018);
  registerProcessorParameter("daughterChargedMass",
			     "Mass [GeV] of daughter decay particles with charge != 0",
			     _daughterChargedMass,
			      daughterChargedMass); 

  std::vector<float> daughterNeutralMass;
  daughterNeutralMass.push_back(0.0);
  registerProcessorParameter("daughterNeutralMass",
			     "Mass [GeV] of daugher decay particles with charge == 0",
			     _daughterNeutralMass,
			      daughterNeutralMass);

  std::string inputParticleCollectionName = "PandoraPFOs";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName" , 
			     "Input Particle Collection Name "  ,
			     _inputParticleCollectionName,
			      inputParticleCollectionName);

  std::string inputTrackCollectionName = "MarlinTrkTracks";
  registerInputCollection( LCIO::TRACK,
				"InputTrackCollectionName" ,
				"Input Track Collection Name " ,
				_inputTrackCollectionName,
				inputTrackCollectionName);
  
  registerInputCollection( LCIO::MCPARTICLE,
				"MCParticleCollection" ,
				"Name of the MCParticle input collection" ,
				_mcParticleCollectionName,
				std::string("MCDecayParticles"));

 
  std::string outputParticleCollectionName = "MassConstraintCandidates";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "OutputParticleCollectionName" , 
			     "Output Particle Collection Name "  ,
			     _outputParticleCollectionName,
			     outputParticleCollectionName);

  std::string outputTrackCollectionName = "MassConstraintCandidates";
  registerOutputCollection( LCIO::TRACK,
			    "OutputTrackCollectionName",
			    "Output Particle Collection Name" ,
			    _outputParticleCollectionName,
			    outputTrackCollectionName);

  registerProcessorParameter( "FitProbabilityCut" , 
			      "Minimum fit probability"  ,
			      _fitProbabilityCut,
			      (double)0.001); 

  registerProcessorParameter( "AllowedMassDeviation" ,
			      "cut on measured mass and true parent mass difference",
			      _allowedMassDeviation,
			      (double) 999.0); 

  registerProcessorParameter( "fitter" ,
                              "0 = OPALFitter, 1 = NewFitter, 2 = NewtonFitter",
                              _ifitter,
                              (int)0);

  registerProcessorParameter("FitAnalysis" ,
			      "0 = No TTree , 1 = make TTree,  generates a TTree of fit and measured particle parameters along with pull distributions",
			      _fitAnalysis,
			      (int)1);
  
  registerProcessorParameter("GeneratorAnalysis",
			     "0 = No TTree , 1 = make TTree, includes generator level information in the TTree if available",
			     _genAnalysis,
			     (int)1);

//secondary mass constraints
  registerProcessorParameter("nMassConstraints",
			     "the number of mass constraints being applied",
			     _nMassConstraints,
			     (int)1);
 
  std::vector<float> secondaryMasses;
  registerProcessorParameter("SecondaryMasses",
			     "vector of the masses that are not the parent mass ",
			     _secondaryMasses,
			     secondaryMasses);

  std::vector<float> secondaryNCharged;
  registerProcessorParameter("SecondaryNCharged",
			     "vector of nCharged for each secondary constraint",
			     _secondaryNCharged,
			     secondaryNCharged);

  std::vector<float> secondaryNNeutral;
  registerProcessorParameter("SecondaryNNeutral",
			     "vector of nNeutral for each secondary constraint",
			     _secondaryNNeutral,
			     secondaryNNeutral);
  
  return;

}

//===================================================================================

void MassConstraintFitter::init() {
  if(_printing>1)printParameters(); 

  evtNo=0;
//  rejects = new TH1D("hrejects","Rejected Events",5,1.0,5.0);
  if(_fitAnalysis)
	rootFile = new TFile(m_rootFile.c_str(),"RECREATE");
        rejects = new TH1D("hrejects","Rejected Events",5,0.5,5.5); 
////////////////////////////////////////////////////////
///fit and measured analysis tree init
//_fitAnalysis includes the general analysis from the reconstructed particles
//it produces a ttree that contains the reconstructed (measured)  particles, fit particles, and pull distributions between
//the fit and measured particles
////////////////////////////////////////////////////////
//initialize tree and variables 
  if(_fitAnalysis){
	tree = new TTree("tree", "tree");
  	tree->SetDirectory(rootFile);

  	//Some useful plots
  	tree->Branch("RecoEnergy", &RecoEnergy);
  	tree->Branch("FitEnergy", &FitEnergy);
        tree->Branch("RecoMass", &RecoMass);
  	tree->Branch("FitProbability", &FitProbability );
	tree->Branch("Chisq", &Chisq);
	tree->Branch("EnergyAsymmetryGen", &EnergyAsymmetryGen);
	tree->Branch("evtNo", &evtNo);
	tree->Branch("rejectNo",&rejectNo);

        tree->Branch("mcParticleVec.", &mcParticleVec);

	tree->Branch("measNeutralVec.", &measNeutralVec);
	tree->Branch("measNeutralParamVec.",&measNeutralParamVec);
        tree->Branch("measChargedVec.", &measChargedVec);
//	tree->Branch("measTrack.", &measTrack);
        tree->Branch("fitNeutralVec.", &fitNeutralVec);
	tree->Branch("fitNeutralParamVec.",&fitNeutralParamVec);
        tree->Branch("fitChargedVec.", &fitChargedVec);
//	tree->Branch("fitTrack.",&fitTrack);

	tree->Branch("measNeutral_err.", &measNeutral_err);
	tree->Branch("measCharged_err.", &measCharged_err);
//	tree->Branch("measTrack_err.",&measTrack_err);
	tree->Branch("fitNeutral_err.", &fitNeutral_err);
	tree->Branch("fitCharged_err.", &fitCharged_err);
//	tree->Branch("fitTrack_err.",&fitTrack_err);

	tree->Branch("measNeutralPdg.", "vector<double>", &measNeutralPdg);
	tree->Branch("measChargedPdg.", "vector<double>", &measChargedPdg);

//	tree->Branch("fitmeas_NeutralPulls.", &fitmeas_NeutralPulls);
//	tree->Branch("fitmeas_ChargedPulls.", &fitmeas_ChargedPulls);
	
	tree->Branch("measParent.", "TLorentzVector", &measParent);
	tree->Branch("fitParent.", "TLorentzVector", &fitParent);
	
	tree->Branch("measParent_err.", "vector<double>", &measParent_err);
	tree->Branch("fitParent_err.", "vector<double>", &fitParent_err);

	tree->Branch("fitmeas_NeutralPulls.",&fitmeas_NeutralPulls);
	tree->Branch("fitmeas_ChargedPulls.",&fitmeas_ChargedPulls);

	
   }
	//////////////mc stuff
	if(_genAnalysis){
//	tree->Branch("genNeutral.", "TLorentzVector", &genNeutral);
//	tree->Branch("genCharged.", "TLorentzVector", &genCharged);
	tree->Branch("mcChargedVec.", &mcChargedVec);
	tree->Branch("mcChargedParamVec.", &mcChargedParamVec);
	tree->Branch("mcNeutralVec.", &mcNeutralVec);
	tree->Branch("mcNeutralParamVec.", &mcNeutralParamVec);
	
	tree->Branch("measgen_NeutralPulls.", &measgen_NeutralPulls);
	tree->Branch("measgen_ChargedPulls.", &measgen_ChargedPulls);
	tree->Branch("fitgen_NeutralPulls.", &fitgen_NeutralPulls);
	tree->Branch("fitgen_ChargedPulls." , &fitgen_ChargedPulls);

	tree->Branch("genParent.", "TLorentzVector", &genParent);
	
	tree->Branch("genNeutralPdg.", "vector<double>", &genNeutralPdg);
	tree->Branch("genChargedPdg.", "vector<double>", &genChargedPdg);
	
    }


   ////////////////////////////////////////////////////////////////////////
  //generator tree init
  //the generator analysis houses the measured particles, the monte carlo particles
  //and pull distributions between the measured/fit particles and monte carlo particles
  //only properly matched particles are added to the tree, match is governed
  //by cuts on the fractional error between the measured and mc particle
  /////////////////////////////////////////////////////////////////////
   //initialize generator tree 
   
    
 
  return;
}

//===================================================================================

void MassConstraintFitter::processRunHeader( LCRunHeader* run) {
  streamlog_out(MESSAGE) << " processRunHeader "  << run->getRunNumber() << std::endl ; 
}

//===================================================================================

void MassConstraintFitter::processEvent( LCEvent * evt ) { 

  // Make a new vector of particles
  LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  // Access PFO collection
  bool found = this->FindPFOs(evt);
  bool found2 = this->FindTracks(evt);
  bool found3 = 0;

 if(_genAnalysis){
   found3 = this->FindMCParticles(evt);
 }

  if((found && found2) || (found && found && _genAnalysis && found3)){
    if(_printing>1)std::cout << "Analysis of resonance with PDG code: " << _parentPdg << std::endl; 
    this->FindMassConstraintCandidates(recparcol);
  } 
  
  // Add new collection to event
  evt->addCollection( recparcol , _outputParticleCollectionName.c_str() );
  
  return;
  
}
//===================================================================================
void MassConstraintFitter::end(){
 if(_fitAnalysis){
// 	rootFile->cd();
rootFile->Write();
//tree->Write();
   //     genTree->Write(); 
 }
  return;
}
//===================================================================================
bool MassConstraintFitter::FindPFOs( LCEvent* evt ) {

  bool tf = false;

  // clear old vector
  _pfovec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){
 // std::cout<<"collection item: "<< *itname << " with " << _inputParticleCollectionName<<std::endl;    
    if(*itname==_inputParticleCollectionName){ 
//	std::cout<< "found matching collection name for photon" <<std::endl;
      LCCollection* col = evt->getCollection(*itname);
      unsigned int nelem = col->getNumberOfElements();
//	std::cout<<" number of elements "<<nelem<<std::endl;
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	_pfovec.push_back(recoPart);
      }
    }
  }

  if(_printing>1)std::cout << "Find PFOs : " << tf << std::endl; 

  return tf;
}
//===================================================================================
bool MassConstraintFitter::FindTracks( LCEvent* evt ) {

  bool tf = false;

  // clear old vector
  _trackvec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){    
    if(*itname==_inputTrackCollectionName){
      LCCollection* col = evt->getCollection(*itname);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	Track* track = dynamic_cast<Track*>(col->getElementAt(i));
	_trackvec.push_back(track);
      }
    }
  }

  if(_printing>1)std::cout << "FindTracks : " << tf << std::endl; 

  return tf;
}
//===================================================================================  
bool MassConstraintFitter::FindMCParticles( LCEvent* evt ){
   
  bool tf = false;

  // clear old vector
  _mcpartvec.clear();

  std::vector<double> PhotonEnergies;

  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){    
    if(*itname==_mcParticleCollectionName){
      LCCollection* col = evt->getCollection(*itname);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	MCParticle* mcPart = dynamic_cast<MCParticle*>(col->getElementAt(i));
	_mcpartvec.push_back(mcPart);
 // Add all filtered mcp's to tree - independent of matching
        TLorentzVector mcp;
        mcp.SetPxPyPzE( mcPart->getMomentum()[0], mcPart->getMomentum()[1], mcPart->getMomentum()[2], mcPart->getEnergy() );
        mcParticleVec.push_back(build4vec(mcp));
        if(mcPart->getPDG() == 22) PhotonEnergies.push_back(mcPart->getEnergy());
      }
    }
  }

  if(_printing>1)std::cout << "Find MCParticles : " << tf << std::endl;
  if(_printing>1)std::cout << "Generator photons found : " << PhotonEnergies.size() << std::endl; 

// Calculate variable that is close to costh*.    (actually asymmetry = beta x costh* )
  EnergyAsymmetryGen = -2.0;

// Only set if there are two only two photons
  if(PhotonEnergies.size() == 2){
     double e1 = PhotonEnergies[0];
     double e2 = PhotonEnergies[1];
     EnergyAsymmetryGen = fabs(e1-e2)/(e1+e2);
  }
  if(_printing>1)std::cout << "Energy asymmetry : " << EnergyAsymmetryGen << std::endl; 

  PhotonEnergies.clear();
 


  return tf;
}

int MassConstraintFitter::getCorrespondingMCParticleIndex(TLorentzVector rec){
	// if(_mcpartvec.size() > 1) return -1 ;
	std::cout<<"size "<<_mcpartvec.size()<<std::endl;
        if(_mcpartvec.size() == 0) return -1;
        int closest_match_index=-1;
        double theta_residual=-1;
        double phi_residual=-1;
        double e_residual=0;
        double tempresidual1=-1;
        double tempresidual2=-1;
        double tempresidual3=0;//was going to be residual for energy but dont include because no distinguishing between neutral and charged particles
       TLorentzVector mc;
      
        for(unsigned int i=0; i<_mcpartvec.size(); i++){
                //compare angles
                mc.SetPxPyPzE(_mcpartvec[i]->getMomentum()[0],_mcpartvec[i]->getMomentum()[1],_mcpartvec[i]->getMomentum()[2],_mcpartvec[i]->getEnergy());
		if((mc.M()-rec.M()) < 1){//mass check, look at difference incase of rounding. Is it within 5%? dont divide because photons
                	if(closest_match_index==-1){
                        	closest_match_index = 0;
                        	theta_residual = fabs(rec.Theta() - mc.Theta());
                        	phi_residual = fabs(rec.Phi() - mc.Phi()); //currently between -pi,pi might need to change for 0,2pi
                    	        e_residual = fabs(rec.E() - mc.E());
                	}
                	tempresidual1 = fabs(rec.Theta() - mc.Theta());
                	tempresidual2 = fabs(rec.Phi() - mc.Phi());//dont forget to change this line also
             	//   tempresidual3 = fabs(rec.E() - mc.E());
                	if((tempresidual1+tempresidual2+tempresidual3) < (theta_residual+phi_residual+e_residual) ){
                	        closest_match_index=i;
                	        theta_residual = tempresidual1;
                	        phi_residual = tempresidual2;
                	        e_residual = tempresidual3;
                	}
		}

        }
	if(closest_match_index != -1){
		mc.SetPxPyPzE(_mcpartvec[closest_match_index]->getMomentum()[0],_mcpartvec[closest_match_index]->getMomentum()[1],_mcpartvec[closest_match_index]->getMomentum()[2],_mcpartvec[closest_match_index]->getEnergy());
        	if(_printing>3){
			std::cout<<"MC Match: "<<std::endl;
			std::cout<<"Reco (E,theta,phi) "<<rec.E()<<", "<<rec.Theta()<<", "<<rec.Phi()<<" "<<std::endl;
			std::cout<<"MC   (E,theta,phi) "<<mc.E()<<", "<<mc.Theta()<<", "<<mc.Phi()<<" " <<std::endl;
		}
	}
	else{
		if(_printing>3){
		std::cout<<"Particle not matched "<<std::endl;
		std::cout<<"Reco (E,theta,phi)"<<rec.E()<<", "<<rec.Theta()<<", "<<rec.Phi()<<" "<<std::endl;
		}
	}
        return closest_match_index;

}

//===================================================================================
//returns error array with parameterization {dE, dTheta, dPhi}
//error for photon energy is sigma/E = 0.18/sqrt(E)
//error on photon angle is sigma/E = 0.001/sqrt(E)
//i need to get errors from the recoparticle error matrix probably not make my own?
std::vector<double> MassConstraintFitter::getNeutralErrors(TLorentzVector pneutral, ReconstructedParticle* pNeutral ){

	std::vector<double> errors;/// = new double[3];
	if(pNeutral->getType() == 22){
	 	errors.push_back(0.18*std::sqrt(pneutral.E()) );
                if(_photonAngularErrorModel == 1){
        	   errors.push_back( _photonAngularError );
		   errors.push_back( _photonAngularError );
                }
                else{
        	   errors.push_back( _photonAngularError/std::sqrt(pneutral.E()) );
		   errors.push_back( _photonAngularError/std::sqrt(pneutral.E()) );
                }
	}
	else{
		errors.push_back(0.55*std::sqrt(pneutral.E()) );
		errors.push_back(0.001/std::sqrt(pneutral.E()) );
		errors.push_back(0.001/std::sqrt(pneutral.E()) );
	}
	return errors;

}
//===================================================================================
//returns error array with parameterization {dk, dTheta, dPhi} with the standard parameterization used in the code
//errors come from track record with (d0,z0,omega,tanLambda,theta, lamda)
std::vector<double> MassConstraintFitter::getChargedParticleErrors(TLorentzVector pcharged, Track* ptrk){
	std::vector<double> errors; 
//	std::cout<<"PRINTING FULL COV"<<std::endl;
//	FloatVec cov = ptrk->getCovMatrix();
//	for(int i=0; i<15; i++){
//		std::cout<<cov[i]<<" ";
//	}
//	std::cout<<std::endl;
//track matrix is lower diagonal
//	errors.push_back(std::sqrt(ptrk->getCovMatrix()[0]));//d0 
//	errors.push_back(std::sqrt(ptrk->getCovMatrix()[2]));//phi
//	errors.push_back(std::sqrt(ptrk->getCovMatrix()[5]));//ome
// errors.push_back(std::sqrt(ptrk->getCovMatrix()[9]));//z0
//	errors.push_back(std::sqrt(ptrk->getCovMatrix()[14]));//tan

	//dk = 1/pt*dOmega/Omega 
	errors.push_back( fabs( 1/pcharged.Perp()*(std::sqrt(ptrk->getCovMatrix()[5])/ptrk->getOmega()) ) );
	//dTheta=cos^2(lamda)*dtanLamda
	errors.push_back( std::sqrt(ptrk->getCovMatrix()[14])/(ptrk->getTanLambda()*ptrk->getTanLambda() + 1) );
	errors.push_back( std::sqrt(ptrk->getCovMatrix()[2]) );//dphi

	return errors;
}
std::vector<double> MassConstraintFitter::getTrackErrors(Track* ptrk){
       std::vector<double> errors;
//      std::cout<<"PRINTING FULL COV"<<std::endl;
//      //      FloatVec cov = ptrk->getCovMatrix();
//      //      for(int i=0; i<15; i++){
//      //              std::cout<<cov[i]<<" ";
//      //      }
//      //      std::cout<<std::endl;
//      //track matrix is lower diagonal
            errors.push_back(std::sqrt(ptrk->getCovMatrix()[0]));//d0 
            errors.push_back(std::sqrt(ptrk->getCovMatrix()[2]));//phi
            errors.push_back(std::sqrt(ptrk->getCovMatrix()[5]));//ome
             errors.push_back(std::sqrt(ptrk->getCovMatrix()[9]));//z0
            errors.push_back(std::sqrt(ptrk->getCovMatrix()[14]));//tan

return errors;
}
//===================================================================================
void MassConstraintFitter::PrintCov(FloatVec cov, int dim){
	int index = 0;
	for(int i=0; i<dim; i++){
		std::cout<< " | " ;
		for(int j=0; j<dim; j++){
				std::cout<<cov[index]<<" ";
				index++;
		}
		std::cout<<" | " << std::endl;
	}
	std::cout<< std::endl;
}
//===================================================================================
void MassConstraintFitter::PrintCov(double* cov, int dim){
	int index = 0;
	for(int i=0; i<dim; i++){
		std::cout<< " | " ;
		for(int j=0; j<dim; j++){
				std::cout<<cov[index]<<" ";
				index++;
		}
		std::cout<<" | " << std::endl;
	}
	std::cout<<std::endl;
}
//===================================================================================
//take in the fit covariance matrix and convert it to global fit error vectors for readability
//*updated now takes in cov, and based on #of particles/neutrals/charged builds the error data structures for fit particles
//it goes charged and then neutrals
void MassConstraintFitter::setFitErrors(double* cov, int dim){
//using the big fit cov matrix, extract the diagonall and put the fit errors onto global vectors
	//diagonals are every 10 indices so add 10
	std::vector<double> temp_err;
	unsigned int chargeCounter=0;
	for(int i=0; i<dim*dim; i = i+1+(_nCharged*_nChargedParams+_nNeutral*_nNeutralParams)){
		temp_err.push_back(std::sqrt(cov[i]));
		if( temp_err.size() == 3){//THIS METHOD IS WRONG  for parameterizations other that k,thetaphi for charged and E,,theta,phi. it assumes both neutral and charged are 3 parameters
			if(chargeCounter < TrackFO.size()){
				fitCharged_err.push_back(temp_err);
				chargeCounter++;
				temp_err.clear();
			}
			else{
				fitNeutral_err.push_back(temp_err);
				temp_err.clear();
			}
		}
	}
 

}
//===================================================================================
//manually constructs the  covariance matrix for parent particle, meant to be transformed into 4x4 4vector form
//***********This is to be used after the global measured errors vectors are  populated**************************
double* MassConstraintFitter::ConstructParentMeasCovMatrix(){
//concat vectors	
	std::vector<double> params;  //THIS may be missing START/END param

	for(unsigned int i=0; i< measCharged_err.size(); i++){
		params.insert(params.end(), measCharged_err.at(i).begin(), measCharged_err.at(i).end());
	}
	for(unsigned int i=0; i< measNeutral_err.size(); i++){
		params.insert(params.end(), measNeutral_err.at(i).begin(), measNeutral_err.at(i).end());
	}
	unsigned int dim = params.size();

	double* meascov = new double[dim*dim];
	int index = 0;
	for(unsigned int i=0; i<dim; i++){
		for(unsigned int j=0; j<dim; j++){
			if(i==j){
				meascov[index]= params[i]*params[j];
			}
			else{
				meascov[index] = 0;
			}
			index ++;
		}
	}
//	std::cout<<"DEBUGGING PRINT ENTIRE MATRIX "<< std::endl;
//	for(int i= 0; i<dim*dim; i++){
//		std::cout<< meascov[i] << " ";
//		if(i+1 % 9 == 0) std::cout<<std::endl; 
//	}

return meascov;	
}
//===================================================================================
//take in the fit and measured covariance matrices (4x4) 4 vector form, and push diagonals onto vectors to be stored in tree
//convention is [dPx, dPy, dPz, dE]
void MassConstraintFitter::setParentErrors(FloatVec meascov, FloatVec fitcov){

	measParent_err.push_back( std::sqrt(meascov[0]) );
	measParent_err.push_back( std::sqrt(meascov[5]) );
	measParent_err.push_back( std::sqrt(meascov[10]) );
	measParent_err.push_back( std::sqrt(meascov[15]) );
	
	fitParent_err.push_back( std::sqrt(fitcov[0]) );
	fitParent_err.push_back( std::sqrt(fitcov[5]) );
	fitParent_err.push_back( std::sqrt(fitcov[10]) );
	fitParent_err.push_back( std::sqrt(fitcov[15]) );
}
//input 6 parameter track stuff// p is 6 params on a vectordouble
//std::vector<double> MassConstraintFitter::ConstructChargedSubMatrix(std::vector<double> trackparams, TLorentzVector p){

std::vector<double> MassConstraintFitter::ConstructChargedSubMatrix(TLorentzVector p){
	std::vector<double> submat;
	//submatrix jacobian is this (derivatives of parent wrt to charged parameter


	/*  
		dpx/dd0 dpx/dphi dpx/dome dpx/dz0 dpx/dtanl dpx/dstart
		dpy/dd0 dpy/dphi dpy/dome dpy/dz0 dpy/dtanl dpy/dstart
		dpz/dd0 dpz/dphi dpz/dome dpz/dz0 dpz/dtanl dpz/dstart
		de/dd0 de/dphi de/dome de/dz0 de/dtanl de/dstart

	*/
	/*
	double lambda = atan(trackparams.at(4));
	double phi = trackparams.at(1);
	double omega = trackparams.at(2);
	double P = ptlv.P();

	const double c = 2.99792458e8; // m*s^-1
  	const double B = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();;          
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = B*c*mm2m*eV2GeV;

	submat.push_back(0); //dpx/dd0
	submat.push_back(P*cos(lambda)*cos(phi)); //dpx/dphi ...
 	submat.push_back(-eB*cos(phi)/(omega*omega)); //dpx/dome
	submat.push_back(0); //dpx/dz0
	submat.push_back(-P*cos(phi)*tan(lambda)/pow(1+(tan(lambda)*tan(lambda)),3.0/2.0)); //dpx/dtanl
	submat.push_back(0);// dpx/dstart
	submat.push_back(0); //dpy/dd0
	submat.push_back(P*cos(lambda)*cos(phi)); //dpy/dphi
	submat.push_back(-eB*sin(phi)/(omega*omega)); //dpy/dome
	submat.push_back(0);//dpy/dz0
	submat.push_back(-P*sin(phi)*tan(lambda)/pow(1+(tan(lambda)*tan(lambda)),3.0/2.0)); //dpy/dtanl
	submat.push_back(0);//dpy/dstart
	submat.push_back(0);//dpz/dd0'
	submat.push_back(0);//dpz/dphi
	submat.push_back(-eB*tan(lambda)/(omega*omega)); //dpz/dome
	submat.push_back(0);//dpz/dz0
	submat.push_back(1/pow(1+(tan(lambda)*tan(lambda)),3.0/2.0));//dpz/dtanl
	submat.push_back(0);//dpz/dstart
	submat.push_back(0);//de/dd0
	submat.push_back(0);//de/dphi
	submat.push_back(0);//de/dome ?? is this right ?? probably not
	submat.push_back(0);//de/dz0
	submat.push_back(0);//de/dtanl
	submat.push_back(0);//de/dstart
	*/
	//submatrix jacobian is this (derivatives of parent wrt to charged parameter
	/* 
	  dPx/dk dPx/dtheta dPx/dphi
	  dPy/dk dPy/dtheta dPy/dphi
	  dPz/dk dPz/dtheta dPz/dphi
	  dE/dk dE/dtheta dE/dphi
			   	*/
	submat.push_back(-p.Perp() *  p.Px());
	submat.push_back(p.Pz()*p.Px()/p.Perp());
	submat.push_back(-p.Py());
	submat.push_back(-p.Perp() * p.Py());
	submat.push_back(p.Pz()*p.Py()/p.Perp());
	submat.push_back(p.Px());
	submat.push_back(-p.Perp()*p.P());
	submat.push_back(-p.Perp());
	submat.push_back(0.0);
	submat.push_back(1.0);
	submat.push_back(0.0);
	submat.push_back(0.0);	

	return submat;
}
std::vector<double> MassConstraintFitter::ConstructNeutralSubMatrix(TLorentzVector p){
	std::vector<double> submat;
	//submatrix jacobian is this (derivatives of parent(Px,Py,Pz,E) wrt to neutral parameter)
	/*
	  dPx/dE dPx/dtheta dPx/dphi
	  dPy/dE dPy/dtheta dPy/dphi
	  dPz/dE dPz/dtheta dPz/dphi
	  dEp/dEc dEp/dtheta dEp/dphi
	*/ 
	submat.push_back(p.Px()/p.P());
	submat.push_back(p.Pz()*p.Px()/p.Perp());
	submat.push_back(-p.Py());
	submat.push_back(p.Py()/p.P());
	submat.push_back(p.Pz()*p.Py()/p.Perp());
	submat.push_back(p.Px());
	submat.push_back(p.Pz()/p.P());
	submat.push_back(-p.Perp());
	submat.push_back(0.0);
	submat.push_back(1.0);
	submat.push_back(0.0);
	submat.push_back(0.0);

	return submat;
	
}
double* MassConstraintFitter::ConcatSubMatrices(std::vector<std::vector<double> > matrices ){
	//this is the correctly arranged jacobian made in a vector and then copied to the desired type double*
	std::vector<double> jacobian;
	std::vector<std::vector<double>::iterator> its;
	for(unsigned int i=0; i<matrices.size(); i++){
		std::vector<double>::iterator it = matrices.at(i).begin();
		its.push_back(it);
	}
/*	std::cout<<"matrices print out "<<std::endl;
	for(int i=0; i<matrices.size(); i++){

		for(int j=0; j<matrices.at(i).size(); j++){			std::cout<<matrices.at(i).at(j)<<" ";
		}
		std::cout<<" | "<<std::endl;
	}
std::cout<<"concat check 1"<<std::endl;	*/
	while(its.at(_nCharged+_nNeutral-1) < matrices.at(_nCharged+_nNeutral-1).end()){

		for(int i=0; i<_nCharged; i++){
			if(its.at(i) < matrices.at(i).end()){
				jacobian.insert(jacobian.end(), its.at(i), its.at(i)+_nChargedParams);
				its.at(i) = its.at(i) + _nChargedParams;
			}
		}
		for(int i=_nCharged; i<_nNeutral + _nCharged; i++){
			if(its.at(i) < matrices.at(i).end()){
				jacobian.insert(jacobian.end(), its.at(i), its.at(i)+_nNeutralParams);
				its.at(i) = its.at(i) + _nNeutralParams;
			}
		}

	}
/*	std::cout<<"jacobian printout "<<std::endl;
	for(int i=0; i<jacobian.size(); i++){
		std::cout<<jacobian.at(i)<<" ";
	}
	std::cout<<std::endl;
	std::cout<<"jacobian size "<<jacobian.size()<<std::endl;

	std::cout<<"concat check 2"<<std::endl;*/
	//copy to double* because thats efficient ;)
	double* jacobian2 = new double[jacobian.size()];
	for(unsigned int i=0; i<jacobian.size(); i++){
		jacobian2[i] = jacobian.at(i);
	}

	return jacobian2;
	
}
//===================================================================================
//FloatVec MassConstraintFitter::ConstructCovMatrix(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, double*  cov){
//FloatVec MassConstraintFitter::ConstructCovMatrix(std::vector<std::vector<double> > trackparams, std::vector<TLorentzVector> charged, std::vector<TLorentzVector> neutral, double* cov){
FloatVec MassConstraintFitter::ConstructCovMatrix(std::vector<TLorentzVector> charged, std::vector<TLorentzVector> neutral, double* cov){

      
            FloatVec covP(16,0.0);
      //  (px,py,pz,E)

      //  Follow for now implementation in FourMomentumCovMat based on TMatrixD

      // Note some cases of cov_dim != 6 have been seen ...
             const int nrows  = (_nChargedParams*_nCharged+_nNeutralParams*_nNeutral);  
             const int ncolumns  =4; 


	std::vector<std::vector<double> > matrices;
	for(unsigned int i =0; i< charged.size(); i++){
//		matrices.push_back(ConstructChargedSubMatrix(trackparams.at(i), charged.at(i)));
		matrices.push_back(ConstructChargedSubMatrix(charged.at(i)));
	}
//	std::cout<<"constructed charged sub "<<std::endl;
	for(unsigned int i= 0; i< neutral.size(); i++){	
		matrices.push_back(ConstructNeutralSubMatrix(neutral.at(i)));
	}
//	std::cout<<"constructed neutral sub "<<std::endl;

	double* jacobian = ConcatSubMatrices(matrices);//jacobian is by columns
//	std::cout<<"subs concatenated"<<std::endl;
	
	
      //  Now set up to calculate the new covariance matrix, namely  V' = D^T V D
             TMatrixD Dmatrix(nrows,ncolumns, jacobian, "F");
             TMatrixD Vmatrix(nrows,nrows, cov, "F");                      // May need to have the filling array explicitly dimensioned ??
             TMatrixD Covmatrix(ncolumns,ncolumns);                        // Hopefully this is good enough for initialization ?

             Covmatrix.Mult( TMatrixD( Dmatrix, TMatrixD::kTransposeMult, Vmatrix) ,Dmatrix);   // Is this what we want ?
	   int index=0;
             for(int ii=0;ii<4;ii++){
                 for(int jj=0;jj<4;jj++){
        //            if(_printing>3)std::cout << "4-vector cov. " << ii << " " << jj << " " << Covmatrix(ii,jj) << " Corr: " << Covmatrix(ii,jj)/std::sqrt(Covmatrix(ii,ii)*Covmatrix(jj,jj)) << std::endl;
        		covP[index] = Covmatrix(ii,jj);
			index++;
                 }                 
	     } 

	return covP;
}
int MassConstraintFitter::getIndexOfMatchingIndices(std::vector<int> indices, int index){
	int matchIndex= -1;
	for(unsigned int i =0; i<indices.size(); i++){
		if(indices.at(i) == index) matchIndex=i;
	}

	return matchIndex;	
}
//compare to arrays, do any of the contents match?
bool MassConstraintFitter::vectorIndicesOverlap(std::vector<int> v1, std::vector<int> v2){
	bool overlap = false;
	for(unsigned int i=0; i<v1.size(); i++){
		for(unsigned int j=0; j<v2.size(); j++){
			if(v1.at(i) == v2.at(j)) overlap = true;
		}
	}

	return overlap;
}
//pass in indices and massconstraint vector, check do these containt overlapping subsets of particles?? if they do this combination is invalid, constraints can not share particles
bool MassConstraintFitter::secondaryConstraintCombinationValid(std::vector<int> indices, std::vector<massconstraint*> constraintVector ){
		//does this set contain same particles???
		bool valid= true;
		std::vector<int> n1,n2;
		std::vector<int> c1,c2;
		//compare indices
		for(unsigned int i=0; i<indices.size(); i++){
			for( unsigned int j=i+1; j<indices.size(); j++){
			//	std::cout<<"constraints "<<indices.at(i)<<" "<<indices.at(j)<<std::endl;
				//compare indices subsets
				n1 = constraintVector.at(i)->neutralIndices;
				n2 = constraintVector.at(j)->neutralIndices;
				c1 = constraintVector.at(i)->chargedIndices;
				c2 = constraintVector.at(j)->chargedIndices;
			/*	for(int i=0; i<n1.size(); i++){
					std::cout<<n1.at(i)<<" ";
				}
				std::cout<<" and ";
				for(int i=0; i<n2.size(); i++){
					std::cout<<n2.at(i)<<" ";
				}	
			*///	std::cout<<std::endl;
				if( vectorIndicesOverlap(n1,n2) ) valid = false;
				if( vectorIndicesOverlap(c1,c2) ) valid = false;	
			//	if(vectorIndicesOverlap(n1,n2) ) std::cout<<"invalid"<<std::endl;
			//	else std::cout<<"valid!!"<<std::endl;
			}
//			if(!valid) std::cout<<"not valid"<<std::endl;
//			if(valid) std::Cout<<"valid"<<std::endl;
		}
	return valid;
}
//===================================================================================
//Current implementation only allows one type of fitter because using the base class is not working correctly
//BaseFitter* MassConstraintFitter::setUpFit(TLorentzVector gamma, TLorentzVector p1, TLorentzVector p2, Track* p1Track, Track* p2Track ){
//OPALFitterGSL* MassConstraintFitter::setUpFit(TLorentzVector gamma, TLorentzVector p1, TLorentzVector p2, Track* p1Track, Track* p2Track){
OPALFitterGSL*  MassConstraintFitter::setUpFit(std::vector<int> neutralIndices, std::vector<int> chargedIndices, std::vector<int> massConstraintIndices, std::vector<massconstraint*> constraintVector, std::vector<TLorentzVector> pneutral, std::vector<TLorentzVector> ptrack,std::vector<ReconstructedParticle*> pNeutralVec, std::vector<Track*> pTrackVec){                                       

	//set up mass constraint
	MassConstraint mc( (double)_parentMass );

	//set up secondary mass constraints
	std::vector<MassConstraint> secondaryMCs;
	for(unsigned int i =0; i<massConstraintIndices.size(); i++){
//			std::cout<<"inited a constraint"<<std::endl;
			MassConstraint secmcs( (double)constraintVector.at(massConstraintIndices.at(i))->mass  );
			secondaryMCs.push_back(secmcs);
	}	

	//construct fit objects for gamma, and 2 charged particles
	std::vector<std::vector<double> > neutralErrors;
	std::vector<std::vector<double> > chargedErrors;//k theta phi
        std::vector<std::vector<double> > chargedTrackErrors;//5 track params

//std::cout<<"in setup "<<neutralIndices.size()<<" "<<chargedIndices.size()<<std::endl;	
//Working here today? implement these get error functions
	for(unsigned int i =0; i< neutralIndices.size(); i++){
		neutralErrors.push_back(getNeutralErrors(pneutral.at( neutralIndices.at(i)),pNeutralVec.at( neutralIndices.at(i)) ) );//should i pass recopart or tlv here?
	}
	for(unsigned int i=0; i < chargedIndices.size(); i++){
		chargedErrors.push_back(getChargedParticleErrors(ptrack.at( chargedIndices.at(i)), pTrackVec.at( chargedIndices.at(i)) ) );
		chargedTrackErrors.push_back(getTrackErrors(pTrackVec.at(chargedIndices.at(i))));
	}
//std::cout<<"did errors "<<std::endl;

	//set up JetFitObjects
	for(unsigned int i=0; i < neutralIndices.size(); i++){
		neutralJets.push_back( new JetFitObject(pneutral.at(neutralIndices.at(i)).E(),pneutral.at(neutralIndices.at(i)).Theta(),pneutral.at(neutralIndices.at(i)).Phi(), neutralErrors.at(i)[0], neutralErrors.at(i)[1], neutralErrors.at(i)[2]));
	}			
//std::cout<<"did jets "<<std::endl;		

	const double B = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();
	for(unsigned int i =0; i < chargedIndices.size(); i++){
		//chargedFO.push_back( new LeptonFitObject( 1/ptrack.at(chargedIndices.at(i)).Perp(), ptrack.at(chargedIndices.at(i)).Theta(), ptrack.at(chargedIndices.at(i)).Phi(), chargedErrors.at(i)[0], chargedErrors.at(i)[1], chargedErrors.at(i)[2]), ptrack.at(chargedIndices.at(i)).M());
		TrackFO.push_back( new LeptonFitObject(pTrackVec.at(chargedIndices.at(i)), B ,ptrack.at(chargedIndices.at(i)).M()  ));
	} 	

	
/*	if(_printing>4)std::cout <<" Fitting Measured Quantities:" <<std::endl;						
	
	if(_printing>4)std::cout <<" Neutral Particles "<<std::endl;
	for(unsigned int i =0; i< neutralIndices.size(); i++){
		if(_printing>4)std::cout << "neutral "<< i <<" E, theta, phi, M, Pdg: " <<  pneutral.at(neutralIndices.at(i)).E()<< " "<< pneutral.at(neutralIndices.at(i)).Theta()<< " "<< pneutral.at(neutralIndices.at(i)).Phi()<< " "<< pneutral.at(neutralIndices.at(i)).M()<< " "<< pNeutralVec.at(neutralIndices.at(i))->getType() << std::endl;
		if(_printing>4)std::cout << "errors "<< i <<" dE, dTheta, dPhi: " << neutralErrors.at(i)[0] << " "<< neutralErrors.at(i)[1] << " "<<neutralErrors.at(i)[2]<< std::endl;
	}	
	if(_printing>4)std::cout <<" Charged Particles "<<std::endl;
	for(unsigned int i=0; i< chargedIndices.size(); i++){
		if(_printing>4)std::cout <<"charged "<< i <<" d0, phi, omega, z0, tanL, M: "<< pTrackVec.at(chargedIndices.at(i))->getD0()<<" "<<pTrackVec.at(chargedIndices.at(i))->getPhi()<<" "<<pTrackVec.at(chargedIndices.at(i))->getOmega()<<" "<<pTrackVec.at(chargedIndices.at(i))->getZ0()<<" "<<pTrackVec.at(chargedIndices.at(i))->getTanLambda()<<" "<<ptrack.at(chargedIndices.at(i)).M() << std::endl;
		if(_printing>4)std::cout<<"errors "<< i << "dd0 dPhi ddome dz0 dtanL ";
		for(unsigned int j=0; j<chargedErrors.at(i).size(); j++){
			if(_printing>4)std::cout<< chargedTrackErrors.at(i).at(j) << " ";
		}
		if(_printing>4) std::cout<<std::endl;
	}*/


	//clean up memory with errors that are not going to be used again in this scope
//	neutralErrors.clear();
//	chargedErrors.clear();
		
//
//	VertexFitObject* vertexFO = new VertexFitObject("commonVtx",0,0,0);
     
 //     vertexFO->setParam(0,0.,false,false); // measured=false, fixed=false
  //    vertexFO->setParam(1,0.,false,false);
  //    vertexFO->setParam(2,0.,false,false);
  //    vertexFO->setError(0,100.);
  //    vertexFO->setError(1,100.);
  //    vertexFO->setError(2,100.);
      
//	for(int i=0; i<TrackFO.size(); i++){
//		vertexFO->addTrack(TrackFO.at(i),false, true);
//	}
		

	
	//add fitobjects to initial mass constraint
	for(unsigned int i= 0; i<TrackFO.size(); i++){
		mc.addToFOList(*TrackFO[i]);
	}
	for(unsigned int i =0; i<neutralJets.size(); i++){
		mc.addToFOList(*neutralJets[i]);
	}


//	std::cout<<"add to mc "<<std::endl;
	//declare fitter
	//include options for other fitters
	OPALFitterGSL *  fitter = new OPALFitterGSL(); //OPALFitterGSL();

//        fitter = new OPALFit

	///////////////////	

	//add constraints and fit objects to fitter
	for(unsigned int i=0; i< TrackFO.size(); i++){
		fitter->addFitObject( TrackFO[i] );
	}
	for(unsigned int i=0; i< neutralJets.size(); i++){
		fitter->addFitObject( neutralJets[i] );
	}

	//add FO to secondary mcs
//	for(int i=0; i< massConstraintIndices.size(); i++){
//		std::cout<<"massConst ind "<<massConstraintIndices.at(i)<<" ";
//	}
//	std::cout<<std::endl;
	//std::cout<<"neut inds"<<std::endl;
//	for(int i=0; i< neutralIndices.size(); i++){
//		std::cout<<neutralIndices.at(i)<<" ";
//	}
	std::cout<<std::endl;
	for(unsigned int i= 0; i<massConstraintIndices.size(); i++){
		for(unsigned int j=0; j<constraintVector.at(massConstraintIndices.at(i))->chargedIndices.size(); j++){
				int trackfoindex;
		//		std::cout<<"here ? c"<<std::endl;
				trackfoindex = getIndexOfMatchingIndices(chargedIndices, constraintVector.at(massConstraintIndices.at(i))->chargedIndices.at(j));
				secondaryMCs.at(i).addToFOList(*TrackFO[trackfoindex]);
						
		}
		for(unsigned int j=0; j<constraintVector.at(massConstraintIndices.at(i))->neutralIndices.size(); j++){
				int neutralfoindex;
		//		std::cout<<"here ? n"<<std::endl;
		//		std::cout<<"looking for "<< constraintVector.at(massConstraintIndices.at(i))->neutralIndices.at(j)<<std::endl;
				neutralfoindex = getIndexOfMatchingIndices(neutralIndices, constraintVector.at(massConstraintIndices.at(i))->neutralIndices.at(j));
		//		std::cout<<"found at "<<neutralfoindex<<std::endl;
				secondaryMCs.at(i).addToFOList(*neutralJets[neutralfoindex]);
		}
		
	}
//	std::cout<<"add secondary mcs?"<<std::endl;
	//
	
//	std::cout<<" add FO "<<std::endl;
//		fitter->addFitObject( vertexFO );

//		vertexFO->addConstraints( *fitter, int(VertexFitObject::VXYZ) );

//		vertexFO->initForFit();
		
	//add the things, what order should i use? charged then neutral?

	fitter->addConstraint(mc);

	//add secondary mcs
	for(unsigned int i=0; i<secondaryMCs.size(); i++){
		fitter->addConstraint(secondaryMCs.at(i) );
	}
//	std::cout<<" add const "<<std::endl;	
	//do the fit (this also returns a fit probability)
	
	double fitp = fitter->fit();
//	std::cout<<" fitted "<<std::endl;
	std::cout<<"Fit Combination [ ";
	for(unsigned int i=0; i<neutralIndices.size(); i++){
		std::cout<<neutralIndices.at(i)<<" ";
	}
	std::cout<<" ] [ ";
	for(unsigned int i=0; i<chargedIndices.size(); i++){
		std::cout<<chargedIndices.at(i)<<" ";
	}
	std::cout<<" ] "<<std::endl;
	if(massConstraintIndices.size()!=0){
		std::cout<<"Secondary constraint Combinations with subset and mass"<<std::endl;
		for(unsigned int i=0; i<massConstraintIndices.size(); i++){
			std::cout<<"secondary constraint "<< i << "[ ";
			for(unsigned int j=0; j<constraintVector.at(massConstraintIndices.at(i))->neutralIndices.size(); j++){
				std::cout<<constraintVector.at(massConstraintIndices.at(i))->neutralIndices.at(j)<<" ";
			}
			std::cout<<" ] [ ";
			for(unsigned int j=0; j<constraintVector.at(massConstraintIndices.at(i))->chargedIndices.size(); j++){
				std::cout<<constraintVector.at(massConstraintIndices.at(i))->chargedIndices.at(j)<<" ";
			}
			std::cout<<" ] "<<constraintVector.at(massConstraintIndices.at(i))->mass<<std::endl;;
		
		}
	}
	std::cout<<"with fit probability "<<fitp<<std::endl;

	//TODO:: look up the function that returns fit prob, dont fit twice if printing	
//	if(_printing>4)std::cout <<" Fit Probability: "<< fitter->fit()<<std::endl;
	neutralErrors.clear();
	chargedErrors.clear();
	chargedTrackErrors.clear();
	
	return fitter;

}
void MassConstraintFitter::generateSubsets(std::vector<int> v, int k, unsigned int start, int currLen, std::vector<bool> used, std::vector<std::vector<int> >& combinations){
	if(currLen == k) {
		std::vector<int> indices;
		for(unsigned int i=0; i<v.size(); i++){
			if(used.at(i) == true) {
				indices.push_back(i);
			}
		}
		combinations.push_back(indices);
		return;
	}
	if (start == v.size()){
		return;
	}
	used[start] = true;
	generateSubsets(v,k,start+1,currLen+1,used, combinations);

	used[start] = false;
	generateSubsets(v,k,start+1,currLen,used, combinations);
}
//n choose k
void MassConstraintFitter::generateIndicesCombinations(int vectorsize, int nparticles, std::vector<std::vector<int> >& combinations){
	int n = vectorsize;
	int k = nparticles;
//	std::cout<<"in gen combs"<<std::endl;
		std::vector<int> master;
		std::vector<bool> used(n);
	//	std::vector<std::vector<int> > combinations;
		for(int i=0; i<n; i++){
			master.push_back(i);
			used.at(i)=false;
		}
//	for(int i=0; i<master.size(); i++){
//		std::cout<<"in gencombs master "<<master.at(i)<<" used "<<used.at(i)<<std::endl;
//	}

	generateSubsets(master,k,0,0,used,combinations);
//	std::cout<<"in combin "<<combinations.size()<<std::endl;
//	for(int i=0; i<combinations.size();i++){
//		for(int j=0; j<combinations.at(i).size();j++){
//		std::cout<<combinations.at(i).at(j)<<" ";
//		}
//		std::cout<<std::endl;
//	}
	
//	generateSubsets(master,k,0,0,used,combinations);
	
//	return combinations;
}
//clone the old RP and update the momentum/Energy from the fit
//these construction methods need the fit errrors too
ReconstructedParticleImpl* MassConstraintFitter::constructFitParticle(TLorentzVector fitp, ReconstructedParticle* measrp){
	ReconstructedParticleImpl* fitrp = new ReconstructedParticleImpl();
	
	fitrp =(ReconstructedParticleImpl*) measrp->clone();
	fitrp->setEnergy(fitp.E());
	double* mom = new double[3];
	mom[0]=fitp.Px();
	mom[1]=fitp.Py();
	mom[2]=fitp.Pz();
	fitrp->setMomentum(mom);

	return fitrp;
	
}
//TODO ADD NEW FIT COV MATRIX
TrackImpl* MassConstraintFitter::constructTrack(TLorentzVector fitp, Track* meast){
	TrackImpl* fitt = new TrackImpl();
	fitt = (TrackImpl*)meast->clone();
	const double c = 2.99792458e8; // m*s^-1
  	const double B = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();;          
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = B*c*mm2m*eV2GeV;
	//update 
	double sinlambda = fitp.Pz()/fitp.P();
	double coslambda = fitp.Px()/(fitp.P()*cos(fitp.Phi()));
	double tanlambda = sinlambda/coslambda;
	double omega = eB/(fitp.P()*coslambda);
	if(meast->getOmega() < 0){
		omega = -omega;
	}
	fitt->setPhi(fitp.Phi());
	fitt->setTanLambda(tanlambda);
	fitt->setOmega(omega);

	return fitt;
}
void MassConstraintFitter::printmassconstraints(std::vector<massconstraint*> cv){
	for(unsigned int i=0; i<cv.size(); i++){
		std::cout<<"Mass "<<cv.at(i)->mass<<" ";
		std::cout<<"N: [ ";
		for(unsigned int j=0; j<cv.at(i)->neutralIndices.size(); j++){
			std::cout<<cv.at(i)->neutralIndices.at(j)<<" ";
		}
		std::cout<<"]  C: [ ";
		for(unsigned int j=0; j<cv.at(i)->chargedIndices.size(); j++){
			std::cout<<cv.at(i)->chargedIndices.at(j)<<" ";
		} 
		std::cout<<"]"<<std::endl;
	}
}
std::vector<massconstraint*> MassConstraintFitter::buildMassConstraint(std::vector<std::vector<int> > neutralIndices, std::vector<std::vector<int> > chargedIndices, double mass, int nNeutral, int nCharged){
	//remeber to implement for use casses no neutral or no charged
	std::vector<massconstraint*> mcs;
	//first case assume both arrays populated
	if(nNeutral > 0 && (nCharged > 0)){
//	std::cout<<"seg here?"<<std::endl;
		for(unsigned int i=0; i<neutralIndices.size(); i++){
			for(unsigned int j=0; j<chargedIndices.size(); j++){
				massconstraint* c  = new massconstraint;
				c->neutralIndices=neutralIndices.at(i);
				c->chargedIndices=chargedIndices.at(j);
				c->mass=mass;
				mcs.push_back(c);
			}		
		}
	}
	//charged only case
	if(nNeutral == 0){
		for(unsigned int i=0; i<chargedIndices.size(); i++){
			massconstraint* c= new massconstraint;
			c->chargedIndices=chargedIndices.at(i);
			c->mass=mass;
			mcs.push_back(c);
		}
	}
	//neutral only case
	if(nCharged == 0){
//		std::cout<<"should be here"<<std::endl;
	//	printCombinations(neutralIndices);
		for(unsigned int i=0; i<neutralIndices.size(); i++){
			massconstraint* c = new massconstraint;
			c->neutralIndices=neutralIndices.at(i);
			c->mass=mass;
			mcs.push_back(c);
		}
	}
//	std::cout<<"made stuff"<<std::endl;
//	printmassconstraints(mcs);	
	return mcs;
}
std::vector<double> MassConstraintFitter::buildTrackVector(Track* t){

	std::vector<double> tvect;
	tvect.push_back(t->getD0());
	tvect.push_back(t->getPhi());
	tvect.push_back(t->getOmega());
	tvect.push_back(t->getZ0());
	tvect.push_back(t->getTanLambda());
	tvect.push_back(0.0);//sStart not measured
	return tvect;
}
std::vector<double> MassConstraintFitter::build4vec(TLorentzVector v){
	std::vector<double> v4;
	v4.push_back(v.Px());
	v4.push_back(v.Py());
	v4.push_back(v.Pz());
	v4.push_back(v.E());

	return v4;
}
std::vector<double> MassConstraintFitter::buildFitTrackVector(TrackParticleFitObject* tfo){
	std::vector<double> tvect;
	// 6 parameters d0 phi omega z0 tanL sStart||sEnd

	for(int i=0; i<6; i++){
		tvect.push_back(tfo->getParam(i));
	}
	return tvect;
}
std::vector<double> MassConstraintFitter::buildLeptonFitVector(LeptonFitObject* lfo){
	std::vector<double> tvect;
	for(int i=0; i<3; i++){
		tvect.push_back(lfo->getParam(i));
	}
	return tvect;

}
std::vector<double> MassConstraintFitter::buildLeptonVector(TLorentzVector v, double q){
	std::vector<double> tvect;
	tvect.push_back(q/v.Perp());
	tvect.push_back(v.Theta());
	tvect.push_back(v.Phi());
	return tvect;
	

}
std::vector<double> MassConstraintFitter::buildNeutralParamVector(TLorentzVector v){
	std::vector<double> tvect;
	tvect.push_back(v.E());
	tvect.push_back(v.Theta());
	tvect.push_back(v.Phi());
	return tvect;
}
void MassConstraintFitter::printCombinations(std::vector<std::vector<int> > combs){
	std::cout<<"combinations "<< combs.size() <<std::endl;
	for(unsigned int i=0; i<combs.size(); i++){
		std::cout<<" [ ";
		for(unsigned int j=0; j<combs.at(i).size(); j++){
			std::cout<< combs.at(i).at(j) << " ";
		}
		std::cout<<" ] "<<std::endl;
	}
}
/*void MassConstraintFitter::printmassconstraints(std::vector<massconstraint*> cv){
	for(int i=0; i<cv.size(); i++){
		std::cout<<"Mass "<<cv.at(i)->mass<<" ";
		std::cout<<"N: [ ";
		for(int j=0; j<cv.at(i)->neutralIndices.size(); j++){
			std::cout<<cv.at(i)->neutralIndices.at(j)<<" ";
		}
		std::cout<<"]  C: [ ";
		for(int j=0; j<cv.at(i)->chargedIndices.size(); j++){
			std::cout<<cv.at(i)->chargedIndices.at(j)<<" ";
		} 
		std::cout<<"]"<<std::endl;
	}
}*/
void MassConstraintFitter::clear(){
	neutralJets.clear();
                TrackFO.clear();

    //            pneutral.clear();
   //     ptrack.clear();
  //      pNeutralVec.clear();
 //       pTrackVec.clear();

        mcParticleVec.clear();

        measNeutral.clear();
        measCharged.clear();
        measTrack.clear();
        fitNeutral.clear();
        fitCharged.clear();
        fitTrack.clear();
        measNeutralVec.clear();
        measNeutralParamVec.clear();
        measChargedVec.clear();
        measTrackVec.clear();
        fitNeutralVec.clear();
        fitNeutralParamVec.clear();
        fitChargedVec.clear();
        fitTrackVec.clear();

        measCharged_err.clear();
        measNeutral_err.clear();
        measParent_err.clear();
        measTrack_err.clear();

        fitCharged_err.clear();
        fitNeutral_err.clear();
        fitTrack_err.clear();

        fitmeas_NeutralPulls.clear();
        fitmeas_ChargedPulls.clear();

        measParent_err.clear();
        fitParent_err.clear();

        measNeutralPdg.clear();
        measChargedPdg.clear();
	 genNeutral.clear();
        genCharged.clear();

        measgen_NeutralPulls.clear();
        measgen_ChargedPulls.clear();
        fitgen_NeutralPulls.clear();
        fitgen_ChargedPulls.clear();

        genNeutralPdg.clear();
        genChargedPdg.clear();

        mcCharged.clear();
        mcNeutral.clear();
       mcNeutralVec.clear();
       mcNeutralParamVec.clear();
       mcChargedVec.clear();
       mcChargedParamVec.clear();



      //DF  delete fitter;
		
}
//===================================================================================
void MassConstraintFitter::FindMassConstraintCandidates(LCCollectionVec * recparcol) {

  rejectNo=0;

  if(_printing>1)std::cout << "FindMassConstraintCandidates : (nPFOs = " << _pfovec.size() << " )" << std::endl; 
  if(_printing>1)std::cout << "FindMassConstraintCandidates : evtNo = " << evtNo << std::endl;
  // Look for candidates
 
  std::vector<TLorentzVector>pneutral;
  std::vector<TLorentzVector>ptrack;

  std::vector<ReconstructedParticle*>pNeutralVec;
  std::vector<Track*>pTrackVec;

//neutral loop, load collection photons into local variables
  for(unsigned int i=0;i<_pfovec.size();i++){

	if(_pfovec[i]->getCharge() == 0 ){
		 if(_printing>1)std::cout << "FindMassConstraintCandidates : " << _pfovec[i]->getType() << std::endl;

	         TLorentzVector tempNeutral(_pfovec[i]->getMomentum()[0], _pfovec[i]->getMomentum()[1], _pfovec[i]->getMomentum()[2], _pfovec[i]->getEnergy() );
		
		if(_printing>1)std::cout << "Neutral candidate (px,py,pz,E) "<<tempNeutral.Px()<<" "<<tempNeutral.Py()<<" "<<tempNeutral.Pz()<<" "<<tempNeutral.E()<<std::endl;
		if(_printing>1)std::cout << "Neutral candidate (E, theta, phi) "<< tempNeutral.E()<<" "<<tempNeutral.Theta()<<" "<<tempNeutral.Phi()<<std::endl;

		/*if(_usePhotonCalibration){//check for ->getType()==22
			calibratePhoton(pGamma);
			if(_printing>1)std::cout <<"g candidate calibration (px,py,pz,E) " <<pGamma.Px()<<" "<<pGamma.Py()<<" "<<pGamma.Pz()<<" "<<pGamma.E()<<std::endl;
			if(_printing>1)std::cout <<"g candidate calibration (E,theta,phi) "<<pGamma.E()<<" "<<pgamma.Theta()<<" "<<pgamma.Phi()<<std::endl;
		}*/

		pneutral.push_back(tempNeutral);
		pNeutralVec.push_back(_pfovec[i]);
	}
  }//end pfovec loop
//start track loop



  const double c = 2.99792458e8; // m*s^-1
  const double B = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();;          
  const double mm2m = 1e-3;
  const double eV2GeV = 1e-9;
  const double eB = B*c*mm2m*eV2GeV;
 //
 for(unsigned int i=0; i<_trackvec.size();i++){
	for(unsigned int j=0; j<_daughterChargedMass.size(); j++){
	
		double cosLambda = 1 / std::sqrt(1 + _trackvec[i]->getTanLambda()*_trackvec[i]->getTanLambda() );
		double P = (eB/fabs(_trackvec[i]->getOmega()))/cosLambda;
		double sinLambda = _trackvec[i]->getTanLambda()*cosLambda;
		double cosPhi = cos(_trackvec[i]->getPhi());
		double sinPhi = sin(_trackvec[i]->getPhi());
		double px = P*cosLambda*cosPhi;
		double py = P*cosLambda*sinPhi;
		double pz = P*sinLambda;
		double E = std::sqrt( _daughterChargedMass[j]*_daughterChargedMass[j] + px*px + py*py + pz*pz );
		

			
			pTrackVec.push_back(_trackvec[i]);

			TLorentzVector trk(px,py,pz,E);

			if(_printing>1)std::cout << "-+ candidate "<< i <<" (px,py,pz,E,M)  "<< px << " " << py << " " << pz << " " << E << " " << trk.M() << std::endl;
			if(_printing>1)std::cout << "-+ candidate "<< i << " (d0,phi,k,z0,tanLambda) " <<_trackvec[i]->getD0()<<" "<<_trackvec[i]->getPhi()<< " "<<_trackvec[i]->getOmega()<< " "<< _trackvec[i]->getZ0()<< " "<<_trackvec[i]->getTanLambda()<<std::endl;
			
			ptrack.push_back(trk);		

	}
	
}  	
  
  if(_printing>1)std::cout << "FindMassConstraintCandidates : (nNeutral = " << pNeutralVec.size() << " " << pneutral.size() << " )" << std::endl;   
  if(_printing>1)std::cout << "FindMassConstraintCandidates : (n+- = "<< pTrackVec.size() << " " << ptrack.size() << " )" <<std::endl;
  // loop over pairs of candidate photons and keep all combinations consistent with the given parent mass assumption

   //size check
    //if( pgamma.size() >= 1 && pplus.size() >= 1 && pminus.size() >= 1){
    if( int(pneutral.size()) >= _nNeutral && int(_trackvec.size()) >=  _nCharged){  //start back up here need to make ijk arrays with indices for both charged and neutral candidates
   std::cout<<"size of pneutral "<<pneutral.size()<<" "<<_nNeutral<<" size of ptrack "<<ptrack.size()<<" "<<_nCharged<<std::endl; 

//store the best fit prob i,j,k then after iterating through all permutations recompute fit and compute analysis vars
	double fitprobmax = -1;
	//double candidate_ijk[3] = {-1,-1,-1};
       std::vector<std::vector<int> > neutralCandidateIndices;
	if(_nNeutral >0){
	 generateIndicesCombinations(pneutral.size(), _nNeutral, neutralCandidateIndices);
	}
	std::vector<std::vector<int> > chargedCandidateIndices; 
	if(_nCharged >0){
	 generateIndicesCombinations(ptrack.size(), _nCharged, chargedCandidateIndices);
	}


	std::vector<massconstraint*> constraintVector;
	std::vector<std::vector<int> > secondaryMassConstraintIndices;
	//build up secondary mass constraints
	/////////////////////////	 
	//loop over mass vector, choose based on params	
//	std::cout<<"n mass constraints "<<_nMassConstraints<<std::endl;
	if(_nMassConstraints > 1){
//	buildmassconstraints(charged indice vector, neutral indice vector, mass)	
		for(unsigned int i =0; i<_secondaryMasses.size(); i++){
//			std::cout<<"in build constraints"<<std::endl;
			std::vector<std::vector<int> > secondaryNeutralIndices;
			std::vector<std::vector<int> > secondaryChargedIndices;
			generateIndicesCombinations(pneutral.size(), _secondaryNNeutral.at(i), secondaryNeutralIndices);
			generateIndicesCombinations(ptrack.size(), _secondaryNCharged.at(i), secondaryChargedIndices);
//			std::cout<<"the combinations inbefore build mc"<<std::endl;
//			printCombinations(secondaryNeutralIndices);
//			printCombinations(secondaryChargedIndices);	
			std::vector<massconstraint*> tempconstraint;
//			std::cout<<"calling build mc"<<std::endl;
			tempconstraint = buildMassConstraint(secondaryNeutralIndices,secondaryChargedIndices, _secondaryMasses.at(i), _secondaryNNeutral.at(i), _secondaryNCharged.at(i));
//			std::cout<<"concating"<<std::endl;
			constraintVector.insert(constraintVector.end(), tempconstraint.begin(), tempconstraint.end() );
//			std::cout<<"concat done"<<std::endl;
			tempconstraint.clear();
			secondaryNeutralIndices.clear();
			secondaryChargedIndices.clear();

		}
//		std::cout<<"size "<<constraintVector.size()<<std::endl;
		//printmassconstraints(constraintVector);
//		std::cout<<"size "<<constraintVector.size()<<std::endl;
//		std::cout<<"now making mc combos"<<std::endl;
	//now generate all combinations of secondary mass constraints
	//	std::vector<std::vector<int> > secondaryMassConstraintIndices;
		generateIndicesCombinations(constraintVector.size(), _nMassConstraints-1, secondaryMassConstraintIndices );	
	}
///	std::cout<<"done building"<<std::endl;
	double fitprob = -1 ;
//test indices
	std::cout<<"neutral Indices"<<std::endl;
	printCombinations(neutralCandidateIndices);
	std::cout<<"charged Indices"<<std::endl;
	printCombinations(chargedCandidateIndices);
	printmassconstraints(constraintVector);
	std::cout<<"Secondary MC Indicies"<<std::endl;
	printCombinations(secondaryMassConstraintIndices);

	std::vector<int> bestCandidatesNeutral;
	std::vector<int> bestCandidatesCharged;
	
	//
	std::vector<int> bestSecondaryConstraints;

	//iterate and fit all particle combinations, save the best one (highest fit probability)
	//assume both indice vectors are non zero size
	double chargeSum= 0.0;
	if( neutralCandidateIndices.size() !=0 && chargedCandidateIndices.size() !=0){

	for(unsigned int i=0; i< neutralCandidateIndices.size(); i++){
		for(unsigned int j=0; j< chargedCandidateIndices.size(); j++){
			//iterate over charged particles and find total charge
			for( unsigned int k=0; k<chargedCandidateIndices.at(j).size(); k++){
				 (pTrackVec[chargedCandidateIndices.at(j).at(k)]->getOmega() < 0 ) ? chargeSum+=-1.0 : chargeSum+=1.0 ;
			}
	//		std::cout<<"here"<<std::endl;
			//if charge is consistent with parent, try fitting
			//TODO check and see if neutral types are consistent?
			if(chargeSum == _parentCharge){
				if(_nMassConstraints > 1){
				//	std::cout<<"in here nmass"<<std::endl;
			//		printCombinations(secondaryMassConstraintIndices);
				//	for(unsigned int l=0; l<_nMassConstraints; l++){
					for(unsigned int l=0; l< secondaryMassConstraintIndices.size(); l++){
						for(unsigned int m=0; m< secondaryMassConstraintIndices.at(l).size(); m++){
						//loop to make sure massconstraint indices ARE A SUBSET of Candidate indices
						bool isSubset=true;
							for(unsigned int n=0; n<constraintVector.at(secondaryMassConstraintIndices.at(l).at(m))->chargedIndices.size() ; n++  ){
								int tempIndex = getIndexOfMatchingIndices(chargedCandidateIndices.at(i), constraintVector.at(secondaryMassConstraintIndices.at(l).at(m))->chargedIndices.at(n));
//								std::cout<<"temp index c"<<tempIndex<<std::endl;
								if(tempIndex == -1) isSubset=false;
							}
							for(unsigned int n=0; n< constraintVector.at(secondaryMassConstraintIndices.at(l).at(m))->neutralIndices.size(); n++  ){
								int tempIndex = getIndexOfMatchingIndices(neutralCandidateIndices.at(i)  , constraintVector.at(secondaryMassConstraintIndices.at(l).at(m))->neutralIndices.at(n));
//								std::cout<<"temp index n"<<tempIndex<<std::endl;
								if(tempIndex == -1) isSubset=false;
							}	
							if(isSubset){// && secondaryConstraintCombinationValid(secondaryMassConstraintIndices.at(l), constraintVector )){
//								std::cout<<"about to fit"<<std::endl;
								OPALFitterGSL*  fitter = setUpFit(neutralCandidateIndices.at(i), chargedCandidateIndices.at(j),secondaryMassConstraintIndices.at(l), constraintVector,  pneutral, ptrack, pNeutralVec, pTrackVec);
								fitprob = fitter->getProbability();
							}
							if(fitprobmax == -1){// && (fitprob > _fitProbabilityCut)){
								fitprobmax = fitprob;
								bestCandidatesNeutral = neutralCandidateIndices.at(i);
								bestCandidatesCharged = chargedCandidateIndices.at(j);
								bestSecondaryConstraints = secondaryMassConstraintIndices.at(l);
							}
							if( (fitprob > fitprobmax )){//{ && (fitprob > _fitProbabilityCut)){
								fitprobmax = fitprob;
								bestCandidatesNeutral = neutralCandidateIndices.at(i);
								bestCandidatesCharged = chargedCandidateIndices.at(j);
								bestSecondaryConstraints = secondaryMassConstraintIndices.at(l);
							}
							neutralJets.clear();
							TrackFO.clear();
//							std::cout<<"about to del "<<std::endl;
							//DF delete fitter;
						}
					}
				}
				else{//there is only 1 constraint
					std::vector<int> noIndices;
					OPALFitterGSL*  fitter = setUpFit(neutralCandidateIndices.at(i), chargedCandidateIndices.at(j),noIndices, constraintVector,  pneutral, ptrack, pNeutralVec, pTrackVec);
					fitprob = fitter->getProbability();

				
			
		//	std::cout<<"probs (max,prob,cut) "<< fitprobmax<< " "<<fitprob<<" "<<_fitProbabilityCut<<std::endl;
					if(fitprobmax == -1){// && (fitprob > _fitProbabilityCut)){
						fitprobmax = fitprob;
						bestCandidatesNeutral = neutralCandidateIndices.at(i);
						bestCandidatesCharged = chargedCandidateIndices.at(j);
			
					}
					if( (fitprob > fitprobmax )){//{ && (fitprob > _fitProbabilityCut)){
						fitprobmax = fitprob;
						bestCandidatesNeutral = neutralCandidateIndices.at(i);
						bestCandidatesCharged = chargedCandidateIndices.at(j);
					}
		
				}
			}//end charge consistency check
			neutralJets.clear();
			TrackFO.clear();
//			std::cout<<"about to del "<<std::endl;
			//DF delete fitter;
//			std::cout<<"after del "<<std::endl;
			chargeSum=0.0;
		}
	}//end condiitonal
	}

///TODO do n constraints for other USE CASES ///////////////////////////////////////////////////////////////////////////
	//if there are no neutrals only loop over charged
	if( neutralCandidateIndices.size() == 0 ){
		for(unsigned int j=0; j< chargedCandidateIndices.size(); j++){
			//iterate over charged particles and find total charge
			
			for( unsigned int k=0; k<chargedCandidateIndices.at(j).size(); k++){
				 (pTrackVec[chargedCandidateIndices.at(j).at(k)]->getOmega() < 0 ) ? chargeSum+=-1.0 : chargeSum+=1.0 ;
			}
	//	
			std::vector<int> noIndices;
			//if charge is consistent with parent, try fitting
			//TODO check and see if neutral types are consistent?
			if(chargeSum == _parentCharge){
				if(_nMassConstraints > 1){
					//for(unsigned int l=0; l<_nMassConstraints; l++){
					for(unsigned int l=0; l<secondaryMassConstraintIndices.size(); l++){
						for(unsigned int m=0; m<secondaryMassConstraintIndices.at(l).size(); m++){
						bool isSubset=true;
					
							for(unsigned int n=0; n<constraintVector.at(secondaryMassConstraintIndices.at(l).at(m))->chargedIndices.size() ; n++  ){
								int tempIndex = getIndexOfMatchingIndices(chargedCandidateIndices.at(j), constraintVector.at(secondaryMassConstraintIndices.at(l).at(m))->chargedIndices.at(n));
								if(tempIndex == -1) isSubset=false;
							}
						//testing if subset has bugs!!! on B+ decay
							if(isSubset){// && secondaryConstraintCombinationValid(secondaryMassConstraintIndices.at(l), constraintVector )){
								OPALFitterGSL*  fitter = setUpFit(noIndices, chargedCandidateIndices.at(j),secondaryMassConstraintIndices.at(l), constraintVector,  pneutral, ptrack, pNeutralVec, pTrackVec);
								fitprob = fitter->getProbability();
							}
							if(fitprobmax == -1){// && (fitprob > _fitProbabilityCut)){
								fitprobmax = fitprob;
								bestCandidatesNeutral = noIndices;
								bestCandidatesCharged = chargedCandidateIndices.at(j);
								bestSecondaryConstraints = secondaryMassConstraintIndices.at(l);
							}	
							if( (fitprob > fitprobmax )){//{ && (fitprob > _fitProbabilityCut)){
								fitprobmax = fitprob;
								bestCandidatesNeutral = noIndices;
								bestCandidatesCharged = chargedCandidateIndices.at(j);
								bestSecondaryConstraints = secondaryMassConstraintIndices.at(l);
							}
							neutralJets.clear();
							TrackFO.clear();
//							std::cout<<"about to del "<<std::endl;
							//DF delete fitter;
						}
					}

				}
				else{
					OPALFitterGSL*  fitter = setUpFit(noIndices, chargedCandidateIndices.at(j),noIndices,constraintVector, pneutral, ptrack, pNeutralVec, pTrackVec);//neutrals at i? what do here
					fitprob = fitter->getProbability();
			
		//	std::cout<<"probs (max,prob,cut) "<< fitprobmax<< " "<<fitprob<<" "<<_fitProbabilityCut<<std::endl;
					if(fitprobmax == -1){// && (fitprob > _fitProbabilityCut)){
						fitprobmax = fitprob;
						bestCandidatesNeutral = noIndices;
						bestCandidatesCharged = chargedCandidateIndices.at(j);
					}
					if( (fitprob > fitprobmax )){//{ && (fitprob > _fitProbabilityCut)){
						fitprobmax = fitprob;
						bestCandidatesNeutral = noIndices;
						bestCandidatesCharged = chargedCandidateIndices.at(j);
					}
				}		

			}//end charge 
				
			neutralJets.clear();
			TrackFO.clear();
//			std::cout<<"about to del "<<std::endl;
			//DF delete fitter;
//			std::cout<<"after del "<<std::endl;
			chargeSum=0.0;
		}
	}	
	//if there are no charged loop over neutrals
	if( chargedCandidateIndices.size() == 0){
			for(unsigned int j=0; j< neutralCandidateIndices.size(); j++){
			//if charge is consistent with parent, try fitting
			//TODO check and see if neutral types are consistent?
			std::vector<int> noIndices;
				if(_nMassConstraints > 1){
				//	for(unsigned int l=0; l<_nMassConstraints; l++){
					for(unsigned int l=0; l< secondaryMassConstraintIndices.size(); l++){
						for(unsigned int m=0; m< secondaryMassConstraintIndices.at(l).size(); m++){
							bool isSubset=true;
							for(unsigned int n=0; n<constraintVector.at(secondaryMassConstraintIndices.at(l).at(m))->chargedIndices.size() ; n++  ){
								int tempIndex = getIndexOfMatchingIndices(chargedCandidateIndices.at(j), constraintVector.at(secondaryMassConstraintIndices.at(l).at(m))->chargedIndices.at(n));
								if(tempIndex == -1) isSubset=false;
							}
						
							if(isSubset){// && secondaryConstraintCombinationValid(secondaryMassConstraintIndices.at(l), constraintVector )){
								OPALFitterGSL*  fitter = setUpFit(neutralCandidateIndices.at(j), noIndices ,secondaryMassConstraintIndices.at(l), constraintVector,  pneutral, ptrack, pNeutralVec, pTrackVec);
								fitprob = fitter->getProbability();
							}
							if(fitprobmax == -1){// && (fitprob > _fitProbabilityCut)){
								fitprobmax = fitprob;
								bestCandidatesNeutral = neutralCandidateIndices.at(j);
								bestCandidatesCharged = noIndices;
								bestSecondaryConstraints = secondaryMassConstraintIndices.at(l);
							}
							if( (fitprob > fitprobmax )){//{ && (fitprob > _fitProbabilityCut)){
								fitprobmax = fitprob;
								bestCandidatesNeutral = neutralCandidateIndices.at(j);
								bestCandidatesCharged = noIndices;
								bestSecondaryConstraints = secondaryMassConstraintIndices.at(l);
							}
							neutralJets.clear();
							TrackFO.clear();
//							std::cout<<"about to del "<<std::endl;
							//DF delete fitter;
						}
					}

				}
				else{	
			//if charge tracks 0 parent is guaranteed to be neutral 
					OPALFitterGSL*  fitter = setUpFit(neutralCandidateIndices.at(j), noIndices,noIndices,constraintVector, pneutral, ptrack, pNeutralVec, pTrackVec);
					fitprob = fitter->getProbability();
		//	}
		//	std::cout<<"probs (max,prob,cut) "<< fitprobmax<< " "<<fitprob<<" "<<_fitProbabilityCut<<std::endl;
					if(fitprobmax == -1){// && (fitprob > _fitProbabilityCut)){
						fitprobmax = fitprob;
						bestCandidatesNeutral = neutralCandidateIndices.at(j);
						bestCandidatesCharged = noIndices;
					}
					if( (fitprob > fitprobmax )){//{ && (fitprob > _fitProbabilityCut)){
						fitprobmax = fitprob;
						bestCandidatesNeutral = neutralCandidateIndices.at(j);
						bestCandidatesCharged = noIndices;
					}
				}
			neutralJets.clear();
			TrackFO.clear();
//			std::cout<<"about to del "<<std::endl;
			//DF delete fitter;
//			std::cout<<"after del "<<std::endl;
			chargeSum=0.0;
			}
	}
     
	std::cout<<"loop completed"<<std::endl;
	//recompute params for best fit probability
//if nothing makes it for this event just return
	if(bestCandidatesNeutral.size()==0 && bestCandidatesCharged.size()==0){
		evtNo++;
		std::cout<<"no particles"<<std::endl; 
		clear();
		rejectNo=5;
		rejects->Fill(rejectNo);
		return;
	}  
		OPALFitterGSL* fitter;
	std::vector<int> noIndices;
	if(_nMassConstraints > 1){
		fitter = setUpFit(bestCandidatesNeutral, bestCandidatesCharged, bestSecondaryConstraints, constraintVector, pneutral, ptrack, pNeutralVec, pTrackVec);
	}
	else{
	//now refit with best candidates
		fitter = setUpFit(bestCandidatesNeutral, bestCandidatesCharged,noIndices,constraintVector, pneutral, ptrack, pNeutralVec, pTrackVec);
	}
	fitprob = fitter->getProbability();	
					
	int cov_dim;//=(_nNeutralParams*_nNeutral + _nChargedParams*_nCharged);
	double * cov = fitter->getGlobalCovarianceMatrix(cov_dim);  //3*(r+k) === 3*(_nNeutral + _nCharged)
	//construct measured cov matrix

//	double * mcov = ConstructParentMeasCovMatrix();
	//if cov matrix dimension is not correct return, fit did not converge
	if(cov_dim == 0){
		evtNo++;
		std::cout<<"no cov"<<std::endl;
		clear();
		rejectNo=2;
		rejects->Fill(rejectNo);
		return;
	}


	
	 
	 //
	 //
	 //
	 //
	 ////////////////////////////
	 
	//build up fit data structures TLVs first
	
	//neutral
	TLorentzVector temp;
	for(unsigned int i=0; i<neutralJets.size() ; i++){
		temp.SetPxPyPzE(neutralJets.at(i)->getPx(),neutralJets.at(i)->getPy(), neutralJets.at(i)->getPz(), neutralJets.at(i)->getE());
//		std::cout<<"fit jet "<<i<<" "<<neutralJets.at(i)->getParam(0)<<" "<<neutralJets.at(i)->getParam(1)<<" "<< neutralJets.at(i)->getParam(2)<<std::endl;
//		std::cout<<"fit jet err "<<i<<" "<<neutralJets.at(i)->getError(0)<<" "<<neutralJets.at(i)->getError(1)<<" "<<neutralJets.at(i)->getError(2)<<std::endl;	
		fitNeutral.push_back(temp);
		fitNeutralVec.push_back(build4vec(temp));
		fitNeutralParamVec.push_back(buildNeutralParamVector(temp));
	}
	for(unsigned int i=0; i<TrackFO.size() ; i++){
		temp.SetPxPyPzE(TrackFO.at(i)->getPx(), TrackFO.at(i)->getPy(), TrackFO.at(i)->getPz(), TrackFO.at(i)->getE());
//	std::cout<<"fit trk "<<i<<" "<<TrackFO.at(i)->getParam(0)<<" "<<TrackFO.at(i)->getParam(1)<<" "<<TrackFO.at(i)->getParam(2)<<" "<<TrackFO.at(i)->getParam(3)<<" "<<TrackFO.at(i)->getParam(4)<<" "<<TrackFO.at(i)->getParam(5)<<std::endl;
//	std::cout<<"fit trk err ";
//	for(int j=0 ;j<6; j++){
//		std::cout<<TrackFO.at(i)->getError(j)<<" ";
//	}
//	std::cout<<std::endl;
//	std::cout<<"fit trk tlv "<<i<<" "<<TrackFO.at(i)->getPx()<<" "<< TrackFO.at(i)->getPy()<<" "<< TrackFO.at(i)->getPz()<< " " << TrackFO.at(i)->getE()<<std::endl;
		fitCharged.push_back(temp);
		fitChargedVec.push_back(build4vec(temp));
		fitTrackVec.push_back(buildLeptonFitVector(TrackFO.at(i)));
		
	}
	fitParent.SetPxPyPzE(0,0,0,0);
	for(unsigned int i=0; i<fitNeutral.size(); i++){
		fitParent += fitNeutral.at(i);
//		std::cout<<"parent e m "<<fitParent.E()<<" "<<fitParent.M()<<std::endl;	
	}
	for(unsigned int i=0; i<fitCharged.size(); i++){
		fitParent += fitCharged.at(i);
//		std::cout<<"parent e m "<<fitParent.E()<<" "<<fitParent.M()<<std::endl;
	}

//std::cout<<"set fit vals "<<std::endl;

	
	//sets an array of charged and neutral errors arrays from cov
//
	setFitErrors(cov, cov_dim);
	std::vector<ReconstructedParticle*> recoNeutraltemp; // be safe and clear these at the end
	//std::vector<Track*> recoTracktemp;
	//put candidates into local LorentzVector for ease of coding and readability
	//do TLVs first then make new RecoParts to get generate errors
	double q;//charge for tlorentzvector for signed 1/pt
	for(unsigned int i=0; i<bestCandidatesCharged.size(); i++){
		measCharged.push_back( ptrack.at(bestCandidatesCharged[i]));
		measChargedVec.push_back(build4vec(ptrack.at(bestCandidatesCharged[i])));
		measTrack.push_back( pTrackVec.at(bestCandidatesCharged[i]));
		//recoTracktemp.push_back(pTrackVec.at(bestCandidatesCharged[i]));
	//	measTrackVec.push_back(buildTrackVector( pTrackVec.at(bestCandidatesCharged[i])));
	        if(pTrackVec.at(bestCandidatesCharged[i])->getOmega() > 0){
			q = 1.0;
		}	
		else{
			q = -1.0;
		}
		measTrackVec.push_back(buildLeptonVector( ptrack.at(bestCandidatesCharged[i]),q));
	}
	for(unsigned int i=0; i<bestCandidatesNeutral.size(); i++){
		measNeutral.push_back( pneutral.at(bestCandidatesNeutral[i]));
		measNeutralVec.push_back(build4vec(pneutral.at(bestCandidatesNeutral[i])));
		measNeutralParamVec.push_back(buildNeutralParamVector(pneutral.at(bestCandidatesNeutral[i])));
		recoNeutraltemp.push_back( pNeutralVec.at(bestCandidatesNeutral[i]));
	}
	//make reco particles
	measParent.SetPxPyPzE(0,0,0,0);
	for(unsigned int i=0; i<measNeutral.size(); i++){
		measParent+= measNeutral.at(i);
//		std::cout<<"parent meas e m "<<measParent.E()<<" "<<measParent.M()<<std::endl;	
	}
	for(unsigned int i=0; i<measCharged.size(); i++){
		measParent+= measCharged.at(i);
//		std::cout<<"parent meas  e m "<<measParent.E()<<" "<<measParent.M()<<std::endl;	
	}

	////////////////////////more quality cuts based on measured particle parameters
	if(fitprob < _fitProbabilityCut){
                evtNo++;
                std::cout<<"fit prob cut not met"<<std::endl;
                clear();
		rejectNo=3;
		rejects->Fill(rejectNo);
                return;
        }
	std::cout<<"mass stuff: "<<fabs(measParent.M()-_parentMass)<<" "<<_allowedMassDeviation<<std::endl;
	std::cout<<"the masses: "<<measParent.M()<<" "<<_parentMass<<std::endl;
	if(fabs(measParent.M()-_parentMass) > _allowedMassDeviation){
		evtNo++;
		std::cout<<"Measured mass deviation too big"<<std::endl;
		clear();
		rejectNo=4;
		rejects->Fill(rejectNo);
		return;	

	}	
	//mass resolution cut would be good
	
	
	////////////////////////////////////////////

//	std::cout<<"set meas vals"<<std::endl;
	//also have error arrays ready locally for readability
	//Photon Errors go {dE,dTheta,dPhi}
	//ChargedParticle errors go {dk, dTheta, dPhi}
	for(unsigned int i=0; i<measCharged.size(); i++){
		measCharged_err.push_back( getChargedParticleErrors(measCharged.at(i), measTrack.at(i)) );
	}
	for(unsigned int i=0; i<measNeutral.size(); i++){
		measNeutral_err.push_back( getNeutralErrors( measNeutral.at(i), recoNeutraltemp.at(i)) );
	}
//	std::cout<<"set meas err "<<std::endl;
	//construct measured cov matrix must be called after measured global vectors are populated
	double * mcov = ConstructParentMeasCovMatrix();

//	std::cout<<"set p mes err matrix "<<std::endl;
	//populate measured  parent errors from covariance matrix
	FloatVec fitParentCov =  ConstructCovMatrix(fitCharged,fitNeutral,cov);
 //	std::cout<<"fist mat constructed "<<std::endl; 
	FloatVec measParentCov = ConstructCovMatrix(measCharged,measNeutral,mcov);
	setParentErrors(measParentCov, fitParentCov);
	
	

 std::cout<<"HERE"<< _printing <<std::endl;
///im to this point
	if(_printing>1){
		std::cout <<"Optimal Candidate Fit Quantities-------------"<<std::endl;
		std::cout << "Constrained fit results  RC: " << fitter->getError() << std::endl;
		std::cout << "Measured mass = " << measParent.M() << std::endl;
		std::cout << "No. of iterations " << fitter->getIterations() << std::endl;
		std::cout << "Fit probability = " << fitprob << std::endl;
		std::cout << "Covariance matrix dimension " << cov_dim << std::endl;
		std::cout << "Candidates Post Fit: "<< std::endl;
		
		for(unsigned int i=0; i< fitCharged.size(); i++){
			
			std::cout <<"Track "<< i <<std::endl;
			std::cout<< "Measured Local Param, M ";
			for(unsigned int j=0; j<measTrackVec.at(i).size(); j++){
				std::cout<<measTrackVec.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			std::cout<< "Fit Local Param, M      ";
			for(unsigned int j=0; j<fitTrackVec.at(i).size(); j++){
				std::cout<<fitTrackVec.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;

			std::cout<< "Measured 4 vec ";
			for(unsigned int j=0; j<measChargedVec.at(i).size(); j++){
				std::cout<<measChargedVec.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			std::cout<< "Fit 4 vec      ";
			for(unsigned int j=0; j<fitChargedVec.at(i).size(); j++){
				std::cout<<fitChargedVec.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"Measured Error  ";
			for(unsigned int j=0; j<measCharged_err.at(i).size(); j++){
				std::cout<<measCharged_err.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"Fit Error       ";
			for(unsigned int j=0; j<fitCharged_err.at(i).size(); j++){
				std::cout<<fitCharged_err.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			
		}
		std::cout<< "measNeutralsize " << measNeutral.size() <<std::endl;
		for(unsigned int i=0; i< measNeutral.size(); i++){
			std::cout <<"Neutral "<< i<<std::endl;
			std::cout <<"Measured Local Param, M ";
			for(unsigned int j=0; j<measNeutralParamVec.at(i).size(); j++){
				std::cout<<measNeutralParamVec.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"Fit Local Param, M       ";
			for(unsigned int j=0; j<fitNeutralParamVec.at(i).size(); j++){
				std::cout<<fitNeutralParamVec.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"Measured 4 vec ";
			for(unsigned int j=0; j<measNeutralVec.at(i).size(); j++){
				std::cout<<measNeutralVec.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"Fit 4 vec      ";
			for(unsigned int j=0; j<fitNeutralVec.at(i).size(); j++){
				std::cout<<fitNeutralVec.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"Measured Error ";
			for(unsigned int j=0; j<measNeutral_err.at(i).size(); j++){
				std::cout<<measNeutral_err.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"Fit Error      ";
			for(unsigned int j=0; j<fitNeutral_err.at(i).size(); j++){
				std::cout<<fitNeutral_err.at(i).at(j)<<" ";
			}
		/*	std::cout << "Neutral "<< i <<" (E,theta,phi,M): "<<std::endl;
			std::cout << " Measured: "<< measNeutral.at(i).E() <<" "<< measNeutral.at(i).Theta() <<" "<<measNeutral.at(i).Phi()<<" "<<measNeutral.at(i).M()<<std::endl;;
			std::cout << " Fit: " << fitNeutral.at(i).E() <<" "<<fitNeutral.at(i).Theta() <<" "<<fitNeutral.at(i).Phi()<<" "<<fitNeutral.at(i).M()<<std::endl;
			std::cout << "Neutral "<< i <<" [px,py,pz,E] "<<std::endl;
			std::cout << " Measured: "<<measNeutral.at(i).Px() <<" "<< measNeutral.at(i).Py() <<" "<<measNeutral.at(i).Pz()<<" "<<measNeutral.at(i).E()<<std::endl;
			std::cout << " Fit: "<<fitNeutral.at(i).Px()<<" "<<fitNeutral.at(i).Py()<<" "<<fitNeutral.at(i).Pz()<<" "<<fitNeutral.at(i).E()<<std::endl;

			std::cout << "Neutral "<< i <<" Errors (dE, dtheta, dphi): "<<std::endl;
			std::cout << "Measured: ";
			for(int j=0; j<measNeutral_err.at(i).size(); j++){
				std::cout<< measNeutral_err.at(i).at(j)<< " ";
			}
			std::cout<<std::endl;
			std::cout <<" Fit: ";
			for(int j=0; j<fitNeutral_err.at(i).size(); j++){
				std::cout<< fitNeutral_err.at(i).at(j)<<" ";
			}
			std::cout<<std::endl;
		*/

		}
//	std::cout << "at parent fit "<<std::endl;
		//parent info
		std::cout << "Parent Measured (px,py,pz,E) "<<measParent.Px()<<" "<<measParent.Py()<<" "<<measParent.Pz()<<" "<<measParent.E()<<std::endl;
		std::cout << "Parent   Fit    (px,py,pz,E) "<<fitParent.Pz()<<" "<<fitParent.Py()<<" "<<fitParent.Pz()<<" "<<fitParent.E()<<std::endl;		
	
		std::cout << "Parent Measured Error ";
                for(unsigned int i=0; i<measParent_err.size(); i++){
                        std::cout<< measParent_err.at(i)<<" ";
                }
                std::cout << std::endl;
                std::cout << "Parent   Fit   Error  ";
                for(unsigned int i=0; i<fitParent_err.size(); i++){
                        std::cout<< fitParent_err.at(i)<<" ";
                }
                std::cout << std::endl; 
		std::cout << "Parent Measured Energy "<< measParent.E() <<" Parent Fit Energy "<< fitParent.E()<<std::endl;
		std::cout <<"----------------------"<<std::endl;

		std::cout<<"Measured Covariance Matrix: "<<std::endl;
		PrintCov(mcov,_nChargedParams*_nCharged+_nNeutralParams*_nNeutral);
		std::cout<<"Fit Covariance Matrix: "<<std::endl;
		PrintCov(cov,_nChargedParams*_nCharged+_nNeutralParams*_nNeutral);
		std::cout<<"Measured Parent 4vector 4x4 Covariance Matrix: "<<std::endl;
		PrintCov(measParentCov,4);
		std::cout<<"Fit Parent 4vector 4x4 Covariance Matrix: "<<std::endl;
		PrintCov(fitParentCov,4);					     
	}//end optimal fit printing
				//////////////IM up to this point, add pull stuff				
			              
	if(cov_dim > 0){//condition required because there is not always convergence in the fitting	
	//	if(fitprob> _fitProbabilityCut){
			//TODO:: add qualifying fits to collection that are less than optimal
			

			if(_parentCharge == 0){
				//parent is neutral store as a reco
				ReconstructedParticleImpl * recoPart = new ReconstructedParticleImpl();
				//adding old particles to fit particle because lazy and dont want to make new RecoParts*
			//	for(int i=0; i<fitNeutral.size(); i++){
					//NOTE: these constiuent particles dont have updated covariance matrices only P,E
			//		recoPart->addParticle(constructFitParticle( fitNeutral.at(i),recoNeutraltemp.at(i)));
			//	}
			//	for(int i=0; i<recoTracktemp.size(); i++){
			//		recoPart->addTrack(constructFitTrack( fitCharged.at(i),rectoTracktemp.at(i)));
			//	}
				recoPart->setEnergy(fitParent.E());
				double mom[3] = {fitParent.Px(), fitParent.Py(), fitParent.Pz() };
				recoPart->setMomentum(mom);
				recoPart->setMass( _parentMass );
				recoPart->setCharge( _parentCharge );
				recoPart->setCovMatrix( fitParentCov );
				recoPart->setType( _parentPdg );
				recoPart->setGoodnessOfPID((float)fitprob);
				if(_printing>1)std::cout << "candidate passes fit probability cut - KEEP" << std::endl;
				recparcol->addElement( recoPart );
				 
			}
			//if the parent is charged store a track, jacobians have not been derived for this yet
			if(_parentCharge != 0){
				//parent is charged store as a track
			//	TrackImpl* recoTrack = TrackImpl();
				//recoTrack = constructFitTrack( fitCharged
				
			}
			
			//add stuff to regular ttree
			if(_fitAnalysis){

 				RecoEnergy=measParent.E();
				FitEnergy=fitParent.E();
				RecoMass=measParent.M();
				FitProbability=fitprob; 
				/*P1_k_FitMeas_pull= (1/p1Fit.Perp() -1/p1meas.Perp())/ std::sqrt( p1meas_err[0]*p1meas_err[0] - cov[0] )  ;
				P1_Theta_FitMeas_pull= (p1Fit.Theta() - p1meas.Theta())/ std::sqrt(p1meas_err[1]*p1meas_err[1] - cov[10]) ;	
				P1_Phi_FitMeas_pull= getPhiResidual(p1Fit.Phi(), p1meas.Phi())/ std::sqrt( p1meas_err[2]*p1meas_err[2] - cov[20]);
		
				P2_k_FitMeas_pull= (1/p2Fit.Perp() - 1/p2meas.Perp())/ std::sqrt( p2meas_err[0]*p2meas_err[0] - cov[30]) ;
				P2_Theta_FitMeas_pull= (p2Fit.Theta() - p2meas.Theta())/ std::sqrt( p2meas_err[1]*p2meas_err[1] - cov[40] );
				P2_Phi_FitMeas_pull= getPhiResidual(p2Fit.Phi(), p2meas.Phi())/ std::sqrt( p2meas_err[2]*p2meas_err[2] - cov[50]);/home/bears/work/research/constrainedFitter/MassConstraintFitter.h
				
		
				Gamma_E_FitMeas_pull= (gammaFit.E() - gammameas.E())/ std::sqrt( gammameas_err[0]*gammameas_err[0] - cov[60] );
				Gamma_Theta_FitMeas_pull= (gammaFit.Theta() - gammameas.Theta())/ std::sqrt( gammameas_err[1]*gammameas_err[1]- cov[70]);
				Gamma_Phi_FitMeas_pull=getPhiResidual(gammaFit.Phi(), gammameas.Phi())/std::sqrt( gammameas_err[2]*gammameas_err[2] - cov[80]) ;*/
	std::cout<<"starint track pull"<<std::endl;
	//these need to change fitmeas_chargePulls.at(i) is unitialized, make a vector of pull vals then push that vector onto pull 2dvector	
				//charged pulls fit and measured
				for(unsigned int i=0; i<measCharged.size(); i++){
					std::vector<double> tempvec;
					for(unsigned int j=0; j<measTrackVec.at(i).size(); j++){
						tempvec.push_back( (measTrackVec.at(i).at(j) - fitTrackVec.at(i).at(j) ) / ( std::sqrt( measCharged_err.at(i).at(j)*measCharged_err.at(i).at(j) - fitCharged_err.at(i).at(j)*fitCharged_err.at(i).at(j)) ) );
					}
					fitmeas_ChargedPulls.push_back(tempvec);
					tempvec.clear();
				}
				//neutral pulls fit andmeasured
				
	std::cout<<"starting neutral pull"<<std::endl;
	std::cout<<"meas neutr size "<< measNeutralVec.size() <<std::endl;
				std::cout<<"meas neutral param"<<std::endl;
				for(unsigned int i=0; i<measNeutralParamVec.size(); i++){
					for(unsigned int j=0; j<measNeutralParamVec.at(i).size(); j++){
					std::cout<<measNeutralParamVec.at(i).at(j)<<" ";
					}
					std::cout<<std::endl;
				}
				std::cout<<"fit neutral param"<<std::endl;
				for(unsigned int i=0; i<fitNeutralParamVec.size(); i++){
					for(unsigned int j=0; j<fitNeutralParamVec.at(i).size(); j++){
					std::cout<<fitNeutralParamVec.at(i).at(j)<<" ";
					}
					std::cout<<std::endl;
				}
				std::cout<<"meas neutral err"<<std::endl;
				for(unsigned int i=0; i<measNeutral_err.size(); i++){
					for(unsigned int j=0; j<measNeutral_err.at(i).size(); j++){
					std::cout<<measNeutral_err.at(i).at(j)<<" ";
					}
					std::cout<<std::endl;
				}
				std::cout<<"fit neutral err"<<std::endl;
				for(unsigned int i=0; i<fitNeutral_err.size(); i++){
					for(unsigned int j=0; j<fitNeutral_err.at(i).size(); j++){
					std::cout<<fitNeutral_err.at(i).at(j)<<" ";
					}
					std::cout<<std::endl;
				}
				for(unsigned int i=0; i<measNeutralParamVec.size(); i++){
				//	for(int j=0; j<measNeutral.at(i).size();j++){
				//		fitmeas_NeutralPulls.push_back( (measNeutral.at(i).at(j)-fitNeutral.at(i).at(j))/ ( std::sqrt( measNeutral_err.at(i).at(j)*measNeutral_err.at(i).at(j) - fitNeutral_err.at(i).at(j)*fitNeutral_err.at(i).at(j)) ) );
				//	}
					std::vector<double>  tempvec;				
			//		std::cout<<"meas err "<<i<<" size "<<measNeutral_err.at(i).size();
			//		std::cout<<"fit err "<<i<<" size "<<fitNeutral_err.at(i).size();				
					for(unsigned int j=0; j<measNeutralParamVec.at(i).size(); j++){
					double number;
				//	tempvec.push_back
					number=((measNeutralParamVec.at(i).at(j) - fitNeutralParamVec.at(i).at(j)) / (std::sqrt(measNeutral_err.at(i).at(j)*measNeutral_err.at(i).at(j)- fitNeutral_err.at(i).at(j)*fitNeutral_err.at(i).at(j)) ) );
 //					std::cout<<"fit meas nuet pull "<<number<<std::endl;
					tempvec.push_back(number);
					}
		//		 	tempvec.push_back( (measNeutral.at(i).E()-fitNeutral.at(i).E())/ (std::sqrt(measNeutral_err.at(i).at(0)*measNeutral_err.at(i).at(0) - fitNeutral_err.at(i).at(0)*fitNeutral_err.at(i).at(0)) ) );
			//		tempvec.push_back( (measNeutral.at(i).Theta()-fitNeutral.at(i).Theta())/ (std::sqrt(measNeutral_err.at(i).at(1)*measNeutral_err.at(i).at(1) -fitNeutral_err.at(i).at(1)*fitNeutral_err.at(i).at(1)) ) );
//					tempvec.push_back( (measNeutral.at(i).Phi()-fitNeutral.at(i).Phi())/ (std::sqrt(measNeutral_err.at(i).at(2)*measNeutral_err.at(i).at(2)- fitNeutral_err.at(i).at(2)*fitNeutral_err.at(i).at(2)) ) );
					
					fitmeas_NeutralPulls.push_back(tempvec);
					tempvec.clear();
					}
				
		//		if(Gamma_Phi_FitMeas_pull > 100){
		//			 std::cout<<"crazy!!"<<std::endl;
		//			std::cout<<gammaFit.Phi()<<" "<<gammameas.Phi()<<" "<<getPhiResidual(gammaFit.Phi(), gammameas.Phi()) <<" "<<gammameas_err[2]<<" "<<cov[20]<<std::endl;
		//		}
		//
//commenting out pulls for now
//	std::cout<<"staring parent pull"<<std::endl;
				fitmeas_ParentPulls.push_back(( fitParent.Px() - measParent.Px())/ std::sqrt(measParent_err[0]*measParent_err[0]-fitParent_err[0]*fitParent_err[0])) ;
				fitmeas_ParentPulls.push_back(( fitParent.Py() - measParent.Py())/ std::sqrt(measParent_err[1]*measParent_err[1]-fitParent_err[1]*fitParent_err[1])) ;
				fitmeas_ParentPulls.push_back(( fitParent.Pz() - measParent.Pz())/ std::sqrt(measParent_err[2]*measParent_err[2]-fitParent_err[2]*fitParent_err[2])) ;
				fitmeas_ParentPulls.push_back(( fitParent.E() - measParent.E())/ std::sqrt(measParent_err[3]*measParent_err[3]-fitParent_err[3]*fitParent_err[3])) ;

				Chisq = fitter->getChi2();
				//tree->Fill();
			}//end fit analysis option
	
			//add stuff to generator ttree
			if(_genAnalysis){

				
			//
			//	std::cout<<"test parent meas vector after tree fill "<<std::endl;
			//	for(int i= 0; i<parentmeas_err.size(); i++){
			//	std::cout<<parentmeas_err[i]<<" ";
			//	}
			//	std::cout<<std::endl;
			//get montecarlo 4 vectors

			//note particle matching is done with mcparticle so only 4vector parameterizations are currently implemented
				//int mcp1index = getCorrespondingMCParticleIndex(p1meas, 1, _daughterPID);
				//int mcp2index = getCorrespondingMCParticleIndex(p2meas, -1, -_daughterPID);
				//int mcpgindex = getCorrespondingMCParticleIndex(gammameas, 0,_neutralPID);
			//these pulls can be made and valid only when each pgetCorrespondingMCParticleIndexarticle has a MC match
				std::vector<int> neutralMCPIndices;
				std::vector<int> chargedMCPIndices;
				std::vector<double> charges;
		//	std::cout<<"test 1"<<std::endl;
				for(unsigned int i=0; i<fitCharged.size(); i++){
					chargedMCPIndices.push_back(getCorrespondingMCParticleIndex(fitCharged.at(i)) );//need signed curvature
				//Note:: this is quick fixed to depend on local parameterization k theta phi
					if(fitTrackVec.at(i).at(0) > 0){
						charges.push_back(1.0);
					}
					else{
						charges.push_back(-1.0);
					}
				}
				for(unsigned int i=0; i<fitNeutral.size(); i++){
					neutralMCPIndices.push_back(getCorrespondingMCParticleIndex(fitNeutral.at(i)) );
				}
		//	std::cout<<"test 2"<<std::endl;
				
				/*if(mcp1index != -1 && mcp2index != -1 && mcpgindex !=-1){
					 mcp1.SetPxPyPzE(_mcpartvec[mcp1index]->getMomentum()[0],
						_mcpartvec[mcp1index]->getMomentum()[1],
						_mcpartvec[mcp1index]->getMomentum()[2],
						_mcpartvec[mcp1index]->getEnergy());
					 mcp2.SetPxPyPzE(_mcpartvec[mcp2index]->getMomentum()[0],
						_mcpartvec[mcp2index]->getMomentum()[1],
						_mcpartvec[mcp2index]->getMomentum()[2],
						_mcpartvec[mcp2index]->getEnergy()); 
					mcgamma.SetPxPyPzE(_mcpartvec[mcpgindex]->getMomentum()[0],
						_mcpartvec[mcpgindex]->getMomentum()[1],
						_mcpartvec[mcpgindex]->getMomentum()[2],
						_mcpartvec[mcpgindex]->getEnergy());
				*/
			//build up mcp 4vetors
				TLorentzVector mcp;
		//		std::cout<<"charged mcp indices "<<std::endl;
		//		for(int i=0; i<chargedMCPIndices.size(); i++){
		//		std::cout<<chargedMCPIndices.at(i);
		//		}
		//		std::cout<<std::endl;
				for(unsigned int i=0; i< chargedMCPIndices.size(); i++){
					//TLorentzvector mcp;
					if( chargedMCPIndices.at(i) != -1){
						
						mcp.SetPxPyPzE(_mcpartvec[chargedMCPIndices.at(i)]->getMomentum()[0],
						_mcpartvec[chargedMCPIndices.at(i)]->getMomentum()[1],
						_mcpartvec[chargedMCPIndices.at(i)]->getMomentum()[2],
						_mcpartvec[chargedMCPIndices.at(i)]->getEnergy());
					}
					else{
						mcp.SetPxPyPzE(0,0,0,0);
					}
					mcCharged.push_back(mcp);

				}
		//	std::cout<<"test 3"<<std::endl;
				for(unsigned int i=0; i< neutralMCPIndices.size(); i++){
				//	TLorentzVector mcp;
					if( neutralMCPIndices.at(i) != -1){
						mcp.SetPxPyPzE(_mcpartvec[neutralMCPIndices.at(i)]->getMomentum()[0],
						_mcpartvec[neutralMCPIndices.at(i)]->getMomentum()[1],
						_mcpartvec[neutralMCPIndices.at(i)]->getMomentum()[2],
						_mcpartvec[neutralMCPIndices.at(i)]->getEnergy());
					}
					else{
						mcp.SetPxPyPzE(0,0,0,0);
					}
					mcNeutral.push_back(mcp);
					
				}	
		//	std::cout<<"test 4"<<std::endl;

			//set up generator data structures
			for(unsigned int i=0; i<mcCharged.size(); i++){
				mcChargedVec.push_back( build4vec(mcCharged.at(i)));
				mcChargedParamVec.push_back( buildLeptonVector(mcCharged.at(i),charges.at(i) ));
			}
			for(unsigned int i=0; i<mcNeutral.size(); i++){
				mcNeutralVec.push_back( build4vec(mcNeutral.at(i)));
				mcNeutralParamVec.push_back( buildNeutralParamVector(mcNeutral.at(i)));
			}
		//	std::cout<<"test 5"<<std::endl;
			//start generator pulls
				for(unsigned int i=0; i<mcNeutralParamVec.size(); i++){	
					std::vector<double>  tempvec;		
					for(unsigned int j=0; j<mcNeutralParamVec.at(i).size(); j++){
					double number;
				//	tempvec.push_back
					number=((measNeutralParamVec.at(i).at(j) - mcNeutralParamVec.at(i).at(j)) / (measNeutral_err.at(i).at(j) ));
 		//			std::cout<<"gen meas nuet pull "<<number<<std::endl;
					tempvec.push_back(number);
					}
		
					
					measgen_NeutralPulls.push_back(tempvec);
					tempvec.clear();
				}
		//	std::cout<<"test 6"<<std::endl;
				
				for(unsigned int i=0; i<mcChargedParamVec.size(); i++){	
					std::vector<double>  tempvec;
						for(unsigned int j=0; j<mcChargedParamVec.at(i).size(); j++){
					double number;
				//	tempvec.push_back
					number=((measTrackVec.at(i).at(j) - mcChargedParamVec.at(i).at(j)) / (measCharged_err.at(i).at(j) ));
 		//			std::cout<<"gen meas nuet pull "<<number<<std::endl;
					tempvec.push_back(number);
					}
		
					
					measgen_ChargedPulls.push_back(tempvec);
					tempvec.clear();
				}

		//	std::cout<<"test 7"<<std::endl;
				//all MCP pulls are going to be 3 parameters E theta phi

				//measgen charged
				//vector<double> temp;
				//for(int i=0; i<measNeutral.size(); i++){
				//	temp.pushback(measNeutral );
	
				//}

				

				/*	P1_k_MeasGen_pull = (1/p1meas.Perp() -  1/mcp1.Perp())/ p1meas_err[0]; 
  					P1_Theta_MeasGen_pull = ( p1meas.Theta() - mcp1.Theta())/ p1meas_err[1];
		  			P1_Phi_MeasGen_pull = getPhiResidual(p1meas.Phi(), mcp1.Phi())/p1meas_err[2] ;
	
  					P2_k_MeasGen_pull = (1/p2meas.Perp() - 1/mcp2.Perp())/ p2meas_err[0];
		  			P2_Theta_MeasGen_pull = ( p2meas.Theta() - mcp2.Theta())/ p2meas_err[1];
  					P2_Phi_MeasGen_pull = getPhiResidual(p2meas.Phi(), mcp2.Phi()) / p2meas_err[2];

					//gama gen err
					gammagen_err = getPhotonErrors(mcgamma);

		  			Gamma_E_MeasGen_pull = ( gammameas.E() - mcgamma.E() )/ gammagen_err[0];
  					Gamma_Theta_MeasGen_pull =  ( gammameas.Theta() - mcgamma.Theta())/ gammagen_err[1];
					Gamma_Phi_MeasGen_pull = getPhiResidual(gammameas.Phi(), mcgamma.Phi()) / gammagen_err[2];
			
                       
			
					P1_k_FitGen_pull = (1/p1Fit.Perp() - 1/mcp1.Perp())/ std::sqrt(cov[0]);
		 	 		P1_Theta_FitGen_pull = (p1Fit.Theta() - mcp1.Theta())/ std::sqrt(cov[10]);
					P1_Phi_FitGen_pull = getPhiResidual(p1Fit.Phi(), mcp1.Phi()) / std::sqrt(cov[20]);			

		  			P2_k_FitGen_pull = (1/p2Fit.Perp() - 1/mcp2.Perp())/ std::sqrt(cov[30]);
  					P2_Theta_FitGen_pull = (p2Fit.Theta() - mcp2.Theta())/ std::sqrt(cov[40]);
					P2_Phi_FitGen_pull = getPhiResidual(p2Fit.Phi(), mcp2.Phi()) / std::sqrt(cov[50]);

  					Gamma_E_FitGen_pull = (gammaFit.E() - mcgamma.E())/ std::sqrt(cov[60]);
		  			Gamma_Theta_FitGen_pull = (gammaFit.Theta() - mcgamma.Theta())/ std::sqrt(cov[70]);
					Gamma_Phi_FitGen_pull = getPhiResidual(gammaFit.Phi(), mcgamma.Phi()) / std::sqrt(cov[80]);	*/		


/*					if(fabs(Gamma_Phi_FitGen_pull)> 100 ||fabs( P2_Phi_FitGen_pull) > 100 || fabs(P1_Phi_FitGen_pull) >100){
						std::cout<<"crazy!! phi from fit gen"<<std::endl;
						std::cout<<"(P1,P2,Gamma) "<<P1_Phi_FitGen_pull<<" "<<P2_Phi_FitGen_pull<<" "<<Gamma_Phi_FitGen_pull<<std::endl;
						std::cout<<"residual (P1,P2,Gamma)"<<getPhiResidual(p1Fit.Phi(), mcp1.Phi())<<" "<<getPhiResidual(p2Fit.Phi(),mcp2.Phi())<<" "<<getPhiResidual(gammaFit.Phi(), mcgamma.Phi())<<std::endl;
						std::cout<<"errors (P1,P2,Gamma)"<<std::sqrt(cov[50])<<" "<<std::sqrt(cov[80])<<" "<<std::sqrt(cov[20])<<std::endl;
					
					}
					if(fabs(Gamma_Phi_MeasGen_pull)> 100 ||fabs(P2_Phi_MeasGen_pull)> 100 || fabs(P1_Phi_MeasGen_pull) >100){
						std::cout<<"crazzy!! phi from meas gen"<<std::endl;
					}*/
					//genTree->Fill();	
                        //	}
 			//	else{
			//		std::cout<<"Could not locate all corresponding MC particles"<<std::endl;
			//	}
			}//end genfit analysis option
              //  }//end prob cut
	}// end cov dim check
     
	tree->Fill();
	clear();
//	std::cout<<"test 8"<<std::endl;
	//memory management
//	neutralCandidateIndices.clear();
//	chargedCandidateIndices.clear();
	//gammagen_err.clear();
/*
		neutralJets.clear();
		TrackFO.clear();

	        pneutral.clear();
        ptrack.clear();
        pNeutralVec.clear();
        pTrackVec.clear();

        measNeutral.clear();
        measCharged.clear();
        measTrack.clear();
        fitNeutral.clear();
        fitCharged.clear();
        fitTrack.clear();
        measNeutralVec.clear();
	measNeutralParamVec.clear();
        measChargedVec.clear();
        measTrackVec.clear();
        fitNeutralVec.clear();
	fitNeutralParamVec.clear();
        fitChargedVec.clear();
	fitTrackVec.clear();

        measCharged_err.clear();
        measNeutral_err.clear();
        measParent_err.clear();
        measTrack_err.clear();

        fitCharged_err.clear();
        fitNeutral_err.clear();
        fitTrack_err.clear();

	fitmeas_NeutralPulls.clear();
	fitmeas_ChargedPulls.clear();

        measParent_err.clear();
        fitParent_err.clear();

        measNeutralPdg.clear();
        measChargedPdg.clear();

        genNeutral.clear();
        genCharged.clear();

        measgen_NeutralPulls.clear();
        measgen_ChargedPulls.clear();
        fitgen_NeutralPulls.clear();
        fitgen_ChargedPulls.clear();

        genNeutralPdg.clear();
        genChargedPdg.clear();

        mcCharged.clear();
        mcNeutral.clear();

	delete fitter;*/
    }//end pvector size check
   else{
	clear();
	rejectNo=1;
	rejects->Fill(rejectNo);
	}
	//track events for each call
//	std::cout<<" test 9a"<<std::endl;
	evtNo++;
// std::cout<<" test 9 "<<std::endl;	
	return;
}


