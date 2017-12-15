#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackImpl.h"
#include "EVENT/Track.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"
#include "TFile.h"
#include "TTree.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TLorentzVector.h"
typedef CLHEP::HepLorentzVector LorentzVector ;
#include "LeptonFitObject.h"
#include "TrackParticleFitObject.h"
#include "JetFitObject.h"
#include "VertexFitObject.h"
#include "OPALFitterGSL.h"
//#include "NewFitterGSL.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"
#include "TH1D.h"
//#pragma link C++ class std::vector<std::vector<double> >+;
//#pragma link C++ class vector<vector<double> >+;
//#pragma link C++ class vector<TLorentzVector>+;
//#pragma link C++ class std::vector<TLorentzVector>+;
//#pragma link C++ class std::vector<Track*>+;
using namespace lcio ;

/** DiTrackGammaCandidateFinder:<br>
 *
 * (modelled after GammaGammaCandidateFinder processor)
 * 
 * @author Justin Anguiano, University of Kansas
 * Code follows convention X -> p1+ p2- gamma
 */
struct massconstraint{
        std::vector<int> chargedIndices{};
        std::vector<int> neutralIndices{};
        double mass{};
 };

class MassConstraintFitter : public marlin::Processor {
  
 public:
  
 virtual marlin::Processor*  newProcessor() { return new MassConstraintFitter ; }

  MassConstraintFitter(const MassConstraintFitter&) = delete ;
  MassConstraintFitter& operator=(const MassConstraintFitter&) = delete ;
  
  MassConstraintFitter() ;

  /** Called at the beginning of the job before anything is read.
   *  Use to initialize the proscessor, e.g. book histograms.
   */
  virtual void init() ;
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 


  /** Called after data processing for clean up.
   */
  virtual void end() ;

private:
  bool FindPFOs( LCEvent* evt );
  bool FindTracks( LCEvent* evt );
  bool FindMCParticles( LCEvent* evt);
  

  int getCorrespondingMCParticleIndex(TLorentzVector rec);
  void FindMassConstraintCandidates( LCCollectionVec* recparcol);
  int getIndexOfMatchingIndices(std::vector<int> indices, int index);
  bool vectorIndicesOverlap(std::vector<int> v1, std::vector<int> v2);
  bool secondaryConstraintCombinationValid(std::vector<int> indices, std::vector<massconstraint*> constraintVector );
  OPALFitterGSL*  setUpFit(std::vector<int> neutralIndices, std::vector<int> chargedIndices, std::vector<int> massConstraintIndices, std::vector<massconstraint*> constraintVector, std::vector<TLorentzVector> pneutral, std::vector<TLorentzVector> ptrack, std::vector<ReconstructedParticle*> pNeutralVec, std::vector<Track*> pTrackVec);
  std::vector<double> getChargedParticleErrors(TLorentzVector pcharged, Track* ptrk);
  std::vector<double> getNeutralErrors(TLorentzVector pneutral, ReconstructedParticle* pNeutral);
  std::vector<double> getTrackErrors(Track* ptrk);
  void PrintCov(FloatVec cov, int dim);
  void PrintCov(double* cov, int dim);
  void setFitErrors(double* cov, int dim);
  double* ConstructParentMeasCovMatrix();
//  std::vector<double> ConstructChargedSubMatrix(std::vector<double> trackparams, TLorentzVector ptlv);
  std::vector<double> ConstructChargedSubMatrix(TLorentzVector ptlv);
  std::vector<double> ConstructNeutralSubMatrix(TLorentzVector p);
  double* ConcatSubMatrices(std::vector<std::vector<double> > matrices);
  void setParentErrors(FloatVec meascov, FloatVec fitcov);
  void generateSubsets(std::vector<int> v, int k, unsigned int start, int currLen, std::vector<bool> used, std::vector<std::vector<int> >& combinations);
//  FloatVec ConstructCovMatrix(std::vector<std::vector<double> > trackparams, std::vector<TLorentzVector> charged, std::vector<TLorentzVector> neutral, double* cov);
  FloatVec ConstructCovMatrix(std::vector<TLorentzVector> charged, std::vector<TLorentzVector> neutral, double* cov);
 void clear(); 
	void generateIndicesCombinations(int vectorsize, int nparticles,std::vector<std::vector<int> >& combinations);

	ReconstructedParticleImpl* constructFitParticle(TLorentzVector fitp, ReconstructedParticle* measrp);

  TrackImpl* constructTrack(TLorentzVector fitp, Track* meast);
  std::vector<double> buildTrackVector(Track* t);
  std::vector<massconstraint*> buildMassConstraint(std::vector<std::vector<int> > neutralIndices, std::vector<std::vector<int> > chargedIndices, double mass, int nNeutral, int nCharged);
  std::vector<double> build4vec(TLorentzVector v);
  std::vector<double> buildFitTrackVector(TrackParticleFitObject* tfo);
  std::vector<double> buildLeptonFitVector(LeptonFitObject* lfo);
  void printCombinations(std::vector<std::vector<int> > combs);
  void printmassconstraints(std::vector<massconstraint*> cv);
  std::vector<double> buildLeptonVector(TLorentzVector v, double q);
  std::vector<double> buildNeutralParamVector(TLorentzVector v);
 // FloatVec ConstructCovMatrix(std::vector<std::vector<double> > trackparams, std::vector<TLorentzVector> charged, std::vector<TLorentzVector> neutral, double* cov);
 
// struct massconstraint{
  //      std::vector<int> chargedIndices;
    //    std::vector<int> neutralIndices;
      //  double mass;
// };

 protected:
//TTree stuff for fit analysis
  TFile *rootFile{};
  //normal analysis tree for comparing fits to measured
  TTree *tree{};

  //tree that includes generator level information
  TTree *genTree{}; 
  
 ///////////////////////////////////////////////////////
 // Vars for regular tree
 /////////////////////////////////////////////////////
  //TTree vars for fit analysis
  //money plots
  double RecoEnergy{};
  double RecoMass{};
  double FitEnergy{};
  double FitProbability{};
  double Chisq{};
  double EnergyAsymmetryGen{};     // Only makes sense for X -> gamma gamma 
  int evtNo{};
  int rejectNo{};
  TH1D* rejects{};

//add variables
  std::vector<TLorentzVector> measNeutral{};
  std::vector<TLorentzVector> measCharged{};
  std::vector<std::vector<double> > mcParticleVec{};
  std::vector<std::vector<double> > measNeutralVec{};
  std::vector<std::vector<double> > measNeutralParamVec{};
  std::vector<std::vector<double> > measChargedVec{};
  std::vector<Track*> measTrack{};
  std::vector<std::vector<double> > measTrackVec{};
  std::vector<TLorentzVector> fitNeutral{};
  std::vector<TLorentzVector> fitCharged{};
  std::vector<std::vector<double> > fitNeutralVec{};//px py pz e
  std::vector<std::vector<double> > fitNeutralParamVec{};//e theta phi
  std::vector<std::vector<double> > fitChargedVec{};//px py pz e
  std::vector<Track*> fitTrack{};
  std::vector<std::vector<double> > fitTrackVec{};//k theta phi
  std::vector<std::vector<double> > measNeutral_err{};
  std::vector<std::vector<double> > measCharged_err{};
  std::vector<std::vector<double> > measTrack_err{};
  std::vector<std::vector<double> > fitTrack_err{};
  std::vector<std::vector<double> > fitNeutral_err{};
  std::vector<std::vector<double> > fitCharged_err{};
  std::vector<double> measNeutralPdg{};
  std::vector<double> measChargedPdg{};
	//vector<double>: (E || k, Theta, Phi); pull order
  std::vector<std::vector<double> > fitmeas_NeutralPulls{};
  std::vector<std::vector<double> > fitmeas_ChargedPulls{};
  TLorentzVector measParent{};
  TLorentzVector fitParent{};
  std::vector<double> measParent_err{};
  std::vector<double> fitParent_err{};
		//(Px, Py, Pz, E) pull order
  std::vector<double> fitmeas_ParentPulls{};
  //monte carlo variables
  std::vector<TLorentzVector> genNeutral{};
  std::vector<TLorentzVector> genCharged{};
	//vector<double>? (E , Theta, Phi); pull order even for charged particles
  std::vector<std::vector<double> > measgen_NeutralPulls{};  
  std::vector<std::vector<double> > measgen_ChargedPulls{};
  std::vector<std::vector<double> > fitgen_NeutralPulls{};
  std::vector<std::vector<double> > fitgen_ChargedPulls{};
  TLorentzVector genParent{};
//probably should update tlvs to genparentvec for plotting 4vec array
  std::vector<double> genNeutralPdg{};
  std::vector<double> genChargedPdg{};
  std::vector<std::vector<double> > mcTrack{};//not implemented yet
  std::vector<TLorentzVector> mcCharged{};
  std::vector<TLorentzVector> mcNeutral{};
  std::vector<std::vector<double> > mcChargedVec{}; //4vectors gen
  std::vector<std::vector<double> > mcChargedParamVec{}; //local param gen
  std::vector<std::vector<double> > mcNeutralVec{}; //4vec
  std::vector<std::vector<double> > mcNeutralParamVec{};// local param gen
  //vectors for particle input information
  int _parentPdg{};
  double _parentMass{};
  double _parentCharge{};
  int _nDaughters{};
  int _nCharged{};
  int _nChargedParams{};
  int _nNeutral{};
  int _nNeutralParams{};
  std::vector<int> _daughterChargedPdgs{};
  std::vector<int> _daughterNeutralPdgs{};
  std::vector<float> _daughterChargedMass{};
  std::vector<float> _daughterNeutralMass{};
  //secondary constraint stuff
  int _nMassConstraints{};
  std::vector<float> _secondaryMasses{};
  std::vector<float> _secondaryNCharged{};
  std::vector<float> _secondaryNNeutral{};

  //global fitobject pointers
  //return Fit objects in the fitter is not retaining the derived class and chopping to base class
// std::vector<std::vector<int> > neutralCandidateIndices;
// std::vector<std::vector<int> > chargedCandidateIndices;
// JetFitObject* gammaJet;
  std::vector<JetFitObject*> neutralJets{};

  //LeptonFitObject* part1;
 // LeptonFitObject* part2;
   std::vector<LeptonFitObject*> TrackFO{}; 
 //   std::vector<TrackParticleFitObject*> TrackFO;
 //  OPALFitterGSL* fitter{};
//    BaseFitter* ftest;

  std::vector<ReconstructedParticle*>_pfovec{};
  std::vector<Track*> _trackvec{};
  std::vector<MCParticle*> _mcpartvec{};
  int   _printing{};
  std::string _inputTrackCollectionName{};
  std::string _inputParticleCollectionName{};
  std::string _mcParticleCollectionName{};
  std::string _outputParticleCollectionName{};
  std::string _outputTrackCollectionName{};
  double _fitProbabilityCut{};
  double _allowedMassDeviation{};
 
  int _ifitter{};
  int _fitAnalysis{};
  int _genAnalysis{};

  double _photonAngularError{};
  int _photonAngularErrorModel{};
 
  std::string m_rootFile{};

} ;
