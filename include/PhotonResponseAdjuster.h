#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"
//#include "TFile.h"
#include "TLorentzVector.h"
#include "IMPL/ParticleIDImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "TRandom3.h"

using namespace lcio ;
	/** PhotonResponseAdjuster:<br>
 *
 * 
 * @author Justin Anguiano, University of Kansas
 * 
 */

 class PhotonResponseAdjuster : public marlin::Processor {

 public:

 virtual marlin::Processor*  newProcessor() { return new PhotonResponseAdjuster ; }

  PhotonResponseAdjuster(const PhotonResponseAdjuster&) = delete ;
  PhotonResponseAdjuster& operator=(const PhotonResponseAdjuster&) = delete ;

  PhotonResponseAdjuster() ;

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

  bool FindMCParticles( LCEvent* evt ); 
  bool FindPFOs( LCEvent* evt );
  int getCorrespondingMCParticleIndex(TLorentzVector rec);
  double* resimulateDirection(TLorentzVector mcgamma);
  double safeAcos(double x);

  private:
  int _nRun{};
  int _nEvt{};
  std::vector<ReconstructedParticle*> _pfovec{};
  std::vector<MCParticle*> _mcpartvec{};
  std::vector<bool> _mcpartflags{};
  int _printing{};

  std::string _inputParticleCollectionName{};
  std::string _mcParticleCollectionName{};
  std::string _outputParticleCollectionName{};
  double _energyScaleFactor{};
  int nrejected{};
  //matching cuts
  TRandom3* rng{};
  int _smearAngles{};
  int _angularSmearingModel{};
  double _dTheta{};
  double _dPhi{};
  double _allowedEnergyDeviation{};
  double _allowedAngularDeviation{};
  std::string m_rootFile{};
};
