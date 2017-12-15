#include "marlin/Processor.h"
#include "EVENT/Track.h"
#include "lcio.h"
#include "TFile.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"
#include "IMPL/ParticleIDImpl.h"
#include "IMPL/TrackImpl.h"

using namespace lcio;

	/** TrackResponseAdjuster:<br>
 *
 * 
 * @author Justin Anguiano, University of Kansas
 * 
 */

 class TrackResponseAdjuster : public marlin::Processor {

 public:

 virtual marlin::Processor*  newProcessor() { return new TrackResponseAdjuster ; }

  TrackResponseAdjuster(const TrackResponseAdjuster&) = delete ;
  TrackResponseAdjuster& operator=(const TrackResponseAdjuster&) = delete ;

  TrackResponseAdjuster() ;

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

  bool FindTracks(LCEvent* evt);
 
  double safeAcos(double x);
 
  private:
  int nEvt{};
  
  std::vector<Track*> _trackvec{};
  int   _printing{};

  double _D0ErrorScaleFactor{};
  double _Z0ErrorScaleFactor{};
  double _OmegaErrorScaleFactor{};
  double _PhiErrorScaleFactor{};
  double _TanLambdaErrorScaleFactor{};
  
// _inputTrackCollectionName 
  std::string _outputTrackCollectionName{};
  std::string _inputTrackCollectionName{};
//  std::string m_rootFile{};
};
