#include <marlin/Processor.h>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCIterator.h"
#include "UTIL/Operators.h"
#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>

// STUFF needed for GEAR
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>
#include "gear/BField.h"

#include <iomanip>

using namespace lcio ;
using namespace marlin ;
using namespace std ;

/** MCParticleFilter processor
 *  Chooses which MCParticles to keep
 * author:  Graham Wilson
*/

class MCParticleFilter : public Processor {
  
 public:
 
  virtual Processor*  newProcessor() { return new MCParticleFilter ; }
   
  MCParticleFilter() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

 protected:

  int nEvt ;

 private:

  int   _printing;
  std::string _mcParticleCollectionName ;
  std::string _outputParticleCollectionName;
  
} ;
