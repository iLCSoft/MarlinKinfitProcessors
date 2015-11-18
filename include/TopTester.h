#ifndef TopTester_h
#define TopTester_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "TopEventILC.h"

using namespace lcio ;
using namespace marlin ;


/**  Test processor using TopEventILC Toy MC
 *   
 *   ... testing a ttbar -> 6jets or 4jets+lnu hypotheses
 *   with energy and momentum conservation
 *   and an equal mass constraint for tops, MW=80.4 GeV for both W's
 * 
 *  <h4>Input - Prerequisites</h4>
 *  nothing
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the ReconstructedParticle collection
 * 
 * @author J. List, DESY
 * @version $Id: TopTester.h,v 1.0 2014/11/26 10:40:46 boehmej Exp $ 
 */

class TopTester : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TopTester ; }
  
  
  TopTester() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
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

  /** Input parameter: center of mass energy, semileptonic?, which fitter.
   */
  float _ecm;
  bool _semileptonic, _traceall, _leptonasjet;
  int _ifitter, _ievttrace, _ntoy; 
 
  float prob, topmass, startmassW1, startmassW2;
  float momentum[3], energy;
           
  int _nRun, _nEvt, nit;
   
  TopEventILC* topevent;
 

  //output
  // TTree *outTree;
         
} ;

#endif



