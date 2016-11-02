#ifndef DijetTester_h
#define DijetTester_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "DijetEventILC.h"

using namespace lcio ;
using namespace marlin ;


/**  Test processor using DijetEventILC Toy MC
 *   
 *   ... testing a e+e- -> 2jets or 2leptons
 *   with energy and momentum conservation
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
 * @version $Id: DijetTester.h,v 1.0 2015/11/18 10:40:46 boehmej Exp $ 
 */

class DijetTester : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DijetTester ; }
  
  
  DijetTester() ;
  
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

  /** Input parameter: center of mass energy, leptonic?, which fitter.
   */
  float _ecm;
  bool _leptonic, _traceall, _leptonasjet;
  int _ifitter, _ievttrace, _ntoy; 
 
  float prob;
  float momentum[3], energy;
           
  int _nRun, _nEvt, nit;
   
  DijetEventILC* dijetevent;
 

  //output
  // TTree *outTree;
         
} ;

#endif



