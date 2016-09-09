#ifndef SimpleVertexProcessor_h
#define SimpleVertexProcessor_h 1
#include "marlin/Processor.h"

using namespace marlin ;

#include "OPALFitterGSL.h"

class SimpleVertexProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SimpleVertexProcessor ; }
  
  SimpleVertexProcessor() ;
  
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

  OPALFitterGSL* _opalFit;

};


#endif



