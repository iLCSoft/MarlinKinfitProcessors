#ifndef ZH5CFit_h
#define ZH5CFit_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "TopEventILC.h"

using namespace lcio ;
using namespace marlin ;


/**  An example processor for a kinematic fit
 *
 *   ... testing a W+W- -> 4jets hypothesis
 *   with energy and momentum conservation
 *   and an equal mass constraint
 *
 *  <h4>Input - Prerequisites</h4>
 *  Needs 4 reconstructed jets
 *
 *  <h4>Output</h4>
 *  A histogram.
 *
 * @param CollectionName Name of the ReconstructedParticle collection
 *
 * @author J. List, DESY
 * @version $Id: ZH5CFit.h,v 1.2 2010/07/19 23:41:46 mbeckman Exp $
 */

class ZH5CFit : public Processor {

 public:

  virtual Processor*  newProcessor() { return new ZH5CFit ; }


  ZH5CFit() ;

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

  double JetEnergyResolution(double E);

 protected:

  /** Input collection name.
   */
  std::string _jetcolName{}, _name{}, _OutputCollection{} ;
  /** Input parameter: center of mass energy.
   */
  float _ecm{}, _isrpzmax{};
  int _fitISR{}, _ifitter{}, _ievttrace{};
  bool _traceall{};
  double _errene{}, _errtheta{}, _errphi{};
  std::string _OutCol{};
  double b{}, ISRPzMaxB{};

  float prob{}, bestprob{}, bestnit{}, bestmassZ{}, bestmassH{}, beststartmassZ{}, beststartmassH{}, bestphotonenergy{}, startmassZ{}, startmassH{}, variable{};
  float chi2best;
  double Zmomentum[3]{}, Hmomentum[3]{}, ISRmomentum[3]{};
  float momentum[3]{}, energy{};

  int _nRun{}, _nEvt{}, nit{};

  int bestperm{}, errorflag{};


  //output
  // TTree *outTree;

} ;

#endif
