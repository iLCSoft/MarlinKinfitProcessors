#ifndef ZHllqq5CFit_h
#define ZHllqq5CFit_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "TopEventILC.h"
#include <vector>

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
 * @version $Id: ZHllqq5CFit.h,v 1.2 2010/07/19 23:41:46 mbeckman Exp $
 */

class ZHllqq5CFit : public Processor {

  public:

    virtual Processor*  newProcessor() { return new ZHllqq5CFit ; }


    ZHllqq5CFit() ;

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
    
    /** Input collections names.
     */
    std::string _jetcolName{}, _colLeptons{};
    /** Input parameters.
     */
    float _ecm{}, _isrpzmax{};
    int _fitISR{}, _ifitter{}, _ievttrace{};
    bool _traceall{};
    double _errene{}, _errtheta{}, _errphi{};
    float _sigmaScaleFactor;
    /** Output collections names.
     */
    std::string _PostFitRecoCol{};
    std::string _OutPullECol{};  
    std::string _OutPullThetaCol{};
    std::string _OutPullPhiCol{};
    
    // helpers forISR and track momentum calculations - needed here?
    double b{}, ISRPzMaxB{};
    float m_Bfield;
    double c;
    double mm2m;
    double eV2GeV;
    double eB;

    double prob{}, startmassZ{}, startmassH{};
    int _nRun{}, _nEvt{}, nit{};

    int trackcounter{};
    int correction{};
    int nCo{};
    //output
    // TTree *outTree;

} ;

#endif
