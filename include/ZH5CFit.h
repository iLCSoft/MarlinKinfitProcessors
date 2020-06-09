#ifndef ZH5CFit_h
#define ZH5CFit_h 1

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

  void compcorrect();
  void SetZero();
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
  float errorcode;
  float Zmomentum[3]{}, Hmomentum[3]{}, ISRmomentum[3]{};
  float Z_Energy{}, H_Energy{};
  float momentum[3]{}, energy{};
  int _nRun{}, _nEvt{}, nit{}, nCo{};

  int bestperm{}, errorflag{};
  TTree *ZHTree{};
  TFile* _fout = nullptr ;

  float Hmass_before_fit{}, Hmass_after_fit{}, Hmass_NoFit{};
  float Error_code{};
  float hpull_jet1_E{}, hpull_jet2_E{}, hpull_jet3_E{}, hpull_jet4_E{};
  float hpull_jet1_th{}, hpull_jet2_th{}, hpull_jet3_th{}, hpull_jet4_th{};
  float hpull_jet1_phi{}, hpull_jet2_phi{}, hpull_jet3_phi{}, hpull_jet4_phi{};
  int jetmatch{}, jetmatchth{}, jetmatchphi{};
  std::string _outfile{};

  int besterr{};
  double bestzvalue{} ;
  double chi2startmassZ{}, chi2startmassH {};

  double Px{}, Px2{}, Py{}, Py2{}, Pz{}, Pz2{}, pT2{}, P{}, P2{};
  double SigPx2{}, SigPxSigPy{}, SigPxSigPz{}, SigPy2{}, SigPySigPz{}, SigPz2{}, SigE2{};
  double dth_dpx{}, dth_dpy{}, dth_dpz{}, dphi_dpx{}, dphi_dpy{}, JetResE{}, JetResTheta{}, JetResPhi{};


  std::string _colMCP{} ;
  std::string _errorflowcollection {};
  std::string _SLDCol {};
  int nSLDB{};
  int nSLDC{};
  int nSLDBC{};
  typedef std::vector<int>		IntVector;
  IntVector B_index{};
  IntVector C_index{};
  double E_lab;
  double Elab;
  int _NuE{};
  int _useErrorFlow{};
  float sigmaScaleFactor{};
  std::string  _NuCorrector{};
  float ENuplus{};
  float ENuminus{};
  float l_px{};
  float l_py{};
  float l_pz{};
  float l_p{};
  float l_theta{};
  float l_phi{};
  float delta_theta[4]{};
  float delta_phi[4]{};
  int bestjet_th{};
  int bestjet_phi{};
  int bestjet{};
  int correction{};
  //output
  // TTree *outTree;

} ;

#endif
