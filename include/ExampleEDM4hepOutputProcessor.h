#ifndef EDM4HEPOUTPUT_EXAMPLEEDM4HEPOUTPUTPROCESSOR_h
#define EDM4HEPOUTPUT_EXAMPLEEDM4HEPOUTPUTPROCESSOR_h 1

//#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "podio/EventStore.h"
#include "podio/ROOTWriter.h"
#include "podio/UserDataCollection.h"
#include "lcio.h"
#include "marlin/Processor.h"

#include <memory>
#include <string>

using namespace lcio;
using namespace marlin;

/**  Example processor for marlin.
 *
 *  If compiled with MARLIN_USE_AIDA
 *  it creates a histogram (cloud) of the MCParticle energies.
 *
 *  <h4>Input - Prerequisites</h4>
 *  Needs a model in onnx that can be loaded via the onnxruntime library
 *
 *  <h4>Output</h4>
 *
 * @param CollectionName Name of the MCParticle collection
 *
 * @author F. Gaede, DESY
 * @version $Id: ExampleEDM4hepOutputProcessor.h,v 1.4 2005-10-11 12:57:39 gaede
 * Exp $
 */

class ExampleEDM4hepOutputProcessor : public Processor {

public:
  virtual Processor *newProcessor() {
    return new ExampleEDM4hepOutputProcessor;
  }

  ExampleEDM4hepOutputProcessor();
  ExampleEDM4hepOutputProcessor(const ExampleEDM4hepOutputProcessor &) = delete;
  ExampleEDM4hepOutputProcessor &
  operator=(const ExampleEDM4hepOutputProcessor &) = delete;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  void init() override;

  /** Called for every run.
   */
  void processRunHeader(LCRunHeader *run) override;

  /** Called for every event - the working horse.
   */
  void processEvent(LCEvent *evt) override;

  void check(LCEvent *evt) override;

  /** Called after data processing for clean up.
   */
  void end() override;

private:
  std::string m_outputFile{};
  std::string m_RecParCollectionName{};
  std::unique_ptr<podio::ROOTWriter> m_writer{nullptr};
  podio::EventStore m_store{};
  podio::UserDataCollection<float> *EnergyPullsmu{nullptr};
  podio::UserDataCollection<float> *ThetaPullsmu{nullptr};
  podio::UserDataCollection<float> *PhiPullsmu{nullptr};
  podio::UserDataCollection<float> *EnergyPullsjet{nullptr};
  podio::UserDataCollection<float> *ThetaPullsjet{nullptr};
  podio::UserDataCollection<float> *PhiPullsjet{nullptr};   
    
  // EDM4hep collections
  //
  // NOTE: This needs a better interface as soon as there are several
  // collections involved.
  // NOTE: Ownership will be taken over by the EventStore
  edm4hep::ReconstructedParticleCollection *m_edm4hepRecParColl{nullptr};
 // edm4hep::Vector3fCollection *m_edm4hepvecColl{nullptr};

};

#endif
