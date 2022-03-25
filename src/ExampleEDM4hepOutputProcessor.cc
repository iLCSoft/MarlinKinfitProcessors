#include "ExampleEDM4hepOutputProcessor.h"
#include <iostream>
#include <string>
#include <array>
#include <math.h>
#include <EVENT/LCCollection.h>
//  #include <EVENT/MCParticle.h>
#include <podio/UserDataCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/LCFloatVec.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/LCIterator.h>
#include <lcio.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include <algorithm>
#include <iterator>

// using namespace lcio;
using namespace marlin;

ExampleEDM4hepOutputProcessor aExampleEDM4hepOutputProcessor;

ExampleEDM4hepOutputProcessor::ExampleEDM4hepOutputProcessor()
    : Processor("ExampleEDM4hepOutputProcessor")
 {

  // modify processor description
  _description = "ExampleEDM4hepOutputProcessor shows how to write EDM4hep "
                 "output files from within a Marlin processor";

  registerProcessorParameter(
      "EDM4hepOutputFile", "Output file in EDM4hep format that will be filled",
      m_outputFile, std::string("edm4hep_outputs.root"));

  registerProcessorParameter(
      "ReconstructedParticleCollection",
      "The Reconstructed Particle collection to store in EDM4hep format",
      m_RecParCollectionName, std::string("FitReco"));



 /* 
  registerProcessorParameter(
      "EnergyPullCollection",
      "The Energy Pull collection to store in EDM4hep format",
      m_RecParCollectionName, std::string("PullEnergy"));
*/
}

void ExampleEDM4hepOutputProcessor::init() {

  streamlog_out(DEBUG) << "   init called  " << std::endl;

  // usually a good idea to
  printParameters();

  // Create the writer and the EDM4hep collections, register them for output
  m_writer = std::make_unique<podio::ROOTWriter>(m_outputFile, &m_store);


// auto& EnergyPulls =m_store.create<podio::UserDataCollection<float>>("EnergyPulls");
// auto& ThetaPulls = m_store.create<podio::UserDataCollection<float>>("ThetaPulls");
// auto& PhiPulls = m_store.create<podio::UserDataCollection<float>>("PhiPulls");
 
  m_edm4hepRecParColl = new edm4hep::ReconstructedParticleCollection();
  m_store.registerCollection(m_RecParCollectionName, m_edm4hepRecParColl);
 EnergyPullsmu = new podio::UserDataCollection<float>();
  m_store.registerCollection("EnergyPullsmu", EnergyPullsmu);
 ThetaPullsmu = new podio::UserDataCollection<float>();
  m_store.registerCollection("ThetaPullsmu", ThetaPullsmu);
 PhiPullsmu = new podio::UserDataCollection<float>();
  m_store.registerCollection("PhiPullsmu", PhiPullsmu);
 EnergyPullsjet = new podio::UserDataCollection<float>();
  m_store.registerCollection("EnergyPullsjet", EnergyPullsjet);
 ThetaPullsjet = new podio::UserDataCollection<float>();
  m_store.registerCollection("ThetaPullsjet", ThetaPullsjet);
 PhiPullsjet = new podio::UserDataCollection<float>();
  m_store.registerCollection("PhiPullsjet", PhiPullsjet);


  m_writer->registerForWrite(m_RecParCollectionName);
  m_writer->registerForWrite("EnergyPullsmu");
  m_writer->registerForWrite("ThetaPullsmu");
  m_writer->registerForWrite("PhiPullsmu");
  m_writer->registerForWrite("EnergyPullsjet");
  m_writer->registerForWrite("ThetaPullsjet");
  m_writer->registerForWrite("PhiPullsjet");
}

void ExampleEDM4hepOutputProcessor::processRunHeader(LCRunHeader * /*run*/) {}

// Helper function for constructing an edm4hep::Vector3f from an array of
// doubles. Necessary because some of the LCIO return types are not the same
// storage type as in edm4hep, the momentum in the MCParticle:
// https://github.com/key4hep/EDM4hep/blob/master/edm4hep.yaml#L128
// vR
// https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/classEVENT_1_1MCParticle.html#aca3a6c449405c2d3492c61dfbefc2b2f
auto toVector3f(const double *vals) {
  // do explicit casts to silence narrowing warnings
  return edm4hep::Vector3f{(float)vals[0], (float)vals[1], (float)vals[2]};
}

void ExampleEDM4hepOutputProcessor::processEvent(LCEvent *evt) {
  // Use an LCIterator to iterate over all the MCs that we have from the LCIO
  // inputs
 const std::vector<std::string>* colnames = evt->getCollectionNames();
 if (std::find(colnames->begin(), colnames->end(), "FitReco") != colnames->end()) {
 
 auto *lcioColl = evt->getCollection(m_RecParCollectionName);
  streamlog_out(MESSAGE) << " number of Reconstructed: "
                         << lcioColl->getNumberOfElements() << '\n';

  // Use an LCIterator for convenience
  auto lcioCollIt = lcio::LCIterator<lcio::ReconstructedParticle>(lcioColl);

  // First fill the data
  while (const auto *lcioRecPar = lcioCollIt.next()) {
  // streamlog_out(MESSAGE) << " test: "<< '\n';

    // Create an edm4hep MCParticle via the create factory method of the
    // MCParticleCollection
    auto edm4hepRecPar = m_edm4hepRecParColl->create();

    // Populate the edm4hep::MCParticle from the LCIO one
//    edm4hepMC.setMass(lcioRecPar->getMass());
//    edm4hepMC.setPDG(lcioRecPar->getPDG());
//    edm4hepMC.setCharge(lcioRecPar->getCharge());
//    edm4hepMC.setTime(lcioRecPar->getTime());
//    edm4hepMC.setGeneratorStatus(lcioRecPar->getGeneratorStatus());
//    edm4hepMC.setSimulatorStatus(lcioRecPar->getSimulatorStatus());

//    edm4hepMC.setVertex(edm4hep::Vector3d(lcioRecPar->getVertex()));
//    edm4hepMC.setEndpoint(edm4hep::Vector3d(lcioRecPar->getEndpoint()));
   

    edm4hep::Vector3f momentum = toVector3f(lcioRecPar->getMomentum());
    edm4hepRecPar.setMomentum(momentum);
    edm4hepRecPar.setEnergy(lcioRecPar->getEnergy());
    edm4hepRecPar.setCharge(lcioRecPar->getCharge());
    edm4hepRecPar.setType(lcioRecPar->getType());
    streamlog_out(MESSAGE)  << "Partcle Type: " << lcioRecPar->getType() <<std::endl;
    
    float energy = lcioRecPar->getEnergy();
   
    float inv_mass = sqrt(pow(energy,2)-(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]));
    streamlog_out(MESSAGE)  << "Partcle Mass: " << inv_mass <<std::endl;
 
    edm4hepRecPar.setMass(inv_mass);
    if (lcioRecPar->getType() == 13){
   
    std::array<float, 10>  CovMatrix;
    FloatVec inputCovMatrix = lcioRecPar->getCovMatrix();
     /*
    for (int i=0; i<10 ; i++){
    CovMatrix[i] = inputCovMatrix[i];
    } */ 
    std::copy(std::begin(inputCovMatrix), std::end(inputCovMatrix), std::begin(CovMatrix));
    edm4hepRecPar.setCovMatrix(CovMatrix);
    
 } 
   
//    edm4hepMC.setMomentumAtEndpoint(
//        toVector3f(lcioRecPar->getMomentumAtEndpoint()));
//    edm4hepMC.setSpin(edm4hep::Vector3f(lcioRecPar->getSpin()));

//    edm4hepMC.setColorFlow(edm4hep::Vector2i(lcioRecPar->getColorFlow()));
  }
/*
  // Now restore parent/daugther relations
  size_t lcioIndex = 0;
  // Making use of the fact that the collection vec actually inherits from an
  // std::vector and we can use a range-based for-loop with that as well as
  // using STL algorithms on that
  const auto *objVec = dynamic_cast<const lcio::LCCollectionVec *>(lcioColl);
  for (const auto *obj : *objVec) {
    const auto *lcioMC = static_cast<const lcio::MCParticle *>(obj);
    auto edm4hepMC = (*m_edm4hepMCColl)[lcioIndex++];

    for (const auto *parent : lcioMC->getParents()) {
      const auto index = std::distance(
          objVec->begin(), std::find(objVec->begin(), objVec->end(), parent));

      auto edm4hepParent = (*m_edm4hepMCColl)[index];
      edm4hepMC.addToParents(edm4hepParent);
      edm4hepParent.addToDaughters(edm4hepMC);
    }
  }
*/
  // Write the event contents and clear the collections for the next event

}
   if (std::find(colnames->begin(), colnames->end(), "PullEnergy") != colnames->end()){ 

  auto *lcioCollE = evt->getCollection("PullEnergy");
  auto *lcioCollTheta = evt->getCollection("PullTheta");
  auto *lcioCollPhi = evt->getCollection("PullPhi");
  streamlog_out(MESSAGE) << "PullEnergy Collection found : "<< '\n';
  
  // Use an LCIterator for convenience
  auto lcioCollItE = lcio::LCIterator<lcio::LCFloatVec>(lcioCollE);
  auto lcioCollItTheta = lcio::LCIterator<lcio::LCFloatVec>(lcioCollTheta);
  auto lcioCollItPhi = lcio::LCIterator<lcio::LCFloatVec>(lcioCollPhi);
    // First fill the data
  const auto lcioPullE = lcioCollItE.next();
  const auto lcioPullTheta = lcioCollItTheta.next();
  const auto lcioPullPhi = lcioCollItPhi.next();
  for (int j=0 ;j <int(lcioPullE[0].size()); j++) {
     if (j <int(lcioPullE[0].size())/2){
      EnergyPullsjet->push_back(lcioPullE[0][j]);
      ThetaPullsjet->push_back(lcioPullTheta[0][j]);
      PhiPullsjet->push_back(lcioPullPhi[0][j]);
  streamlog_out(MESSAGE) << "PullEnergy found : "<<lcioPullE[0][j] << '\n';
   }
     else{ 
      EnergyPullsmu->push_back(lcioPullE[0][j]);
      ThetaPullsmu->push_back(lcioPullTheta[0][j]);
      PhiPullsmu->push_back(lcioPullPhi[0][j]);
    }

  }



    // streamlog_out(MESSAGE) << " test: "<< '\n'; 
    // Create an edm4hep MCParticle via the create factory method of the
    //auto edm4hepPullE = m_edm4hepRecParColl->create();
 
  }
  m_writer->writeEvent();
  m_store.clearCollections();

}


void ExampleEDM4hepOutputProcessor::check(LCEvent * /*evt*/) {
  // nothing to check here - could be used to fill checkplots in
  // reconstruction processor
}

void ExampleEDM4hepOutputProcessor::end() {
  streamlog_out(MESSAGE) << "Finalizing writer\n";
  m_writer->finish();
}
