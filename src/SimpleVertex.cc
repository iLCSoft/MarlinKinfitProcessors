#include "SimpleVertex.h"
#include <iostream>
#include <map>
#include <cassert>

#include <marlin/Global.h>
#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/TrackState.h>

#include "TrackParticleFitObject.h"
#include "TwoVector.h"
#include "BaseHardConstraint.h"
#include "VertexFitObject.h"
#include "TextTracer.h"

using namespace lcio ;
using namespace marlin ;
using std::cout;
using std::endl;

SimpleVertexProcessor aSimpleVertexProcessor ;

SimpleVertexProcessor::SimpleVertexProcessor() : Processor("SimpleVertexProcessor") {
  // modify processor description
  _description = "SimpleVertexProcessor fits 2 charged tracks to a single vertex" ;
  _opalFit = new OPALFitterGSL();
}


void SimpleVertexProcessor::init() {
  cout << "hello from SimpleVertexProcessor::init" << endl;
  return;
}

void SimpleVertexProcessor::processRunHeader( LCRunHeader* run) {
  cout << "hello from SimpleVertexProcessor::processRunHeader" << endl;
}

void SimpleVertexProcessor::processEvent( LCEvent * evt ) {

  TextTracer tracer (std::cout);
  
  _opalFit->reset();
  _opalFit->setTracer (tracer);

  // designed to work with single K0short events
  // fit tracks from kshort decays to a vertex

  const double mass_pi = 0.14; // assumed mass of tracks
  
  try {
    LCCollection* col_trk=evt->getCollection("MarlinTrkTracks");
    if ( col_trk->getNumberOfElements()==2 ) { // only use events with exactly 2 reconstructed tracks

      const Track* trk[2];
      const TrackState* ts[2];
      TrackParticleFitObject* tfo[2];

      for (int i=0; i<2; i++) {
        trk[i] = dynamic_cast < Track* > ( col_trk->getElementAt(i) );
        assert( trk[i] );

        ts[i] = trk[i]->getTrackState( TrackState::AtFirstHit ); // use track state at first hit  - this should be closest to the vertex

        if ( ts[i] ) {
          tfo[i] = new TrackParticleFitObject(ts[i], mass_pi);
        } else {
          cout << "WARNING: no track state found @ first hit, using default track parameters" << endl;
          tfo[i] = new TrackParticleFitObject(trk[i], mass_pi);
        }
      }

      // name the fit objects
      tfo[0]->setName("track0");
      tfo[1]->setName("track1");

      // define vertex object (initially at 0,0,0)
      VertexFitObject* vertexFO = new VertexFitObject("commonVtx",0,0,0);
      vertexFO->setParam(0,0.,false,false); // measured=false, fixed=false
      vertexFO->setParam(1,0.,false,false);
      vertexFO->setParam(2,0.,false,false);
      vertexFO->setError(0,100.); 
      vertexFO->setError(1,100.);
      vertexFO->setError(2,100.);

      // add tracks to vertex      
      vertexFO->addTrack( tfo[0], false, true ); // inbound=false, measured=true
      vertexFO->addTrack( tfo[1], false, true );

      // add tracks and vertex to fitter
      for (int i=0; i<2; i++) {
	_opalFit->addFitObject(tfo[i]);
      }
      _opalFit->addFitObject( vertexFO );

      // add a 3d vertex constraint
      vertexFO->addConstraints( *_opalFit, int(VertexFitObject::VXYZ) );

      // prepare for fit
      vertexFO->initForFit();
      
      // do the fit       
      _opalFit->fit();
      
      // get the results
      int    ierr =  _opalFit->getError();
      double fprob = _opalFit->getProbability();
      double chi2 =  _opalFit->getChi2();
      int    dof =   _opalFit->getDoF();
      int    nit =   _opalFit->getIterations();
      
      cout << "done : ierr " << ierr << " fitProb " << fprob << " chisq " << chi2 << " nDOF " << dof << " nIter " << nit << endl;

// Vertex FO error matrix
      for (unsigned int i=0; i<3; i++){
          for (unsigned int j=i; j<3; j++){
               cout << "Cov " << i << " " << j << " " << vertexFO->getCov(i,j) << endl;
          }
      }

      ThreeVector vertex;
      vertexFO->getVertexEx(vertex);
      cout << "Vertex coordinates " << vertex.getX() << " " << vertex.getY() << " " << vertex.getZ() << endl;
      cout << "           errors  " << vertexFO->getError(0) << " " << vertexFO->getError(1) << " " << vertexFO->getError(2) << endl;



    }

  }
  catch(DataNotAvailableException &e) {};

  return;
}



void SimpleVertexProcessor::check( LCEvent * evt ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void SimpleVertexProcessor::end(){
  std::cout << "SimpleVertexProcessor::end()  " << name()
            << std::endl ;
}

