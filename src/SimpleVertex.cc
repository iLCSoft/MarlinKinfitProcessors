#include "SimpleVertex.h"
#include <iostream>
#include <map>
#include <cassert>

#include <marlin/Global.h>
#include "lcio.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>
#endif

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/TrackState.h>
#include <EVENT/MCParticle.h>

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
  
  // register processor parameters
  registerProcessorParameter( "traceall" ,
                              "set true if every event should be traced",
                              _traceall,
                              (bool)false);
                              
  registerProcessorParameter( "ievttrace" ,
                              "number of individual event to be traced",
                              _ievttrace,
                              (int)0);
                              
  _opalFit = new OPALFitterGSL();
}


void SimpleVertexProcessor::init() {
  message<MESSAGE>( log()  << "hello from SimpleVertexProcessor::init");
  printParameters() ;
  return;
}

void SimpleVertexProcessor::processRunHeader( LCRunHeader* ) {
  message<MESSAGE>( log()  << "hello from SimpleVertexProcessor::processRunHeader");
}

void SimpleVertexProcessor::processEvent( LCEvent * evt ) {

   message<MESSAGE>( log() 
		      << " processing event " << evt->getEventNumber() 
		      << "  in run "          << evt->getRunNumber() 
		      ) ;
                      
#ifdef MARLIN_USE_AIDA
  static AIDA::IHistogram1D* hFitProb;
  static AIDA::IHistogram1D* hNIt;
  static AIDA::IHistogram1D* hFitError;
  static AIDA::IHistogram1D* hErrVtxX; 
  static AIDA::IHistogram1D* hErrVtxY; 
  static AIDA::IHistogram1D* hErrVtxZ; 
  static AIDA::IHistogram1D* hPullVtxX; 
  static AIDA::IHistogram1D* hPullVtxY; 
  static AIDA::IHistogram1D* hPullVtxZ; 
  
  if( isFirstEvent() ) { 
    hFitProb = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hFitProb", "fit probability", 100, 0., 1. ) ; 
    hNIt = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hNIt", "number of iterations", 201, -0.5, 200.5 ) ; 
    hFitError = AIDAProcessor::histogramFactory(this)->createHistogram1D( "hFitError", "Error flag", 10, -0.5, 9.5 ) ; 

    hErrVtxX =  AIDAProcessor::histogramFactory(this)->createHistogram1D( "hErrVtxX", "error of fitted vertex x / mm", 100, 0, 3. ) ;    
    hErrVtxY =  AIDAProcessor::histogramFactory(this)->createHistogram1D( "hErrVtxY", "error of fitted vertex y / mm", 100, 0, 3. ) ;    
    hErrVtxZ =  AIDAProcessor::histogramFactory(this)->createHistogram1D( "hErrVtxZ", "error of fitted vertex z / mm", 100, 0, 3. ) ;    
   
    hPullVtxX =  AIDAProcessor::histogramFactory(this)->createHistogram1D( "hPullVtxX", "pull wrt truth of fitted vertex x", 100, -5., 5. ) ;    
    hPullVtxY =  AIDAProcessor::histogramFactory(this)->createHistogram1D( "hPullVtxY", "pull wrt truth of fitted vertex y", 100, -5., 5. ) ;    
    hPullVtxZ =  AIDAProcessor::histogramFactory(this)->createHistogram1D( "hPullVtxZ", "pull wrt truth of fitted vertex z", 100, -5., 5. ) ;    
  }  
#endif   

  _opalFit->reset();
     
  TextTracer tracer (std::cout);
  if (evt->getEventNumber() == _ievttrace || _traceall) {
    _opalFit->setTracer (tracer);
  }
  else {
    _opalFit->setTracer (0);
  }  
  
  // designed to work with single K0short events
  // fit tracks from kshort decays to a vertex

  const double mass_pi = 0.14; // assumed mass of tracks
  
  try {
    LCCollection* col_trk=evt->getCollection("MarlinTrkTracks");
    LCCollection* col_MCP=evt->getCollection("MCParticlesSkimmed");
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
          message<WARNING>( log()  << "no track state found @ first hit, using default track parameters" );
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
      
      message<MESSAGE>( log() << "done : ierr " << ierr << " fitProb " << fprob << " chisq " << chi2 << " nDOF " << dof << " nIter " << nit );

      hFitError->fill( ierr ) ;
      if (ierr == 0) {
        hFitProb->fill( fprob ) ;
        hNIt->fill( nit ) ;
       
// Vertex FO error matrix
        for (unsigned int i=0; i<3; i++){
          for (unsigned int j=i; j<3; j++){
               message<MESSAGE>( log() << "Cov " << i << " " << j << " " << vertexFO->getCov(i,j) );
          }
        }

        ThreeVector vertex;
        vertexFO->getVertexEx(vertex);
        message<MESSAGE>( log() << "Vertex coordinates " << vertex.getX() << " " << vertex.getY() << " " << vertex.getZ() );
        message<MESSAGE>( log() << "           errors  " << vertexFO->getError(0) << " " << vertexFO->getError(1) << " " << vertexFO->getError(2) );
        hErrVtxX->fill (vertexFO->getError(0));
        hErrVtxY->fill (vertexFO->getError(1)); 
        hErrVtxZ->fill (vertexFO->getError(2));
        
        double trueVtx[3] = {0.,0.,0.};
        for (int imcp = 0; imcp < col_MCP->getNumberOfElements(); imcp++) {
          MCParticle* mcp = dynamic_cast < MCParticle* > ( col_MCP->getElementAt(imcp) );
          if (abs(mcp->getPDG()) == 310) {
             trueVtx[0] = mcp->getEndpoint()[0];
             trueVtx[1] = mcp->getEndpoint()[1];
             trueVtx[2] = mcp->getEndpoint()[2];
             message<MESSAGE>( log() << "found true Kaon decay vertex with coordinates (x,y,z) = (" << trueVtx[0] << ", " << trueVtx[1] << ", " << trueVtx[2] << ")" );
          } 
        }
        hPullVtxX->fill ((vertex.getX()-trueVtx[0])/vertexFO->getError(0));
        hPullVtxY->fill ((vertex.getY()-trueVtx[1])/vertexFO->getError(1)); 
        hPullVtxZ->fill ((vertex.getZ()-trueVtx[2])/vertexFO->getError(2));
        
      }  // ierr == 0

    }

  }
  catch(DataNotAvailableException &e) {};

  return;
}



void SimpleVertexProcessor::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void SimpleVertexProcessor::end(){
   message<MESSAGE>( log() << "SimpleVertexProcessor::end()  " << name()) ;
}

