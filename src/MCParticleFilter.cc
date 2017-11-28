#include "MCParticleFilter.h"

#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>


MCParticleFilter aMCParticleFilter ;


MCParticleFilter::MCParticleFilter() : Processor("MCParticleFilter") {
  
  // modify processor description
  _description = "MCParticleFilter does whatever it does ..." ;

  
  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "Printing" , 
                              "Print certain messages"  ,
                              _printing,
                               (int)1 ) ;

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle input collection"  ,
			   _mcParticleCollectionName ,
			   std::string("MCParticle") ) ;

  std::string outputParticleCollectionName = "MCDecayParticles";
  registerOutputCollection( LCIO::MCPARTICLE,
                             "OutputParticleCollectionName" , 
                             "Output Particle Collection Name "  ,
                             _outputParticleCollectionName,
                             outputParticleCollectionName);
}


void MCParticleFilter::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;

  // usually a good idea to
  printParameters() ;
  nEvt = 0;

}

void MCParticleFilter::processRunHeader( LCRunHeader* run) { 
    
}
 

void MCParticleFilter::processEvent( LCEvent * evt ) {

  // Make a new vector of particles
  LCCollectionVec * mcparcol = new LCCollectionVec(LCIO::MCPARTICLE);
  mcparcol->setSubset(true);      

  streamlog_out(MESSAGE) << " start processing event " << std::endl;
  cout << "======================================== event " << nEvt << std::endl ;

  streamlog_out(DEBUG) << " mcp iterator and navigator " << std::endl;
  LCIterator<MCParticle> mcpIt( evt, _mcParticleCollectionName ) ;
  streamlog_out(DEBUG) << " got mcpIt with length of " << mcpIt.size() << std::endl;
  
  //----------------------------------------------------------------------------------------------------------------------------
  // loop over MCParticles
  //----------------------------------------------------------------------------------------------------------------------------
  
  int nmcp = 0;
 
  while( MCParticle* mcp = mcpIt.next()  ) {
    
    if (!mcp) continue;

    streamlog_out(DEBUG) << " MCParticle " << nmcp << " pdgid = " << mcp->getPDG() << " E = " << mcp->getEnergy() << " genStatus = " << mcp->getGeneratorStatus() << std::endl; 
    streamlog_out(DEBUG) << " Vertex   = " << mcp->getVertex()[0] << " " << mcp->getVertex()[1] << " " << mcp->getVertex()[2] << std::endl;
    streamlog_out(DEBUG) << " Momentum = " << mcp->getMomentum()[0] << " " << mcp->getMomentum()[1] << " " << mcp->getMomentum()[2] << std::endl; 

    double rsq = mcp->getVertex()[0]*mcp->getVertex()[0] + mcp->getVertex()[1]*mcp->getVertex()[1] + mcp->getVertex()[2]*mcp->getVertex()[2];    

    nmcp++;

    // Put conditional statements here to positively select MCParticles

    // Add this MCParticle to the output collection
    // Quick hack - require stable particles for the output collection 
    if(rsq<1.0 && mcp->getGeneratorStatus()==1)mcparcol->addElement(mcp);
       
  } // end loop over MCParticles
  

  nEvt++;

  // Add new collection to event
  evt->addCollection( mcparcol , _outputParticleCollectionName.c_str() );
  
//  cout << "======================================== event " << nEvt << std::endl ;
  
}

void MCParticleFilter::check( LCEvent * evt ) { 

}

void MCParticleFilter::end(){ 

}
