#include "TrackResponseAdjuster.h"

TrackResponseAdjuster aTrackResponseAdjuster;


TrackResponseAdjuster::TrackResponseAdjuster() : Processor("TrackResponseAdjuster") {

  // modify processor description
  _description = "Smears Charged Particles from filtered Monte Carlo Particle collection" ;


  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "Printing" ,
                              "Print certain messages"  ,
                              _printing,
                               (int)5 ) ;

  registerProcessorParameter("D0ErrorScaleFactor",
			      "D0 Error Scale Factor" ,
			      _D0ErrorScaleFactor,
			      (double)1.0);

   registerProcessorParameter("Z0ErrorScaleFactor",
			      "Z0 Error Scale Factor",
			      _Z0ErrorScaleFactor,
			      (double)1.0);

  registerProcessorParameter("OmegaErrorScaleFactor",
			      "Omega Error Scale Factor",
			      _OmegaErrorScaleFactor,
			      (double)1.0);

  registerProcessorParameter("PhiErrorScaleFactor",
			      "Phi Error Scale Factor",
			      _PhiErrorScaleFactor,
			      (double)1.0);

  registerProcessorParameter("TanLambdaErrorScaleFactor",
			      "TanLambda Error Scale Factor",
			      _TanLambdaErrorScaleFactor,
			      (double)1.0);

   std::string inputTrackCollectionName = "MarlinTrkTracks";
  registerInputCollection( LCIO::TRACK,
                                "InputTrackCollectionName" ,
                                "Input Track Collection Name " ,
                                _inputTrackCollectionName,
                                inputTrackCollectionName);  

  std::string outputTrackCollectionName = "CalibratedTracks";
  registerOutputCollection( LCIO::TRACK,
			    "OutputTrackCollectionName",
			    "Output Track Collection Name",
			    _outputTrackCollectionName,
			    outputTrackCollectionName);
}


void TrackResponseAdjuster::init() {

  streamlog_out(DEBUG) << "   init called  "
                       << std::endl ;

 //rng = new TRandom3();
  // usually a good idea to
  printParameters() ;
  nEvt = 0;

//  gROOT->ProcessLine("#include <vector>");
}

void TrackResponseAdjuster::processRunHeader( LCRunHeader* run) {

}

bool TrackResponseAdjuster::FindTracks( LCEvent* evt ) {

  bool tf = false;

  // clear old vector
  _trackvec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){
    if(*name==_inputTrackCollectionName){
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
        Track* track = dynamic_cast<Track*>(col->getElementAt(i));
        _trackvec.push_back(track);
      }
    }
  }

  if(_printing>1)std::cout << "FindTracks : " << tf << std::endl;

  return tf;
}


double TrackResponseAdjuster::safeAcos(double x){
	if (x < -1.0) x = -1.0 ;
	else if (x > 1.0) x = 1.0 ;
	return acos (x) ;
 }

void TrackResponseAdjuster::processEvent( LCEvent * evt ) {
 //FindMCParticles(evt);
// = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
 FindTracks(evt);
  LCCollectionVec * calrectrk = new LCCollectionVec(LCIO::TRACK);


  streamlog_out(MESSAGE) << " start processing event " << std::endl;
 

 std::cout<<"track vector size "<<_trackvec.size()<<std::endl;


	for(int i =0; i<_trackvec.size(); i++){
		
			TrackImpl* CalTrack = new TrackImpl();

			//FloatVec oldcov = _trackvec.at(i)->getCovMatrix();
			
			std::vector<float> newcov;
			//deep copy
			for(int j=0; j<_trackvec.at(i)->getCovMatrix().size(); j++){
				newcov.push_back((double)_trackvec.at(i)->getCovMatrix().at(j));
			}
			
			//d0 d0 
			newcov.at(0) = pow(_D0ErrorScaleFactor*std::sqrt(newcov.at(0)),2);
			//d0 phi
			newcov.at(1) = _D0ErrorScaleFactor*_PhiErrorScaleFactor*newcov.at(1);
			//phi phi
			newcov.at(2) = pow(_PhiErrorScaleFactor*std::sqrt(newcov.at(2)),2);
			//d0 ome
			newcov.at(3) = _D0ErrorScaleFactor*_OmegaErrorScaleFactor*newcov.at(3);
			//phi ome
			newcov.at(4) = _OmegaErrorScaleFactor*_PhiErrorScaleFactor*newcov.at(4);
			//ome ome
			newcov.at(5) = pow(_OmegaErrorScaleFactor*std::sqrt(newcov.at(5)),2);
			//d0 z0
			newcov.at(6) = _D0ErrorScaleFactor*_Z0ErrorScaleFactor*newcov.at(6);
			//phi z0
			newcov.at(7) = _Z0ErrorScaleFactor*_PhiErrorScaleFactor*newcov.at(7);
			//ome z0
			newcov.at(8) = _Z0ErrorScaleFactor*_OmegaErrorScaleFactor*newcov.at(8);
			//z0 z0
			newcov.at(9) = pow(_Z0ErrorScaleFactor*std::sqrt(newcov.at(9)),2);
			//tl d0
			newcov.at(10) = _TanLambdaErrorScaleFactor*_D0ErrorScaleFactor*newcov.at(10);
			//tl phi
			newcov.at(11) = _TanLambdaErrorScaleFactor*_PhiErrorScaleFactor*newcov.at(11);
			//tl ome
			newcov.at(12) = _TanLambdaErrorScaleFactor*_OmegaErrorScaleFactor*newcov.at(12);
			//tl z0
			newcov.at(13) = _TanLambdaErrorScaleFactor*_Z0ErrorScaleFactor*newcov.at(13);
			//tl tl
			newcov.at(14) = pow(_TanLambdaErrorScaleFactor*std::sqrt(newcov.at(14)),2);
			
		
			CalTrack->setD0(_trackvec.at(i)->getD0()); //Impact parameter in r-phi
			CalTrack->setPhi(_trackvec.at(i)->getPhi()); //phi of track at reference point (primary vertex)
			CalTrack->setOmega(_trackvec.at(i)->getOmega());// signed curvature in 1/mm 
			CalTrack->setZ0(_trackvec.at(i)->getZ0()); //Impact parameeter in r-z
			CalTrack->setTanLambda(_trackvec.at(i)->getTanLambda());// dip of the track in r-z at primary vertex
			CalTrack->setCovMatrix(newcov);

		    //    fastreccol->addElement( fastRecoPart );
			calrectrk->addElement( CalTrack );

			if(_printing>1){
				
				std::cout<<"Event No. :"<< nEvt <<std::endl;
				std::cout<<"old matrix"<<std::endl;
				for(int j=0; j<_trackvec.at(i)->getCovMatrix().size(); j++){
					std::cout<<_trackvec.at(i)->getCovMatrix().at(j)<<" ";
				}
				std::cout<<std::endl;
				std::cout<<"new matrix"<<std::endl;
				for(int j=0; j<newcov.size(); j++){
					std::cout<<newcov.at(j)<<" ";
				}
				std::cout<<std::endl;
				
			}
	
			
 		
	}

  nEvt++;

  // Add new collection to event

  evt->addCollection(calrectrk , _outputTrackCollectionName.c_str() ); 
 std::cout << "======================================== event " << nEvt << std::endl ;
//delete fastreccol;
//delete fastrectrk;
}


//void TrackResponseAdjuster::check( LCEvent * evt ) {

//}


void TrackResponseAdjuster::end(){

}

