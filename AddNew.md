# Adding new sensitve detectors and corresponding hits

CaTs is easily extendable. To add a new Detector (newDet) and corresponding Hits (newDetHit) the following new files have to be provided:
  * newDetSD.hh(cc)  : describes the new sensitive Detector class.
  * newDetHit.hh(cc) : describes the Hits created by the new sensitive detectors.

one should also provide a file:

  * readnewDetHit.cc : main executable that demonstrates how to access the new Hits and makes a few histograms. 

In addition the following files have to be modified:

  * DetectorConstruction.cc  : check for new keyword in gdml file. Attach newDetSD to a corresponding logical volume when keyword is provided.
  * EventAction.cc: make sure that Hits are written out at the end of the Event.
  * CaTSClasses.hh : needed for Hits root persistency.
  * selection.xml  : needed for Hits root persistency.
  * CMakeLists.txt : make sure the new executable readnewDetHit is build and installed in the bin area. 
