#pragma once

#include <vector>
#include <iostream>

#include "PETU4SensorIdentifier.h"
#include "G4PVPlacement.hh"

using namespace std;

struct MyU4SensorIdentifier : public U4SensorIdentifier 
{
 //   int getGlobalIdentity(const G4VPhysicalVolume* pv ) const ; 
    int getInstanceIdentity(const G4VPhysicalVolume* instance_outer_pv ) const ; 
    static void FindSD_r( std::vector<const G4VPhysicalVolume*>& sdpv , const G4VPhysicalVolume* pv, int depth );  
}; 


//inline int MyU4SensorIdentifier::getGlobalIdentity( const G4VPhysicalVolume* ) const 
//{
    
//    return -1;
//}

/**
MyU4SensorIdentifier::getInstanceIdentity
---------------------------------------------------

Canonically used from U4Tree::identifySensitiveInstances

**/

inline int MyU4SensorIdentifier::getInstanceIdentity( const G4VPhysicalVolume* instance_outer_pv ) const 
{
//    const G4PVPlacement* pvp = dynamic_cast<const G4PVPlacement*>(instance_outer_pv) ;
//    int copyno = pvp ? pvp->GetCopyNo() : -1 ;

//    std::vector<const G4VPhysicalVolume*> sdpv ; 
//    FindSD_r(sdpv, instance_outer_pv, 0 );  

//    unsigned num_sd = sdpv.size() ; 
//    int sensor_id = num_sd == 0 ? -1 : copyno ; 

    //bool dump = copyno < 10 ; 
//    bool dump = false ; 
//    if(dump) std::cout 
//        << "MyU4SensorIdentifier::getIdentity" 
//        << " copyno " << copyno
//        << " num_sd " << num_sd
//        << " sensor_id " << sensor_id 
//        << std::endl 
//        ;      


//    return sensor_id ; 

    int sensor_id = -1;
    const G4PVPlacement* pvp = dynamic_cast<const G4PVPlacement*>(instance_outer_pv) ;
    int copyno = pvp ? pvp->GetCopyNo() : -1 ;
  
    
    G4LogicalVolume *lv = instance_outer_pv->GetLogicalVolume();
    G4VSolid *solid = lv->GetSolid();
    G4String name = solid->GetName();
    char sdname[256];
    if(name.find("lysoBox") != string::npos) {
	 sensor_id =  copyno ;
    } 

    return sensor_id ;
}

inline void MyU4SensorIdentifier::FindSD_r( std::vector<const G4VPhysicalVolume*>& sdpv , const G4VPhysicalVolume* pv, int depth )
{
    const G4LogicalVolume* lv = pv->GetLogicalVolume() ;
    G4VSensitiveDetector* sd = lv->GetSensitiveDetector() ;
    if(sd) sdpv.push_back(pv); 
    for (size_t i=0 ; i < size_t(lv->GetNoDaughters()) ; i++ ) FindSD_r( sdpv, lv->GetDaughter(i), depth+1 );
}


