#include "G4Event.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"

#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"

#include "G4UImanager.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4VHitsCollection.hh"

#include "HitSegment.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "EventActionMessenger.hh"

// This class provides user hooks at the beginning and end of an event.

// Constructor
//
EventAction::EventAction(RunAction* run) 
             : printModulo(100), 
               drawFlag("all"),
               nSelected(0),
               iSelected(0),
               debugStarted(false)
{
    runAction = run;
    eventMessenger = new EventActionMessenger(this);

    UI = G4UImanager::GetUIpointer();
    selectedEvents.clear();

}

// Destructor
//
EventAction::~EventAction() {}

// Called at beginning of each event.
//
void EventAction::BeginOfEventAction(const G4Event* aEvent) {

    runAction->ClearVariable();

    // Get Geant4 event number for this event
    //
    G4int eventNumber = aEvent->GetEventID();

    // Print out event number every printModulo events
    //
    if (eventNumber % printModulo == 0)
       G4cout << "\n---> Begin of event: " << eventNumber << G4endl;

    if ( nSelected > 0 ) {
       for (G4int i=iSelected; i<nSelected; i++) {
         if (eventNumber == selectedEvents[i]) {
            iSelected = i++;
            debugStarted = true;
//            G4RunManager::GetRunManager()->rndmSaveThisEvent();
            UI->ApplyCommand("/random/saveThisEvent");
            UI->ApplyCommand("/tracking/verbose 1");
            break;
         }
       }
    }
          
}

G4bool EventAction::SortPredicate(const HitSegment* a,
                                  const HitSegment* b) {

    return a->GetParentID() > b->GetParentID();

}

// Called at end of each event.
//
void EventAction::EndOfEventAction(const G4Event* aEvent) {

    // Get the hit collection for this event.
    G4HCofThisEvent* HCofEvent = aEvent->GetHCofThisEvent();

    // There is no collected information for this event, so exit.
    if(!HCofEvent) return;

    // Get a pointer to the sensitive detector manager
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    // Get the table of hit collections (one collection for each
    // sensitive detector) from the sensitive detector manager
    G4HCtable *HCT = SDman->GetHCtable();

    // Loop through the entries in the table
    for (G4int i=0; i<HCT->entries(); ++i) {

        // Get the name of the ith sensitive detector
        G4String SDname = HCT->GetSDname(i);

        // Get the name of the ith sensitive detector's hit
        // collection
        G4String HCname = HCT->GetHCname(i);

	G4int volumeID[3];

        if (SDname == "crystalSD") { volumeID[0] = 5;
        } else if (SDname == "scintSD") { volumeID[0] = 1; }
	else if (SDname =="wireChamberSD") { volumeID[0] = 2;}
	else if (SDname =="zirconiumSD") {volumeID[0] = 3;}   // DvB addition

	//        G4cout << SDname << " " << volumeID[0] << G4endl;

        // Get the unique id number of the ith hit collection
        G4int fHCID = SDman->GetCollectionID(SDname+"/"+HCname);

        // Get the hit collection corresponding to the unique id above
        // from the collection of hit collections in this event
        HitSegmentCollection* fhits =
                 dynamic_cast<HitSegmentCollection*>(HCofEvent->GetHC(fHCID));

        // If there is no hits collection...continue
        if (!fhits) continue;

        // If the hits collection has no entries...continue
        if (fhits->GetSize()<1) continue;

        // Allocate hit pointers repeatedly used within the loops.
        HitSegment* fhit = NULL;

        // Get a hit out of the hit collection and cast it as a
        // HitSegment. If the cast is not sucessful we have
        // a hit collection of hits that are not of type
        // HitSegment.

        fhit = dynamic_cast<HitSegment*>(fhits->GetHit(0));

        // If there are no HitSegment hits...continue.
        // The collection was mislabeled.
        if (!fhit) continue;

        // Ask how many hits there are
        const G4int numberOfHits = fhits->GetSize();

        std::vector<HitSegment*>* vhits = fhits->GetVector();
        std::stable_sort(vhits->begin(),
                         vhits->end(),
                         &EventAction::SortPredicate);

        G4int numberOfMergeHits = 0;

        HitSegment* jhit = NULL;
        HitSegment* khit = NULL;

        // Loop through the hits in the G4VHitsCollection and
        // merge delta rays
        for (G4int j = 0; j < numberOfHits; j++) {

            jhit = dynamic_cast<HitSegment*>(fhits->GetHit(j));

            for (G4int k = (j + 1); k < numberOfHits; k++) {

                khit = dynamic_cast<HitSegment*>(fhits->GetHit(k));

                if (jhit->GetVolumeID() != khit->GetVolumeID()) continue;

                if (jhit->CheckMergeCondition(khit)) {
                    *khit += *jhit;
                    jhit->SetMergedFlag(true);
                    break;
                }
            }
        }

        // Rearrange the array to place merged hits in front and
        // remove the remaining hits.
        for (G4int j = 0; j < numberOfHits; j++) {

            jhit = dynamic_cast<HitSegment*>(fhits->GetHit(j));

            if (!jhit->GetMergedFlag()) {

                khit = dynamic_cast<HitSegment*>
                                          (fhits->GetHit(numberOfMergeHits));

                *khit = *jhit;

//                G4cout << "numberOfMergeHits:         " << numberOfMergeHits
//                << " TrackID:   " << khit->GetTrackID()
//                << " ParentID: " << khit->GetParentID()
//                << " " << khit->GetVolumeID().GetVolume(0)->GetName()
//                << " " << khit->GetVolumeID().GetReplicaNumber(0)
//                << " " << khit->GetVolumeID().GetVolume(1)->GetName()
//                << " " << khit->GetVolumeID().GetReplicaNumber(1)
//                << " " << khit->GetVolumeID().GetVolume(2)->GetName()
//                << " " << khit->GetVolumeID().GetReplicaNumber(2) << G4endl;

                numberOfMergeHits++;
            }
        }

        // Remove excess hits in the array
        if (numberOfMergeHits < numberOfHits) {
            for (G4int j = numberOfMergeHits; j < numberOfHits; j++){

                delete (*vhits)[j];

            }

            vhits->erase(vhits->begin()+numberOfMergeHits,vhits->end());
        }

        // Fill the Ntuple

//        G4cout << "numberOfMergeHits: " << numberOfMergeHits << G4endl; 

        for (G4int j = 0; j < numberOfMergeHits; j++) {
            jhit = dynamic_cast<HitSegment*>(fhits->GetHit(j));
            if (SDname == "crystalSD") {

	      //	        G4cout << (jhit->GetVolumeID().GetVolume())->GetName()<< G4endl;
               if ((jhit->GetVolumeID().GetVolume())->GetName() == "NaI") {
		 volumeID[1] = jhit->GetVolumeID().GetReplicaNumber(0);
		 volumeID[2] = 100;
               } 
	       
	       else if((jhit->GetVolumeID().GetVolume())->GetName() == "front"){//chloe's addition
		 volumeID[1] = 175;//chloe's addition
		 volumeID[2] = 175;//chloe's addition
		 // G4cout << (jhit->GetVolumeID().GetVolume())->GetName()<< G4endl;
		 //	 G4cout << volumeID[0] << "  " << volumeID[1] << "  " << volumeID[2] << G4endl;
	       }
	       
	       else {
		 volumeID[1] = jhit->GetVolumeID().GetReplicaNumber(0);
		 volumeID[2] = jhit->GetVolumeID().GetReplicaNumber(1);
               }
	       //chloe's addition
	    }


	    else if (SDname == "scintSD") {
	      volumeID[1] = jhit->GetVolumeID().GetReplicaNumber(0);
	      volumeID[2] = 200;
            } else if (SDname == "wireChamberSD"){
	      
	      //    G4cout << (jhit->GetVolumeID().GetVolume())->GetName()<< G4endl;
	      
	      //  volumeID[1] = jhit->GetVolumeID().GetReplicaNumber(0);
	      // volumeID[2] = jhit->GetVolumeID().GetReplicaNumber(1);   
	      
	      if((jhit->GetVolumeID().GetVolume())->GetName() == "sistripSistrip"){
	      volumeID[1] = jhit->GetVolumeID().GetReplicaNumber(0);
	      volumeID[2] = jhit->GetVolumeID().GetReplicaNumber(3)*10+1+jhit->GetVolumeID().GetReplicaNumber(1);}
	      else if  ((jhit->GetVolumeID().GetVolume())->GetName() == "Cellb" ) 
		{
		  volumeID[1] = jhit->GetVolumeID().GetReplicaNumber(0);
		  volumeID[2] = 301+jhit->GetVolumeID().GetReplicaNumber(1);
		}
	      else if((jhit->GetVolumeID().GetVolume())->GetName() == "Cell" )
		{
		  volumeID[1] = jhit->GetVolumeID().GetReplicaNumber(0);
		  volumeID[2] = 401+jhit->GetVolumeID().GetReplicaNumber(1);
		}
	      
	    }
	    

	    
//            G4cout << "EventID: " << aEvent->GetEventID() << G4endl;
//            G4cout << "volumeID[0]: " << volumeID[0] << G4endl;
//            G4cout << "volumeID[1]: " << volumeID[1] << G4endl;
//            G4cout << "volumeID[2]: " << volumeID[2] << G4endl;
//            G4cout << "Energy:      " << jhit->GetEnergyDeposit()/MeV << "MeV" << G4endl;
//            G4cout << "GetStartX:   " << jhit->GetStartX()/cm << "cm" << G4endl;
//            G4cout << "GetStartY:   " << jhit->GetStartY()/cm << "cm" << G4endl;
//            G4cout << "GetStartZ:   " << jhit->GetStartZ()/cm << "cm" << G4endl;
//            G4cout << "GetStartT:   " << jhit->GetStartT()/ns << "ns" << G4endl;
//            G4cout << "GetStopX:    " << jhit->GetStopX()/cm << "cm" << G4endl;
//            G4cout << "GetStopY:    " << jhit->GetStopY()/cm << "cm" << G4endl;
//            G4cout << "GetStopZ:    " << jhit->GetStopZ()/cm << "cm" << G4endl;
//            G4cout << "GetStopT:    " << jhit->GetStopT()/ns << "ns" << G4endl;
//            G4cout << "Ebirk:       " << jhit->GetEbirk()/MeV << "MeV" << G4endl;
//            G4cout << "PDG:         " << jhit->GetPDG() << G4end;
//            G4cout << "CreatorFlag: " << jhit->GetCreatorFlag() << G4endl;

            runAction->FillTuple(aEvent->GetEventID(),
                                 volumeID[0], volumeID[1], volumeID[2],
                                 jhit->GetEnergyDeposit(),
                                 jhit->GetStartX(),
                                 jhit->GetStartY(),
                                 jhit->GetStartZ(),
                                 jhit->GetStartT(),
                                 jhit->GetStopX(),
                                 jhit->GetStopY(),
                                 jhit->GetStopZ(),
                                 jhit->GetStopT(),
                                 jhit->GetEbirk(),
                                 jhit->GetPDG(),
                                 jhit->GetCreatorFlag());
#ifdef G4VIS_USE
            jhit->Draw();
#endif
        }
#ifdef G4VIS_USE
//        if (fhits) fhits->DrawAllHits();
#endif
    }

    runAction->SetEventID(aEvent->GetEventID());

    runAction->FillTree();
    runAction->FillUntriggeredHistos();
    runAction->FillTargetStopHistos();

    if (debugStarted) {
       debugStarted = false;
       UI->ApplyCommand("/tracking/verbose 0");
    }

#ifdef G4VIS_USE
    G4TrajectoryContainer* trajectoryContainer =
                                             aEvent->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
    for (G4int i=0; i<n_trajectories; i++) {
        G4Trajectory* trj = (G4Trajectory *)
                                   ((*(aEvent->GetTrajectoryContainer()))[i]);
        if (drawFlag == "all") trj->DrawTrajectory(1000);
        else if ((drawFlag == "charged") && (trj->GetCharge() != 0.))
                                                    trj->DrawTrajectory(1000);
        else if ((drawFlag == "neutral") && (trj->GetCharge() == 0.))
                                                    trj->DrawTrajectory(1000);
        else if ((drawFlag == "charged+n") &&
              ((trj->GetCharge() != 0.) ||
               (trj->GetCharge() == 0. && trj->GetParticleName()=="neutron")))
                                                    trj->DrawTrajectory(1000);
    }
#endif
}
