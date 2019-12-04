#include "G4TouchableHandle.hh"

#include "VolumeID.hh"

VolumeID::VolumeID() { }

VolumeID::VolumeID(const G4TouchableHandle& handle) {

    if(fVolumes.size()>0)fVolumes.clear();

    G4int historyDepth = handle->GetHistoryDepth();

    for (G4int i=0; i< historyDepth; ++i) {
        AddVolume(handle->GetVolume(i),handle->GetReplicaNumber(i));
    }

}

VolumeID::~VolumeID() { }

VolumeID& VolumeID::operator = (const G4TouchableHandle& handle) {

    if(fVolumes.size()>0)fVolumes.clear();

    G4int historyDepth = handle->GetHistoryDepth();

    for (G4int i=0; i< historyDepth; ++i) {
        AddVolume(handle->GetVolume(i),handle->GetReplicaNumber(i));
    }
    
    return *this;
}

VolumeID& VolumeID::operator = (const VolumeID& id) {

    if(fVolumes.size()>0)fVolumes.clear();

    for (VolumeID::Volumes::const_iterator i = id.fVolumes.begin();
         i != id.fVolumes.end();
         ++i) {
        AddVolume(i->fPhysVol, i->fReplica);
    }
    
    return *this;
}

void VolumeID::AddVolume(G4VPhysicalVolume* physvol, G4int replica) {

    Volume vol;
    vol.fPhysVol = physvol;
    vol.fReplica = replica;
    fVolumes.push_back(vol);

}

G4bool operator == (const VolumeID& x, const VolumeID& y) {

    if (x.fVolumes.size() != y.fVolumes.size()) return false;

    VolumeID::Volumes::const_iterator a;
    VolumeID::Volumes::const_iterator b;

    for (a = x.fVolumes.begin(), b = y.fVolumes.begin();
         a != x.fVolumes.end();
         ++a, ++b) {
        if (a->fPhysVol != b->fPhysVol) return false;
        if (a->fReplica != b->fReplica) return false;
    }
    return true;
}

G4bool operator == (const VolumeID& x, const G4TouchableHandle& y) {

    unsigned int historyDepth = y->GetHistoryDepth();

    if (x.fVolumes.size() != historyDepth) return false;

    VolumeID::Volumes::const_iterator a;
    unsigned int i;
    for (a = x.fVolumes.begin(), i = 0;
         a != x.fVolumes.end() && i < historyDepth;
         ++a, ++i) {

        if (a->fPhysVol != y->GetVolume(i)) {
            return false;
        }
        if (a->fReplica != y->GetReplicaNumber(i)) {
            return false;
        }
    }
    return true;
}

G4bool operator == (const G4TouchableHandle& x, const VolumeID& y) {
    return (y==x);
}

G4bool operator != (const VolumeID& x, const VolumeID& y) {
    return !(x==y);
}

G4bool operator != (const VolumeID& x, const G4TouchableHandle& y) {
    return !(x==y);
}

G4bool operator != (const G4TouchableHandle& x, const VolumeID& y) {
    return !(y==x);
}
/*
G4bool operator < (const VolumeID& x, const VolumeID& y) {

    if (x.fVolumes.size() < y.fVolumes.size()) return true;
    if (x.fVolumes.size() > y.fVolumes.size()) return false;

    VolumeID::Volumes::const_iterator a;
    VolumeID::Volumes::const_iterator b;
    for (a = x.fVolumes.begin(), b = y.fVolumes.begin();
         a != x.fVolumes.end();
         ++a, ++b) {
        if (a->fPhysVol < b->fPhysVol) return true;
        if (a->fPhysVol > b->fPhysVol) return false;
        if (a->fPhysVol == b->fPhysVol) {
            if (a->fReplica < b->fReplica) return true;
            if (a->fReplica > b->fReplica) return false;
        }
    }
    return false;
}
*/
