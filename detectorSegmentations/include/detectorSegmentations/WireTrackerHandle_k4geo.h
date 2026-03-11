#ifndef DETECTORSEGMENTATIONS_WIRETRACKERHANDLE_K4GEO_H
#define DETECTORSEGMENTATIONS_WIRETRACKERHANDLE_K4GEO_H

// FCCSW
#include "detectorSegmentations/WireTracker_k4geo.h"

// DD4hep
#include "DD4hep/Segmentations.h"
#include "DD4hep/detail/SegmentationsInterna.h"

namespace dd4hep {

// Forward declarations
class Segmentation;
template <typename T>
class SegmentationWrapper;

typedef Handle<SegmentationWrapper<DDSegmentation::WireTracker_k4geo>> WireTrackerHandle_k4geo;

class WireTracker_k4geo : public WireTrackerHandle_k4geo {
public:
    /// Defintiion of the basic handled object
    typedef WireTrackerHandle_k4geo::Object Object;
public:
    /// Default constructor
    WireTracker_k4geo() = default;
    /// Copy constructor
    WireTracker_k4geo(const WireTracker_k4geo& e) = default;
    /// Copy Constructor from segmentation base object
    WireTracker_k4geo(const Segmentation& e) : Handle<Object>(e) {}
    /// Copy constructor from handle
    WireTracker_k4geo(const Handle<Object>& e) : Handle<Object>(e) {}
    /// Copy constructor from other polymorph/equivalent handle
    template <typename Q>
    WireTracker_k4geo(const Handle<Q>& e) : Handle<Object>(e) {}
    /// Assignment operator
    WireTracker_k4geo& operator=(const WireTracker_k4geo& seg) = default;
    /// Equality operator
    bool operator==(const WireTracker_k4geo& seg) const { return m_element == seg.m_element; }


    /// determine the position based on the cell ID
    inline Position position(const CellID id) const { return Position(access()->implementation->position(id)); }

    /// determine the cell ID based on the position
    inline dd4hep::CellID cellID(const Position& local, const Position& global, const VolumeID volID) const {
        return access()->implementation->cellID(local, global, volID);
    }

    inline CellID setCellID(int System, int Superlayer, int Layer, int Nphi, int StereoSign) const {
        return access()->implementation->setCellID(System, Superlayer, Layer, Nphi, StereoSign);
    }

    inline double phiFromXY(const Position& aposition) const {
        return access()->implementation->phiFromXY(aposition);
    }

    inline void setDCHinfo(dd4hep::rec::DCH_info* dch_info) {
        access()->implementation->setDCHinfo(dch_info);
    }

    inline const std::string& fieldNameSystem() const {
        return access()->implementation->fieldNameSystem();
    }
    inline const std::string& fieldNameSuperlayer() const {
        return access()->implementation->fieldNameSuperlayer();
    }
    inline const std::string& fieldNameLayer() const {
        return access()->implementation->fieldNameLayer();
    }
    inline const std::string& fieldNameNphi() const {
        return access()->implementation->fieldNameNphi();
    }
    inline const std::string& fieldNameStereosign() const {
        return access()->implementation->fieldNameStereosign();
    }

    inline void setFieldNameSystem(const std::string& fieldName) {
        access()->implementation->setFieldNameSystem(fieldName);
    }
    inline void setFieldNameSuperlayer(const std::string& fieldName) {
        access()->implementation->setFieldNameSuperlayer(fieldName);
    }
    inline void setFieldNameLayer(const std::string& fieldName) {
        access()->implementation->setFieldNameLayer(fieldName);
    }
    inline void setFieldNameNphi(const std::string& fieldName) {
        access()->implementation->setFieldNameNphi(fieldName);
    }
    inline void setFieldNameStereosign(const std::string& fieldName) {
        access()->implementation->setFieldNameStereosign(fieldName);
    }

};

}

#endif  // DETECTORSEGMENTATIONS_WIRETRACKERHANDLE_K4GEO_H