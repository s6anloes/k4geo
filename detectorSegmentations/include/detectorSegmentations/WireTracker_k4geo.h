#ifndef DETSEGMENTATION_WIRETRACKER_H
#define DETSEGMENTATION_WIRETRACKER_H

#include "detectorSegmentations/DRparamBarrel_k4geo.h"
#include "detectorSegmentations/DRparamEndcap_k4geo.h"

#include "DDSegmentation/Segmentation.h"

// DDRec
#include "DDRec/DCH_info.h"

namespace dd4hep {
namespace DDSegmentation {



    class WireTracker_k4geo : public Segmentation {
    public:
        /// default constructor using an arbitrary type
        WireTracker_k4geo(const std::string& aCellEncoding);
        /// Default constructor used by derived classes passing an existing decoder
        WireTracker_k4geo(const BitFieldCoder* decoder);
        /// destructor
        virtual ~WireTracker_k4geo() = default;

        /**  Determine the global position based on the cell ID.
         *   @param[in] aCellId ID of a cell.
         */
        virtual Vector3D position(const CellID& aCellID) const override;

        /**  Determine the cell ID based on the position.
         *   @param[in] aLocalPosition
         *   @param[in] aGlobalPosition
         *   @param[in] aVolumeId ID of the Geant4 volume
         *   return Cell ID.
         */
        virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                            const VolumeID& aVolumeID) const override;

        void setDCHinfo(dd4hep::rec::DCH_info* dch_info) { m_dch_info = dch_info; }


        CellID setCellID(int System, int Superlayer, int Layer, int Nphi, int StereoSign) const;

        inline double phiFromXY(const Vector3D& aposition) const { return std::atan2(aposition.Y, aposition.X) + M_PI; }

        inline const std::string& fieldNameSystem() const { return m_systemId; }
        inline const std::string& fieldNameSuperlayer() const { return m_superlayerId; }
        inline const std::string& fieldNameLayer() const { return m_layerId; }
        inline const std::string& fieldNameNphi() const { return m_nphiId; }
        inline const std::string& fieldNameStereosign() const { return m_stereosignId; }

        inline void setFieldNameSystem(const std::string& fieldName) { m_systemId = fieldName; }
        inline void setFieldNameSuperlayer(const std::string& fieldName) { m_superlayerId = fieldName; }
        inline void setFieldNameLayer(const std::string& fieldName) { m_layerId = fieldName; }
        inline void setFieldNameNphi(const std::string& fieldName) { m_nphiId = fieldName; }
        inline void setFieldNameStereosign(const std::string& fieldName) { m_stereosignId = fieldName; }

    private:
        /// Initialization common to all ctors.
        void commonSetup();

        std::string m_systemId;
        std::string m_superlayerId;
        std::string m_layerId;
        std::string m_nphiId;
        std::string m_stereosignId;

        int m_systemIndex = -1;
        int m_superlayerIndex = -1;
        int m_layerIndex = -1;
        int m_nphiIndex = -1;
        int m_stereosignIndex = -1;

        /// Drift chamber info extension for geometry calculations
        dd4hep::rec::DCH_info* m_dch_info{nullptr};

    };



} // namespace DDSegmentation
} // namespace dd4hep

#endif
