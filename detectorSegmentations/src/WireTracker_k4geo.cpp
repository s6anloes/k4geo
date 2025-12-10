#include "detectorSegmentations/WireTracker_k4geo.h"
#include "DD4hep/Printout.h"

namespace dd4hep {
namespace DDSegmentation {

    /// default constructor using an encoding string
    WireTracker_k4geo::WireTracker_k4geo(const std::string& cellEncoding) : Segmentation(cellEncoding) { commonSetup(); }

    WireTracker_k4geo::WireTracker_k4geo(const BitFieldCoder* decoder) : Segmentation(decoder) { commonSetup(); }

    // WireTracker_k4geo::~WireTracker_k4geo() {
    // }

    CellID WireTracker_k4geo::cellID(const Vector3D& /*localPosition*/, const Vector3D& globalPosition, const VolumeID& vID) const {

        CellID cID = vID;

        unsigned int superlayer = decoder()->get(vID, m_superlayerIndex);
        unsigned int layer = decoder()->get(vID, m_layerIndex);

        // global layer number
        unsigned int ilayer = m_dch_info->CalculateILayerFromCellIDFields(layer, superlayer);

        // get coordinates of first wire (wire 0, nphi = 0) centre (z=0) and wire direction in this layer
        auto wire0_pos = m_dch_info->Calculate_wire_z0_point(ilayer, 0);
        auto wire0_dir = m_dch_info->Calculate_wire_vector_ez(ilayer, 0);

        // get the coordinates of wire 0 at the z position of the hit
        double dz = globalPosition.Z;
        auto wire0_at_zhit = wire0_pos + dz * wire0_dir;

        // get the phi of the hit and of the wire 0 at the hit z position
        double phi_hit = phiFromXY(globalPosition);
        double phi_wire0 = phiFromXY(wire0_at_zhit);

        double dphi = phi_hit - phi_wire0;
        if (dphi < 0) {
            dphi += 2 * M_PI;
        }

        double cell_phi_width = m_dch_info->Get_phi_width(ilayer);

        double offset = cell_phi_width * 0.25 * (ilayer % 2); // staggering of layers in phi

        unsigned int nphi = positionToBin(dphi, cell_phi_width, offset);

        decoder()->set(cID, m_nphiIndex, static_cast<VolumeID>(nphi));

        int stereosign = m_dch_info->database.at(ilayer).StereoSign();
        decoder()->set(cID, m_stereosignIndex, static_cast<VolumeID>(stereosign));

        if (ilayer%2 != 0) {
            std::cout << "phi_hit: " << phi_hit << "\n"
                      << "phi_wire0: " << phi_wire0 << "\n"
                      << "dphi: " << dphi << "\n"
                      << "nphi: " << nphi << "\n";
        }

        return cID;
    }

    Vector3D WireTracker_k4geo::position(const CellID& cellId) const {
        return Vector3D(0, 0, 0);
    }

    /// Initialization common to all ctors.
    void WireTracker_k4geo::commonSetup() {
        // define type and description
        _type = "WireTracker_k4geo";
        _description = "WireTracker segmentation based on a cylindrical tracking volume with wires in z direction arranged in radially outgoing layers, potentially with stereo angle";

        // register all necessary parameters
        registerIdentifier("identifier_system", "Cell ID identifier for System", m_systemId, "system");
        registerIdentifier("identifier_superlayer", "Cell ID identifier for Superlayer", m_superlayerId, "superlayer");
        registerIdentifier("identifier_layer", "Cell ID identifier for Layer", m_layerId, "layer");
        registerIdentifier("identifier_nphi", "Cell ID identifier for Nphi", m_nphiId, "nphi");
        registerIdentifier("identifier_stereosign", "Cell ID identifier for StereoSign", m_stereosignId, "stereosign");

        m_systemIndex = decoder()->index(m_systemId);
        m_superlayerIndex = decoder()->index(m_superlayerId);
        m_layerIndex = decoder()->index(m_layerId);
        m_nphiIndex = decoder()->index(m_nphiId);
        m_stereosignIndex = decoder()->index(m_stereosignId);
    }


    CellID WireTracker_k4geo::setCellID(int System, int Superlayer, int Layer, int Nphi, int StereoSign) const {
        
        VolumeID systemId = static_cast<VolumeID>(System);
        VolumeID superlayerId = static_cast<VolumeID>(Superlayer);
        VolumeID layerId = static_cast<VolumeID>(Layer);
        VolumeID nphiId = static_cast<VolumeID>(Nphi);
        VolumeID stereosignId = static_cast<VolumeID>(StereoSign);
        
        VolumeID vID = 0;
        decoder()->set(vID, m_systemIndex, systemId);
        decoder()->set(vID, m_superlayerIndex, superlayerId);
        decoder()->set(vID, m_layerIndex, layerId);
        decoder()->set(vID, m_nphiIndex, nphiId);
        decoder()->set(vID, m_stereosignIndex, stereosignId);

        return vID;
  }

} // namespace DDSegmentation
} // namespace dd4hep
