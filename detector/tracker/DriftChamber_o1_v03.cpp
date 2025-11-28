//----------------------------------
//         DCH detector v3
//----------------------------------

/*!
 *  \brief     Detector constructor of DCH v3
 *  \details   This code creates full geometry of IDEA DCH subdetector
 *  \author    Andreas Loeschcke Centeno      andreas.loeschcke.centeno@cern.ch
 *  \author    Brieuc Francois                  brieuc.francois@cern.ch
 *  \version   3
 *  \date      2025
 *  \pre       DD4hep compiled with Geant4+Qt
 */

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Shapes.h"
#include "XML/Utilities.h"

#include "DDRec/DCH_info.h"

namespace DCH_v3 {

using DCH_length_t = dd4hep::rec::DCH_info_struct::DCH_length_t;
using DCH_angle_t = dd4hep::rec::DCH_info_struct::DCH_angle_t;
using DCH_layer = dd4hep::rec::DCH_info_struct::DCH_layer;


/// Function to build DCH
static dd4hep::Ref_t create_DCH_o1_v03(dd4hep::Detector& desc, dd4hep::xml::Handle_t handle,
                                       dd4hep::SensitiveDetector sens) {
  dd4hep::xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  dd4hep::DetElement det(detName, detID);
  dd4hep::xml::setDetectorTypeFlag(detElem, det);
  sens.setType("tracker");

  // initialize empty DCH_info object
  // data extension mechanism requires it to be a raw pointer
  dd4hep::rec::DCH_info* DCH_i = new dd4hep::rec::DCH_info();
  // fill DCH_i with information from the XML file
  {
    // DCH outer geometry dimensions
    DCH_i->Set_rin(desc.constantAsDouble("DCH_gas_inner_cyl_R"));
    DCH_i->Set_rout(desc.constantAsDouble("DCH_gas_outer_cyl_R"));
    DCH_i->Set_lhalf(desc.constantAsDouble("DCH_gas_Lhalf"));

    // guard wires position, fix position
    DCH_i->Set_guard_rin_at_z0(desc.constantAsDouble("DCH_guard_inner_r_at_z0"));
    DCH_i->Set_guard_rout_at_zL2(desc.constantAsDouble("DCH_guard_outer_r_at_zL2"));

    DCH_angle_t dch_alpha = desc.constantAsDouble("DCH_alpha");
    DCH_i->Set_twist_angle(2 * dch_alpha);

    DCH_i->Set_nsuperlayers(desc.constantAsLong("DCH_nsuperlayers"));
    DCH_i->Set_nlayersPerSuperlayer(desc.constantAsLong("DCH_nlayersPerSuperlayer"));

    DCH_i->Set_ncell0(desc.constantAsLong("DCH_ncell"));
    DCH_i->Set_ncell_increment(desc.constantAsLong("DCH_ncell_increment"));
    DCH_i->Set_ncell_per_sector(desc.constantAsLong("DCH_ncell_per_sector"));

    DCH_i->Set_first_width(desc.constantAsDouble("DCH_first_width"));
    DCH_i->Set_first_sense_r(desc.constantAsDouble("DCH_first_sense_r"));

    bool buildLayers = detElem.attr<bool>(_Unicode(buildLayers));
    if (buildLayers) {
      DCH_i->BuildLayerDatabase();
      // safety check just in case something went wrong...
      if (DCH_i->IsDatabaseEmpty())
        throw std::runtime_error("Empty database");
    }

    bool printExcelTable = detElem.attr<bool>(_Unicode(printExcelTable));
    if (printExcelTable)
      DCH_i->Show_DCH_info_database(std::cout);
  }

  bool debugGeometry = detElem.hasChild(_Unicode(debugGeometry));
  bool useG4TT = detElem.hasChild(_Unicode(useG4TT));
  auto gasElem = detElem.child("gas");
  auto gasvolMat = desc.material(gasElem.attr<std::string>(_Unicode(material)));
  auto gasvolVis = desc.visAttributes(gasElem.attr<std::string>(_Unicode(vis)));

  auto vesselElem = detElem.child("vessel");
  auto vesselSkinVis = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(visSkin)));
  auto vesselBulkVis = desc.visAttributes(vesselElem.attr<std::string>(_Unicode(visBulk)));
  auto vessel_mainMaterial = desc.material(vesselElem.attr<std::string>(_Unicode(mainMaterial)));
  auto vessel_fillmaterial_outerR = desc.material(vesselElem.attr<std::string>(_Unicode(fillmaterial_outerR)));
  auto vessel_fillmaterial_endcap = desc.material(vesselElem.attr<std::string>(_Unicode(fillmaterial_endcap)));
  DCH_length_t vessel_fillmaterial_fraction_outerR = vesselElem.attr<double>(_Unicode(fillmaterial_fraction_outerR));
  DCH_length_t vessel_fillmaterial_fraction_endcap = vesselElem.attr<double>(_Unicode(fillmaterial_fraction_endcap));

  if (0 > vessel_fillmaterial_fraction_outerR || 1 < vessel_fillmaterial_fraction_outerR)
    throw std::runtime_error("vessel_fillmaterial_fraction_outerR must be between 0 and 1");
  if (0 > vessel_fillmaterial_fraction_endcap || 1 < vessel_fillmaterial_fraction_endcap)
    throw std::runtime_error("vessel_fillmaterial_fraction_z must be between 0 and 1");

  auto wiresElem = detElem.child("wires");
  auto wiresVis = desc.visAttributes(wiresElem.attr<std::string>(_Unicode(vis)));
  bool buildSenseWires = wiresElem.attr<bool>(_Unicode(buildSenseWires));
  bool buildFieldWires = wiresElem.attr<bool>(_Unicode(buildFieldWires));

  DCH_length_t dch_SWire_thickness = wiresElem.attr<double>(_Unicode(SWire_thickness));
  DCH_length_t dch_FSideWire_thickness = wiresElem.attr<double>(_Unicode(FSideWire_thickness));
  DCH_length_t dch_FCentralWire_thickness = wiresElem.attr<double>(_Unicode(FCentralWire_thickness));

  auto dch_SWire_material = desc.material(wiresElem.attr<std::string>(_Unicode(SWire_material)));
  auto dch_FSideWire_material = desc.material(wiresElem.attr<std::string>(_Unicode(FSideWire_material)));
  auto dch_FCentralWire_material = desc.material(wiresElem.attr<std::string>(_Unicode(FCentralWire_material)));

  /* Geometry tree:
   * Gas (tube) -> Layer_1 (hyp) -> cell_1 (twisted tube)
   *                             -> cell_... (twisted tube)
   *            -> Layer_... (hyp) -> cell_1 (twisted tube)
   *                               -> cell_... (twisted tube)
   *            -> Inner radius vessel wall
   *            -> Outer radius vessel wall -> fill made of foam
   *            -> Endcap disk  vessel wall -> fill made of kapton
   *
   * Layers represent a segmentation in radius
   * Sectors represent a segmentation in phi
   * Each cell corresponds to a Detector Element
   * Vessel wall has to be defined as 3 volumes to account for independent thickness and materials
   */

  DCH_length_t safety_r_interspace = 1 * dd4hep::nm;
  DCH_length_t safety_z_interspace = 1 * dd4hep::nm;
  DCH_length_t safety_phi_interspace = 1e-6 * dd4hep::rad;

  DCH_length_t vessel_thickness_innerR = desc.constantAsDouble("DCH_vessel_thickness_innerR");
  DCH_length_t vessel_thickness_outerR = desc.constantAsDouble("DCH_vessel_thickness_outerR");
  DCH_length_t vessel_endcapdisk_zmin = desc.constantAsDouble("DCH_vessel_disk_zmin");
  DCH_length_t vessel_endcapdisk_zmax = desc.constantAsDouble("DCH_vessel_disk_zmax");

  // if( 0 > vessel_thickness_z )
  // throw std::runtime_error("vessel_thickness_z must be positive");
  if (0 > vessel_thickness_innerR)
    throw std::runtime_error("vessel_thickness_innerR must be positive");
  if (0 > vessel_thickness_outerR)
    throw std::runtime_error("vessel_thickness_outerR must be positive");

  // // // // // // // // // // // // // // //
  // // // // // MAIN VOLUME // // // // // //
  // // // // // // // // // // // // // // //
  dd4hep::Tube gas_s(DCH_i->rin - vessel_thickness_innerR, DCH_i->rout + vessel_thickness_outerR,
                     vessel_endcapdisk_zmax);
  dd4hep::Volume gas_v(detName + "_gas", gas_s, gasvolMat);
  gas_v.setVisAttributes(gasvolVis);
  gas_v.setRegion(desc, detElem.regionStr());
  gas_v.setLimitSet(desc, detElem.limitsStr());

  DCH_length_t vessel_innerR_start = DCH_i->rin - vessel_thickness_innerR + safety_r_interspace;
  DCH_length_t vessel_innerR_end = DCH_i->rin;
  DCH_length_t vessel_outerR_start = DCH_i->rout;
  DCH_length_t vessel_outerR_end = DCH_i->rout + vessel_thickness_outerR - safety_r_interspace;
  DCH_length_t vessel_R_zhalf = vessel_endcapdisk_zmax - safety_z_interspace;
  // // // // // // // // // // // // // // //
  // // // // //  INNER R WALL  // // // // //
  // // // // // // // // // // // // // // //
  dd4hep::Tube vessel_innerR_s(vessel_innerR_start, vessel_innerR_end, vessel_R_zhalf);
  dd4hep::Volume vessel_innerR_v(detName + "_vessel_innerR", vessel_innerR_s, vessel_mainMaterial);
  vessel_innerR_v.setVisAttributes(vesselSkinVis);
  gas_v.placeVolume(vessel_innerR_v);
  // // // // // // // // // // // // // // //
  // // // // //  OUTER R WALL  // // // // //
  // // // // // // // // // // // // // // //
  dd4hep::Tube vessel_outerR_s(vessel_outerR_start, vessel_outerR_end, vessel_R_zhalf);
  dd4hep::Volume vessel_outerR_v(detName + "_vessel_outerR", vessel_outerR_s, vessel_mainMaterial);
  vessel_outerR_v.setVisAttributes(vesselSkinVis);

  // if thickness fraction of bulk material is defined, build the bulk material
  if (0 < vessel_fillmaterial_fraction_outerR) {
    double f = vessel_fillmaterial_fraction_outerR;
    DCH_length_t fillmaterial_thickness = f * (vessel_outerR_end - vessel_outerR_start);
    DCH_length_t rstart = vessel_outerR_start + 0.5 * (1 - f) * fillmaterial_thickness;
    DCH_length_t rend = vessel_outerR_end - 0.5 * (1 - f) * fillmaterial_thickness;
    dd4hep::Tube vessel_fillmat_outerR_s(rstart, rend, vessel_R_zhalf - safety_z_interspace);
    dd4hep::Volume vessel_fillmat_outerR_v(detName + "_vessel_fillmat_outerR", vessel_fillmat_outerR_s,
                                           vessel_fillmaterial_outerR);
    vessel_fillmat_outerR_v.setVisAttributes(vesselBulkVis);
    vessel_outerR_v.placeVolume(vessel_fillmat_outerR_v);
  }
  gas_v.placeVolume(vessel_outerR_v);

  // // // // // // // // // // // // // // //
  // // // // //  ENDCAP WALL   // // // // //
  // // // // // // // // // // // // // // //
  DCH_length_t vessel_endcap_thickness = vessel_endcapdisk_zmax - vessel_endcapdisk_zmin - 2 * safety_z_interspace;
  DCH_length_t vessel_endcap_zpos = 0.5 * (vessel_endcapdisk_zmax + vessel_endcapdisk_zmin);
  DCH_length_t vessel_endcap_rstart = vessel_innerR_end + safety_r_interspace;
  DCH_length_t vessel_endcap_rend = vessel_outerR_start - safety_r_interspace;
  dd4hep::Tube vessel_endcap_s(vessel_endcap_rstart, vessel_endcap_rend, 0.5 * vessel_endcap_thickness);
  dd4hep::Volume vessel_endcap_v(detName + "_vessel_endcap", vessel_endcap_s, vessel_mainMaterial);
  vessel_endcap_v.setVisAttributes(vesselSkinVis);

  // if thickness fraction of bulk material is defined, build the bulk material
  if (0 < vessel_fillmaterial_fraction_endcap) {
    double f = vessel_fillmaterial_fraction_endcap;
    DCH_length_t fillmaterial_thickness = f * vessel_endcap_thickness;
    dd4hep::Tube vessel_fillmat_endcap_s(vessel_endcap_rstart + safety_r_interspace,
                                         vessel_endcap_rend - safety_r_interspace, 0.5 * fillmaterial_thickness);
    dd4hep::Volume vessel_fillmat_endcap_v(detName + "_vessel_fillmat_endcap", vessel_fillmat_endcap_s,
                                           vessel_fillmaterial_endcap);
    vessel_fillmat_endcap_v.setVisAttributes(vesselBulkVis);
    vessel_endcap_v.placeVolume(vessel_fillmat_endcap_v);
  }
  // place endcap wall at +/- z
  gas_v.placeVolume(vessel_endcap_v, dd4hep::Position(0, 0, vessel_endcap_zpos));
  gas_v.placeVolume(vessel_endcap_v, dd4hep::Position(0, 0, -vessel_endcap_zpos));

  // // // // // // // // // // // // // // //
  // // // // //  DCH layers    // // // // //
  // // // // // // // // // // // // // // //
  for (const auto& [ilayer, l] : DCH_i->database) {

    // // // // // // // // // // // // // // // // // // // // /
    // // // // // INITIALIZATION OF THE LAYER_TRIPLET // // // // // //
    // // // // // // // // // // // // // // // // // // // //
    // Layer triplet is a hyperboloid that contains three hyperboloid layers:
    // - thin hyperboloid for inner field wires (thickness only the wire thickness)
    // - thick hyperboloid for the sense wires (thickness of cell drift region)
    // - thin hyperboloid for outer field wires (thickness only the wire thickness)
    
    // Hyperboloid parameters:
    /// inner radius at z=0
    DCH_length_t triplet_inner_r = l.radius_fdw_z0 + safety_r_interspace;
    /// inner stereoangle, calculated from rin(z=0)
    DCH_angle_t  triplet_inner_stereo = DCH_i->stereoangle_z0(triplet_inner_r);
    /// outer radius at z=0
    DCH_length_t triplet_outer_r = l.radius_fuw_z0 - safety_r_interspace;
    /// outer stereoangle, calculated from rout(z=0)
    DCH_angle_t  triplet_outer_stereo = DCH_i->stereoangle_z0(triplet_outer_r);
    /// half-length
    DCH_length_t triplet_half_length = DCH_i->Lhalf + safety_z_interspace;

    dd4hep::Hyperboloid triplet_solid(triplet_inner_r, triplet_inner_stereo, triplet_outer_r, triplet_outer_stereo, triplet_half_length);

    std::string triplet_name = detName + "_triplet" + std::to_string(ilayer);
    dd4hep::Volume triplet_volume(triplet_name, triplet_solid, gasvolMat);
    triplet_volume.setVisAttributes(desc.visAttributes(Form("dch_layer_vis%d", ilayer % 2)));
    if (ilayer < 156){
      auto triplet_placed = gas_v.placeVolume(triplet_volume, ilayer);

      int ilayerWithinSuperlayer = (ilayer - 1) % DCH_i->nlayersPerSuperlayer;
      int nsuperlayer_minus_1 = DCH_i->Get_nsuperlayer_minus_1(ilayer);
      triplet_placed.addPhysVolID("superlayer", nsuperlayer_minus_1);
    }

    // ilayer is a counter that runs from 1 to 112 (nsuperlayers * nlayersPerSuperlayer)
    // it seems more convenient to store the layer number within the superlayer
    // ilayerWithinSuperlayer runs from 0 to 7 (nlayersPerSuperlayer-1)

    // dd4hep::DetElement layer_DE(det, layer_name + "DE", ilayer);
    // layer_DE.setPlacement(layer_pv);

    // // // // // // // // // // // // // // // // // // // //
    // // // // // SEGMENTATION OF THE TRIPLET  // // // // // //
    // // // // // INTO LAYERS (F + S + F WIRES) // // // // // //
    // // // // // // // // // // // // // // // // // // // //

    DCH_length_t sense_layer_thickness_z0 = triplet_outer_r - triplet_inner_r - 2*(dch_FSideWire_thickness + 2*safety_r_interspace);
    DCH_length_t field_layer_thickness_z0 = dch_FSideWire_thickness + 2*safety_r_interspace;


    // inner radius of the inner field layer
    DCH_length_t inner_field_layer_inner_r = triplet_inner_r + safety_r_interspace;

    // outer radius of the inner field layer
    DCH_length_t inner_field_layer_outer_r = inner_field_layer_inner_r + field_layer_thickness_z0;

    dd4hep::Hyperboloid inner_field_layer_solid(
        inner_field_layer_inner_r,
        DCH_i->stereoangle_z0(inner_field_layer_inner_r),
        inner_field_layer_outer_r,
        DCH_i->stereoangle_z0(inner_field_layer_outer_r),
        triplet_half_length);

    dd4hep::Volume inner_field_layer_volume(
        triplet_name + "_inner_field_layer",
        inner_field_layer_solid,
        gasvolMat);

    inner_field_layer_volume.setVisAttributes(desc.visAttributes(Form("dch_innerouter_vis%d", ilayer % 2)));


    ////// OUTER FIELD LAYER //////
    DCH_length_t outer_field_layer_outer_r = triplet_outer_r - safety_r_interspace;
    DCH_length_t outer_field_layer_inner_r = outer_field_layer_outer_r - field_layer_thickness_z0;

    dd4hep::Hyperboloid outer_field_layer_solid(
        outer_field_layer_inner_r,
        DCH_i->stereoangle_z0(outer_field_layer_inner_r),
        outer_field_layer_outer_r,
        DCH_i->stereoangle_z0(outer_field_layer_outer_r),
        triplet_half_length);

    dd4hep::Volume outer_field_layer_volume(
        triplet_name + "_outer_field_layer",
        outer_field_layer_solid,
        gasvolMat);
    outer_field_layer_volume.setVisAttributes(desc.visAttributes(Form("dch_innerouter_vis%d", ilayer % 2)));


    ///// SENSE LAYER //////

    DCH_length_t sense_layer_inner_r = inner_field_layer_outer_r + safety_r_interspace;
    DCH_length_t sense_layer_outer_r = outer_field_layer_inner_r - safety_r_interspace;

    dd4hep::Hyperboloid sense_layer_solid(
        sense_layer_inner_r,
        DCH_i->stereoangle_z0(sense_layer_inner_r),
        sense_layer_outer_r,
        DCH_i->stereoangle_z0(sense_layer_outer_r),
        triplet_half_length);

    dd4hep::Volume sense_layer_volume(
        triplet_name + "_sense_layer",
        sense_layer_solid,
        gasvolMat);

    sense_layer_volume.setVisAttributes(desc.visAttributes(Form("dch_sense_vis%d", ilayer % 2)));
    

    auto inner_field_layer_placed = triplet_volume.placeVolume(inner_field_layer_volume);
    auto sense_layer_placed = triplet_volume.placeVolume(sense_layer_volume);
    auto outer_field_layer_placed = triplet_volume.placeVolume(outer_field_layer_volume);


    /// WIRES ///
    // Radius and length of the wires
    DCH_length_t sense_wire_radius = dch_SWire_thickness / 2.0;
    DCH_length_t sense_wire_placement_radius = (sense_layer_inner_r + sense_layer_outer_r) / 2.0;
    // DCH_length_t sense_wire_length = DCH_i->Lhalf ;// / cos(DCH_i->stereoangle_z0(sense_wire_placement_radius)) - tan(DCH_i->stereoangle_z0(sense_wire_placement_radius)) * sense_wire_radius - safety_z_interspace;
    DCH_length_t sense_wire_length = 0.5 * DCH_i->WireLength(ilayer, sense_wire_placement_radius) -
                              sense_wire_radius * cos(DCH_i->stereoangle_z0(sense_wire_placement_radius)) - safety_z_interspace;

    if (sense_wire_placement_radius - l.radius_sw_z0 > 1e-6*dd4hep::mm) {
      std::cout << "Warning: sense wire placement radius " << sense_wire_placement_radius/dd4hep::mm
                << " mm is larger than database value " << l.radius_sw_z0/dd4hep::mm << " mm for layer " << ilayer << " \n by " << sense_wire_placement_radius - l.radius_sw_z0/dd4hep::mm << " mm"
                << ". Check calculation." << std::endl;
    }

    if (sense_wire_placement_radius - (0.5 * (l.radius_fdw_z0 + l.radius_fuw_z0)) > 1e-6*dd4hep::mm) {
      std::cout << "Warning: sense wire placement radius " << sense_wire_placement_radius/dd4hep::mm
                << " mm is larger than average of field wire radii "
                << 0.5 * (l.radius_fdw_z0 + l.radius_fuw_z0) / dd4hep::mm << " mm for layer " << ilayer << " \n by " << sense_wire_placement_radius - (0.5 * (l.radius_fdw_z0 + l.radius_fuw_z0))/dd4hep::mm << " mm"
                << ". Check calculation." << std::endl;
    }


    dd4hep::Tube    sense_wire_solid(0., sense_wire_radius, sense_wire_length);
    dd4hep::Volume  sense_wire_volume(triplet_name + "_sense_wire", sense_wire_solid, dch_SWire_material);

    sense_wire_volume.setVisAttributes(wiresVis);

    
    DCH_length_t side_field_wire_radius = dch_FSideWire_thickness / 2.0;

    // DCH_length_t inner_field_wire_length = DCH_i->Lhalf/cos(DCH_i->stereoangle_z0(l.radius_fdw_z0)) - tan(DCH_i->stereoangle_z0(l.radius_fdw_z0))*side_field_wire_radius - safety_z_interspace;
    DCH_length_t inner_field_wire_placement_radius = (inner_field_layer_inner_r + inner_field_layer_outer_r) / 2.0;
    DCH_length_t inner_field_wire_length = 0.5 * DCH_i->WireLength(ilayer, inner_field_wire_placement_radius) - side_field_wire_radius * cos(DCH_i->stereoangle_z0(inner_field_wire_placement_radius)) - safety_z_interspace;
    // DCH_length_t inner_field_wire_length = 0.5 * DCH_i->WireLength(ilayer, inner_field_wire_placement_radius) -
    //                                   side_field_wire_radius * cos(DCH_i->stereoangle_z0(inner_field_wire_placement_radius)) - safety_z_interspace;
    dd4hep::Tube    inner_field_wire_solid(0., side_field_wire_radius, inner_field_wire_length);
    dd4hep::Volume  inner_field_wire_volume(triplet_name + "_inner_field_wire", inner_field_wire_solid, dch_FSideWire_material);
    inner_field_wire_volume.setVisAttributes(wiresVis);

    
    DCH_length_t outer_field_wire_placement_radius = (outer_field_layer_inner_r + outer_field_layer_outer_r) / 2.0; 
    DCH_length_t outer_field_wire_length = DCH_i->Lhalf/cos(DCH_i->stereoangle_z0(outer_field_wire_placement_radius)) - tan(DCH_i->stereoangle_z0(outer_field_wire_placement_radius))*side_field_wire_radius - safety_z_interspace;
    // DCH_length_t outer_field_wire_length = 0.5 * DCH_i->WireLength(ilayer, outer_field_wire_placement_radius) - side_field_wire_radius * cos(DCH_i->stereoangle_z0(outer_field_wire_placement_radius)) - safety_z_interspace;
    // DCH_length_t outer_field_wire_length = 0.5 * DCH_i->WireLength(ilayer, outer_field_wire_placement_radius) -
    //                                   side_field_wire_radius * cos(DCH_i->stereoangle_z0(outer_field_wire_placement_radius)) - safety_z_interspace;
    dd4hep::Tube    outer_field_wire_solid(0., side_field_wire_radius, outer_field_wire_length);
    dd4hep::Volume  outer_field_wire_volume(triplet_name + "_outer_field_wire", outer_field_wire_solid, dch_FSideWire_material);
    outer_field_wire_volume.setVisAttributes(wiresVis);


    DCH_length_t central_field_wire_radius = dch_FCentralWire_thickness / 2.0;
    // DCH_length_t central_field_wire_length = sense_wire_length;
    
    DCH_length_t central_field_wire_length = 0.5 * DCH_i->WireLength(ilayer, sense_wire_placement_radius) - central_field_wire_radius * cos(DCH_i->stereoangle_z0(sense_wire_placement_radius)) - safety_z_interspace;
    // DCH_length_t central_field_wire_length = 0.5 * DCH_i->WireLength(ilayer, sense_wire_placement_radius) -
    //                                   central_field_wire_radius * cos(DCH_i->stereoangle_z0(sense_wire_placement_radius)) - safety_z_interspace;
    dd4hep::Tube    field_central_wire_solid(0., central_field_wire_radius, central_field_wire_length);
    dd4hep::Volume  field_central_wire_volume(triplet_name + "_central_field_wire", field_central_wire_solid, dch_FCentralWire_material);
    field_central_wire_volume.setVisAttributes(wiresVis);



    // Loop over phi to place the wires in all the layers
    DCH_angle_t phi_step = (TMath::TwoPi() / l.nwires) * dd4hep::rad;

    // width of one cell, used for staggering the cells in phi
    DCH_angle_t cell_phi_width = 2*phi_step - safety_phi_interspace;

    dd4hep::RotationX sense_stereo_rot((-1.) * l.StereoSign() * DCH_i->stereoangle_z0(sense_wire_placement_radius));
    dd4hep::Transform3D sense_transform(sense_stereo_rot * dd4hep::Translation3D(sense_wire_placement_radius, 0., 0.));

    dd4hep::RotationX inner_field_stereo_rot((-1.) * l.StereoSign() * DCH_i->stereoangle_z0(inner_field_wire_placement_radius));
    dd4hep::Transform3D inner_field_transform(inner_field_stereo_rot * dd4hep::Translation3D(inner_field_wire_placement_radius, 0., 0.));

    dd4hep::RotationX outer_field_stereo_rot((-1.) * l.StereoSign() * DCH_i->stereoangle_z0(outer_field_wire_placement_radius));
    dd4hep::Transform3D outer_field_transform(outer_field_stereo_rot * dd4hep::Translation3D(outer_field_wire_placement_radius, 0., 0.));

    std::cout << "Rotation Angle Sense: " << (-1.) * l.StereoSign() * DCH_i->stereoangle_z0(sense_wire_placement_radius)/dd4hep::deg << " deg" << std::endl;
    std::cout << "Rotation Angle Inner Field: " << (-1.) * l.StereoSign() * DCH_i->stereoangle_z0(inner_field_wire_placement_radius)/dd4hep::deg << " deg" << std::endl;
    std::cout << "Rotation Angle Outer Field: " << (-1.) * l.StereoSign() * DCH_i->stereoangle_z0(outer_field_wire_placement_radius)/dd4hep::deg << " deg" << std::endl;

    for (int nphi=0; nphi<1; ++nphi) 
    {
      DCH_angle_t wire_phi_angle = phi_step * nphi + 0.25 * cell_phi_width * (ilayer % 2);

      dd4hep::RotationZ wire_phi_rot(wire_phi_angle);

      // Inner field wire
      dd4hep::Transform3D inner_field_wire_final_transform(wire_phi_rot * inner_field_transform);
      auto inner_field_placed = inner_field_layer_volume.placeVolume(inner_field_wire_volume, inner_field_wire_final_transform);

      // Outer field wire
      dd4hep::Transform3D outer_field_wire_final_transform(wire_phi_rot * outer_field_transform);
      auto outer_field_placed = outer_field_layer_volume.placeVolume(outer_field_wire_volume, outer_field_wire_final_transform);

      // Alternating sense wire and central field wire
      dd4hep::Transform3D sense_or_central_final_transform(wire_phi_rot * sense_transform);
      dd4hep::Volume* wire_to_be_placed = (nphi % 2 == 0) ? &sense_wire_volume : &field_central_wire_volume;
      auto sense_or_central_placed = sense_layer_volume.placeVolume(*wire_to_be_placed, sense_or_central_final_transform);
    }

    
  }

  // Place our mother volume in the world
  dd4hep::Volume wVol = desc.pickMotherVolume(det);
  dd4hep::PlacedVolume vessel_pv = wVol.placeVolume(gas_v);
  // Associate the wall to the detector element.
  det.setPlacement(vessel_pv);
  // Assign the system ID to our mother volume
  vessel_pv.addPhysVolID("system", detID);
  det.addExtension<dd4hep::rec::DCH_info>(DCH_i);

  return det;
}

}; // namespace DCH_v3

DECLARE_DETELEMENT(DriftChamber_o1_v03_T, DCH_v3::create_DCH_o1_v03)
