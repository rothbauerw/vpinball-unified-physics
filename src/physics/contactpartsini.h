// license:GPLv3+
//
// Settings of the rebuilt gates and spinners (docs/GateSpinnerDesign.txt section 10) from the unified vpx_physics.ini: [Gates] / [Spinners] for all, then
// [Gate:<name>] / [Spinner:<name>] for one part. Values not stated by the user are UNVERIFIED defaults (see the design document and NewPhysicsReference).

#pragma once

#include <filesystem>
#include <fstream>
#include <string>

#include "physics/physicsini.h"
#include "physics/contactgate.h"
#include "physics/contactspinner.h"
#include "physics/contactbumper.h"

namespace PartsIni
{
   // Calls f(key, value, lowerValue) for every key of the sections `plain` (e.g. "gates") and then `kind:name` (e.g. "gate:gate1"), in file order within each pass
   template <class F> inline void ReadSections(const std::filesystem::path &file, const std::string &plain, const std::string &kindName, F f)
   {
      for (int pass = 0; pass < 2; pass++)
      {
         std::ifstream in(file);
         std::string line;
         bool on = false;
         const std::string want = PhysicsIni::Lower(pass == 0 ? plain : kindName);
         while (std::getline(in, line))
         {
            const size_t cmt = line.find_first_of("#;");
            if (cmt != std::string::npos)
               line.resize(cmt);
            line = PhysicsIni::Trim(line);
            if (line.empty())
               continue;
            if (line.front() == '[' && line.back() == ']')
            {
               on = PhysicsIni::Lower(PhysicsIni::Trim(line.substr(1, line.size() - 2))) == want;
               continue;
            }
            const size_t eq = line.find('=');
            if (!on || eq == std::string::npos)
               continue;
            const std::string key = PhysicsIni::Lower(PhysicsIni::Trim(line.substr(0, eq)));
            const std::string val = PhysicsIni::Trim(line.substr(eq + 1));
            f(key, atof(val.c_str()), PhysicsIni::Lower(val));
         }
      }
   }

   // Gate settings over `s` (the part's own values already filled in). Unknown keys are reported in the log.
   inline void LoadGate(GateModel::Settings &s, const std::filesystem::path &file, const std::string &name)
   {
      ReadSections(file, "gates", "gate:" + name, [&](const std::string &k, double v, const std::string &lv) {
         if (k == "pivot_height_mm") s.pivotHeightMm = v;
         else if (k == "hit_offset_mm") s.hitOffsetMm = v;
         else if (k == "depth_mm") s.depthMm = v;
         else if (k == "wire_d_mm") s.wireDMm = v;
         else if (k == "w_apex") s.wApex = v;
         else if (k == "bottom_width_frac") s.bottomFrac = v;
         else if (k == "plate_top_mm") s.plateTopMm = v;
         else if (k == "plate_bottom_mm") s.plateBottomMm = v;
         else if (k == "plate_len_mm") s.plateLenMm = v;
         else if (k == "long_top_mm") s.longTopMm = v;
         else if (k == "long_bottom_mm") s.longBottomMm = v;
         else if (k == "long_len_mm") s.longLenMm = v;
         else if (k == "plate_thickness_mm") s.plateThickMm = v;
         else if (k == "lip_mm") s.lipMm = v;
         else if (k == "lip_bend_deg") s.lipBendDeg = v;
         else if (k == "rest_angle_deg") s.restAngleDeg = v;
         else if (k == "open_stop_deg") s.openStopDeg = v;
         else if (k == "stop_arm_mm") s.stopArmMm = v;
         else if (k == "stop_e") s.stopE = v;
         else if (k == "bracket_edge_mm") s.bracketEdgeMm = v;
         else if (k == "mu_pivot") s.muPivot = v;
         else if (k == "density_kg_m3") s.density = v;
         else if (k == "switch") s.sw = lv == "none" ? GateModel::Switch::None : GateModel::Switch::Leaf;
         else if (k == "switch_lead_deg") s.switchLeadDeg = v;
         else if (k == "switch_hyst_deg") s.switchHystDeg = v;
         else if (k == "switch_force_n") s.switchForceN = v;
         else if (k == "hold_angle_deg") s.holdAngleDeg = v;
         else if (k == "hold_speed_dps") s.holdSpeedDps = v;
         else if (k == "steel_e_gpa") s.steelEGpa = v;
         else if (k == "release_gap_mm") s.releaseGapMm = v;
         else if (k == "hit_window_deg") s.hitWindowDeg = v;
         else if (k == "rubber") {}   // the opt-out key (HitRubber::HingeWanted)
         else PLOGW << "Gate " << name << ": unknown key '" << k << "' ignored";
      });
   }

   // Jet bumper settings over `p` (docs/JetBumperDesign.txt section 5): [Bumpers] for all, then [Bumper:<name>]. Lengths in mm, masses in g, spring rates in N/mm unless the key says otherwise.
   inline void LoadBumper(JetBumper::Params &p, const std::filesystem::path &file, const std::string &name)
   {
      using JetBumper::DriveMode;
      ReadSections(file, "bumpers", "bumper:" + name, [&](const std::string &k, double v, const std::string &lv) {
         if (k == "body_dia_mm") p.bodyDiaMm = v;
         else if (k == "body_height_mm") p.bodyTopMm = v;
         else if (k == "skirt_dia_mm") p.skirtOdMm = v;
         else if (k == "skirt_inner_dia_mm") p.skirtInnerDiaMm = v;
         else if (k == "skirt_t_inner_mm") p.skirtTInnerMm = v;
         else if (k == "skirt_t_edge_mm") p.skirtTEdgeMm = v;
         else if (k == "skirt_clearance_mm") p.skirtClearanceMm = v;
         else if (k == "lip_dia_mm") p.lipDiaMm = v;
         else if (k == "stem_length_mm") p.stemLenMm = v;
         else if (k == "stem_dia_mm") p.stemDiaMm = v;
         else if (k == "stem_tip_radius_mm") p.tipRadiusMm = v;
         else if (k == "spoon_cone_radius_mm") p.coneRadiusMm = v;
         else if (k == "spoon_cone_depth_mm") p.coneDepthMm = v;
         else if (k == "ring_od_mm") p.ringOdMm = v;
         else if (k == "ring_id_mm") p.ringIdMm = v;
         else if (k == "ring_thickness_mm") p.ringSheetMm = v;
         else if (k == "ring_rise_mm") p.ringRiseMm = v;
         else if (k == "ring_mid_rest_mm") p.ringMidRestMm = v;
         else if (k == "ring_stroke_mm") p.ringStrokeMm = v;
         else if (k == "skirt_mass_g") p.skirtMass = v * 1e-3;
         else if (k == "ring_assembly_mass_g") p.ringMass = v * 1e-3;
         else if (k == "spoon_mass_g") p.spoonMass = v * 1e-3;
         else if (k == "skirt_spring_n_per_mm") p.skirtKNPerM = v * 1e3;
         else if (k == "skirt_spring_preload_n") p.skirtPreloadN = v;
         else if (k == "skirt_spring_radius_mm") p.springRadiusMm = v;
         else if (k == "skirt_damping_zeta") p.skirtZeta = v;
         else if (k == "spoon_spring_n_per_mm") p.spoonKNPerM = v * 1e3;
         else if (k == "spoon_spring_preload_n") p.spoonPreloadN = v;
         else if (k == "spoon_damping_zeta") p.spoonZeta = v;
         else if (k == "spoon_travel_mm") p.spoonTravelMm = v;
         else if (k == "return_spring_n_per_mm") p.returnKNPerM = v * 1e3;
         else if (k == "return_spring_preload_n") p.returnPreloadN = v;
         else if (k == "ring_friction_ns_per_m") p.ringViscousNsPerM = v;
         else if (k == "ring_rest_stop_e") p.ringRestE = v;
         else if (k == "ring_coil_stop_e") p.ringCoilE = v;
         else if (k == "ring_rest_stop_e_mpa") p.ringRestEMpa = v;
         else if (k == "switch_close_mm") p.switchCloseMm = v;
         else if (k == "switch_hyst_mm") p.switchHystMm = v;
         else if (k == "coil_drive") p.mode = lv == "held" ? DriveMode::Held : lv == "external" ? DriveMode::External : DriveMode::Pulse;
         else if (k == "pulse_ms") p.pulseMs = v;
         else if (k == "coil") { if (lv != "ae-23-800") PLOGW << "Bumper " << name << ": coil '" << lv << "' is not in the table (only AE-23-800); set coil_r_ohm, coil_l_air_mh, coil_l_arm_mh, coil_bobbin_mm, coil_x0 instead"; }
         else if (k == "coil_v_supply") p.coil.vSupply = v;
         else if (k == "coil_r_ohm") p.coil.rOhm = v;
         else if (k == "coil_l_air_mh") p.coil.lAirMH = v;
         else if (k == "coil_l_arm_mh") p.coil.lArmMH = v;
         else if (k == "coil_bobbin_mm") p.coil.bobbinMm = v;
         else if (k == "coil_x0") p.coil.x0 = v;
         else if (k == "coil_plunger_closed_frac") p.plungerClosedFrac = v;
         else if (k == "plastic_e_mpa") p.plasticEMpa = v;
         else if (k == "plastic_restitution") p.plasticE = v;
         else if (k == "plastic_mu") p.plasticMu = v;
         else if (k == "steel_restitution") p.steelRingE = v;
         else if (k == "steel_mu") p.steelMu = v;
         else if (k == "physical") {}   // physical = 0 keeps the stock bumper (HitRubber::BumperWanted)
         else PLOGW << "Bumper " << name << ": unknown key '" << k << "' ignored";
      });
   }

   inline void LoadSpinner(SpinnerModel::Settings &s, const std::filesystem::path &file, const std::string &name)
   {
      ReadSections(file, "spinners", "spinner:" + name, [&](const std::string &k, double v, const std::string &lv) {
         if (k == "blade_height_mm") s.heightMm = v;
         else if (k == "blade_thickness_mm") s.thickMm = v;
         else if (k == "offset_mm") s.offsetMm = v;
         else if (k == "pivot_height_mm") s.pivotHeightMm = v;
         else if (k == "wire_d_mm") s.wireDMm = v;
         else if (k == "density_kg_m3") s.density = v;
         else if (k == "mu_pivot") s.muPivot = v;
         else if (k == "damping") s.damping = v;
         else if (k == "switch") s.sw = lv == "none" ? SpinnerModel::Switch::None : lv == "opto" ? SpinnerModel::Switch::Opto : lv == "micro" ? SpinnerModel::Switch::Micro : SpinnerModel::Switch::Leaf;
         else if (k == "switch_angle_deg") s.switchAngleDeg = v;
         else if (k == "opto_half_width_deg") s.optoHalfWidthDeg = v;
         else if (k == "cam_r_mm") s.camRMm = v;
         else if (k == "contact_half_deg") s.contactHalfDeg = v;
         else if (k == "closed_half_deg") s.closedHalfDeg = v;
         else if (k == "hyst_deg") s.hystDeg = v;
         else if (k == "leaf_force_n") s.leafForceN = v;
         else if (k == "micro_force_n") s.microForceN = v;
         else if (k == "micro_k_n_per_m") s.microKNPerM = v;
         else if (k == "switch_loss_ns_per_m") s.switchLossNsPerM = v;
         else if (k == "stop_arm_mm") s.stopArmMm = v;
         else if (k == "stop_e") s.stopE = v;
         else if (k == "bracket_edge_mm") s.bracketEdgeMm = v;
         else if (k == "steel_e_gpa") s.steelEGpa = v;
         else if (k == "release_gap_mm") s.releaseGapMm = v;
         else if (k == "hit_window_deg") s.hitWindowDeg = v;
         else if (k == "rubber") {}
         else PLOGW << "Spinner " << name << ": unknown key '" << k << "' ignored";
      });
   }
}
