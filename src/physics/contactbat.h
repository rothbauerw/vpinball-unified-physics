#pragma once
// A bat as a body of the contact system (docs/UnifiedContactDesign.txt, "Bat"): one rotation about the pivot, driven by the plunger and its coils through the pin-slop
// link, the return spring and the rest stop. Everything here is SI and knows nothing about flippers as table parts, balls, or the engine: the owner feeds it the
// solenoid state and the external torque of the substep (the balls' contact forces) and reads the angle back.
//
// The mechanics (settings, link kinematics, force laws) are in coilmodel.h; this class owns the state and the substep (Step), the settings file, and the pivot nudge.
// The model is a port of rothbauerw's real-data flipper model: see the comment in coilmodel.h and Step.

#include <filesystem>
#include <string>

#include "physics/coilmodel.h"
#include "physics/rubbercontact.h"

// Settings helpers shared with the sling plate
namespace FlipperModel
{
   RubberContact::RubberParams MakeRubberParams(const FlipperSettings &s);   // the bat rubber's contact parameters
}

class ContactBat
{
public:
   // ---- settings ----
   FlipperSettings set;
   FlipperDerived der;
   std::filesystem::path iniPath;
   std::filesystem::file_time_type iniTime{};
   int iniCheckCounter = 0;
   // (Re)builds the settings from the ini file (the table's vpx_physics.ini; none = the defaults); the state is left untouched. In contactbatini.cpp (VPX side).
   void LoadSettings(const std::string &who, const std::filesystem::path &tableFile, bool force);

   // ---- state ----
   bool rigid = false;            // rigid one-body mode (coil off and link taut); two-body otherwise
   bool wasPreEos = true;         // series wiring: the EOS switch was closed on the previous step
   bool coilOffPrev = true;       // the coil was off on the previous step
   bool init = false;
   double xp = 0.0;               // plunger position, mm along the coil axis
   double vp = 0.0;               // plunger velocity, mm/s
   double phi = 0.0;              // flipper angle from rest toward the end stop, rad
   double omega = 0.0;            // flipper angular velocity, rad/s
   double iHigh = 0.0;            // power coil current (fraction of steady state)
   double iLow = 0.0;             // hold coil current (fraction of steady state)
   double lastAlpha = 0.0;        // last flipper angular acceleration, rad/s^2
   bool plungerHitStop = false;   // the plunger reached the coil stopper during the last Step
   double coilWorkJ = 0.0;        // work the coil has done on the plunger since this was last read (energy books, ModelStandards 2h); a report, not a force
   double hitL = 0.0;             // angular momentum at the last plunger hit on the coil stop, kg m^2/s
   // Dormancy: a bat at rest (coil unchanged, nothing has touched it, the state did not move over a whole physics step) would only repeat the same substeps every
   // millisecond. idle is set after such a step; the substeps are skipped until any part of the state differs from the snapshot taken then.
   bool idle = false;
   bool idleSol = false, idleRigid = false;
   double idlePhi = 0.0, idleOmega = 0.0, idleXp = 0.0, idleVp = 0.0, idleIHigh = 0.0, idleILow = 0.0;

   // The bat at angle phi0 (rad from rest) with the coil in state sol, at rest, plunger slaved
   void Reset(double phi0, bool sol);
   // One SI substep of h seconds. extTorqueNm: external torque on the bat in the closing direction (a ball's contact), 0 = none. solenoid: the coil is energised.
   void Step(double h, double extTorqueNm, bool solenoid);
   // Moment of inertia the balls feel (kg m^2): the bat alone while the link is loose, bat + the plunger's reflected mass while the link is loaded or rigid
   double InertiaSI() const;
   // Snapshot / compare for dormancy
   bool IdleMatches(bool sol) const;
   void IdleCapture(bool sol);

   // ---- pivot nudge ----
   // The striking bat's coil stop hit starts a one-shot pulse on the target bat's pivot: x displacement away from the striker, half-cosine rise then half-cosine return.
   bool nzActive = false;
   double nzT = 0.0, nzA = 0.0, nzRiseS = 0.001, nzFallS = 0.008, nzBaseX = 0.0;
   float nzVelU = 0.f;            // pivot x velocity, VPX units per VPT
   void NudgeEval(double t, double &dispM, double &velMs) const;   // pivot displacement (m) and velocity (m/s) t seconds into the pulse
   void NudgeStart(double baseXU, double dispM, double riseS, double fallS);   // baseXU: the pivot's x now (VPX units)
   void NudgeStep(double dtS, float &pivotXU);                      // advances the pulse; pivotXU is set to where the pivot is at the end of dtS
};
