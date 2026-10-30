// license:GPLv3+

#pragma once

#include "physics/coilmodel.h"

// Slingshot kicker plate: the plate (pivoted lever arm), the link to the coil's plunger, the return spring, the stops, the single-winding coil and the
// driver that makes the coil pulse. No VPX objects here: HitRubber (hitrubber.cpp) gives it the band displacement at the switches and takes the plate angle.
//
// The mechanics are those of the model flipper (coilmodel.h): the plate is "the flipper" (its inertia i_bat + i_crank about the pivot, angle phi from
// rest toward the coil stop), the coil pulls a plunger that drives the plate through a slider-crank link with pin slop. The circuit is one winding with the
// inductance of the plunger position, V across it while the driver pulses, a decay through its resistance otherwise.

struct SlingSettings
{
   FlipperSettings mech;        // plate, link, coil, spring and stops; the ini keys sling_<flipper key> (sling_r_crank_mm, sling_stop_deg, sling_r_power, ...)
   double pulseMs = 100.0;     // the driver's pulse once the switch has closed, ms (a guess: to be measured)
   double coilDelayMs = 0.0;   // from the switch closing to the coil's pulse, ms: each sling its own (user, 2026-10-30; two videos: about 20 and 33 - 39 ms to the plate seen moving)
   double minClosedMs = 2.0;   // a switch counts as closed for at least this long, ms (a guess)
   double switchGapMm = 1.5875; // band displacement at a switch that closes it, mm (1/16 in)
   double plateContactZeta = 0.1; // damping ratio of the plate face pressing on the band (a guess: steel on rubber)
   SlingSettings();
};

class SlingPlate
{
public:
   // State, SI: the plate angle from rest toward the coil stop, the plunger along the coil axis, the circuit current as a fraction of V/R
   struct State
   {
      double phi = 0.0, omega = 0.0;      // rad, rad/s
      double xp = 0.0, vp = 0.0;          // mm, mm/s
      double iHigh = 0.0;                 // fraction of V/R
      bool rigid = true;                  // one rigid body (coil off and the link taut) or two bodies
      bool operator==(const State &o) const { return phi == o.phi && omega == o.omega && xp == o.xp && vp == o.vp && iHigh == o.iHigh && rigid == o.rigid; }
   };

   // (Re)derive from the settings; the state is kept (a hot reload) unless `reset`.
   void Configure(const SlingSettings &s, bool reset);
   // One substep of h seconds. coilOn: the driver's pulse. extTorqueNm: torque on the plate from outside (the band), in the direction of increasing phi.
   void Step(double h, bool coilOn, double extTorqueNm);

   double PhiDeg() const { return st.phi / FlipperModel::DEG; }
   double CurrentA() const { return st.iHigh * set.mech.vSupply / set.mech.rP; }
   bool HitStop() const { return hitStop; }     // the plunger reached the coil stopper during the last Step
   double Alpha() const { return alpha; }
   State st;
   SlingSettings set;
   FlipperDerived der;

private:
   bool hitStop = false;
   double alpha = 0.0;
};

// The coil driver: a switch closing triggers one pulse of fixed length, whatever the switch does afterwards; the next pulse needs the pulse to be over and the
// switch to have reopened.
class SlingDriver
{
public:
   void SetPulse(double pulseS) { pulse = pulseS; }
   void SetDelay(double delayS) { delay = delayS; }
   // once per cycle, with whether any switch is closed
   void Update(bool switchClosed)
   {
      if (remaining <= 0.0 && pending <= 0.0 && !switchClosed)
         armed = true;
      if (armed && switchClosed && remaining <= 0.0 && pending <= 0.0)
      {
         armed = false;
         if (delay > 0.0)   // on its way: it fires whatever the switch does in the meantime
            pending = delay;
         else
         {
            remaining = pulse;
            fired++;
         }
      }
   }
   // script trigger: a pulse starts now unless one is running (the switches do not matter)
   bool Trigger()
   {
      if (remaining > 0.0)
         return false;
      pending = 0.0;
      remaining = pulse;
      armed = false;
      fired++;
      return true;
   }
   void Cancel()   // disabled: a running pulse stops, and nothing fires until the switches have been open again
   {
      remaining = 0.0;
      pending = 0.0;
      armed = false;
   }
   bool Coil() const { return remaining > 0.0; }
   void Tick(double h)
   {
      if (remaining > 0.0)
         remaining = remaining > h ? remaining - h : 0.0;
      if (pending > 0.0)
      {
         pending -= h;
         if (pending <= 0.0)
         {
            pending = 0.0;
            remaining = pulse;
            fired++;
         }
      }
   }
   double remaining = 0.0;   // s of pulse left
   double pending = 0.0;     // s until a pulse on its way starts
   int fired = 0;            // pulses started
   bool armed = true;        // the switches have been open since the last pulse (state: saved with a band part's snapshot)
private:
   double pulse = 0.03;
   double delay = 0.0;
};
