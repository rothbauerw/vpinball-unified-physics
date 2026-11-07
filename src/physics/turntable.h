#pragma once
// Turntables: a motor-driven disc set into the playfield (or a lower / upper playfield). Configured in [Turntable:<name>] sections of vpx_physics.ini;
// the table script drives the motor with Turntable "<name>", volts (or TurntableDuty "<name>", duty for a PWM drive) and reads TurntableSpeed("<name>") (rpm) and
// TurntableAngle("<name>") (degrees).
//
// Physics (SI inside):
//  - disc: rotational inertia I (default a solid disc of mass_g), viscous bearing drag b from ramp_down_s (time from max_rpm to 5% with the motor off).
//  - motor: a DC motor with a linear torque-speed curve. Driven with V (signed, volts) it gives T = Ts V / Vrated - k w, k = Ts / w_run, w_run = max_rpm: the speed the disc settles at
//    with supply_volts and no load, so the speed is proportional to the voltage. Ts is chosen so the disc reaches 95% of its speed in ramp_up_s (or given as stall_torque_nm).
//    V = 0 is an unpowered motor: the disc coasts down through the bearing drag alone. A PWM drive is its average voltage (duty x supply_volts): the motor inductance and the
//    disc inertia filter the switching, so the mean is what the disc responds to.
//  - balls: the disc is a moving part of the unified contact model on the group step (DiscModel::Engine, contactdiscpart.h; docs/unified-physics/archive/UnifiedContactDesign.txt 12.1). The disc is a flat surface whose
//    points move at w x r; the friction of the contact law (Coulomb, mu of the textured surface) carries a ball along, and the equal and opposite force and torque of every ball
//    act on the disc in the same substep (a ball dropped on the rim slows it, a heavy load stalls it). Throw-off, slip and spin-up of the ball come from the contact law: nothing
//    about them is set. The motor and the disc are advanced by the contact system with the balls, and alone (ContactPartHost::PartAdvance) when no ball is near.
//  - hole: the playfield (plane or mesh) has a round hole of radius + gap_mm (1/16 in) at the disc; a ball over the gap rests on the two rims.
// Needs [Features] playfield_contact = 1 (the disc is a surface of the playfield's force model); without it the sections are ignored with a warning.
// The disc never changes shape, so its collision geometry is fixed for the session; the motor and surface keys are read when the table loads.
// Direction: the sign of the voltage's effect is the ini key direction (cw / ccw, as seen from above). Speed and angle are read in one fixed frame whatever the key:
// TurntableSpeed is positive clockwise and negative counter-clockwise, TurntableAngle (0 .. 360) is measured clockwise from screen right (plus ref_angle_deg).

#include <vector>
#include <memory>
#include "physics/rubbercontact.h"

class PhysicsEngine;
class PinTable;
class MoverObject;
class Renderer;
class ContactHost;   // (global: an elaborated name inside the namespace below would declare Turntables::ContactHost when nothing has declared it yet)

namespace Turntables
{
   void Setup(PhysicsEngine *physics, PinTable *table, std::vector<MoverObject *> &movers, std::vector<class ContactHost *> &hosts);   // builds the discs, their colliders and the holes (called when the physics engine is built)
   bool Active();                                            // at least one turntable exists
   bool SetVoltage(const char *name, float volts);           // false = no such turntable
   bool SetDuty(const char *name, float duty);               // PWM: mean voltage = duty x rated volts (-1..1)
   bool GetSpeedRpm(const char *name, float &rpm);
   bool GetAngleDeg(const char *name, float &deg);
   // contact side
   void Holes(const RubberContact::CarvedHole *&holes, int &n);
   void Render(Renderer *renderer);                          // draws the discs that have visible = 1 (and their gap rings)
   void ReleaseRender();
}
