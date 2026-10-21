// license:GPLv3+

// Unified-engine additions to HitBall (members declared in hitball.h): resting-contact spin loss and the external velocity increment.

#include "core/stdafx.h"
#include "physics/hitball.h"
#include "parts/ball.h"
#include "parts/pintable.h"
#include "physics/contactsettings.h"
#include "physics/magnet.h"
#include "physics/PhysicsEngine.h"

void HitBall::ApplyRestingSpinLoss(const Vertex3Ds& hitnormal, const float fricCoeff, const float normalImpulse)
{
   const ContactSettings& cs = ContactSettings::Get();
   if (!cs.enable || !(normalImpulse > 0.0f))
      return;

   // Angular momentum is what the impulses act on (units: impulse * length), so each loss is capped at the spin it removes: it stops the spin, never reverses it.
   const float Ln = m_angularmomentum.Dot(hitnormal);
   const Vertex3Ds Lt = m_angularmomentum - Ln * hitnormal;

   // Rolling resistance (the pressure centre of a rolling ball sits ahead of the contact point): torque = rollingResistance * R * N against the spin in the tangent plane.
   const float ltLen = Lt.Length();
   if (cs.rollingResistance > 0.0 && ltLen > 1e-12f)
   {
      const float dL = std::min(ltLen, (float)cs.rollingResistance * m_d.m_radius * normalImpulse);
      m_angularmomentum -= Lt * (dL / ltLen);
   }

   // Twist (drilling) friction: the contact patch resists spin about the normal with torque (3 pi / 16) * mu_dynamic * a * N (the patch is slipping).
   if (cs.twist && cs.twistPatchMm > 0.0 && fabsf(Ln) > 1e-12f)
   {
      constexpr float kMmPerUnit = 0.53975f;
      const float a = (float)cs.twistPatchMm / kMmPerUnit;
      const float dL = std::min(fabsf(Ln), (float)(3.0 * M_PI / 16.0) * (float)cs.muDynamicRatio * fricCoeff * a * normalImpulse);
      m_angularmomentum -= hitnormal * (Ln > 0.0f ? dL : -dL);
   }
}

Vertex3Ds HitBall::ExternalVelocityIncrement(const bool withMagnets) const
{
   Vertex3Ds dv = (float)PHYS_FACTOR * m_physics->GetGravity() /* * m_d.m_mass / m_d.m_mass */; // Gravity F = m.G

   // Table velocity due to nudge (fictitious force due to change of reference frame, therefore mass is not applied)
   const float slope = ANGTORAD(m_pBall->GetPTable()->GetPlayfieldSlope()); // nudge acceleration is in the horizontal cabinet plane, reference frame is the playfield
   const Vertex2D cabinetAcceleration = m_physics->GetCabinetAcceleration();
   dv.x -= (float)PHYS_FACTOR * MS2TOVPUVPT2(cabinetAcceleration.x);
   dv.y -= (float)PHYS_FACTOR * MS2TOVPUVPT2(cabinetAcceleration.y) * cosf(slope);
   dv.z -= (float)PHYS_FACTOR * MS2TOVPUVPT2(cabinetAcceleration.y) * sinf(slope);

   if (withMagnets && m_magnetic && Magnets::Active())
   {
      const float p[3] = { m_d.m_pos.x, m_d.m_pos.y, m_d.m_pos.z };
      float mdv[3] = { m_magPair[0], m_magPair[1], m_magPair[2] };
      Magnets::AddVelocityIncrement(p, m_d.m_radius, m_d.m_mass, mdv);
      dv.x += mdv[0]; dv.y += mdv[1]; dv.z += mdv[2];
   }
   return dv;
}

// A ball the force-model contact holds on a surface (the playfield, a disc, a rubber) must not be shoved down into it: the surface is stiff (a 40 um overlap is kilonewtons),
// so a push of a few tenths of a unit from a ball resting on top of it comes back as a launch. The push may still lift it or move it sideways.
void HitBall::Shove(Vertex3Ds delta)
{
   if (delta.z < 0.f && HitRubberContactHolds(this))
      delta.z = 0.f;
   m_d.m_pos += delta;
}
