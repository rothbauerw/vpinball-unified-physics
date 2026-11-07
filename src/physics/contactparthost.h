// license:GPLv3+

#pragma once
// The host of a moving part of the contact system (a gate, a spinner, a jet bumper, a turntable disc; docs/unified-physics/archive/UnifiedContactDesign.txt 12.5). It knows only the moving-part
// interface (contactmovingpart.h): the part probes, touches, steps and saves itself; the host registers it with the contact system, advances it alone when no island stepped
// it, and hands its pose and events to the table part's adapter. The adapter is the table part's own (gate, spinner, bumper, turntable files): the host knows no kind of part.

#include <memory>
#include <vector>

#include "physics/contacthost.h"
#include "physics/contactmovingpart.h"

// What a table part does with its moving part's state (VPX side). One per kind of table part, written in that part's own files.
class PartAdapter
{
public:
   virtual ~PartAdapter() = default;
   // Once per physics interval, after the part advanced: write its pose where the renderer and the script read it, and raise its detectors (rule 8)
   virtual void Publish() = 0;
   // The part's events of the interval, in the order they happened. Called after the step (ContactSystem::DeferEvent): a script may move or delete balls.
   virtual void Fire(const std::vector<ContactMovingPart::Event> &events) = 0;
   // The part's contact settings come from the table load only: a change of the ini while running does not reload them
   virtual bool FixedSettings() const { return false; }
};

class ContactPartHost final : public ContactHost
{
public:
   // lo, hi: the lowest and highest corner of what the part can sweep, VPX units (the engagement and the stock tree see this box)
   ContactPartHost(IEditable *editable, std::shared_ptr<ContactMovingPart> part, std::unique_ptr<PartAdapter> adapter, const RubberObjectSettings &settings,
                   const std::filesystem::path &iniPath, const std::string &name, const char *kind, const float lo[3], const float hi[3]);
   // The gravity a part feels (m/s^2, table frame): gravity minus the cabinet's nudge acceleration, as the balls feel it
   static RubberContact::V3 PartGravitySI();
   // Once per physics interval (the mover): advance the part alone for the time no island stepped it, publish its pose, queue its events
   void PartAdvance(float dtime) override;
   // The host's contact settings (the release gap, the material the contact system's step rule reads) when the part's own settings change (a flipper's rubber keys)
   void SetContactSettings(const RubberObjectSettings &s) { m_set = s; }

   // ContactSurface: the part is the surface
   ContactMovingPart *MovingPart() const override { return m_part.get(); }
   const RubberObjectSettings *Settings() const override { return nullptr; }   // a moving part sets its own step (BallStep), not substep_us
   void Refresh(SurfaceContact &it, double ballR, double reachM, const RubberContact::V3 &c) const override;
   RubberContact::TubeProbe Probe(const RubberContact::BallSI &bs, const SurfaceContact &it, const RubberContact::V3 &c) const override;
   void VelocityAt(const RubberContact::TubeProbe &pr, RubberContact::V3 &v, RubberContact::V3 &w) const override;
   RubberContact::ContactResult EvaluateOwn(SurfaceContact &it, const RubberContact::BallSI &bs, const RubberContact::V3 &c, const RubberContact::V3 &v, const RubberContact::V3 &w, double h) const override;
   double SurfaceSpeedBoundMS(double cycleS) const override;
   bool HitEventWanted(const RubberContact::V3 &) const override { return false; }   // the part fires its own Hit (its events)

   // HitObject: no collider of its own (HitTest finds nothing: the contact system engages the balls near the part)
   int GetType() const override { return e3DLine; }
   void Collide(const CollisionEvent &) override {}
   void CalcHitBBox() override;
   MoverObject *GetMoverObject() override { return &m_mover; }
   void DrawUI(std::function<Vertex2D(Vertex3Ds)>, ImDrawList *, bool) const override {}

private:
   void Reload() const override;

   class Mover final : public MoverObject
   {
   public:
      explicit Mover(ContactPartHost *owner) : m_owner(owner) {}
      bool AddToList() const override { return true; }
      void UpdateDisplacements(const float) override {}
      void UpdateVelocities() override {}
      void PreUpdateDisplacements(const float dtime) override { m_owner->PartAdvance(dtime); }

   private:
      ContactPartHost *m_owner;
   };

   std::shared_ptr<ContactMovingPart> m_part;
   std::unique_ptr<PartAdapter> m_adapter;
   float m_lo[3], m_hi[3];
   Mover m_mover{ this };
   bool m_gravityLogged = false;
};
