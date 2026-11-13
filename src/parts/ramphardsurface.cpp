// license:GPLv3+

// Unified engine: the ramp's floor and side walls as a hard surface of the contact model (Ramp::AddHardSurface, called from Ramp::PhysicSetup in ramp.cpp).

#include "core/stdafx.h"
#include "ramp.h"
#include "physics/PhysicsEngine.h"
#include "physics/contactsurfacehost.h"
#include "physics/contactsystem.h"

void Ramp::AddHardSurface(PhysicsEngine* physics, const bool isUI, const Vertex2D* const rgvLocal, const float* const rgheight1, const int cvertex, const float wallheightright, const float wallheightleft, const size_t collidersBefore)
{
   // unified engine: the ramp's floor and side walls are a hard surface of the contact model (sheets: the ball side is taken from where the ball is)
   std::shared_ptr<RubberSurface> rubberFaces;
   RubberObjectSettings rubberSettings;
   const std::filesystem::path rubberIni = ContactSurfaceHost::IniPathFor(m_ptable->m_filename);
   if (!isUI && ContactSystem::Unified() && ContactSurfaceHost::HardSurfaceWanted(rubberIni, "ramp", GetName()))
      rubberFaces = std::make_shared<RubberSurface>();
   if (rubberFaces)
   {
      // a sheet face; a degenerate one (the ramp is 0 wide at a vertex) is skipped
      auto addFace = [&](const Vertex3Ds &a, const Vertex3Ds &b, const Vertex3Ds &c) {
         if (CrossProduct(b - a, c - a).LengthSquared() > 1e-6f)
            rubberFaces->Add(a, b, c);
      };
      for (int i = 0; i < (cvertex - 1); i++)
      {
         // floor, as the quad layout above
         const Vertex2D &pv1 = rgvLocal[i];
         const Vertex2D &pv2 = rgvLocal[cvertex * 2 - i - 1];
         const Vertex2D &pv3 = rgvLocal[cvertex * 2 - i - 2];
         const Vertex2D &pv4 = rgvLocal[i + 1];
         const Vertex3Ds p1(pv1.x, pv1.y, rgheight1[i]), p2(pv2.x, pv2.y, rgheight1[i]), p3(pv3.x, pv3.y, rgheight1[i + 1]), p4(pv4.x, pv4.y, rgheight1[i + 1]);
         addFace(p2, p1, p3);
         addFace(p3, p1, p4);
         // side walls: a sloped quad from the floor to the wall height at both ends
         if (wallheightright > 0.f)
         {
            const Vertex3Ds t1(pv1.x, pv1.y, rgheight1[i] + wallheightright), t4(pv4.x, pv4.y, rgheight1[i + 1] + wallheightright);
            addFace(p1, p4, t4);
            addFace(p1, t4, t1);
         }
         if (wallheightleft > 0.f)
         {
            const Vertex3Ds t2(pv2.x, pv2.y, rgheight1[i] + wallheightleft), t3(pv3.x, pv3.y, rgheight1[i + 1] + wallheightleft);
            addFace(p2, p3, t3);
            addFace(p2, t3, t2);
         }
      }
      if (!rubberFaces->tris.empty())
      {
         const HitObject *const first = collidersBefore < m_vhoCollidable.size() ? m_vhoCollidable[collidersBefore] : nullptr;
         ContactSurfaceHost::LoadHardSurfaceSettings(rubberSettings, rubberIni, "ramp", GetName(), ContactSurfaceHost::HardSurfaceE(first), first != nullptr ? first->m_friction : 0.3);
         ContactSurfaceHost *const host = new ContactSurfaceHost(this, rubberFaces, rubberSettings, rubberIni, GetName(), "ramp");
         if (first != nullptr)
            host->SetHardSurface(ContactSurfaceHost::HardSurfaceE(first), first->m_friction);
         host->SetEventCarrier(collidersBefore < m_vhoCollidable.size() ? m_vhoCollidable[collidersBefore] : nullptr);
         for (size_t i = collidersBefore; i < m_vhoCollidable.size(); i++)
            m_vhoCollidable[i]->m_surfaceHost = host;   // the ramp's own colliders hand a hit to the host
         physics->KeepHostedCopiesOutOfTrees(m_vhoCollidable, collidersBefore, rubberIni);   // shared_finder P9: the copies out of stock's trees (PhysicsEngine.cpp)
         SetupHitObject(physics, host, isUI);           // the host is enabled and disabled with them (Collidable)
         PLOGI << "Ramp " << GetName() << ": hard surface on (" << rubberFaces->tris.size() << " faces, mu " << rubberSettings.mu << ", alpha " << rubberSettings.alphaSPerM << ")";
      }
   }
}
