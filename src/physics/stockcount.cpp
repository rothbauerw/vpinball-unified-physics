// license:GPLv3+

#include "core/stdafx.h"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <typeinfo>
#include <unordered_map>
#include <vector>
#include "physics/stockcount.h"
#include "physics/contacthost.h"
#include "physics/contactsystem.h"
#include "physics/hitball.h"
#include "physics/contactsurfacehost.h"
#include "physics/testdriver.h"
#include "parts/ball.h"
#include "parts/trigger.h"
#include "core/player.h"
#include "physics/PhysicsEngine.h"
#include "ui/live/LiveUI.h"

namespace TestDriver { uint32_t ReportT(); }

namespace StockCount
{
   namespace
   {
      struct Row
      {
         std::string name, cls, kind;   // the part, the stock collider class, what kind of object
         const char *allowed = nullptr;   // why it may use stock (an unconverted part); null = converted: any use fails the run
         bool stockOnly = false;          // a collider with no contact model behind it: reaching its Collide is stock, whatever it did to the ball
         unsigned long long collide = 0, collideActs = 0, wall = 0, contact = 0, contactActs = 0, statics = 0;
         uint32_t firstT = 0;
         std::string firstBall;
         bool warned = false;
         bool handIn = false;             // a collider with a contact model behind it: its Collide or Contact hands the ball in (the engagement may seat it, that is the model's); stock here is the two laws by name
         bool Used() const { return wall != 0 || statics != 0 || (!handIn && (collideActs != 0 || contactActs != 0)) || (stockOnly && collide != 0); }
      };

      const void *s_owner = nullptr;
      const HitObject *s_playfield = nullptr, *s_glass = nullptr;
      std::vector<Row *> s_rows;
      std::unordered_map<const void *, Row *> s_byObject;
      Row *s_cur = nullptr, s_outside;
      HitBall *s_ball = nullptr, *s_other = nullptr;
      Vertex3Ds s_pos, s_vel, s_L, s_oPos, s_oVel, s_oL;
      bool s_contact = false;

      // What may use stock: the parts not converted yet, and nothing else. A part leaves this list the day it is converted.
      const char *AllowedStock(const ItemTypeEnum t)
      {
         switch (t)
         {
         case eItemPlunger: return "plunger: not converted (map step 8)";
         case eItemHitTarget: return "drop or stand-up target: not converted (map step 8)";
         case eItemTrigger: return "trigger: a sensor, not converted (map step 8)";
         default: return nullptr;
         }
      }

      Row *MakeRow(const void *key, const HitObject *obj)
      {
         auto r = new Row();
         if (obj != nullptr)
         {
            const char *n = typeid(*obj).name();
            r->cls = std::strncmp(n, "class ", 6) == 0 ? n + 6 : n;
         }
         const int type = obj != nullptr ? obj->GetType() : -1;
         r->stockOnly = type == eGate || type == eSpinner || type == ePlunger;
         r->handIn = obj != nullptr && (obj->m_surfaceHost != nullptr || dynamic_cast<const ContactHost *>(obj) != nullptr);
         if (obj == s_playfield)
         {
            r->name = "playfield";
            r->kind = "playfield";
         }
         else if (obj == s_glass)
         {
            r->name = "top glass";
            r->kind = "top glass";
            r->allowed = "top glass: no host yet (map step 8)";
         }
         else if (type == eBall)
         {
            r->name = "ball pair";
            r->kind = "ball pair";
         }
         else if (obj != nullptr && obj->m_editable != nullptr)
         {
            r->name = obj->m_editable->GetName();
            const ItemTypeEnum t = obj->m_editable->GetItemType();
            r->kind = t == eItemSurface ? "wall" : t == eItemRamp ? "ramp" : t == eItemPrimitive ? "primitive" : t == eItemRubber ? "rubber" : t == eItemFlipper ? "flipper" : t == eItemBumper ? "bumper"
               : t == eItemGate ? "gate" : t == eItemSpinner ? "spinner" : t == eItemKicker ? "kicker" : t == eItemTrigger ? "trigger" : t == eItemPlunger ? "plunger" : t == eItemHitTarget ? "target"
               : t == eItemTable ? "table" : "part";
            r->allowed = AllowedStock(t);
            if (t == eItemTrigger && dynamic_cast<const Trigger *>(obj->m_editable) != nullptr)
            {
               // converted, no stock allowed: a wire trigger on the wire switch part, a Star or Button (rollover caps of the group step) and a None (the engine's sensing zone, trigger.cpp:812-906).
               // Stock stays only for an Inder, and for a wire trigger with [Features] triggers = 0
               const Trigger *const tr = dynamic_cast<const Trigger *>(obj->m_editable);
               const TriggerShape sh = tr->m_d.m_shape;
               if (tr->OnWirePart() || sh == TriggerStar || sh == TriggerButton || sh == TriggerNone)
                  r->allowed = nullptr;
            }
            if (t == eItemTable)   // the table's own edges (PhysicsEngine): no part of the table to convert yet
            {
               r->name = "table edge";
               r->allowed = "table edge: no host yet, not in the map";
            }
         }
         else
         {
            r->name = "(no part)";
            r->kind = "?";
         }
         if (type == eLineSegSlingshot)
         {
            r->kind = "wall slingshot";
            r->allowed = "stock wall slingshot: decide convert or drop (map step 8)";
         }
         s_byObject[key] = r;
         s_rows.push_back(r);
         return r;
      }

      void Fresh()
      {
         if (s_owner == g_pplayer)
            return;
         s_owner = g_pplayer;
         for (Row *const r : s_rows)
            delete r;
         s_rows.clear();
         s_byObject.clear();
         s_cur = nullptr;
      }

      void Warn(Row &r)
      {
         if (r.warned || r.allowed != nullptr)
            return;
         r.warned = true;
         const std::string msg = "Stock used by a converted part: " + r.kind + " " + r.name + " (" + r.cls + "), t=" + std::to_string(ReportTime()) + ": a route to remove";
         PLOGW << msg;
         if (g_pplayer->m_liveUI != nullptr)
            g_pplayer->m_liveUI->PushNotification(msg, 8000);
      }
   }

   uint32_t ReportTime() { return TestDriver::ReportT(); }

   void SetPlanes(const HitObject *const playfield, const HitObject *const topGlass)
   {
      s_playfield = playfield;
      s_glass = topGlass;
   }

   bool Off()
   {
      static const bool off = []() { const char *const e = getenv("VPX_CHECKS"); return e != nullptr && *e == '0'; }();
      return off;
   }

   void Begin(HitObject *const obj, HitBall *const ball, const bool contact)
   {
      if (Off() || !ContactSystem::Unified())   // stock engine: nothing is converted, nothing to count
         return;
      Fresh();
      const void *key = obj;
      const HitBall *const other = obj->GetType() == eBall ? static_cast<const HitBall *>(obj) : nullptr;
      const auto it = s_byObject.find(key);
      Row *r = it != s_byObject.end() ? it->second : nullptr;
      if (r == nullptr)
      {
         r = MakeRow(key, obj);
      }
      s_cur = r;
      s_ball = ball;
      s_other = other != nullptr ? const_cast<HitBall *>(other) : nullptr;
      s_contact = contact;
      s_pos = ball->m_d.m_pos;
      s_vel = ball->m_d.m_vel;
      s_L = ball->m_angularmomentum;
      if (s_other != nullptr)
      {
         s_oPos = s_other->m_d.m_pos;
         s_oVel = s_other->m_d.m_vel;
         s_oL = s_other->m_angularmomentum;
      }
      if (contact)
         r->contact++;
      else
         r->collide++;
   }

   void End(const bool ballAlive)
   {
      Row *const r = s_cur;
      s_cur = nullptr;
      if (r == nullptr || !ballAlive)
         return;
      const HitBall &b = *s_ball;
      bool acted = std::memcmp(&s_pos, &b.m_d.m_pos, sizeof s_pos) != 0 || std::memcmp(&s_vel, &b.m_d.m_vel, sizeof s_vel) != 0 || std::memcmp(&s_L, &b.m_angularmomentum, sizeof s_L) != 0;
      if (s_other != nullptr)
      {
         const HitBall &o = *s_other;
         if (std::memcmp(&s_oPos, &o.m_d.m_pos, sizeof s_oPos) != 0 || std::memcmp(&s_oVel, &o.m_d.m_vel, sizeof s_oVel) != 0 || std::memcmp(&s_oL, &o.m_angularmomentum, sizeof s_oL) != 0)
            acted = true;
      }
      if (acted)
         (s_contact ? r->contactActs : r->collideActs)++;
      if (acted)   // test tooling (VPX_TEST_IMPULSES=1): the change this stock call made to the ball's momentum and angular momentum, by its object's name (no-op otherwise)
      {
         const float v0[3] = { s_vel.x, s_vel.y, s_vel.z }, L0[3] = { s_L.x, s_L.y, s_L.z };
         TestDriver::AddBallChange(s_ball, r->name.empty() ? std::string("stock") : r->name, v0, L0);
         if (s_other != nullptr)
         {
            const float ov0[3] = { s_oVel.x, s_oVel.y, s_oVel.z }, oL0[3] = { s_oL.x, s_oL.y, s_oL.z };
            const char *const n = TestDriver::BallName(s_ball);
            TestDriver::AddBallChange(s_other, std::string("ball") + (n != nullptr ? n : "?"), ov0, oL0);
         }
      }
      if (r->firstBall.empty() && r->Used())
      {
         r->firstT = ReportTime();
         const char *const n = TestDriver::BallName(s_ball);
         r->firstBall = n != nullptr ? n : "?";
      }
      if (r->Used())
         Warn(*r);
   }

   void Wall()
   {
      if (Off() || !ContactSystem::Unified())
         return;
      (s_cur != nullptr ? *s_cur : s_outside).wall++;
   }

   void Static()
   {
      if (Off() || !ContactSystem::Unified())
         return;
      (s_cur != nullptr ? *s_cur : s_outside).statics++;
   }

   bool ConvertedUsed()
   {
      if (s_outside.wall != 0 || s_outside.statics != 0)
         return true;
      for (const Row *const r : s_rows)
         if (r->allowed == nullptr && r->Used())
            return true;
      return false;
   }

   void Report(std::ostream &o, const std::string &csvPath)
   {
      if (Off())
      {
         o << "\nstock counter: off (VPX_CHECKS=0)\n";
         return;
      }
      // one line per object (a wall has many colliders): name, kind and class together
      std::vector<Row> merged;
      for (const Row *const r : s_rows)
      {
         auto m = std::find_if(merged.begin(), merged.end(), [&](const Row &x) { return x.name == r->name && x.kind == r->kind && x.cls == r->cls; });
         if (m == merged.end())
         {
            merged.push_back(*r);
            continue;
         }
         m->collide += r->collide; m->collideActs += r->collideActs; m->wall += r->wall; m->contact += r->contact; m->contactActs += r->contactActs; m->statics += r->statics;
         if (m->firstBall.empty() || (!r->firstBall.empty() && r->firstT < m->firstT))
         {
            m->firstT = r->firstT;
            m->firstBall = r->firstBall;
         }
      }
      std::stable_sort(merged.begin(), merged.end(), [](const Row &a, const Row &b) { return (a.allowed == nullptr) != (b.allowed == nullptr) ? a.allowed == nullptr : a.name < b.name; });
      int bad = s_outside.wall != 0 || s_outside.statics != 0 ? 1 : 0, allowed = 0, seen = 0;
      for (const Row &r : merged)
         if (r.Used())
            (r.allowed == nullptr ? bad : allowed)++;
         else
            seen++;
      o << "\nstock counter: " << (bad == 0 ? "PASS" : "FAIL") << ", " << bad << " converted object(s) used stock, " << allowed << " allowed unconverted part(s) used it, " << seen << " object(s) reached with no stock act (sensor or hand-in)\n";
      char line[400];
      if (s_outside.wall != 0 || s_outside.statics != 0)
         o << "   FAIL: stock law run outside the dispatch: Collide3DWall " << s_outside.wall << ", HandleStaticContact " << s_outside.statics << "\n";
      if (!merged.empty())
         o << "   object | kind | collider | collide calls / changed the ball / wall law | contact calls / changed the ball / static law | first t, ball | verdict\n";
      for (const Row &r : merged)
      {
         std::snprintf(line, sizeof line, " | %llu / %llu / %llu | %llu / %llu / %llu | %u, %s | ", r.collide, r.collideActs, r.wall, r.contact, r.contactActs, r.statics, r.firstT, r.firstBall.c_str());
         o << "   " << r.name << " | " << r.kind << " | " << r.cls << line << (!r.Used() ? (r.handIn ? "hand-in only" : "no stock act") : r.allowed != nullptr ? std::string("allowed, ") + r.allowed : std::string("FAIL: converted, route to remove")) << '\n';
      }
      if (FILE *const f = fopen(csvPath.c_str(), "w"))
      {
         fputs("object,kind,collider,collide_calls,collide_acted,wall_law,contact_calls,contact_acted,static_law,first_t_ms,first_ball,verdict\n", f);
         for (const Row &r : merged)
            fprintf(f, "\"%s\",%s,\"%s\",%llu,%llu,%llu,%llu,%llu,%llu,%u,%s,%s\n", r.name.c_str(), r.kind.c_str(), r.cls.c_str(), r.collide, r.collideActs, r.wall, r.contact, r.contactActs, r.statics, r.firstT, r.firstBall.c_str(),
                    !r.Used() ? "no_act" : r.allowed != nullptr ? "allowed" : "FAIL");
         fclose(f);
      }
   }
}
