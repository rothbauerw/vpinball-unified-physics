// license:GPLv3+

#include "core/stdafx.h"
#include "physics/turntable.h"
#include "physics/meshholes.h"
#include "physics/contactparthost.h"
#include "physics/contactdiscpart.h"
#include "physics/physicsini.h"
#include "physics/contactsystem.h"
#include "physics/physconst.h"
#include "physics/PhysicsEngine.h"
#include "physics/collideex.h"
#include "parts/pintable.h"
#include "core/player.h"
#include "renderer/Renderer.h"
#include "renderer/RenderDevice.h"
#include "renderer/MeshBuffer.h"
#include "renderer/VertexBuffer.h"
#include "renderer/IndexBuffer.h"
#include "renderer/Shader.h"
#include "ui/live/LiveUI.h"
#include <cmath>
#include <map>

namespace Turntables
{
   namespace
   {
      constexpr double PI = 3.14159265358979323846;
      constexpr double UM = 0.53975e-3;   // metres per table unit
      constexpr int SEG = 96;             // sides of the disc polygon

      struct Turntable
      {
         std::string name;
         // ini
         double x = 0, y = 0;                  // centre, table units
         double z = 0;                         // height of the disc top, table units; default = the playfield surface (flush)
         double playfieldZ = 0;                // height of the playfield the hole is cut in; default = the playfield surface at the centre
         bool zSet = false, playfieldZSet = false;
         double radiusMm = 76.2;               // disc radius (3 in)
         double gapMm = 1.5875;                // gap between the disc and the playfield (1/16 in)
         bool clockwise = true;                // direction a positive voltage turns the disc, seen from above
         double massG = 800;                   // disc and hub: inertia = 1/2 m r^2 unless inertia_kgm2 is given
         double inertia = 0;                   // kg m^2, 0 = from the mass
         double freeRpm = 250;                 // max_rpm: the speed at supply_volts with no load (the Whirlwind discs top out at 200 - 300 rpm)
         double ratedVolts = 12;               // supply_volts
         double spinupS = 1.0;                 // ramp_up_s: time from rest to 95% of max_rpm at supply_volts, no load
         double rampDownS = 4.0;               // ramp_down_s: time from max_rpm to 5% of it with the motor off, no load (bearing drag)
         double frictionNm = 0.02;             // friction_nm: constant (Coulomb) torque of the bearing and gear train: it holds the disc still until a ball's push exceeds it
         double mu = 0.6;                      // friction coefficient of the textured surface against the steel ball
         double rollingResistance = 0.03;      // rolling_resistance: coefficient of a ball rolling on the textured (rubber-like) surface; torque = coefficient * R * N. The playfield's is [Contact] rolling_resistance (0.003); a negative value here = use that one
         bool visible = false;                 // visible: the disc drawn as the physical disc (radius_mm at z), no gap ring: the playfield's hole, gap included, is drawn by the hole drawing
         std::string marker = "spokes";        // visible: spokes | dot | stripe | none
         double refDeg = 0;                    // angle reference, degrees clockwise from the +x direction (screen right)
         MeshHoles::HoleKeys hole;             // cut_hole, draw, wall_depth_mm, host: the disc's hole in the mesh under it and the mesh whose top it is flush with (meshholes.h HoleKeys)
         // the disc as a moving part of the contact system, on the group step (contactdiscpart.h, discgroup.h): motor, bearing, speed and angle (clockwise positive), the script's voltage
         std::shared_ptr<DiscModel::Engine> part;
         std::shared_ptr<MeshBuffer> discMesh, markMesh;
      };

      // The table side of the disc's moving part (UnifiedContactDesign 12.1): the renderer and the script read the part itself, it raises no events (it has no table object to raise them on), and its settings are
      // read when the table loads only
      class DiscAdapter final : public PartAdapter
      {
      public:
         DiscAdapter(std::shared_ptr<DiscModel::Engine> model, const std::string &name, const std::string &note) : m_model(std::move(model)), m_name(name), m_note(note) {}
         void Publish() override
         {
            m_model->bookMs = (int)g_pplayer->m_time_msec;
            if (!m_note.empty() && g_pplayer->m_liveUI != nullptr)   // (what the table load found, shown once the player can show it)
            {
               g_pplayer->m_liveUI->PushNotification(m_note, 6000);
               m_note.clear();
            }
            // detector (rule 8): a step of the disc alone was taken on energy, or could not be solved
            if (m_model->energySteps > m_esSeen)
            {
               m_esSeen = m_model->energySteps;
               if (m_esLines < 5)
               {
                  m_esLines++;
                  PLOGW << "Group step: energy step alone, part " << m_name << ", balls 0, t=" << g_pplayer->m_time_msec << " ms, " << m_esSeen << " so far; " << m_model->StateNote() << ". A net, not expected to fire: keep vpinball.log and report it";
               }
            }
            if (m_model->failedSteps > 0 && !m_failShown)
            {
               m_failShown = true;
               const std::string msg = m_name + ": a step of the disc alone could not be solved: the rest of that cycle was lost";
               PLOGW << msg;
               if (g_pplayer->m_liveUI != nullptr)
                  g_pplayer->m_liveUI->PushNotification(msg, 6000);
            }
         }
         void Fire(const std::vector<ContactMovingPart::Event> &) override {}
         bool ReloadSettings() override { return true; }   // (the disc keeps the settings it was built with: its keys are read when the table loads)

      private:
         std::shared_ptr<DiscModel::Engine> m_model;
         std::string m_name, m_note;
         long m_esSeen = 0;
         int m_esLines = 0;
         bool m_failShown = false;
      };

      std::vector<Turntable> g_tt;
      std::map<std::string, size_t> g_byName;
      std::map<std::string, std::pair<bool, float>> g_pending;   // script calls made before the physics started: (isDuty, value)

      void Derive(Turntable &t)
      {
         DiscGroup::Motor &m = t.part->gp.motor;
         const double r = t.radiusMm * 1e-3;
         m.I = t.inertia > 0 ? t.inertia : 0.5 * t.massG * 1e-3 * r * r;
         m.freeRpm = t.freeRpm;
         m.ratedVolts = t.ratedVolts;
         m.spinupS = t.spinupS;
         m.rampDownS = t.rampDownS;
         m.frictionNm = t.frictionNm;
         m.dir = t.clockwise ? 1.0 : -1.0;
         m.dynRatio = ContactSettings::Get().muDynamicRatio;
         m.radiusM = r;
         m.Derive();
         PLOGI << "Turntable " << t.name << ": inertia " << m.I << " kg m^2 from ramp_down_s " << t.rampDownS << " s and the sliding torque " << m.tSlide << " N m (" << m.dynRatio << " x friction_nm " << t.frictionNm << "), drag " << m.b << " N m s"
               << (m.floorApplied ? "; floor applied: the friction alone would stop the disc's own inertia (1/8 in steel plate) too soon, the drag makes up the rest" : "; no floor");
      }

      // The playfield's surface around a disc: the height the disc is flush with. In each of 24 directions, the nearest playfield face outward from the disc's edge, taking the highest face at a point; a direction that finds none has no say (past the table's edge).
      // So it holds for a mesh cut under the disc's centre, cut wider than the disc and its gap, or cut off-centre. The height is the median of what the directions found, so a uniform slope gives the height at the disc's centre and a ramp's faces in less than half
      // the directions have no say. farU: how far out from the disc's edge the furthest direction had to go to find a face, and which way (degrees clockwise from screen right).
      struct AroundDisc
      {
         double z = 0.0, lo = 0.0, hi = 0.0, farU = 0.0, farDeg = 0.0;
         int points = 0;
      };
      bool SurfaceAroundDisc(const double x, const double y, const double rU, AroundDisc &o)
      {
         constexpr int kDirs = 24;
         constexpr double kFirstU = 0.5, kLastU = 8192.0, kTolU = 0.1;
         const auto height = [](const double px, const double py, double &z) {
            bool any = false;
            for (ContactHost *const h : ContactHost::AllHosts())
            {
               double zh;
               if (h->IsPlayfieldKind() && h->SurfaceHeightAt(px, py, zh) && (!any || zh > z))
               {
                  z = zh;
                  any = true;
               }
            }
            return any;
         };
         std::vector<double> found;
         for (int k = 0; k < kDirs; k++)
         {
            const double a = 2.0 * PI * k / kDirs, cx = std::cos(a), cy = std::sin(a);
            double miss = 0.0, hit = kFirstU, z = 0.0;
            while (hit < kLastU && !height(x + (rU + hit) * cx, y + (rU + hit) * cy, z))
            {
               miss = hit;
               hit *= 2.0;
            }
            if (hit >= kLastU)
               continue;
            while (hit - miss > kTolU)
            {
               const double mid = 0.5 * (hit + miss);
               double zm = 0.0;
               if (height(x + (rU + mid) * cx, y + (rU + mid) * cy, zm))
               {
                  hit = mid;
                  z = zm;
               }
               else
                  miss = mid;
            }
            found.push_back(z);
            if (hit > o.farU)
            {
               o.farU = hit;
               o.farDeg = 360.0 * k / kDirs;
            }
         }
         if (found.empty())
            return false;
         std::sort(found.begin(), found.end());
         const size_t n = found.size();
         o.z = 0.5 * (found[(n - 1) / 2] + found[n / 2]);
         o.lo = found.front();
         o.hi = found.back();
         o.points = (int)n;
         return true;
      }

      void Parse(const std::filesystem::path &file)
      {
         std::ifstream in(file);
         std::string line;
         Turntable *cur = nullptr;
         while (std::getline(in, line))
         {
            const size_t c = line.find_first_of("#;");
            if (c != std::string::npos)
               line.resize(c);
            line = PhysicsIni::Trim(line);
            if (line.empty())
               continue;
            if (line.front() == '[' && line.back() == ']')
            {
               cur = nullptr;
               const std::string head = PhysicsIni::Trim(line.substr(1, line.size() - 2));
               if (PhysicsIni::Lower(head).rfind("turntable:", 0) == 0)
               {
                  g_tt.emplace_back();
                  cur = &g_tt.back();
                  cur->name = PhysicsIni::Trim(head.substr(10));
               }
               continue;
            }
            const size_t eq = line.find('=');
            if (cur == nullptr || eq == std::string::npos)
               continue;
            const std::string key = PhysicsIni::Lower(PhysicsIni::Trim(line.substr(0, eq)));
            const std::string raw = PhysicsIni::Trim(line.substr(eq + 1)), sval = PhysicsIni::Lower(raw);
            const double v = atof(sval.c_str());
            if (key == "x") cur->x = v;
            else if (key == "y") cur->y = v;
            else if (key == "z") { cur->z = v; cur->zSet = true; }
            else if (key == "playfield_z") { cur->playfieldZ = v; cur->playfieldZSet = true; }
            else if (key == "radius_mm") cur->radiusMm = v;
            else if (key == "gap_mm") cur->gapMm = v;
            else if (key == "direction") cur->clockwise = !(sval == "ccw" || sval == "counterclockwise" || sval == "counter-clockwise");
            else if (key == "mass_g") cur->massG = v;
            else if (key == "inertia_kgm2") cur->inertia = v;
            else if (key == "max_rpm" || key == "free_rpm") cur->freeRpm = v;
            else if (key == "supply_volts" || key == "rated_volts") cur->ratedVolts = v;
            else if (key == "ramp_up_s" || key == "spinup_s") cur->spinupS = v;
            else if (key == "ramp_down_s") cur->rampDownS = v;
            else if (key == "coast_s") cur->rampDownS = 3.0 * v;   // old key: the time constant of the slow-down
            else if (key == "friction_nm") cur->frictionNm = v;
            else if (key == "mu") cur->mu = v;
            else if (key == "rolling_resistance") cur->rollingResistance = v;
            else if (key == "visible") cur->visible = v != 0;
            else if (key == "marker") cur->marker = sval;
            else if (key == "ref_angle_deg") cur->refDeg = v;
            else if (MeshHoles::HoleKey(cur->hole, key, v, raw, "Turntable " + cur->name)) {}   // cut_hole, draw, wall_depth_mm, host (the disc sits on that mesh's top; host wins over playfield_z) (meshholes.h HoleKeys)
         }
      }

   }

   void Setup(PinTable *const table, std::vector<MoverObject *> &movers, std::vector<ContactHost *> &hosts)
   {
      g_tt.clear();
      g_byName.clear();
      const std::filesystem::path file = PhysicsIni::FindUnified(table->m_filename);
      if (file.empty())
         return;
      const PhysicsIni::Features feat = PhysicsIni::Active(file);
      if (!feat.turntables)
         return;
      g_tt.reserve(16);
      Parse(file);
      if (g_tt.empty())
         return;
      if (!ContactSystem::Unified())
      {
         // the disc is a part of the unified contact model only (the contact finder takes it from the host registry): stock physics has no collider for it
         const std::string msg = "Turntables: " + std::to_string(g_tt.size()) + " [Turntable:<name>] section(s) in " + file.filename().string() + " ignored: turntables need unified physics ([Features] unified = 1)";
         PLOGW << msg;
         if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
            g_pplayer->m_liveUI->PushNotification(msg, 10000);
         g_tt.clear();
         return;
      }
      if (!feat.playfieldContact)
      {
         // the disc is a surface of the playfield's force model and the playfield is cut round it: with the stock playfield contact the plane would still carry the ball
         PLOGW << "Turntables: " << g_tt.size() << " [Turntable:<name>] section(s) ignored: turntables need [Features] playfield_contact = 1";
         g_tt.clear();
         return;
      }
      for (size_t i = 0; i < g_tt.size(); i++)
      {
         Turntable &t = g_tt[i];
         t.part = std::make_shared<DiscModel::Engine>();
         Derive(t);
         g_byName[PhysicsIni::Lower(t.name)] = i;
         const auto pend = g_pending.find(PhysicsIni::Lower(t.name));
         if (pend != g_pending.end())
            t.part->SetVolts(pend->second.first ? pend->second.second * t.ratedVolts : pend->second.second);
         std::string note;
         bool onHost = false, cutHole = true;   // (cutHole false: the named host refused)
         {
            AroundDisc around;
            const double gapU = t.gapMm * 1e-3 / UM;
            const bool found = SurfaceAroundDisc(t.x, t.y, t.radiusMm * 1e-3 / UM, around);
            if (found)
               PLOGI << "Turntable " << t.name << ": the playfield around it is at z " << around.z << " (" << around.points << " directions, " << around.lo << " .. " << around.hi << "), furthest " << around.farU * UM * 1e3 << " mm from its edge";
            if (!t.hole.host.empty())   // host = <static mesh>: the disc is flush with that mesh's top at its centre (host wins over playfield_z); refused (PartOnHost says why): as without it, no hole
            {
               double z = 0.0;
               cutHole = onHost = MeshHoles::PartOnHost("Turntable " + t.name, t.hole.host, t.x, t.y, z);
               if (onHost && t.playfieldZSet)
               {
                  const std::string msg = "Turntable " + t.name + ": playfield_z " + std::to_string(t.playfieldZ) + " ignored: host " + t.hole.host + " sets its height (" + std::to_string(z) + "). Remove one of them from its ini section";
                  PLOGW << msg;
                  if (g_pplayer != nullptr && g_pplayer->m_liveUI != nullptr)
                     g_pplayer->m_liveUI->PushNotification(msg, 10000);
               }
               if (onHost)
                  t.playfieldZ = z;
            }
            if (!onHost && !t.playfieldZSet)
            {
               t.playfieldZ = found ? around.z : 0.0;
               if (!found)
                  note = "Turntable " + t.name + ": no playfield surface found around it: the disc is at z 0. Set z in its ini section, or check the playfield";
            }
            if (found && around.farU > gapU + 0.1)   // the table's own geometry: an open gap a ball can drop into. Said, not filled
               note = "Turntable " + t.name + ": the playfield's own hole reaches " + std::to_string((around.farU - gapU) * UM * 1e3) + " mm past the disc and its gap, at " + std::to_string((int)around.farDeg) + " deg clockwise from screen right: a ball can drop in there";
            if (!note.empty())
               PLOGW << note;
         }
         if (!t.zSet)
            t.z = t.playfieldZ;
         const double R = t.radiusMm * 1e-3 / UM;   // table units
         const double hole = (t.radiusMm + t.gapMm) * 1e-3;
         if (cutHole)
            MeshHoles::Add("Turntable " + t.name, MeshHoles::Round(t.x * UM, t.y * UM, hole), t.playfieldZ * UM, t.hole);   // the disc and its gap cut into the mesh under it at load (meshholes.h)

         // the disc: a flat polygon at the top height, plus a skirt down to the playfield when it stands proud of it: the host's box (the contact finder takes the disc from the host registry; no stock collider)
         std::vector<Vertex3Ds> top;
         for (int k = 0; k < SEG; k++)
         {
            const double a = 2.0 * PI * k / SEG;
            top.emplace_back((float)(t.x + R * std::cos(a)), (float)(t.y + R * std::sin(a)), (float)t.z);
         }
         const Vertex3Ds ctr((float)t.x, (float)t.y, (float)t.z);
         const bool proud = t.z > t.playfieldZ + 0.5;
         const float zBottom = proud ? (float)t.playfieldZ - 2.f : (float)t.z;
         struct Tri { Vertex3Ds v[3]; };
         std::vector<Tri> tris;
         for (int k = 0; k < SEG; k++)
         {
            const Vertex3Ds &p0 = top[k], &p1 = top[(k + 1) % SEG];
            tris.push_back({ { ctr, p1, p0 } });   // the face looks up
            if (proud)
            {
               const Vertex3Ds q0(p0.x, p0.y, zBottom), q1(p1.x, p1.y, zBottom);
               tris.push_back({ { p0, p1, q1 } });
               tris.push_back({ { p0, q1, q0 } });
            }
         }
         float lo[3] = { tris[0].v[0].x, tris[0].v[0].y, tris[0].v[0].z }, hi[3] = { lo[0], lo[1], lo[2] };
         for (const Tri &tr : tris)
            for (const Vertex3Ds &v : tr.v)
               for (int k = 0; k < 3; k++)
               {
                  lo[k] = std::min(lo[k], (&v.x)[k]);
                  hi[k] = std::max(hi[k], (&v.x)[k]);
               }

         // the disc is the host's moving part: it probes its own shape (a flat disc, and its side when it stands proud), turns with the balls in one solve (the group step);
         // it takes the host's contact law (material, layer, edge rounding)
         const RubberObjectSettings ps = ContactHost::TurntableSettings(t.mu, t.rollingResistance);
         ContactPartHost *const host = new ContactPartHost(table, t.part, std::make_unique<DiscAdapter>(t.part, t.name, note), ps, file, t.name, "turntable", lo, hi);
         t.part->gp.Build(t.part->gp.motor, GroupStep::V3(t.x * UM, t.y * UM, t.z * UM), t.radiusMm * 1e-3, (double)zBottom * UM, ps.edgeRadiusMm * 1e-3);
         t.part->gp.dynRatio = ContactSettings::Get().muDynamicRatio;
         t.part->rp = ps.ToParams();
         t.part->layerM = host->LayerM();
         hosts.push_back(host);   // the physics engine owns it and deletes it with itself (the host is in the process-wide registry: a leaked one is tested against every ball in every later session)
         host->CalcHitBBox();
         if (MoverObject *const pmo = host->GetMoverObject(); pmo != nullptr && pmo->AddToList())
            movers.push_back(pmo);
         PLOGI << "Turntable " << t.name << ": centre " << t.x << ", " << t.y << ", top z " << t.z << ", radius " << t.radiusMm << " mm, gap " << t.gapMm << " mm, "
               << (t.clockwise ? "clockwise" : "counter-clockwise") << ", " << t.freeRpm << " rpm at " << t.ratedVolts << " V, friction " << t.frictionNm << " N m, ramp up " << t.spinupS << " s, ramp down " << t.rampDownS << " s, I " << t.part->gp.motor.I << " kg m^2, stall torque " << t.part->gp.motor.Ts
               << " N m, drag " << t.part->gp.motor.b << " N m s, mu " << t.mu << ", rolling resistance " << (t.rollingResistance >= 0.0 ? t.rollingResistance : ContactSettings::Get().rollingResistance) << (proud ? ", disc stands proud of the playfield" : ", disc flush with the playfield");
      }
      g_pending.clear();
   }

   bool Active() { return !g_tt.empty(); }

   static Turntable *Find(const char *name)
   {
      const auto it = g_byName.find(PhysicsIni::Lower(name));
      return it == g_byName.end() ? nullptr : &g_tt[it->second];
   }

   bool SetVoltage(const char *name, float volts)
   {
      if (g_tt.empty())
      {
         g_pending[PhysicsIni::Lower(name)] = { false, volts };
         return true;
      }
      Turntable *const t = Find(name);
      if (t == nullptr)
         return false;
      t->part->SetVolts(volts);
      return true;
   }

   bool SetDuty(const char *name, float duty)
   {
      if (g_tt.empty())
      {
         g_pending[PhysicsIni::Lower(name)] = { true, duty };
         return true;
      }
      Turntable *const t = Find(name);
      if (t == nullptr)
         return false;
      t->part->SetVolts(std::max(-1.0, std::min(1.0, (double)duty)) * t->ratedVolts);
      return true;
   }

   bool GetSpeedRpm(const char *name, float &rpm)
   {
      rpm = 0.f;
      const Turntable *const t = Find(name);
      if (t == nullptr)
         return false;
      rpm = (float)(t->part->gp.rate * 60.0 / (2.0 * PI));   // fixed frame: clockwise seen from above is positive, whatever the direction key
      return true;
   }

   bool GetAngleDeg(const char *name, float &deg)
   {
      deg = 0.f;
      const Turntable *const t = Find(name);
      if (t == nullptr)
         return false;
      double a = std::fmod((t->part->gp.Angle() * 180.0 / PI) + t->refDeg, 360.0);   // degrees clockwise from screen right, plus the reference
      if (a < 0)
         a += 360.0;
      deg = (float)a;
      return true;
   }

   namespace
   {
      void FillVertex(Vertex3D_NoTex2 &q, double x, double y, double z)
      {
         q = {};
         q.x = (float)x;
         q.y = (float)y;
         q.z = (float)z;
         q.nz = 1.f;
      }
   }

   // The disc (as the physical disc: radius_mm at its top z) and (when marker is not none) the pattern that makes the turning visible. Only the pattern is rewritten every frame; the plain disc looks the same at every angle.
   void Render(Renderer *const renderer)
   {
      for (Turntable &t : g_tt)
      {
         if (!t.visible)
            continue;
         RenderDevice *const rd = renderer->m_renderDevice;
         const double R = t.radiusMm * 1e-3 / UM;
         const double zTop = t.z, zMark = t.z + 0.3;   // the disc at its physical top; the turning pattern just above it
         auto build = [&](std::vector<Vertex3D_NoTex2> &v, std::vector<unsigned int> &idx, const double ang, const bool markers)
         {
            v.clear();
            idx.clear();
            auto sector = [&](double r0, double r1, double a0, double a1, double z, int steps)
            {
               const unsigned int base = (unsigned int)v.size();
               for (int k = 0; k <= steps; k++)
               {
                  const double a = a0 + (a1 - a0) * k / steps;
                  v.emplace_back();
                  FillVertex(v.back(), t.x + r0 * std::cos(a), t.y + r0 * std::sin(a), z);
                  v.emplace_back();
                  FillVertex(v.back(), t.x + r1 * std::cos(a), t.y + r1 * std::sin(a), z);
               }
               for (int k = 0; k < steps; k++)
               {
                  const unsigned int a = base + 2 * k;
                  idx.insert(idx.end(), { a, a + 1, a + 2, a + 1, a + 3, a + 2 });
               }
            };
            if (!markers)
            {
               sector(0.0, R, 0.0, 2.0 * PI, zTop, SEG);
               return;
            }
            if (t.marker == "spokes")   // six wedges of 30 degrees
               for (int k = 0; k < 6; k++)
                  sector(0.0, 0.92 * R, ang + k * PI / 3.0, ang + k * PI / 3.0 + PI / 6.0, zMark, 6);
            else if (t.marker == "stripe")   // one radial bar
               sector(0.0, 0.92 * R, ang - 0.06, ang + 0.06, zMark, 2);
            else if (t.marker == "dot")      // one spot between the hub and the rim
            {
               const double cx = 0.6 * R * std::cos(ang), cy = 0.6 * R * std::sin(ang), rd2 = 0.12 * R;
               const unsigned int base = (unsigned int)v.size();
               v.emplace_back();
               FillVertex(v.back(), t.x + cx, t.y + cy, zMark);
               for (int k = 0; k <= 24; k++)
               {
                  v.emplace_back();
                  FillVertex(v.back(), t.x + cx + rd2 * std::cos(2.0 * PI * k / 24), t.y + cy + rd2 * std::sin(2.0 * PI * k / 24), zMark);
               }
               for (int k = 0; k < 24; k++)
                  idx.insert(idx.end(), { base, base + 1 + k, base + 2 + k });
            }
         };
         std::vector<Vertex3D_NoTex2> v;
         std::vector<unsigned int> idx;
         if (!t.discMesh)
         {
            build(v, idx, 0.0, false);
            t.discMesh = std::make_shared<MeshBuffer>("Turntable disc"s, std::make_shared<VertexBuffer>(rd, (unsigned int)v.size(), (const float *)v.data()), std::make_shared<IndexBuffer>(rd, idx), true);
            PLOGI << "Turntable " << t.name << " drawn: the disc only, radius " << t.radiusMm << " mm at its top z " << zTop << " VPU (the physical disc), " << SEG
                  << " segments; no gap ring (the playfield's hole, gap included, is drawn by the hole drawing)";
         }
         if (t.marker != "none")
         {
            build(v, idx, t.part->gp.Angle() + t.refDeg * PI / 180.0, true);
            if (!v.empty())
            {
               if (!t.markMesh || t.markMesh->m_vb->m_count != (unsigned int)v.size())
                  t.markMesh = std::make_shared<MeshBuffer>("Turntable marker"s, std::make_shared<VertexBuffer>(rd, (unsigned int)v.size(), (const float *)v.data(), true), std::make_shared<IndexBuffer>(rd, idx), true);
               else
               {
                  Vertex3D_NoTex2 *buf = nullptr;
                  t.markMesh->m_vb->Lock(buf);
                  memcpy(buf, v.data(), v.size() * sizeof(Vertex3D_NoTex2));
                  t.markMesh->m_vb->Unlock();
               }
            }
         }
         rd->ResetRenderState();
         rd->SetRenderState(RenderState::ZENABLE, RenderState::RS_TRUE);
         rd->SetRenderState(RenderState::ZWRITEENABLE, RenderState::RS_TRUE);
         rd->SetRenderState(RenderState::CULLMODE, RenderState::CULL_NONE);
         rd->SetRenderState(RenderState::ALPHABLENDENABLE, RenderState::RS_FALSE);
         rd->m_basicShader->SetTechnique(ShaderTechnique::unshaded_without_texture);
         auto draw = [&](const std::shared_ptr<MeshBuffer> &m, float r, float g, float b)
         {
            if (!m)
               return;
            rd->m_basicShader->SetVector(ShaderUniform::staticColor_Alpha, r, g, b, 1.f);
            rd->DrawMesh(rd->m_basicShader, false, Vertex3Ds((float)t.x, (float)t.y, (float)t.z), 0.f, m, RenderDevice::TRIANGLELIST, 0, m->m_ib->m_count);
         };
         draw(t.discMesh, 0.35f, 0.36f, 0.40f);
         if (t.marker != "none")
            draw(t.markMesh, 0.95f, 0.72f, 0.10f);
         rd->m_basicShader->SetVector(ShaderUniform::staticColor_Alpha, 1.f, 1.f, 1.f, 1.f);
      }
   }

   void ReleaseRender()
   {
      for (Turntable &t : g_tt)
      {
         t.discMesh.reset();
         t.markMesh.reset();
      }
   }
}
