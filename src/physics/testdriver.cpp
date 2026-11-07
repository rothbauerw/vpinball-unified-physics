// license:GPLv3+

#include "core/stdafx.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include "physics/testdriver.h"
#include "parts/ball.h"
#include "parts/flipper.h"
#include "parts/rubber.h"
#include "parts/pintable.h"
#include "core/player.h"
#include "ui/live/LiveUI.h"
#include "physics/magnet.h"
#include "physics/turntable.h"
#include "physics/stockcount.h"
#include "physics/contactsystem.h"

// Test driver (docs/unified-physics/TestMatrixDesign.txt section 2; keys and file format in docs/unified-physics/archive/NewPhysicsReference.txt). Harness only: off unless VPX_TEST_INPUT is set.
// Input file: one action per line, "<t_ms> <action> <args>", '#' starts a comment; actions with the same time run in file order. Time 0 is the first cycle after the
// table's start-up (section 3 point 10g: the first frame that took less than half the anti-hung limit), so the start-up's slow frames never skip driver time.
//    ball NAME x y z vx vy vz [radius mass]   create a ball: x, y in VPX units, z = height of the ball's bottom above the playfield (U), velocity in m/s
//    destroy NAME                             remove the ball
//    flip PART 1|0                            flipper solenoid on / off (as the script's RotateToEnd / RotateToStart)
//    enabled PART 1|0                         a flipper's Enabled (as the script's Flipper.Enabled)
//    sling PART                               fire the sling coil of a rubber band (as Rubber.SlingTrigger)
//    magnet NAME V                            supply voltage of a [Magnet:NAME] in the physics ini (as the script's Magnet NAME, V; 0 = off)
//    turntable NAME V                         motor voltage of a [Turntable:NAME] in the physics ini (as the script's Turntable NAME, V; 0 = off)
//    respawn Y                                from now on a ball whose centre passes y = Y (U) toward the drain (y > Y) is removed and created again at its own start state
//    respawn_z Z                              from now on a ball whose centre drops below z = Z (U) is removed and created again at its own start state (off by default;
//                                             for tables with pockets or subways, where a ball that is not kicked out would stay)
//    speed_limit V                            launch limit in m/s (default 8, TestMatrixDesign 1 limit 3)
//    floor Z                                  a ball centre below z = Z (U) is through the playfield (default -200: pockets and subways go below 0)
//    record N                                 record the balls every N cycles (default 1; 0 = off, for cost runs)
//    timing detail|total                      cycle timing with every split (default) or the cycle total only (section 3 point 10h: the cost of the splits themselves)
//    end                                      write report.txt and done.txt, stop
namespace TestDriver
{
   namespace
   {
      constexpr double kMps = 18.53;   // VPX velocity units (U per 10 ms) per m/s
      constexpr uint32_t kStartFrameUs = 100000;   // start-up is over at the first frame shorter than half the anti-hung limit (200 ms in PhysicsEngine::UpdatePhysics)

      struct Action
      {
         uint32_t t;
         int line;
         string verb;
         vector<string> args;
      };

      struct TBall
      {
         string name;
         Ball *ball = nullptr;   // null when removed
         float x, y, z, vx, vy, vz, radius, mass;   // start state (VPX units; z = centre)
         float wx = 0.f, wy = 0.f, wz = 0.f;        // start spin, rad/s
         bool respawnDue = false;
      };

      struct Count
      {
         uint64_t n = 0;
         uint32_t firstT = 0;
         string firstBall;
         double worst = 0.0;
         uint32_t worstT = 0;
         string worstBall;
      };

      struct CycleTime
      {
         uint32_t t;
         float cycle, script;   // us; cycle = the whole loop pass without the script time inside it
         float part[kParts];    // us, each without the script inside it
         ContactCost cc;        // the contact system's split (ticks) and its slowest island
      };

      struct FrameTime
      {
         uint32_t t, cycles;
         uint32_t physics, script, prepare, misc, sleep, frame;   // us, from VPX's frame profiler (each section without the ones nested in it)
      };

      int64_t NowNs() { return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

      struct State
      {
         bool loaded = false, on = false, ended = false, started = false, startupDone = false;
         string input, outDir;
         int runIndex = 0;
         vector<Action> actions;
         size_t next = 0;
         vector<TBall> balls;
         uint32_t t0 = 0, lastAbs = 0, lastT = 0;   // t0 = engine time of the driver's time 0; lastT = the driver's time of this cycle
         uint64_t cycles = 0, skips = 0, respawns = 0, respawnsZ = 0, recordEvery = 1;
         double speedLimit = 8.0, floorZ = -200.0, respawnY = 0.0, respawnZ = 0.0;
         bool respawnOn = false, respawnZOn = false;
         Count fast, below, outside, gone;
         vector<string> errors;
         FILE *csv = nullptr;
         // timing (TestMatrixDesign section 3 points 9 and 10)
         bool detail = true;
         vector<CycleTime> cycleTimes;
         vector<FrameTime> frameTimes;
         bool loopOpen = false;
         int64_t loopT0 = 0;
         uint32_t script0 = 0, cyclesInFrame = 0;
         float partUs[kParts] = {};
         ContactCost cc;
         uint64_t tick0 = 0;   // the tick counter against the clock at the start: ticks per microsecond at the end
         int64_t ns0 = 0;
      } s;

      void Warn(const string &msg)
      {
         PLOGW << "Test driver: " << msg;
         if (g_pplayer && g_pplayer->m_liveUI)
            g_pplayer->m_liveUI->PushNotification("Test driver: " + msg, 8000);
         s.errors.push_back(msg);
      }

      void Note(Count &c, uint32_t t, const string &ball, double value)
      {
         if (c.n++ == 0)
         {
            c.firstT = t;
            c.firstBall = ball;
         }
         if (value > c.worst)
         {
            c.worst = value;
            c.worstT = t;
            c.worstBall = ball;
         }
      }

      string OutPath(const string &file)
      {
         // a second table run in the same process writes its own files (balls_2.csv, ...), so nothing is overwritten
         if (s.runIndex <= 1)
            return s.outDir + file;
         const size_t dot = file.rfind('.');
         return s.outDir + file.substr(0, dot) + '_' + std::to_string(s.runIndex) + file.substr(dot);
      }

      void Load()
      {
         s.loaded = true;
         const char *const in = getenv("VPX_TEST_INPUT");
         if (in == nullptr || *in == 0)
            return;
         s.input = in;
         const char *const out = getenv("VPX_TEST_OUT");
         s.outDir = out != nullptr ? out : "";
         if (!s.outDir.empty() && s.outDir.back() != '\\' && s.outDir.back() != '/')
            s.outDir += '\\';
         s.on = true;
      }

      void ReadInput()
      {
         s.actions.clear();
         std::ifstream f(s.input);
         if (!f)
         {
            Warn("input file not found: " + s.input);
            return;
         }
         string line;
         int n = 0;
         while (std::getline(f, line))
         {
            n++;
            const size_t hash = line.find('#');
            if (hash != string::npos)
               line.resize(hash);
            std::istringstream ls(line);
            Action a;
            a.line = n;
            if (!(ls >> a.t))
               continue;   // empty or comment line
            if (!(ls >> a.verb))
            {
               Warn("line " + std::to_string(n) + ": no action");
               continue;
            }
            string w;
            while (ls >> w)
               a.args.push_back(w);
            s.actions.push_back(std::move(a));
         }
         std::stable_sort(s.actions.begin(), s.actions.end(), [](const Action &a, const Action &b) { return a.t < b.t; });
      }

      void Start(uint32_t tMs)
      {
         s.runIndex++;
         s.started = true;
         s.ended = false;
         s.next = 0;
         s.balls.clear();   // a new table: the old ball pointers are gone with the old player
         s.t0 = tMs;
         s.lastT = 0;
         s.cycles = s.skips = s.respawns = s.respawnsZ = 0;
         s.recordEvery = 1;
         s.speedLimit = 8.0;
         s.floorZ = -200.0;
         s.respawnOn = s.respawnZOn = false;
         s.fast = s.below = s.outside = s.gone = Count();
         s.errors.clear();
         s.detail = true;
         s.cycleTimes.clear();
         s.cycleTimes.reserve(65536);
         s.frameTimes.clear();
         s.frameTimes.reserve(16384);
         s.cyclesInFrame = 0;
         s.tick0 = Ticks();
         s.ns0 = NowNs();
         if (s.csv != nullptr)
            fclose(s.csv);
         ReadInput();
         // the random generator also serves sounds, lights and stock scatter: start it from its fixed seed (core/def.cpp) so the table load does not shift it
         mwc64x_state = 4077358422479273989ull;
         s.csv = fopen(OutPath("balls.csv").c_str(), "w");
         if (s.csv != nullptr)
            fputs("t_ms,ball,x,y,z,vx,vy,vz,wx,wy,wz\n", s.csv);   // VPX units: position U, velocity U per 10 ms, angular momentum as stored
         else
            Warn("cannot write " + OutPath("balls.csv"));
         PLOGI << "Test driver: " << s.actions.size() << " actions from " << s.input << ", output " << s.outDir << ", time 0 = engine time " << tMs << " ms (start-up done)";
      }

      bool Alive(const Ball *b)
      {
         return b != nullptr && std::find(g_pplayer->m_vball.begin(), g_pplayer->m_vball.end(), b) != g_pplayer->m_vball.end();
      }

      TBall *FindBall(const string &name)
      {
         for (TBall &b : s.balls)
            if (b.name == name)
               return &b;
         return nullptr;
      }

      IEditable *FindPart(const string &name, ItemTypeEnum type)
      {
         for (IEditable *const e : g_pplayer->m_ptable->GetParts())
            if (e->GetItemType() == type && e->GetName() == name)
               return e;
         return nullptr;
      }

      void Create(TBall &b)
      {
         b.ball = g_pplayer->CreateBall(b.x, b.y, b.z - b.radius, b.vx, b.vy, b.vz, b.radius, b.mass);   // CreateBall adds the radius to z
         if (b.ball != nullptr && (b.wx != 0.f || b.wy != 0.f || b.wz != 0.f))
            b.ball->m_hitBall.m_angularmomentum = Vertex3Ds(b.wx, b.wy, b.wz) * (b.ball->m_hitBall.Inertia() * (float)DEFAULT_STEPTIME_S);
         b.respawnDue = false;
      }

      void Remove(TBall &b)
      {
         if (Alive(b.ball))
            g_pplayer->DestroyBall(b.ball);
         b.ball = nullptr;
      }

      bool Num(const Action &a, size_t i, float &v)
      {
         if (i >= a.args.size())
            return false;
         char *end = nullptr;
         v = strtof(a.args[i].c_str(), &end);
         return end != nullptr && *end == 0;
      }

      void Report(uint32_t t);

      void Apply(const Action &a, uint32_t t)
      {
         const string at = "line " + std::to_string(a.line) + " (" + a.verb + "): ";
         if (a.verb == "ball")
         {
            float v[11] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, DEFAULT_BALL_SIZE, 1.f, 0.f, 0.f, 0.f };   // DEFAULT_BALL_SIZE is the radius (25 U), as Kicker.CreateBall uses it
            const size_t n = a.args.size();
            bool ok = n == 7 || n == 9 || n == 12;
            for (size_t i = 1; ok && i < n; i++)
               ok = Num(a, i, v[i - 1]);
            if (!ok)
            {
               Warn(at + "expected NAME x y z vx vy vz [radius mass [wx wy wz]]");
               return;
            }
            TBall *b = FindBall(a.args[0]);
            if (b == nullptr)
            {
               s.balls.push_back(TBall());
               b = &s.balls.back();
               b->name = a.args[0];
            }
            else
               Remove(*b);
            b->x = v[0];
            b->y = v[1];
            b->radius = v[6];
            b->z = v[2] + b->radius;
            b->vx = (float)(v[3] * kMps);
            b->vy = (float)(v[4] * kMps);
            b->vz = (float)(v[5] * kMps);
            b->mass = v[7];
            b->wx = v[8];
            b->wy = v[9];
            b->wz = v[10];
            Create(*b);
         }
         else if (a.verb == "roll")
         {
            TBall *const r = a.args.size() == 1 ? FindBall(a.args[0]) : nullptr;
            if (r != nullptr && Alive(r->ball))
            {
               HitBall &hb = r->ball->m_hitBall;
               const float L = hb.Inertia() * (float)DEFAULT_STEPTIME_S / (hb.m_d.m_radius * 0.53975e-3f) / (float)kMps;   // (spin = z x v / R; angular momentum = spin x the engine's inertia x its time unit; v in U per 10 ms over kMps)
               hb.m_angularmomentum = Vertex3Ds(-hb.m_d.m_vel.y * L, hb.m_d.m_vel.x * L, 0.f);
            }
            else
               Warn(at + "no such ball");
         }
         else if (a.verb == "destroy")
         {
            TBall *const b = a.args.size() == 1 ? FindBall(a.args[0]) : nullptr;
            if (b == nullptr)
               Warn(at + "no such ball");
            else
               Remove(*b);
         }
         else if (a.verb == "flip")
         {
            Flipper *const f = a.args.size() == 2 ? static_cast<Flipper *>(FindPart(a.args[0], eItemFlipper)) : nullptr;
            if (f == nullptr || (a.args[1] != "1" && a.args[1] != "0"))
               Warn(at + "expected a flipper name and 1 or 0");
            else if (a.args[1] == "1")
               f->RotateToEnd();
            else
               f->RotateToStart();
         }
         else if (a.verb == "enabled")
         {
            Flipper *const f = a.args.size() == 2 ? static_cast<Flipper *>(FindPart(a.args[0], eItemFlipper)) : nullptr;
            if (f == nullptr || (a.args[1] != "1" && a.args[1] != "0"))
               Warn(at + "expected a flipper name and 1 or 0");
            else
               f->put_Enabled(a.args[1] == "1" ? VARIANT_TRUE : VARIANT_FALSE);
         }
         else if (a.verb == "sling")
         {
            Rubber *const r = a.args.size() == 1 ? static_cast<Rubber *>(FindPart(a.args[0], eItemRubber)) : nullptr;
            if (r == nullptr)
               Warn(at + "no such rubber");
            else
               r->SlingTrigger();
         }
         else if (a.verb == "magnet" || a.verb == "turntable")
         {
            float v;
            if (a.args.size() != 2 || !Num(a, 1, v))
               Warn(at + "expected a name and a voltage");
            else if (!(a.verb == "magnet" ? Magnets::SetVoltage(a.args[0].c_str(), v) : Turntables::SetVoltage(a.args[0].c_str(), v)))
               Warn(at + "no such " + a.verb + " in the physics ini");
         }
         else if (a.verb == "timing")
         {
            if (a.args.size() != 1 || (a.args[0] != "detail" && a.args[0] != "total"))
               Warn(at + "expected detail or total");
            else
               s.detail = a.args[0] == "detail";
         }
         else if (a.verb == "respawn" || a.verb == "respawn_z" || a.verb == "speed_limit" || a.verb == "floor" || a.verb == "record")
         {
            float v;
            if (a.args.size() != 1 || !Num(a, 0, v))
            {
               Warn(at + "expected one number");
               return;
            }
            if (a.verb == "respawn")
            {
               s.respawnY = v;
               s.respawnOn = true;
            }
            else if (a.verb == "respawn_z")
            {
               s.respawnZ = v;
               s.respawnZOn = true;
            }
            else if (a.verb == "speed_limit")
               s.speedLimit = v;
            else if (a.verb == "floor")
               s.floorZ = v;
            else
               s.recordEvery = v < 0.5f ? 0 : (uint64_t)std::lround(v);
         }
         else if (a.verb == "end")
            Report(t);
         else
            Warn(at + "unknown action");
      }

      void Line(std::ostream &o, const char *what, const Count &c, const char *unit)
      {
         o << what << ": " << (c.n == 0 ? "PASS" : "FAIL") << ", " << c.n << " ball-cycles";
         if (c.n != 0)
            o << ", first at t=" << c.firstT << " ms ball " << c.firstBall << ", worst " << c.worst << ' ' << unit << " at t=" << c.worstT << " ms ball " << c.worstBall;
         o << '\n';
      }

      // ---- cost (TestMatrixDesign section 3 points 9 and 10) ----

      // The columns of the cost split, in report and csv order. Values are us per cycle, each without the script run inside it.
      enum Col
      {
         cCycle, cStep, cContact, cSurfaces, cPairs, cBats, cBands, cMovingParts, cContactShared, cStockSearch, cMoveParts, cStepOther,
         cMagnets, cTurntable, cInput, cTimers, cDriver, cNudge, cVelocities, cWatchdog, cGuides, cTrails, cOther, cScript, cCols
      };
      const char *const kColName[cCols] = { "cycle", "step", "contact", "surfaces", "pairs", "bats", "bands", "moving_parts", "contact_shared", "stock_search", "move_parts",
         "step_other", "magnets", "turntable", "input", "timers", "driver", "nudge", "velocities", "watchdog", "guides", "trails", "other", "script" };
      const char *const kColLabel[cCols] = { "cycle", "physics step", "  contact system", "    surfaces (playfield, walls, posts, rubbers)", "    ball-ball pairs",
         "    flippers", "    bands and slings", "    gates, spinners, bumpers", "    shared (islands, step choice, integration, write-back)", "  stock hit search",
         "  movers' displacement", "  rest of the step (stock collisions, events)", "magnets", "turntable", "input", "script timers (without the script)",
         "test driver", "nudge and plunger", "movers' velocities", "fast-ball watchdog", "magnet guide tubes", "ball trails", "other (loop overhead)", "script" };
      const char *const kWhyName[ContactCost::kWhys] = { "end of the cycle", "impact arriving", "surface stiffness", "ball-ball stiffness", "moving part's step", "ini substep_us",
         "flipper step", "band step", "smallest step allowed", "taken again (rewind)" };

      double s_ticksPerUs = 1.0;

      void Cols(const CycleTime &c, double v[cCols])
      {
         const double k = 1.0 / s_ticksPerUs;
         double kinds = 0.0;
         v[cCycle] = c.cycle;
         v[cStep] = c.part[kStep];
         v[cContact] = c.part[kContact];
         v[cSurfaces] = c.cc.kindTicks[ContactCost::kSurfaces] * k;
         v[cPairs] = c.cc.kindTicks[ContactCost::kPairs] * k;
         v[cBats] = c.cc.kindTicks[ContactCost::kBats] * k;
         v[cBands] = c.cc.kindTicks[ContactCost::kBands] * k;
         v[cMovingParts] = c.cc.kindTicks[ContactCost::kMovingParts] * k;
         for (int i = cSurfaces; i <= cMovingParts; i++)
            kinds += v[i];
         v[cContactShared] = v[cContact] - kinds;
         v[cStockSearch] = c.part[kStockSearch];
         v[cMoveParts] = c.part[kMoveParts];
         v[cStepOther] = v[cStep] - v[cContact] - v[cStockSearch] - v[cMoveParts];
         v[cMagnets] = c.part[kMagnets];
         v[cTurntable] = c.part[kTurntable];
         v[cInput] = c.part[kInput];
         v[cTimers] = c.part[kTimers];
         v[cDriver] = c.part[kDriver];
         v[cNudge] = c.part[kNudge];
         v[cVelocities] = c.part[kVelocities];
         v[cWatchdog] = c.part[kWatchdog];
         v[cGuides] = c.part[kGuides];
         v[cTrails] = c.part[kTrails];
         double top = v[cStep] + v[cMagnets] + v[cTurntable];
         for (int i = cInput; i <= cTrails; i++)
            top += v[i];
         v[cOther] = v[cCycle] - top;
         v[cScript] = c.script;
      }

      void WriteTimes()
      {
         if (FILE *const f = fopen(OutPath("cycles.csv").c_str(), "w"))
         {
            fputs("t_ms", f);
            for (int i = 0; i < cCols; i++)
               fprintf(f, ",%s_us", kColName[i]);
            fputs(",islands,top_steps,top_us,top_min_step_us,why_end,why_arrival,why_surface_stiffness,why_pair_stiffness,why_part_step,why_setting,why_bat,why_band,why_minimum,why_retake,top_members\n", f);
            double v[cCols];
            for (const CycleTime &c : s.cycleTimes)
            {
               Cols(c, v);
               fprintf(f, "%u", c.t);
               for (int i = 0; i < cCols; i++)
                  fprintf(f, ",%.2f", v[i]);
               fprintf(f, ",%d,%d,%.2f,%.3f", c.cc.islands, c.cc.topSteps, c.cc.topTicks / s_ticksPerUs, c.cc.topMinStepUs);
               for (int i = 0; i < ContactCost::kWhys; i++)
                  fprintf(f, ",%d", c.cc.topWhy[i]);
               fprintf(f, ",\"%s\"\n", c.cc.topMembers.c_str());
            }
            fclose(f);
         }
         if (FILE *const f = fopen(OutPath("frames.csv").c_str(), "w"))
         {
            fputs("t_ms,cycles,physics_us,script_us,prepare_us,misc_us,sleep_us,frame_us,work_us\n", f);
            for (const FrameTime &r : s.frameTimes)
               fprintf(f, "%u,%u,%u,%u,%u,%u,%u,%u,%u\n", r.t, r.cycles, r.physics, r.script, r.prepare, r.misc, r.sleep, r.frame, r.physics + r.script + r.prepare);
            fclose(f);
         }
      }

      // The split of a run of cycles (a window, or one cycle) and its slowest island
      void Split(std::ostream &o, const vector<const CycleTime *> &c, size_t first, size_t n)
      {
         double sum[cCols] = {}, v[cCols];
         int why[ContactCost::kWhys] = {};
         const CycleTime *top = c[first];
         for (size_t k = first; k < first + n; k++)
         {
            Cols(*c[k], v);
            for (int i = 0; i < cCols; i++)
               sum[i] += v[i];
            for (int i = 0; i < ContactCost::kWhys; i++)
               why[i] += c[k]->cc.topWhy[i];
            if (c[k]->cc.topSteps > top->cc.topSteps)
               top = c[k];
         }
         for (int i = 0; i < cCols; i++)
            if (i == cCycle || std::fabs(sum[i]) >= 0.5)
               o << "         " << kColLabel[i] << ": " << sum[i] << '\n';
         if (top->cc.topSteps > 0)
         {
            o << "         slowest island (t=" << top->t << " ms): " << top->cc.topSteps << " substeps, " << top->cc.topTicks / s_ticksPerUs << " us, smallest step "
              << top->cc.topMinStepUs << " us; " << top->cc.topMembers << '\n';
            o << "         what set the substeps of each cycle's slowest island:";
            for (int i = 0; i < ContactCost::kWhys; i++)
               if (why[i] != 0)
                  o << ' ' << kWhyName[i] << ' ' << why[i] << ';';
            o << '\n';
         }
      }

      // Cost part of the report. Returns false only when a focused run is over limit 6.
      bool Cost(std::ostream &o)
      {
         constexpr size_t kWindow = 7;
         constexpr double kWindowLimit = 3000.0, kCycleLimit = 1500.0, kFrameWork = 6900.0;
         s_ticksPerUs = std::max(1.0, (double)(Ticks() - s.tick0) / std::max(1.0, (NowNs() - s.ns0) * 1e-3));
         vector<const CycleTime *> c;
         for (const CycleTime &x : s.cycleTimes)
            c.push_back(&x);
         o << "\nCost (TestMatrixDesign section 3 points 9 and 10): every cycle from time 0 (after start-up); a cycle = the whole 1 ms loop pass without its script time.\n"
           << "Timing: " << (s.detail ? "detail (every split)" : "total only (point 10h)") << ". Each part leaves out the script run inside it.\n";
         if (c.size() < kWindow)
         {
            o << "   not enough cycles timed (" << c.size() << ")\n";
            return true;
         }
         // mean, 99th percentile and maximum of every column
         {
            vector<vector<float>> col(cCols);
            double v[cCols];
            for (const CycleTime *x : c)
            {
               Cols(*x, v);
               for (int i = 0; i < cCols; i++)
                  col[i].push_back((float)v[i]);
            }
            o << "   per cycle (us): mean / 99th percentile / max\n";
            for (int i = 0; i < cCols; i++)
            {
               vector<float> &a = col[i];
               double sum = 0.0;
               for (const float x : a)
                  sum += x;
               std::sort(a.begin(), a.end());
               if (i != cCycle && a.back() < 0.5f && sum / a.size() < 0.05)
                  continue;
               o << "      " << kColLabel[i] << ": " << sum / a.size() << " / " << a[std::min(a.size() - 1, (size_t)std::ceil(0.99 * a.size()) - 1)] << " / " << a.back() << '\n';
            }
         }
         const CycleTime *worst = c[0];
         size_t worstI = 0;
         for (size_t i = 0; i < c.size(); i++)
            if (c[i]->cycle > worst->cycle)
            {
               worst = c[i];
               worstI = i;
            }
         vector<std::pair<double, size_t>> win;   // sum of 7 consecutive cycles, first index
         double w = 0.0;
         for (size_t i = 0; i < c.size(); i++)
         {
            w += c[i]->cycle;
            if (i >= kWindow)
               w -= c[i - kWindow]->cycle;
            if (i + 1 >= kWindow)
               win.emplace_back(w, i + 1 - kWindow);
         }
         std::sort(win.begin(), win.end(), [](const auto &a, const auto &b) { return a.first > b.first; });
         const bool focused = !Player::PlayUnfocused();
         const bool pass = win[0].first <= kWindowLimit && worst->cycle <= kCycleLimit;
         o << "   worst single cycle: " << worst->cycle << " us at t=" << worst->t << " ms (limit " << kCycleLimit << ")\n";
         Split(o, c, worstI, 1);
         o << "   worst 7 consecutive cycles: " << win[0].first << " us from t=" << c[win[0].second]->t << " ms (limit " << kWindowLimit << ")\n";
         o << "   5 worst 7-cycle windows, split (us, summed over the 7 cycles; parts under 0.5 us left out):\n";
         vector<size_t> taken;
         for (const auto &[ws, i] : win)
         {
            if (taken.size() == 5)
               break;
            if (std::any_of(taken.begin(), taken.end(), [i = i](size_t j) { return (i > j ? i - j : j - i) < kWindow; }))
               continue;
            taken.push_back(i);
            o << "      from t=" << c[i]->t << " ms:\n";
            Split(o, c, i, kWindow);
         }
         size_t nf = 0, over = 0, cyclesSum = 0;
         uint32_t cyclesMax = 0, tFirst = 0, tLast = 0;
         for (size_t k = 1; k < s.frameTimes.size(); k++)   // (the first frame began before time 0)
         {
            const FrameTime &f = s.frameTimes[k];
            if (nf++ == 0)
               tFirst = f.t;
            tLast = f.t;
            cyclesSum += f.cycles;
            cyclesMax = std::max(cyclesMax, f.cycles);
            over += (f.physics + f.script + f.prepare > kFrameWork) ? 1 : 0;
         }
         if (nf > 1)
            o << "   frames: " << nf << ", rate " << (nf - 1) * 1000.0 / std::max(1u, tLast - tFirst) << " per s, cycles per frame mean " << (double)cyclesSum / nf << " max "
              << cyclesMax << "; frames with work (physics + script + preparation) over " << kFrameWork << " us: " << over << '\n';
         else
            o << "   frames: too few timed (" << nf << ")\n";
         o << "   limit 6: " << (focused ? (pass ? "PASS" : "FAIL") : "not judged (unfocused run: VPX_PLAY_UNFOCUSED)") << '\n';
         return pass || !focused;
      }

      void Report(uint32_t t)
      {
         if (s.ended)
            return;
         s.ended = true;
         if (s.csv != nullptr)
         {
            fclose(s.csv);
            s.csv = nullptr;
         }
         std::ostringstream cost;
         const bool costPass = Cost(cost);   // (sets the tick rate the csv uses)
         WriteTimes();
         const bool stockPass = !ContactSystem::Unified() || !StockCount::ConvertedUsed();
         const bool pass = s.skips == 0 && s.fast.n == 0 && s.below.n == 0 && s.outside.n == 0 && s.gone.n == 0 && s.errors.empty() && costPass && stockPass;
         std::ofstream o(OutPath("report.txt"));
         o << "Test driver report (docs/unified-physics/TestMatrixDesign.txt)\n" << (Player::BenchMode() ? "BENCH RUN: every cost below is bench cost (no frames, no GPU, virtual frames of 7 ms); limit 6 is judged on a focused VPX run\n" : "") << "input: " << s.input << "\ncycles: " << s.cycles << " (t 0 to " << t << " ms; time 0 = engine time " << s.t0
           << " ms, after start-up), balls: " << s.balls.size() << ", respawns: " << s.respawns;
         if (s.respawnZOn)
            o << " (of them " << s.respawnsZ << " for a ball centre below z = " << s.respawnZ << " U)";
         o << "\n\n";
         o << "time skips (cycles not simulated; repeatability): " << (s.skips == 0 ? "PASS" : "FAIL") << ", " << s.skips << "\n";
         Line(o, "limit 3 launch: speed over the limit", s.fast, "m/s");
         o << "   (limit " << s.speedLimit << " m/s)\n";
         Line(o, "limit 4 ball below the playfield", s.below, "U below the floor");
         o << "   (floor z = " << s.floorZ << " U, ball centre)\n";
         Line(o, "limit 4 ball outside the table", s.outside, "U outside");
         Line(o, "ball removed by the table, not by the driver", s.gone, "");
         o << "input errors: " << (s.errors.empty() ? "none" : "FAIL") << '\n';
         for (const string &e : s.errors)
            o << "   " << e << '\n';
         o << cost.str();
         if (ContactSystem::Unified())
            StockCount::Report(o, OutPath("stock.csv"));
         o << "\nNOT CHECKED by the driver: limit 1 energy (check not built, to-do 1.6), limit 2 overlap depth, limit 4 ball inside a solid, limit 5 safety nets\n"
              "(look for WARN lines in vpinball.log), limit 6 per object (the report gives the time per kind of object, not per object).\n";
         o << "\nRESULT of the checks above: " << (pass ? "PASS" : "FAIL") << '\n';
         o.close();
         std::ofstream(OutPath("done.txt")) << (pass ? "PASS" : "FAIL") << '\n';
         PLOGI << "Test driver: end at t=" << t << " ms, " << (pass ? "PASS" : "FAIL") << ", report " << OutPath("report.txt");
      }
   }

   void BeginCycle(uint32_t tMs)
   {
      if (!s.loaded)
         Load();
      if (!s.on || g_pplayer == nullptr)
         return;
      if (s.runIndex > 0 && tMs < s.lastAbs)   // the table was restarted (time starts again): wait for its start-up again
         s.started = s.startupDone = false;
      s.lastAbs = tMs;
      s.cyclesInFrame++;
      if (!s.started)
      {
         if (!s.startupDone)
            return;
         Start(tMs);
      }
      const uint32_t t = tMs - s.t0;
      if (t != 0 && t != s.lastT + 1 && !s.ended)
         s.skips += t - s.lastT - 1;
      s.lastT = t;
      if (s.ended)
         return;
      for (TBall &b : s.balls)
         if (b.respawnDue)
         {
            Remove(b);
            Create(b);
            s.respawns++;
         }
      while (s.next < s.actions.size() && s.actions[s.next].t <= t && !s.ended)
         Apply(s.actions[s.next++], t);
   }

   void EndCycle(uint32_t tMs)
   {
      if (!s.on || !s.started || s.ended || g_pplayer == nullptr)
         return;
      const uint32_t t = tMs - s.t0;
      s.cycles++;
      const PinTable &tb = *g_pplayer->m_ptable;
      const bool rec = s.csv != nullptr && s.recordEvery != 0 && s.cycles % s.recordEvery == 0;
      for (TBall &b : s.balls)
      {
         if (b.ball == nullptr)
            continue;
         if (!Alive(b.ball))
         {
            Note(s.gone, t, b.name, 0.0);
            b.ball = nullptr;
            continue;
         }
         const BallS &d = b.ball->m_hitBall.m_d;
         const Vertex3Ds &w = b.ball->m_hitBall.m_angularmomentum;
         if (rec)
            fprintf(s.csv, "%u,%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n", t, b.name.c_str(), d.m_pos.x, d.m_pos.y, d.m_pos.z, d.m_vel.x, d.m_vel.y, d.m_vel.z, w.x, w.y, w.z);
         const double speed = std::sqrt((double)d.m_vel.x * d.m_vel.x + (double)d.m_vel.y * d.m_vel.y + (double)d.m_vel.z * d.m_vel.z) / kMps;
         if (speed > s.speedLimit)
            Note(s.fast, t, b.name, speed);
         if (d.m_pos.z < s.floorZ)
            Note(s.below, t, b.name, s.floorZ - d.m_pos.z);
         const double out = std::max({ (double)tb.m_left - d.m_pos.x, (double)d.m_pos.x - tb.m_right, (double)tb.m_top - d.m_pos.y, (double)d.m_pos.y - tb.m_bottom });
         if (out > 0.0)
            Note(s.outside, t, b.name, out);
         if (s.respawnZOn && d.m_pos.z < s.respawnZ && !b.respawnDue)
         {
            b.respawnDue = true;
            s.respawnsZ++;
         }
         if (s.respawnOn && d.m_pos.y > s.respawnY)
            b.respawnDue = true;   // taken at the start of the next cycle, so this cycle's checks and record see where the ball went
      }
   }

   void LoopStart()
   {
      s.loopOpen = s.on && s.started && !s.ended && g_pplayer != nullptr;   // a pass is timed from the cycle after the driver's start
      if (!s.loopOpen)
         return;
      s.loopT0 = NowNs();
      s.script0 = g_pplayer->m_logicProfiler.Get(FrameProfiler::PROFILE_SCRIPT);   // script time so far in this frame; the loop pass never spans two frames
      std::fill(std::begin(s.partUs), std::end(s.partUs), 0.f);
      s.cc = ContactCost();
   }

   void LoopEnd()
   {
      if (!s.loopOpen)
         return;
      s.loopOpen = false;
      if (!s.on || s.ended)
         return;
      const float total = (float)((NowNs() - s.loopT0) * 1e-3);
      const float script = (float)(g_pplayer->m_logicProfiler.Get(FrameProfiler::PROFILE_SCRIPT) - s.script0);
      CycleTime c;
      c.t = s.lastT;
      c.cycle = total - script;
      c.script = script;
      std::copy(std::begin(s.partUs), std::end(s.partUs), c.part);
      c.cc = std::move(s.cc);
      s.cycleTimes.push_back(std::move(c));
   }

   void FrameStart()
   {
      if (!s.on || g_pplayer == nullptr)
         return;
      const FrameProfiler &p = g_pplayer->m_logicProfiler;
      if (!s.started)
      {
         if (s.cyclesInFrame > 0 && p.GetPrev(FrameProfiler::PROFILE_FRAME) < kStartFrameUs)
            s.startupDone = true;   // physics ran in a frame of normal length: the start-up is over
         s.cyclesInFrame = 0;
         return;
      }
      if (s.ended)
         return;
      s.frameTimes.push_back({ g_pplayer->m_time_msec - s.t0, s.cyclesInFrame, p.GetPrev(FrameProfiler::PROFILE_PHYSICS), p.GetPrev(FrameProfiler::PROFILE_SCRIPT),
         p.GetPrev(FrameProfiler::PROFILE_PREPARE_FRAME), p.GetPrev(FrameProfiler::PROFILE_MISC), p.GetPrev(FrameProfiler::PROFILE_SLEEP), p.GetPrev(FrameProfiler::PROFILE_FRAME) });
      s.cyclesInFrame = 0;
   }

   bool Ended() { return s.ended; }

   ContactCost *Contact() { return s.loopOpen && s.detail ? &s.cc : nullptr; }

   uint32_t ReportT() { return s.on ? s.lastT : (g_pplayer != nullptr ? g_pplayer->m_time_msec : 0); }

   const char *BallName(const void *hitBall)
   {
      for (const TBall &b : s.balls)
         if (b.ball != nullptr && &b.ball->m_hitBall == hitBall)
            return b.name.c_str();
      return nullptr;
   }

   Timed::Timed(Part p) : part(p), t0(0), script0(0)
   {
      if (!s.loopOpen || !s.detail)
         return;
      t0 = NowNs();
      script0 = g_pplayer->m_logicProfiler.Get(FrameProfiler::PROFILE_SCRIPT);
   }

   Timed::~Timed()
   {
      if (t0 == 0 || !s.loopOpen)
         return;
      const float script = (float)(g_pplayer->m_logicProfiler.Get(FrameProfiler::PROFILE_SCRIPT) - script0);
      s.partUs[part] += (float)((NowNs() - t0) * 1e-3) - script;
   }
}
