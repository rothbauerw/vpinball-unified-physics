#pragma once
// Unified physics settings file: <table>.physics.ini next to the table (testing fallback: vpx_physics.ini in the folder), or the file named by VPX_PHYSICS_INI. [Features] unified = 1 switches the table to the unified engine; no file or no key = stock VPX physics.
// A unified file is recognised by its [Features] section. Sections:
//   [Features]        on/off switch for each feature
//   [Contact]         resting-contact friction of the ball (table wide)
//   [Flippers]        every flipper;            [Flipper:<name>]  one flipper
//   [Rubbers]         rings and bands;          [Rubber:<name>]   one rubber
//   [Wall:<name>], [Primitive:<name>]   a wall / primitive flagged  rubber = 1  uses the rubber force model; other keys override [Rubbers] for it
//   [Surfaces]       unified engine: defaults for every wall and primitive that is not a rubber (hard surfaces: the same keys as [Rubbers]); [Wall:<name>] / [Primitive:<name>] override, rubber = 1 makes one a rubber surface
//   [Turntable:<name>]   a motor-driven disc set into the playfield; the script drives it with Turntable "<name>", volts   (see turntable.h)
//   [Magnet:<name>]   an electromagnet (position, coil, core); the script drives it with Magnet "<name>", volts   (see magnet.h)
//   [Slings]         every sling (keys without the sling_ prefix);  [Sling:<rubber name>]  one sling (plate, switches, overrides)
// With unified = 1 every model is on by default (the file may name no other key); a feature is switched off with  <key> = 0  in [Features] (a retired key = 0 is ignored and told: RETIRED below). On the stock engine (no unified = 1)
// every model is off, whatever the file says (Active below): no mixing of the new models with the stock engine.
// vpx_physics.ini is the only settings file (the old flippers.ini, rubbers.ini, contact.ini and their VPX_*_INI variables were removed 2026-10-20).

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <vector>

namespace ContactSystem
{
   bool Unified();   // contactsystem.cpp
   void NoticeRetired(const std::string &msg);   // contactsystem.cpp
}

namespace PhysicsIni
{
   struct Features
   {
      bool unified = false;       // master switch: the table runs on the unified engine (off = stock VPX physics); read by ContactSystem::SetTable
      bool surfaceRubbers = true; // walls and primitives flagged rubber = 1 in [Wall:<name>] / [Primitive:<name>]; off = they are hard surfaces of the contact model like the rest
      bool slings = true;         // slings (plate, coil, switches) on bands; off = a plain band
      bool contact = true;        // resting-contact friction; off = plain Coulomb friction with the object's friction
      float hitSpeedMin = 0.5f;   // [Features] hit_speed_min: impacts slower than this (ball velocity units, 18.53 = 1 m/s) do not fire the playfield Hit event
      bool magnets = true;        // electromagnets ([Magnet:<name>] sections); off = none
      bool turntables = true;     // turntables ([Turntable:<name>] sections); off = none
      bool meshOverlapRepair = true;   // mesh_overlap_repair (on while tested, rule 6; to-do N3 piece 5): 1 = at load, in every static collidable mesh (playfield, wall, ramp, primitive), faces of one
                                       // host on one plane (normals parallel either way, plane distance <= 1 um) that overlap cover the surface once: a face wholly covered by smaller ones is dropped,
                                       // one partly covered is clipped by them (meshholes.h RepairOverlaps; logged per host, never on screen); 0 = the faces as the table has them
      bool aircoreInStep = true;  // aircore_in_step (on while tested, rule 6; runs/2026-11-13-aircore-ff/design.md): 1 = an air-core coil's pull on a ball in the contact system's step is read at the
                                  // step's midpoint in every evaluation (Magnets::AircoreForce, its slope in Newton's Jacobian); 0 = held at the cycle-start place for the whole cycle (Ball::fExt), as today
      bool seatEditorBalls = true; // seat_editor_balls (on while tested, rule 6; runs/2026-11-13-editor-balls): 1 = at play start each ball placed in the table editor is moved out of any surface it overlaps
                                  // (playfield, walls, ramps, primitives, rubbers, moving parts at rest; not other balls) by the shortest move, the no-past seating (ContactSystem::PushOutAtStart); 0 = left where it was placed, as today
      std::vector<std::string> retiredOff;   // the retired keys (RETIRED) this file sets to 0, in file order: told once per table start by ContactSystem::SetTable (WARN and a notice), never acted on
      bool hitSpeed = true;       // PROTOTYPE: Ball.HitSpeed (normal speed before an impact) and the table's Hit event for the playfield; off = neither
      float holeRoundMm = 0.2f;   // hole_round_mm: a round hole cut in the mesh (Button, kicker, turntable, [PlayfieldCut] circle) is a polygon whose edges touch the true circle at their midpoints;
                                  // each vertex stands at most this far outside it, mm (MeshHoles::Round: N = max(24, smallest N with r (1 / cos(pi / N) - 1) <= this); smaller = more edges, a rounder
                                  // rim). TUNED by feel (user 2026-11-17; 0.2 mm from paper section 4); 0.01 .. 1 (outside: a notice and the nearest)
   };

   inline std::string Trim(const std::string &s)
   {
      const size_t a = s.find_first_not_of(" \t\r\n");
      if (a == std::string::npos)
         return std::string();
      return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
   }

   inline std::string Lower(std::string s)
   {
      for (char &c : s)
         c = (char)tolower((unsigned char)c);
      return s;
   }

   // the unified file for this table, or an empty path if there is none
   inline std::filesystem::path FindUnified(const std::filesystem::path &tableFile)
   {
      if (const char *const e = getenv("VPX_PHYSICS_INI"))
         return e;   // (the test harness)
      if (tableFile.empty())
         return std::filesystem::path();
      std::error_code ec;
      std::filesystem::path p = tableFile;   // <table>.physics.ini next to the table; vpx_physics.ini (shared by the folder) is a testing fallback
      p.replace_extension(".physics.ini");
      if (std::filesystem::exists(p, ec))
         return p;
      p = tableFile.parent_path() / "vpx_physics.ini";
      return std::filesystem::exists(p, ec) ? p : std::filesystem::path();
   }

   // A section header nothing reads (a misspelt kind such as [Primitve:Post1]) is reported once.
   inline void WarnUnknownSection(const std::string &lowerName)
   {
      // every section a loader reads (PartsIni::ReadSections' plain and kind names, ContactHost::ApplySection's, the magnet/turntable/playfield-cut parsers)
      static const char *const plain[] = { "features", "contact", "flippers", "rubbers", "slings", "playfield", "surfaces", "gates", "spinners", "bumpers", "triggers", "plungers", "targets", "tableedge", "glass" };
      static const char *const kinds[] = { "flipper:", "rubber:", "sling:", "wall:", "primitive:", "magnet:", "turntable:", "playfieldcut:", "gate:", "spinner:", "ramp:", "bumper:", "kicker:", "trigger:", "plunger:", "target:" };
      for (const char *p : plain)
         if (lowerName == p)
            return;
      for (const char *k : kinds)
         if (lowerName.rfind(k, 0) == 0)
            return;
      static std::set<std::string> warned;
      if (warned.insert(lowerName).second)
         PLOGW << "Physics settings: section [" << lowerName << "] is not a known kind (features, contact, flippers, rubbers, slings, playfield, surfaces, gates, spinners, bumpers, triggers, plungers, targets, tableedge, glass, or a part's kind: flipper, rubber, sling, wall, primitive, magnet, turntable, playfieldcut, gate, spinner, ramp, bumper, kicker, trigger, plunger, target); it is ignored";
   }

   // [Features] keys whose 0 once sent an object or a law to stock, or asked for a path a part does not have (flippers, flipper_rubber, bumpers). On the unified engine each object
   // always takes its unified branch (DevelopmentRules 6: unified = 0 runs stock only, 1 bypasses stock entirely). Kept only so a file that still says  <key> = 0  is told it is ignored (rule 8); = 1 is silent
   inline constexpr const char *RETIRED[] = { "hard_surfaces", "playfield_contact", "table_bounds", "triggers", "hinges", "plunger", "rubbers", "bands", "ball_gravity_in_step",
                                              "flippers", "flipper_rubber", "bumpers" };

   // true if the file has a [Features] section (it is a unified file); `f` gets its switches (a missing key = on)
   inline bool ReadFeatures(const std::filesystem::path &file, Features &f)
   {
      std::ifstream in(file);
      if (!in)
         return false;
      std::string line, cur;
      bool found = false;
      while (std::getline(in, line))
      {
         const size_t c = line.find_first_of("#;");
         if (c != std::string::npos)
            line.resize(c);
         line = Trim(line);
         if (line.empty())
            continue;
         if (line.front() == '[' && line.back() == ']')
         {
            cur = Lower(Trim(line.substr(1, line.size() - 2)));
            found = found || cur == "features";
            WarnUnknownSection(cur);
            continue;
         }
         if (cur != "features")
            continue;
         const size_t eq = line.find('=');
         if (eq == std::string::npos)
            continue;
         const std::string key = Lower(Trim(line.substr(0, eq))), v = Lower(Trim(line.substr(eq + 1)));
         const bool on = (v == "1" || v == "true" || v == "yes" || v == "on");
         if (key == "unified")
            f.unified = on;
         else if (std::find(std::begin(RETIRED), std::end(RETIRED), key) != std::end(RETIRED))
         {
            if (!on && std::find(f.retiredOff.begin(), f.retiredOff.end(), key) == f.retiredOff.end())
               f.retiredOff.push_back(key);   // (= 1 is silent: it is what runs)
         }
         else if (key == "surface_rubbers")
            f.surfaceRubbers = on;
         else if (key == "slings")
            f.slings = on;
         else if (key == "contact")
            f.contact = on;
         else if (key == "magnets")
            f.magnets = on;
         else if (key == "turntables")
            f.turntables = on;
         else if (key == "mesh_overlap_repair")
            f.meshOverlapRepair = on;
         else if (key == "aircore_in_step")
            f.aircoreInStep = on;
         else if (key == "hit_speed")
            f.hitSpeed = on;
         else if (key == "seat_editor_balls")
            f.seatEditorBalls = on;
         else if (key == "hit_speed_min")
            f.hitSpeedMin = (float)atof(v.c_str());   // (a key not listed here is ignored)
         else if (key == "hole_round_mm")
            f.holeRoundMm = (float)atof(v.c_str());
      }
      return found;
   }

   // The features in force for a table: on the stock engine (no unified = 1) all off; on the unified engine every one on unless its [Features] key is 0 in the
   // unified file (file may be empty or missing: all on).
   inline Features Active(const std::filesystem::path &file)
   {
      Features f;
      if (!ContactSystem::Unified())
      {
         f.surfaceRubbers = f.slings = f.contact = f.magnets = f.turntables = f.meshOverlapRepair = f.hitSpeed = f.aircoreInStep = f.seatEditorBalls = false;
         return f;
      }
      if (!file.empty())
         ReadFeatures(file, f);
      return f;
   }

}
