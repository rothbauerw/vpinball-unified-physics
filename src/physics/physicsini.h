#pragma once
// Unified physics settings file: <table>.physics.ini next to the table (testing fallback: vpx_physics.ini in the folder), or the file named by VPX_PHYSICS_INI. [Features] unified = 1 switches the table to the unified engine; no file or no key = stock VPX physics.
// A unified file is recognised by its [Features] section. Sections:
//   [Features]        on/off switch for each feature
//   [Contact]         resting-contact friction of the ball (table wide)
//   [Flippers]        every flipper;            [Flipper:<name>]  one flipper
//   [Rubbers]         rings and bands;          [Rubber:<name>]   one rubber
//   [Wall:<name>], [Primitive:<name>]   a wall / primitive flagged  rubber = 1  uses the rubber force model; other keys override [Rubbers] for it
//   [Surfaces]       unified engine: defaults for every wall and primitive that is not a rubber (hard surfaces: the same keys as [Rubbers]); [Wall:<name>] / [Primitive:<name>] override, rubber = 0 opts out; [Features] hard_surfaces = 0 switches it off
//   [Turntable:<name>]   a motor-driven disc set into the playfield; the script drives it with Turntable "<name>", volts   (see turntable.h)
//   [Magnet:<name>]   an electromagnet (position, coil, core); the script drives it with Magnet "<name>", volts   (see magnet.h)
//   [Slings]         every sling (keys without the sling_ prefix);  [Sling:<rubber name>]  one sling (plate, switches, overrides)
// With unified = 1 every model is on by default (the file may name no other key); a feature is switched off with  <key> = 0  in [Features]. On the stock engine (no unified = 1)
// every model is off, whatever the file says (Active below): no mixing of the new models with the stock engine.
// vpx_physics.ini is the only settings file (the old flippers.ini, rubbers.ini, contact.ini and their VPX_*_INI variables were removed 2026-10-20).

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>

namespace ContactSystem
{
   bool Unified();   // contactsystem.cpp
}

namespace PhysicsIni
{
   struct Features
   {
      bool unified = false;       // master switch: the table runs on the unified engine (off = stock VPX physics); read by ContactSystem::SetTable
      bool flippers = true;       // read so a table that sets it to 0 is told it is ignored (a bat has no stock path, FlipperMoverObject::MakePart)
      bool flipperRubber = true;  // the same: the bat has no impulse path
      bool rubbers = true;        // force model for ring rubbers (posts, wall rubbers); off = stock collisions
      bool playfieldContact = true; // the playfield (plane) on the force model: one material ([Playfield]); off = the stock plane contact
      bool surfaceRubbers = true; // walls and primitives flagged rubber = 1 in [Wall:<name>] / [Primitive:<name>]; off = stock collisions
      bool hardSurfaces = true;   // unified engine: every other wall and primitive is a hard surface of the contact model (material from the part's elasticity and friction, [Surfaces] and [Wall:<name>] / [Primitive:<name>] override); off = stock collisions
      bool proximityEngage = true; // unified engine: a free ball is taken into the contact model when it is about to touch a surface (reach 2 v dt + the release gap), before the stock hit search can find it; the part hit events come from the contact
      bool triggers = true;       // unified engine: wire triggers (WireA-D) are rollover wire switches of the contact system (design-wire.md); off = stock sensors. On by default (user 2026-11-04)
      bool hinges = true;         // unified engine: gates and spinners are hinged plates of the contact model (mass, torque from the balls touching them); off = the stock kick of the plate by the ball
      bool bumpers = true;        // read only so bumpers = 0 can be told as ignored (bumper.cpp): a unified table has no other bumper
      bool bands = true;          // band model (rubber stretched between posts); off = stock collisions
      bool slings = true;         // slings (plate, coil, switches) on bands; off = a plain band
      bool contact = true;        // resting-contact friction; off = plain Coulomb friction with the object's friction
      float hitSpeedMin = 0.5f;   // [Features] hit_speed_min: impacts slower than this (ball velocity units, 18.53 = 1 m/s) do not fire the playfield Hit event
      bool magnets = true;        // electromagnets ([Magnet:<name>] sections); off = none
      bool turntables = true;     // turntables ([Turntable:<name>] sections); off = none
      bool sharedFinder = true;   // shared_finder (on while tested, rule 6; to-do 1.34, runs/2026-11-11-finder/design.md): 1 = the shared finder's index over every static mesh's faces (a box tree per mesh object over its faces' boxes, faces joined through shared edges into mesh pieces, one box tree over the objects) is built once before the first cycle; in this step (P1) only check mode (VPX_FINDER_CHECK=1) reads it, so no result changes; 0 = no index
      bool playfieldCuts = true;  // playfield_cuts: holes cut in the playfield ([PlayfieldCut:<name>] sections, playfieldcut.h; the slot under a wire switch); off = every cut ignored, the playfield solid (on while tested, rule 6)
      bool freeFlight = true;     // free_flight (on while tested, rule 6; runs/2026-11-12-freeflight/log.md): 1 = every ball the contact system does not hold is still stepped by it each cycle, as a body with no contacts (its gravity, nudge and magnets in the substeps), so no ball is ever moved by the stock flight (pos += time x end velocity); 0 = such a ball is moved by stock as before
      bool kickers = true;        // kickers (on while tested, rule 6; KickerDesign.txt): 1 = a kicker of type physical ([Kicker:<name>] type; default virtual) is a cup of the
                                  // contact model under a hole carved in the playfield, with Hit / Unhit from its held test; a virtual kicker is a sensing zone that parks the ball, a held ball of infinite mass (kickerarm.h); 0 = every kicker on today's path, unchanged
      bool hitSpeed = true;       // PROTOTYPE: Ball.HitSpeed (normal speed before an impact) and the table's Hit event for the playfield; off = neither
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
      static const char *const plain[] = { "features", "contact", "flippers", "rubbers", "slings", "playfield", "surfaces", "gates", "spinners", "bumpers" };
      static const char *const kinds[] = { "flipper:", "rubber:", "sling:", "wall:", "primitive:", "magnet:", "turntable:", "playfieldcut:", "gate:", "spinner:", "ramp:", "bumper:", "kicker:" };
      for (const char *p : plain)
         if (lowerName == p)
            return;
      for (const char *k : kinds)
         if (lowerName.rfind(k, 0) == 0)
            return;
      static std::set<std::string> warned;
      if (warned.insert(lowerName).second)
         PLOGW << "Physics settings: section [" << lowerName << "] is not a known kind (features, contact, flippers, rubbers, slings, flipper:, rubber:, sling:, wall:, primitive:); it is ignored";
   }

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
         else if (key == "flippers")
            f.flippers = on;
         else if (key == "flipper_rubber")
            f.flipperRubber = on;
         else if (key == "rubbers")
            f.rubbers = on;
         else if (key == "playfield_contact")
            f.playfieldContact = on;
         else if (key == "surface_rubbers")
            f.surfaceRubbers = on;
         else if (key == "hard_surfaces")
            f.hardSurfaces = on;
         else if (key == "proximity_engage")
            f.proximityEngage = on;
         else if (key == "hinges")
            f.hinges = on;
         else if (key == "triggers")
            f.triggers = on;
         else if (key == "bumpers")
            f.bumpers = on;
         else if (key == "bands")
            f.bands = on;
         else if (key == "slings")
            f.slings = on;
         else if (key == "contact")
            f.contact = on;
         else if (key == "magnets")
            f.magnets = on;
         else if (key == "turntables")
            f.turntables = on;
         else if (key == "playfield_cuts")
            f.playfieldCuts = on;
         else if (key == "shared_finder")
            f.sharedFinder = on;
         else if (key == "free_flight")
            f.freeFlight = on;
         else if (key == "hit_speed")
            f.hitSpeed = on;
         else if (key == "kickers")
            f.kickers = on;
         else if (key == "hit_speed_min")
            f.hitSpeedMin = (float)atof(v.c_str());
         else
            PLOGW << "Physics settings: [Features] ignored line '" << line << "'";
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
         f.flippers = f.flipperRubber = f.rubbers = f.playfieldContact = f.surfaceRubbers = f.hardSurfaces = f.proximityEngage = f.hinges = f.triggers = false;
         f.bands = f.slings = f.contact = f.magnets = f.turntables = f.playfieldCuts = f.hitSpeed = f.sharedFinder = f.freeFlight = f.kickers = false;
         return f;
      }
      if (!file.empty())
         ReadFeatures(file, f);
      return f;
   }

}
