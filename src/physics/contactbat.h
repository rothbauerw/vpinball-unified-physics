#pragma once
// The bat's settings file: what is read from the table's physics ini and what is derived from it. The bat itself, its state and its step are BatGroup::Part (batgroup.h).

#include <filesystem>
#include <string>

#include "physics/coilmodel.h"
#include "physics/rubbercontact.h"

// Settings helpers shared with the sling plate
namespace FlipperModel
{
   RubberContact::RubberParams MakeRubberParams(const FlipperSettings &s);   // the bat rubber's contact parameters
}

class ContactBat
{
public:
   // ---- settings ----
   FlipperSettings set;
   FlipperDerived der;
   std::filesystem::path iniPath;
   // (Re)builds the settings from the ini file (the table's vpx_physics.ini; none = the defaults): at table load, and again when the file changed (ContactSystem::ReloadChangedIni,
   // on playfield focus gain, through the flipper's PartAdapter::ReloadSettings). In contactbatini.cpp (VPX side).
   void LoadSettings(const std::string &who, const std::filesystem::path &tableFile);
};
