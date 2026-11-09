// license:GPLv3+
//
// A rollover wire switch of the contact system: the part (wireswitchgroup.h) and its engine side (contacthingepart.h), as a gate. The table side is in trigger.cpp.

#pragma once

#include "physics/contacthingepart.h"
#include "physics/wireswitchgroup.h"

namespace WireSwitchModel
{
   using WireSwitchGroup::Settings;

   class WireSwitch : public HingeModel::Engine
   {
   public:
      WireSwitch() : HingeModel::Engine(part) {}
      WireSwitchGroup::Part part;

      // chain: the traced centreline in metres, world frame (wireshape.h Build); surfaceZ: the playfield under the wire, m
      void Build(const Settings &s, const double surfaceZ, const std::vector<WireShape::Capsule> &chain, const RubberContact::BallSI &b)
      {
         ball = b;
         part.Build(s, surfaceZ, chain, b, V3(0.0, 0.0, -9.81));
      }
   };
}
