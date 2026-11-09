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

      // body: the placed wire in metres, world frame (wireshape.h Place); surfaceZ: the playfield top under the wire, m
      void Build(const Settings &s, const double surfaceZ, const WireShape::WireBody &body, const RubberContact::BallSI &b)
      {
         ball = b;
         part.Build(s, surfaceZ, body, b, V3(0.0, 0.0, -9.81));
      }
   };
}
