// license:GPLv3+
//
// A gate of the contact system: the part (gategroup.h) and its engine side (contacthingepart.h).

#pragma once

#include "physics/contacthingepart.h"
#include "physics/gategroup.h"

namespace GateModel
{
   using GateGroup::Settings;
   using GateGroup::Type;

   class Gate : public HingeModel::Engine
   {
   public:
      Gate() : HingeModel::Engine(part) {}
      GateGroup::Part part;

      void Build(const Settings &s, const V3 &surfacePoint, const V3 &front, const RubberContact::BallSI &b)
      {
         ball = b;
         part.Build(s, surfacePoint, front, b, V3(0.0, 0.0, -9.81));
      }
      bool CanPass() const { return part.CanPass(); }
      double PassAngleDeg() const { return part.PassAngleDeg(); }
      // The script's Open / Close
      void Hold(const bool open)
      {
         part.Hold(open);
         Wake();
      }
   };
}
