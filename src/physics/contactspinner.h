// license:GPLv3+
//
// A spinner of the contact system: the part (spinnergroup.h) and its engine side (contacthingepart.h).

#pragma once

#include "physics/contacthingepart.h"
#include "physics/spinnergroup.h"

namespace SpinnerModel
{
   using SpinnerGroup::Settings;
   using SpinnerGroup::Switch;

   class Spinner : public HingeModel::Engine
   {
   public:
      Spinner() : HingeModel::Engine(part) {}
      SpinnerGroup::Part part;

      void Build(const Settings &s, const V3 &surfacePoint, const V3 &front, const RubberContact::BallSI &b)
      {
         ball = b;
         part.Build(s, surfacePoint, front, b, V3(0.0, 0.0, -9.81));
      }
   };
}
