// license:GPLv3+

// The bat rubber's contact parameters from its settings (the settings file is read in contactbatini.cpp). See contactbat.h.

#include "contactbat.h"
#include "physics/rubbercontact.h"
#include "physics/physconst.h"


// ---------------------------------------------------------------------------
namespace FlipperModel
{
   RubberContact::RubberParams MakeRubberParams(const FlipperSettings& s)
   {
      RubberContact::RubberParams p;
      p.eMpa = s.rubberEMpa;
      p.thicknessMm = s.rubberThicknessMm;
      p.alphaSPerM = s.rubberAlphaSPerM;
      p.betaLayer = s.rubberBetaLayer;
      p.mu = s.rubberMu;
      p.shearScale = s.rubberShearScale;
      p.twist = s.rubberTwist;
      p.tanDelta = s.rubberTanDelta;
      p.shearTan = s.rubberTanDelta;   // the same material
      p.zLoMm = s.rubberZLoMm;
      p.heightMm = s.rubberHeightMm;
      return p;
   }
}
