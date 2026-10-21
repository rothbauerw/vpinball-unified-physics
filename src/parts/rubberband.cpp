// license:GPLv3+

// Unified engine additions to Rubber: the band and sling state for scripts (Rubber.Band*, Rubber.Sling*) and the animated band mesh (ApplyBandDeformation).
// Members of Rubber (declared in rubber.h), kept apart from the stock rubber.cpp.

#include "core/stdafx.h"
#include "rubber.h"
#include "physics/hitrubber.h"


// ---- band state for scripts (rubber-contact branch) ----

#define BAND_GETTER(NAME, TYPE, EXPR) \
   STDMETHODIMP Rubber::get_Band##NAME(TYPE *pVal) \
   { \
      const BandInfo bi = m_pBand ? m_pBand->GetBandInfo() : BandInfo(); \
      *pVal = (EXPR); \
      return S_OK; \
   }

BAND_GETTER(Active, VARIANT_BOOL, FTOVB(bi.active))
BAND_GETTER(Span, int, bi.span)
BAND_GETTER(Deflection, float, bi.deflectionMm)
BAND_GETTER(PeakDeflection, float, bi.peakMm)
BAND_GETTER(Position, float, bi.posSpan)
BAND_GETTER(RubberPosition, float, bi.posRubber)
BAND_GETTER(NodeX, float, bi.nodeX)
BAND_GETTER(NodeY, float, bi.nodeY)
BAND_GETTER(NodeZ, float, bi.nodeZ)
BAND_GETTER(HitTime, int, bi.hitTimeMs)
BAND_GETTER(SpeedIn, float, bi.speedIn)
BAND_GETTER(SpeedOut, float, bi.speedOut)
BAND_GETTER(HitSpan, int, bi.hitSpan)
BAND_GETTER(HitPosition, float, bi.hitPosSpan)
BAND_GETTER(HitRubberPosition, float, bi.hitPosRubber)
#undef BAND_GETTER

#define SLING_GETTER(NAME, TYPE, EXPR) \
   STDMETHODIMP Rubber::get_Sling##NAME(TYPE *pVal) \
   { \
      const SlingInfo si = m_pBand ? m_pBand->GetSlingInfo() : SlingInfo(); \
      *pVal = (EXPR); \
      return S_OK; \
   }

SLING_GETTER(Enabled, VARIANT_BOOL, FTOVB(si.enabled))
SLING_GETTER(PlateAngle, float, si.plateAngleDeg)
SLING_GETTER(SwitchClosed, VARIANT_BOOL, FTOVB(si.switchClosed))
SLING_GETTER(CoilOn, VARIANT_BOOL, FTOVB(si.coilOn))
SLING_GETTER(CoilCurrent, float, si.coilCurrentA)
SLING_GETTER(PulseCount, int, si.pulses)
SLING_GETTER(FireTime, int, si.fireTimeMs)
#undef SLING_GETTER

STDMETHODIMP Rubber::get_BandSwitchCount(int *pVal)
{
   *pVal = m_pBand ? m_pBand->BandSwitchCount() : 0;
   return S_OK;
}

STDMETHODIMP Rubber::get_BandSwitchClosed(int index, VARIANT_BOOL *pVal)
{
   *pVal = FTOVB(m_pBand != nullptr && m_pBand->BandSwitchClosedAt(index));
   return S_OK;
}

STDMETHODIMP Rubber::put_SlingEnabled(VARIANT_BOOL newVal)
{
   if (m_pBand)
      m_pBand->SetSlingEnabled(VBTOb(newVal));
   return S_OK;
}

STDMETHODIMP Rubber::SlingTrigger()
{
   if (m_pBand)
      m_pBand->SlingTrigger();
   return S_OK;
}

// The rings of the mesh follow the deflected band (HitRubber::BandDeformAt): each ring is shifted and turned about its centre, normals turn with it.
void Rubber::ApplyBandDeformation(Vertex3D_NoTex2 *buf)
{
   const int S = m_numSegments;
   const bool animating = m_pBand != nullptr && m_pBand->BandAnimating();
   bool any = false;
   if (S > 0 && animating)
   {
      const int rings = m_numVertices / S;
      for (int ring = 0; ring < rings; ring++)
      {
         Vertex3D_NoTex2 *const v = buf + ring * S;
         Vertex3Ds centre(0.f, 0.f, 0.f);
         for (int j = 0; j < S; j++)
            centre += Vertex3Ds(v[j].x, v[j].y, v[j].z);
         centre /= (float)S;
         Vertex3Ds shift;
         float R[9];
         if (!m_pBand->BandDeformAt(centre, shift, R))
            continue;
         any = true;
         m_bandPeakShift = std::max(m_bandPeakShift, shift.Length());
         static const char *const meshLog = getenv("VPX_BAND_MESH_LOG");   // diagnostics: the ring centres before and after the deformation (table units)
         static FILE *const mf = (meshLog && *meshLog) ? fopen(meshLog, "w") : nullptr;
         if (mf)
         {
            fprintf(mf, "%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f%c", m_bandFrames, ring, centre.x, centre.y, centre.z, centre.x + shift.x, centre.y + shift.y, centre.z + shift.z, 10);
            fflush(mf);
         }
         for (int j = 0; j < S; j++)
         {
            const float px = v[j].x - centre.x, py = v[j].y - centre.y, pz = v[j].z - centre.z;
            v[j].x = centre.x + shift.x + R[0] * px + R[1] * py + R[2] * pz;
            v[j].y = centre.y + shift.y + R[3] * px + R[4] * py + R[5] * pz;
            v[j].z = centre.z + shift.z + R[6] * px + R[7] * py + R[8] * pz;
            const float nx = v[j].nx, ny = v[j].ny, nz = v[j].nz;
            v[j].nx = R[0] * nx + R[1] * ny + R[2] * nz;
            v[j].ny = R[3] * nx + R[4] * ny + R[5] * nz;
            v[j].nz = R[6] * nx + R[7] * ny + R[8] * nz;
         }
      }
   }
   if (any)
      m_bandFrames++;
   if (any != m_bandDeformed)
   {
      if (any)
      {
         m_bandPeakShift = 0.f;
         m_bandFrames = 0;
      }
      else
         PLOGI << "Rubber " << GetName() << ": band mesh animation ended after " << m_bandFrames << " frames, largest ring shift " << m_bandPeakShift * 0.53975f << " mm";
   }
   m_bandDeformed = any;
}
