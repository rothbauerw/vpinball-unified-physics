// license:GPLv3+

// Definition of the Kicker class

#pragma once

#include "core/resourceid.h"
#include "parts/pintable.h"
#include "physics/collide.h"
#include "physics/hitable.h"
#include "renderer/Renderable.h"
#include "utils/eventproxy.h"

class MeshBuffer;
class Texture;

class KickerData final : public BaseProperty
{
public:
   Vertex2D m_vCenter;
   float m_radius;
   string m_szSurface;
   KickerType m_kickertype;
   float m_hitAccuracy; // kicker hit grabbing object height ... default ballsize*0.7
   float m_orientation;
   float m_hit_height;
   bool m_enabled;
   bool m_fallThrough;
   bool m_legacyMode;
};

class KickerHitCircle;

class Kicker :
   public IDispatchImpl<IKicker, &IID_IKicker, &LIBID_VPinballLib>,
   //public ISupportErrorInfo,
   public CComObjectRoot,
   public CComCoClass<Kicker, &CLSID_Kicker>,
   public EventProxy<Kicker, &DIID_IKickerEvents>,
   public IConnectionPointContainerImpl<Kicker>,
   public IProvideClassInfo2Impl<&CLSID_Kicker, &DIID_IKickerEvents, &LIBID_VPinballLib>,
   public IEditable,
   public IHitable,
   public IRenderable,
   public IScriptable,
   public IFireEvents,
   public IPerPropertyBrowsing // Ability to fill in dropdown in property browser
{
public:
#ifdef VPX_MANUAL_SCRIPT_DISPATCH
   STDMETHOD(GetIDsOfNames)(REFIID /*riid*/, LPOLESTR* rgszNames, UINT cNames, LCID lcid,DISPID* rgDispId);
   STDMETHOD(Invoke)(DISPID dispIdMember, REFIID /*riid*/, LCID lcid, WORD wFlags, DISPPARAMS* pDispParams, VARIANT* pVarResult, EXCEPINFO* pExcepInfo, UINT* puArgErr);
   STDMETHOD(GetDocumentation)(MEMBERID index, BSTR *pBstrName, BSTR *pBstrDocString, DWORD *pdwHelpContext, BSTR *pBstrHelpFile);
   HRESULT FireDispID(const DISPID dispid, DISPPARAMS * const pdispparams) final;
#endif
   Kicker() { }
   virtual ~Kicker();

   BEGIN_COM_MAP(Kicker)
      COM_INTERFACE_ENTRY(IDispatch)
      COM_INTERFACE_ENTRY(IKicker)
      //COM_INTERFACE_ENTRY(ISupportErrorInfo)
      COM_INTERFACE_ENTRY_IMPL(IConnectionPointContainer)
      COM_INTERFACE_ENTRY(IPerPropertyBrowsing)
      COM_INTERFACE_ENTRY(IProvideClassInfo)
      COM_INTERFACE_ENTRY(IProvideClassInfo2)
   END_COM_MAP()
   //DECLARE_NOT_AGGREGATABLE(Kicker)
   // Remove the comment from the line above if you don't want your object to
   // support aggregation.

   BEGIN_CONNECTION_POINT_MAP(Kicker)
      CONNECTION_POINT_ENTRY(DIID_IKickerEvents)
   END_CONNECTION_POINT_MAP()

   STANDARD_EDITABLE_DECLARES(Kicker, eItemKicker, KICKER)

   DECLARE_REGISTRY_RESOURCEID(IDR_KICKER)
   // ISupportsErrorInfo
   STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);

   void Translate(const Vertex2D &offset) final;
   // Multi-object manipulation
   Vertex2D GetCenter() const final;

   void SetDefaultPhysics(const bool fromMouseClick) final;
   void ExportMesh(ObjLoader& loader) final;

   void WriteRegDefaults() final;

   KickerData m_d;

   vector<Vertex3Ds> m_hitMesh;

   // Physical kicker (unified engine, [Features] kickers = 1, [Kicker:<name>] type = physical; KickerDesign.txt): the held test (ContactSystem::NoneZones) calls these.
   void HeldChange(HitBall *ball, bool in);   // during the step: the ball's held state changed; the kicker's held ball follows it
   void HeldEvent(bool hit);                  // after the step: Hit / Unhit for the script (its ActiveBall set by the caller)
   bool IsPhysical() const { return m_physical; }
   bool IsConverted() const { return m_physical || m_physicalRefused || m_virtual; }   // on the unified kicker path (no stock kicker code; stock counter, design 7f)
   // Virtual kicker (unified engine, [Features] kickers = 1, type virtual, the default; KickerDesign.txt 1g-1i, 4b; Requirements F 10): its sensing zone (ContactSystem::NoneZones) calls these.
   bool Park(HitBall *ball);   // a ball's centre entered the zone: an empty, enabled kicker parks it at its centre at stock's hold height (true: listed, Hit); else false (not caught)
   bool Parked(const HitBall *ball) const { return m_virtual && ball != nullptr && m_heldBall == ball; }   // the parked ball: held, no Unhit until it is kicked
   void Capture(HitBall *ball);                      // Park: the capture (design 1i): placed at the centre at stock's hold height, speed and spin to 0, the energy change booked; no surface contact
   void CheckParked();                               // KickerHold::CycleEnd: the parked ball still at its place, at rest (nothing moves a held ball): a fault said once if not
   bool IsVirtual() const { return m_virtual; }
   HitBall *HeldBall() const { return m_heldBall; }

private:
   bool PhysicalSetup(PhysicsEngine *physics);   // PhysicSetup on the unified engine: false = today's kicker (kickers = 0)
   void VirtualSetup(PhysicsEngine *physics, bool fromIni);   // the virtual kicker: its sensing zone, no KickerHitCircle, no playfield hole
   bool m_virtual = false;               // this kicker is a virtual kicker on the unified path (a sensing zone and its parked, held ball)
   Vertex3Ds m_parkPlace;                // virtual: where the parked ball is held (its centre where it crossed into the zone, or where CreateBall made it), VPX units
   bool m_parkMovedNoted = false;        // virtual: CheckParked said once that something moved the parked ball
   bool m_physical = false;              // this kicker is a physical kicker (no KickerHitCircle; its cup host, carved hole and held test instead)
   bool m_physicalRefused = false;       // a physical kicker that could not be built (not on the playfield): it does nothing (no stock path, design 6a)
   HitBall *m_heldBall = nullptr;        // physical: the ball the held test holds; virtual: the parked ball (DestroyBall, Kick; null = empty)
   HitBall *m_lastHeldBall = nullptr;    // physical and virtual: LastCapturedBall
   bool m_offsetNoted = false;           // physical: "KickXYZ / KickZ offsets ignored" said once (design 4a)
   bool m_enabledNoted = false;          // physical: "Enabled is ignored" said once (Requirements F 9)
   bool m_createNoted = false;           // physical: "CreateBall refused, the kicker is occupied" said once (design 6c)
   bool m_inclNoted = false;             // physical: "inclination outside the arm's 0..90 deg" said once (design 4c)
   void GenerateMesh(Vertex3D_NoTex2 *const buf) const;

   Renderer* m_renderer = nullptr;
   std::shared_ptr<MeshBuffer> m_plateMeshBuffer;
   std::shared_ptr<MeshBuffer> m_meshBuffer;
   std::shared_ptr<MeshBuffer> m_plateMeshEdgeBuffer;
   std::shared_ptr<MeshBuffer> m_meshEdgeBuffer;
   unsigned int m_numVertices = 0;
   unsigned int m_numIndices = 0;
   std::unique_ptr<Texture> m_texture;

   float m_baseHeight = 0.f;

   KickerHitCircle *m_phitkickercircle = nullptr;

public:
   // IKicker
   STDMETHOD(get_Material)(/*[out, retval]*/ BSTR *pVal);
   STDMETHOD(put_Material)(/*[in]*/ BSTR newVal);
   STDMETHOD(get_DrawStyle)(/*[out, retval]*/ KickerType *pVal);
   STDMETHOD(put_DrawStyle)(/*[in]*/ KickerType newVal);
   STDMETHOD(get_Enabled)(/*[out, retval]*/ VARIANT_BOOL *pVal);
   STDMETHOD(put_Enabled)(/*[in]*/ VARIANT_BOOL newVal);
   STDMETHOD(get_Surface)(/*[out, retval]*/ BSTR *pVal);
   STDMETHOD(put_Surface)(/*[in]*/ BSTR newVal);
   STDMETHOD(get_Y)(/*[out, retval]*/ float *pVal);
   STDMETHOD(put_Y)(/*[in]*/ float newVal);
   STDMETHOD(get_X)(/*[out, retval]*/ float *pVal);
   STDMETHOD(put_X)(/*[in]*/ float newVal);
   STDMETHOD(Kick)(float angle, float speed, float inclination);
   STDMETHOD(KickZ)(float angle, float speed, float inclination, float heightz);
   STDMETHOD(KickXYZ)(float angle, float speed, float inclination, float x, float y, float z);
   //STDMETHOD(DestroyBall)();
   void BallGone(const HitBall *ball);   // a ball destroyed by anything but this kicker: a kicker that no longer holds a ball is empty
   STDMETHOD(DestroyBall)(/*[out, retval]*/ int *pVal);
   STDMETHOD(CreateBall)(IBall **pBall);
   STDMETHOD(CreateSizedBall)(/*[in]*/float radius, /*out, retval]*/ IBall **pBall);
   STDMETHOD(CreateSizedBallWithMass)(/*[in]*/float radius, /*[in]*/float mass, /*out, retval]*/ IBall **pBall);
   STDMETHOD(BallCntOver)(/*[out, retval]*/ int *pVal);

   STDMETHOD(get_Scatter)(/*[out, retval]*/ float *pVal);
   STDMETHOD(put_Scatter)(/*[in]*/ float newVal);
   STDMETHOD(get_HitAccuracy)(/*[out, retval]*/ float *pVal);
   STDMETHOD(put_HitAccuracy)(/*[in]*/ float newVal);
   STDMETHOD(get_HitHeight)(/*[out, retval]*/ float *pVal);
   STDMETHOD(put_HitHeight)(/*[in]*/ float newVal);
   STDMETHOD(get_Orientation)(/*[out, retval]*/ float *pVal);
   STDMETHOD(put_Orientation)(/*[in]*/ float newVal);
   STDMETHOD(get_Radius)(/*[out, retval]*/ float *pVal);
   STDMETHOD(put_Radius)(/*[in]*/ float newVal);
   STDMETHOD(get_FallThrough)(/*[out, retval]*/ VARIANT_BOOL *pVal);
   STDMETHOD(put_FallThrough)(/*[in]*/ VARIANT_BOOL newVal);
   STDMETHOD(get_Legacy)(/*[out, retval]*/ VARIANT_BOOL *pVal);
   STDMETHOD(put_Legacy)(/*[in]*/ VARIANT_BOOL newVal);
   STDMETHOD(get_LastCapturedBall)(/*[out, retval]*/ IBall **pVal);
};

class KickerHitCircle : public HitCircle
{
public:
   KickerHitCircle(Kicker * const kicker, const Vertex2D& c, const float r, const float zlow, const float zhigh)
      : HitCircle(kicker, c,r,zlow,zhigh)
   {
   }

   float HitTest(const BallS& ball, const float dtime, CollisionEvent& coll) const final;
   int GetType() const final { return eTrigger; }
   void Collide(const CollisionEvent& coll) final { DoCollide(coll.m_ball, coll.m_hitnormal, coll.m_hitflag, false); }

   void DoChangeBallVelocity(HitBall *const pball, const Vertex3Ds &hitnormal) const;
   void DoCollide(HitBall *const pball, const Vertex3Ds &hitnormal, const bool hitbit, const bool newBall);

   Kicker *m_pkicker = nullptr;
   HitBall *m_pHitBall = nullptr; // The ball inside this kicker
   HitBall *m_lastCapturedBall = nullptr; // same as m_pHitBall but this one won't be nulled only overwritten from another captured ball
};
