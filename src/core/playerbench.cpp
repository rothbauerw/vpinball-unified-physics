// license:GPLv3+

// The real bench: VPX itself with no frames drawn, the physics stepped as fast as the CPU goes (test harness only, VPX_BENCH=1).
// One definition of the physics: this file only replaces the game loop. Each pass makes the same calls in the same order as
// UpdateGameLogic + PrepareFrame in a 144 Hz game, without input and without Render3DScene.
//   VPX_BENCH_PREROLL  physics ms before the test driver's time 0 (default 3000): the table's parts settle, as in the start-up of a VPX run
//   VPX_BENCH_MS       stop at this engine time in ms, for a run without an input file (default: run until the driver's end)
// One run per process: it ends here, with no teardown (the render thread waits for a frame that never comes; the VPX regression covers teardown).

#include "core/stdafx.h"
#include "core/player.h"
#include "core/ieditable.h"
#include "parts/pintable.h"
#include "physics/PhysicsEngine.h"
#include "physics/testdriver.h"
#include "renderer/Renderer.h"
#include "renderer/RenderDevice.h"

bool Player::BenchMode()
{
   static const bool on = getenv("VPX_BENCH") != nullptr;
   return on;
}

void Player::BenchLoop()
{
   constexpr uint32_t frameMs = 7; // virtual frame, as 144 Hz. Only what VPX does per frame sees it: animations, frame timers
   const auto envMs = [](const char *const name, const uint32_t def) { const char *const e = getenv(name); return e != nullptr && *e ? (uint32_t)atoi(e) : def; };
   const uint32_t prerollMs = envMs("VPX_BENCH_PREROLL", 3000), stopMs = envMs("VPX_BENCH_MS", 0);
   PLOGI << "Bench: no frames, pre-roll " << prerollMs << " ms, graphics backend " << m_renderer->m_renderDevice->m_driver_name;
   while (!TestDriver::Ended() && (stopMs == 0 || m_time_msec < stopMs))
   {
      m_logicProfiler.NewFrame(m_time_msec);
      if (m_time_msec >= prerollMs)
         TestDriver::FrameStart();
      m_logicProfiler.EnterProfileSection(FrameProfiler::PROFILE_PREPARE_FRAME);
      m_overall_frames++;
      const float diff_time_msec = (float)(m_time_msec - m_last_frame_time_msec);
      m_last_frame_time_msec = m_time_msec;
      if (diff_time_msec > 0.f)
         for (IEditable *editable : m_ptable->GetParts())
            if (auto ph = editable->GetIRenderable(); ph)
               ph->UpdateAnimation(diff_time_msec);
      FireTimers(-1);
      m_physics->ResetPerFrameStats();
      m_logicProfiler.ExitProfileSection();
      m_physics->UpdatePhysics(m_physics->GetCurrentTime() + frameMs * PHYSICS_STEPTIME + 1); // the +1: UpdatePhysics runs the steps strictly before its target
      FireTimers(-2);
      m_pluginManager.ProcessAsyncCallbacks();
   }
   PLOGI << "Bench: end at engine time " << m_time_msec << " ms";
#ifdef _MSC_VER
   TerminateProcess(GetCurrentProcess(), 0);
#else
   _Exit(0);
#endif
}
