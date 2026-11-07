// license:GPLv3+
#pragma once
// The stock counter (docs/unified-physics/NewPhysicsReference.txt section 11): every stock collision and stock contact, by object, in every driver run. PhysicsEngine wraps each dispatch of a stock Collide or
// Contact in Begin / End; what counts is what happened to the ball (its position, velocity or spin changed), plus the two stock laws by name (Collide3DWall, HandleStaticContact reached)
// and a stock-only part's Collide. So a route nobody listed is still counted. A converted object that uses stock fails the run; the only exceptions are the unconverted parts in AllowedStock.
// Reads the balls, never writes them. Costs a hash look-up and a copy of nine floats per stock dispatch (eighteen for a ball pair), nothing otherwise.
// Not counted: what is not a dispatch: free flight, the hit search, a script's or a kicker's move of a ball.

#include <cstdint>
#include <ostream>
#include <string>

class HitObject;
class HitBall;

namespace StockCount
{
   void SetPlanes(const HitObject *playfield, const HitObject *topGlass);
   void Begin(HitObject *obj, HitBall *ball, bool contact);
   void End(bool ballAlive);   // ballAlive false: a script run by the collision removed the ball
   void Wall();                // HitBall::Collide3DWall reached
   void Static();              // HitBall::HandleStaticContact reached
   bool Off();                 // VPX_CHECKS=0: off (a diagnostic, for its cost)
   bool ConvertedUsed();
   uint32_t ReportTime();      // the time a report gives: the driver's (balls.csv), else the engine's (testdriver.cpp)
   void Report(std::ostream &o, const std::string &csvPath);
}
