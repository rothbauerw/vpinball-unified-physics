// license:GPLv3+

#include "core/stdafx.h"

#include "physics/carvedholes.h"

static std::vector<RubberContact::CarvedHole> g_carvedHoles;   // the turntables' and the physical kickers' holes, metres

void CarvedHoles::Add(const RubberContact::CarvedHole &h)
{
   g_carvedHoles.push_back(h);
}

void CarvedHoles::Get(const RubberContact::CarvedHole *&holes, int &n)
{
   holes = g_carvedHoles.data();
   n = (int)g_carvedHoles.size();
}

void CarvedHoles::Clear()
{
   g_carvedHoles.clear();   // (the kickers add theirs in PhysicSetup, the turntables in Turntables::Setup, both after this)
}
