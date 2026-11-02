// license:GPLv3+

#pragma once
// The table side of a gate's or spinner's moving part (HingeModel::Engine, contacthingepart.h): the host is built here, and the adapter writes the part's angle where the
// stock renderer and the script read it and turns the part's events into the table part's script events (docs/GateSpinnerDesign.txt, UnifiedContactDesign.txt 12.5).

#include <filesystem>
#include <memory>
#include <string>

#include "physics/contactparthost.h"

namespace HingeModel { class Engine; }
class Gate;
class Spinner;

class HingePartAdapter final : public PartAdapter
{
public:
   // gate or spinner: the table part (the other is null); angle, angleSpeed: the stock mover's fields the renderer and the script read (rad, rad per VPT)
   HingePartAdapter(std::shared_ptr<HingeModel::Engine> model, Gate *gate, Spinner *spinner, float *angle, float *angleSpeed, const std::string &name);
   void Publish() override;
   void Fire(const std::vector<ContactMovingPart::Event> &events) override;

   // The host of a gate or spinner of the unified engine: the part model is the surface and the body; the host engages the balls near it and its mover advances
   // the part when no ball is near. kind: "gate" or "spinner" (the ini section prefix).
   static ContactPartHost *MakeHost(IEditable *part, const char *kind, const std::string &name, const std::filesystem::path &ini, std::shared_ptr<HingeModel::Engine> model, Gate *gate,
                                    Spinner *spinner, float *angle, float *angleSpeed, double elasticity, double friction);

private:
   std::shared_ptr<HingeModel::Engine> m_model;
   Gate *m_gate;
   Spinner *m_spinner;
   float *m_angle, *m_angleSpeed;
   std::string m_name;
   long m_esSeen = 0;
   int m_esLines = 0;
   bool m_failShown = false;
};
