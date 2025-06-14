// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Event> : std::false_type { };
  template<> struct DefaultConstructible<G4Event> : std::false_type { };
}

// Class generating the wrapper for type G4Event
// signature to use in the veto file: G4Event
struct JlG4Event: public Wrapper {

  JlG4Event(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Event (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:58:7
    jlcxx::TypeWrapper<G4Event>  t = jlModule.add_type<G4Event>("G4Event");
    jlcxx::stl::apply_stl<G4Event*>(jlModule);
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Event>>(new jlcxx::TypeWrapper<G4Event>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes    );


    DEBUG_MSG("Adding wrapper for void G4Event::G4Event(G4int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:62:14
    t.constructor<G4int>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("evID")    );
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4Event::operator==(const G4Event &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Event::operator==(const G4Event &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:71:12
    t.method("==", [](G4Event const& a, const G4Event & arg0)->G4bool { return a.operator==(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));
    t.method("==", [](G4Event const* a, const G4Event & arg0)->G4bool { return a->operator==(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));

    DEBUG_MSG("Adding wrapper for G4bool G4Event::operator!=(const G4Event &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Event::operator!=(const G4Event &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:72:12
    t.method("!=", [](G4Event const& a, const G4Event & arg0)->G4bool { return a.operator!=(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));
    t.method("!=", [](G4Event const* a, const G4Event & arg0)->G4bool { return a->operator!=(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for void G4Event::Print() (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::Print()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:74:10
    t.method("Print", [](G4Event const& a)->void { a.Print(); }, jlcxx::arg("this"));
    t.method("Print", [](G4Event const* a)->void { a->Print(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::Draw() (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::Draw()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:76:10
    t.method("Draw", [](G4Event const& a)->void { a.Draw(); }, jlcxx::arg("this"));
    t.method("Draw", [](G4Event const* a)->void { a->Draw(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::SetEventID(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::SetEventID(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:81:17
    t.method("SetEventID", [](G4Event& a, G4int arg0)->void { a.SetEventID(arg0); }, jlcxx::arg("this"), jlcxx::arg("i"));
    t.method("SetEventID", [](G4Event* a, G4int arg0)->void { a->SetEventID(arg0); }, jlcxx::arg("this"), jlcxx::arg("i"));

    DEBUG_MSG("Adding wrapper for void G4Event::SetHCofThisEvent(G4HCofThisEvent *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::SetHCofThisEvent(G4HCofThisEvent *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:83:17
    t.method("SetHCofThisEvent", [](G4Event& a, G4HCofThisEvent * arg0)->void { a.SetHCofThisEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("value"));
    t.method("SetHCofThisEvent", [](G4Event* a, G4HCofThisEvent * arg0)->void { a->SetHCofThisEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("value"));

    DEBUG_MSG("Adding wrapper for void G4Event::SetDCofThisEvent(G4DCofThisEvent *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::SetDCofThisEvent(G4DCofThisEvent *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:85:17
    t.method("SetDCofThisEvent", [](G4Event& a, G4DCofThisEvent * arg0)->void { a.SetDCofThisEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("value"));
    t.method("SetDCofThisEvent", [](G4Event* a, G4DCofThisEvent * arg0)->void { a->SetDCofThisEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("value"));

    DEBUG_MSG("Adding wrapper for void G4Event::SetTrajectoryContainer(G4TrajectoryContainer *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::SetTrajectoryContainer(G4TrajectoryContainer *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:87:17
    t.method("SetTrajectoryContainer", [](G4Event& a, G4TrajectoryContainer * arg0)->void { a.SetTrajectoryContainer(arg0); }, jlcxx::arg("this"), jlcxx::arg("value"));
    t.method("SetTrajectoryContainer", [](G4Event* a, G4TrajectoryContainer * arg0)->void { a->SetTrajectoryContainer(arg0); }, jlcxx::arg("this"), jlcxx::arg("value"));

    DEBUG_MSG("Adding wrapper for void G4Event::SetEventAborted() (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::SetEventAborted()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:89:17
    t.method("SetEventAborted", [](G4Event& a)->void { a.SetEventAborted(); }, jlcxx::arg("this"));
    t.method("SetEventAborted", [](G4Event* a)->void { a->SetEventAborted(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::SetRandomNumberStatus(G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::SetRandomNumberStatus(G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:91:17
    t.method("SetRandomNumberStatus", [](G4Event& a, G4String & arg0)->void { a.SetRandomNumberStatus(arg0); }, jlcxx::arg("this"), jlcxx::arg("st"));
    t.method("SetRandomNumberStatus", [](G4Event* a, G4String & arg0)->void { a->SetRandomNumberStatus(arg0); }, jlcxx::arg("this"), jlcxx::arg("st"));

    DEBUG_MSG("Adding wrapper for void G4Event::SetRandomNumberStatusForProcessing(G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::SetRandomNumberStatusForProcessing(G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:96:17
    t.method("SetRandomNumberStatusForProcessing", [](G4Event& a, G4String & arg0)->void { a.SetRandomNumberStatusForProcessing(arg0); }, jlcxx::arg("this"), jlcxx::arg("st"));
    t.method("SetRandomNumberStatusForProcessing", [](G4Event* a, G4String & arg0)->void { a->SetRandomNumberStatusForProcessing(arg0); }, jlcxx::arg("this"), jlcxx::arg("st"));

    DEBUG_MSG("Adding wrapper for void G4Event::KeepTheEvent(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::KeepTheEvent(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:101:17
    t.method("KeepTheEvent", [](G4Event const& a)->void { a.KeepTheEvent(); }, jlcxx::arg("this"), jlcxx::arg("vl"));
    t.method("KeepTheEvent", [](G4Event const& a, G4bool arg0)->void { a.KeepTheEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("vl"));
    t.method("KeepTheEvent", [](G4Event const* a)->void { a->KeepTheEvent(); }, jlcxx::arg("this"), jlcxx::arg("vl"));
    t.method("KeepTheEvent", [](G4Event const* a, G4bool arg0)->void { a->KeepTheEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("vl"));

    DEBUG_MSG("Adding wrapper for G4bool G4Event::KeepTheEventFlag() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Event::KeepTheEventFlag()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:103:19
    t.method("KeepTheEventFlag", [](G4Event const& a)->G4bool { return a.KeepTheEventFlag(); }, jlcxx::arg("this"));
    t.method("KeepTheEventFlag", [](G4Event const* a)->G4bool { return a->KeepTheEventFlag(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4bool G4Event::ToBeKept() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Event::ToBeKept()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:105:19
    t.method("ToBeKept", [](G4Event const& a)->G4bool { return a.ToBeKept(); }, jlcxx::arg("this"));
    t.method("ToBeKept", [](G4Event const* a)->G4bool { return a->ToBeKept(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::KeepForPostProcessing() (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::KeepForPostProcessing()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:112:17
    t.method("KeepForPostProcessing", [](G4Event const& a)->void { a.KeepForPostProcessing(); }, jlcxx::arg("this"));
    t.method("KeepForPostProcessing", [](G4Event const* a)->void { a->KeepForPostProcessing(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::PostProcessingFinished() (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::PostProcessingFinished()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:114:17
    t.method("PostProcessingFinished", [](G4Event const& a)->void { a.PostProcessingFinished(); }, jlcxx::arg("this"));
    t.method("PostProcessingFinished", [](G4Event const* a)->void { a->PostProcessingFinished(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4int G4Event::GetNumberOfGrips() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Event::GetNumberOfGrips()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:123:18
    t.method("GetNumberOfGrips", [](G4Event const& a)->G4int { return a.GetNumberOfGrips(); }, jlcxx::arg("this"));
    t.method("GetNumberOfGrips", [](G4Event const* a)->G4int { return a->GetNumberOfGrips(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4int G4Event::GetEventID() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Event::GetEventID()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:126:18
    t.method("GetEventID", [](G4Event const& a)->G4int { return a.GetEventID(); }, jlcxx::arg("this"));
    t.method("GetEventID", [](G4Event const* a)->G4int { return a->GetEventID(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::AddPrimaryVertex(G4PrimaryVertex *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::AddPrimaryVertex(G4PrimaryVertex *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:129:17
    t.method("AddPrimaryVertex", [](G4Event& a, G4PrimaryVertex * arg0)->void { a.AddPrimaryVertex(arg0); }, jlcxx::arg("this"), jlcxx::arg("aPrimaryVertex"));
    t.method("AddPrimaryVertex", [](G4Event* a, G4PrimaryVertex * arg0)->void { a->AddPrimaryVertex(arg0); }, jlcxx::arg("this"), jlcxx::arg("aPrimaryVertex"));

    DEBUG_MSG("Adding wrapper for G4int G4Event::GetNumberOfPrimaryVertex() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Event::GetNumberOfPrimaryVertex()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:141:18
    t.method("GetNumberOfPrimaryVertex", [](G4Event const& a)->G4int { return a.GetNumberOfPrimaryVertex(); }, jlcxx::arg("this"));
    t.method("GetNumberOfPrimaryVertex", [](G4Event const* a)->G4int { return a->GetNumberOfPrimaryVertex(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4PrimaryVertex * G4Event::GetPrimaryVertex(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4PrimaryVertex * G4Event::GetPrimaryVertex(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:145:29
    t.method("GetPrimaryVertex", [](G4Event const& a)->G4PrimaryVertex * { return a.GetPrimaryVertex(); }, jlcxx::arg("this"));
    t.method("GetPrimaryVertex", [](G4Event const& a, G4int arg0)->G4PrimaryVertex * { return a.GetPrimaryVertex(arg0); }, jlcxx::arg("this"), jlcxx::arg("i"));
    t.method("GetPrimaryVertex", [](G4Event const* a)->G4PrimaryVertex * { return a->GetPrimaryVertex(); }, jlcxx::arg("this"));
    t.method("GetPrimaryVertex", [](G4Event const* a, G4int arg0)->G4PrimaryVertex * { return a->GetPrimaryVertex(arg0); }, jlcxx::arg("this"), jlcxx::arg("i"));

    DEBUG_MSG("Adding wrapper for G4HCofThisEvent * G4Event::GetHCofThisEvent() (" __HERE__ ")");
    // signature to use in the veto list: G4HCofThisEvent * G4Event::GetHCofThisEvent()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:164:29
    t.method("GetHCofThisEvent", [](G4Event const& a)->G4HCofThisEvent * { return a.GetHCofThisEvent(); }, jlcxx::arg("this"));
    t.method("GetHCofThisEvent", [](G4Event const* a)->G4HCofThisEvent * { return a->GetHCofThisEvent(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4DCofThisEvent * G4Event::GetDCofThisEvent() (" __HERE__ ")");
    // signature to use in the veto list: G4DCofThisEvent * G4Event::GetDCofThisEvent()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:166:29
    t.method("GetDCofThisEvent", [](G4Event const& a)->G4DCofThisEvent * { return a.GetDCofThisEvent(); }, jlcxx::arg("this"));
    t.method("GetDCofThisEvent", [](G4Event const* a)->G4DCofThisEvent * { return a->GetDCofThisEvent(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4TrajectoryContainer * G4Event::GetTrajectoryContainer() (" __HERE__ ")");
    // signature to use in the veto list: G4TrajectoryContainer * G4Event::GetTrajectoryContainer()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:168:35
    t.method("GetTrajectoryContainer", [](G4Event const& a)->G4TrajectoryContainer * { return a.GetTrajectoryContainer(); }, jlcxx::arg("this"));
    t.method("GetTrajectoryContainer", [](G4Event const* a)->G4TrajectoryContainer * { return a->GetTrajectoryContainer(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4bool G4Event::IsAborted() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Event::IsAborted()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:175:19
    t.method("IsAborted", [](G4Event const& a)->G4bool { return a.IsAborted(); }, jlcxx::arg("this"));
    t.method("IsAborted", [](G4Event const* a)->G4bool { return a->IsAborted(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::SetUserInformation(G4VUserEventInformation *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::SetUserInformation(G4VUserEventInformation *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:179:17
    t.method("SetUserInformation", [](G4Event& a, G4VUserEventInformation * arg0)->void { a.SetUserInformation(arg0); }, jlcxx::arg("this"), jlcxx::arg("anInfo"));
    t.method("SetUserInformation", [](G4Event* a, G4VUserEventInformation * arg0)->void { a->SetUserInformation(arg0); }, jlcxx::arg("this"), jlcxx::arg("anInfo"));

    DEBUG_MSG("Adding wrapper for G4VUserEventInformation * G4Event::GetUserInformation() (" __HERE__ ")");
    // signature to use in the veto list: G4VUserEventInformation * G4Event::GetUserInformation()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:181:37
    t.method("GetUserInformation", [](G4Event const& a)->G4VUserEventInformation * { return a.GetUserInformation(); }, jlcxx::arg("this"));
    t.method("GetUserInformation", [](G4Event const* a)->G4VUserEventInformation * { return a->GetUserInformation(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for const G4String & G4Event::GetRandomNumberStatus() (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4Event::GetRandomNumberStatus()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:185:28
    t.method("GetRandomNumberStatus", [](G4Event const& a)->const G4String & { return a.GetRandomNumberStatus(); }, jlcxx::arg("this"));
    t.method("GetRandomNumberStatus", [](G4Event const* a)->const G4String & { return a->GetRandomNumberStatus(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for const G4String & G4Event::GetRandomNumberStatusForProcessing() (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4Event::GetRandomNumberStatusForProcessing()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:193:28
    t.method("GetRandomNumberStatusForProcessing", [](G4Event const& a)->const G4String & { return a.GetRandomNumberStatusForProcessing(); }, jlcxx::arg("this"));
    t.method("GetRandomNumberStatusForProcessing", [](G4Event const* a)->const G4String & { return a->GetRandomNumberStatusForProcessing(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4SubEvent * G4Event::PopSubEvent(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4SubEvent * G4Event::PopSubEvent(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:244:17
    t.method("PopSubEvent", [](G4Event& a, G4int arg0)->G4SubEvent * { return a.PopSubEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));
    t.method("PopSubEvent", [](G4Event* a, G4int arg0)->G4SubEvent * { return a->PopSubEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));

    DEBUG_MSG("Adding wrapper for G4int G4Event::TerminateSubEvent(G4SubEvent *) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Event::TerminateSubEvent(G4SubEvent *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:247:11
    t.method("TerminateSubEvent", [](G4Event& a, G4SubEvent * arg0)->G4int { return a.TerminateSubEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));
    t.method("TerminateSubEvent", [](G4Event* a, G4SubEvent * arg0)->G4int { return a->TerminateSubEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));

    DEBUG_MSG("Adding wrapper for G4int G4Event::StoreSubEvent(G4int, G4SubEvent *) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Event::StoreSubEvent(G4int, G4SubEvent *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:251:11
    t.method("StoreSubEvent", [](G4Event& a, G4int arg0, G4SubEvent * arg1)->G4int { return a.StoreSubEvent(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("arg0"), jlcxx::arg("arg1"));
    t.method("StoreSubEvent", [](G4Event* a, G4int arg0, G4SubEvent * arg1)->G4int { return a->StoreSubEvent(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("arg0"), jlcxx::arg("arg1"));

    DEBUG_MSG("Adding wrapper for G4int G4Event::SpawnSubEvent(G4SubEvent *) (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Event::SpawnSubEvent(G4SubEvent *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:255:11
    t.method("SpawnSubEvent", [](G4Event& a, G4SubEvent * arg0)->G4int { return a.SpawnSubEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));
    t.method("SpawnSubEvent", [](G4Event* a, G4SubEvent * arg0)->G4int { return a->SpawnSubEvent(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));

    DEBUG_MSG("Adding wrapper for G4int G4Event::GetNumberOfRemainingSubEvents() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Event::GetNumberOfRemainingSubEvents()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:257:11
    t.method("GetNumberOfRemainingSubEvents", [](G4Event const& a)->G4int { return a.GetNumberOfRemainingSubEvents(); }, jlcxx::arg("this"));
    t.method("GetNumberOfRemainingSubEvents", [](G4Event const* a)->G4int { return a->GetNumberOfRemainingSubEvents(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4int G4Event::GetNumberOfCompletedSubEvent() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Event::GetNumberOfCompletedSubEvent()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:260:18
    t.method("GetNumberOfCompletedSubEvent", [](G4Event const& a)->G4int { return a.GetNumberOfCompletedSubEvent(); }, jlcxx::arg("this"));
    t.method("GetNumberOfCompletedSubEvent", [](G4Event const* a)->G4int { return a->GetNumberOfCompletedSubEvent(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::MergeSubEventResults(const G4Event *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::MergeSubEventResults(const G4Event *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:263:10
    t.method("MergeSubEventResults", [](G4Event& a, const G4Event * arg0)->void { a.MergeSubEventResults(arg0); }, jlcxx::arg("this"), jlcxx::arg("se"));
    t.method("MergeSubEventResults", [](G4Event* a, const G4Event * arg0)->void { a->MergeSubEventResults(arg0); }, jlcxx::arg("this"), jlcxx::arg("se"));

    DEBUG_MSG("Adding wrapper for void G4Event::FlagAsSubEvent(G4Event *, G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::FlagAsSubEvent(G4Event *, G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:287:10
    t.method("FlagAsSubEvent", [](G4Event& a, G4Event * arg0, G4int arg1)->void { a.FlagAsSubEvent(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("me"), jlcxx::arg("ty"));
    t.method("FlagAsSubEvent", [](G4Event* a, G4Event * arg0, G4int arg1)->void { a->FlagAsSubEvent(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("me"), jlcxx::arg("ty"));

    DEBUG_MSG("Adding wrapper for G4Event * G4Event::GetMotherEvent() (" __HERE__ ")");
    // signature to use in the veto list: G4Event * G4Event::GetMotherEvent()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:292:21
    t.method("GetMotherEvent", [](G4Event const& a)->G4Event * { return a.GetMotherEvent(); }, jlcxx::arg("this"));
    t.method("GetMotherEvent", [](G4Event const* a)->G4Event * { return a->GetMotherEvent(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4int G4Event::GetSubEventType() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4Event::GetSubEventType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:294:18
    t.method("GetSubEventType", [](G4Event const& a)->G4int { return a.GetSubEventType(); }, jlcxx::arg("this"));
    t.method("GetSubEventType", [](G4Event const* a)->G4int { return a->GetSubEventType(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::ScoresRecorded() (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::ScoresRecorded()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:301:10
    t.method("ScoresRecorded", [](G4Event const& a)->void { a.ScoresRecorded(); }, jlcxx::arg("this"));
    t.method("ScoresRecorded", [](G4Event const* a)->void { a->ScoresRecorded(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4bool G4Event::ScoresAlreadyRecorded() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Event::ScoresAlreadyRecorded()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:302:12
    t.method("ScoresAlreadyRecorded", [](G4Event const& a)->G4bool { return a.ScoresAlreadyRecorded(); }, jlcxx::arg("this"));
    t.method("ScoresAlreadyRecorded", [](G4Event const* a)->G4bool { return a->ScoresAlreadyRecorded(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Event::EventCompleted() (" __HERE__ ")");
    // signature to use in the veto list: void G4Event::EventCompleted()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:303:10
    t.method("EventCompleted", [](G4Event const& a)->void { a.EventCompleted(); }, jlcxx::arg("this"));
    t.method("EventCompleted", [](G4Event const* a)->void { a->EventCompleted(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4bool G4Event::IsEventCompleted() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Event::IsEventCompleted()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Event.hh:304:12
    t.method("IsEventCompleted", [](G4Event const& a)->G4bool { return a.IsEventCompleted(); }, jlcxx::arg("this"));
    t.method("IsEventCompleted", [](G4Event const* a)->G4bool { return a->IsEventCompleted(); }, jlcxx::arg("this"));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Event>> type_;
};
std::shared_ptr<Wrapper> newJlG4Event(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Event(module));
}
