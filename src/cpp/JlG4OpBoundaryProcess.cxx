// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4OpBoundaryProcess> : std::false_type { };
  template<> struct DefaultConstructible<G4OpBoundaryProcess> : std::false_type { };
}

// Class generating the wrapper for type G4OpBoundaryProcess
// signature to use in the veto file: G4OpBoundaryProcess
struct JlG4OpBoundaryProcess: public Wrapper {

  JlG4OpBoundaryProcess(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4OpBoundaryProcess (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:119:7
    jlcxx::TypeWrapper<G4OpBoundaryProcess>  t = jlModule.add_type<G4OpBoundaryProcess>("G4OpBoundaryProcess");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4OpBoundaryProcess>>(new jlcxx::TypeWrapper<G4OpBoundaryProcess>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes    );


    DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::G4OpBoundaryProcess(const G4String &, G4ProcessType) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:122:12
    t.constructor<const G4String &>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("processName")    );
    t.constructor<const G4String &, G4ProcessType>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("processName"), jlcxx::arg("type")    );

    DEBUG_MSG("Adding wrapper for G4bool G4OpBoundaryProcess::IsApplicable(const G4ParticleDefinition &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4OpBoundaryProcess::IsApplicable(const G4ParticleDefinition &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:126:18
    t.method("IsApplicable", [](G4OpBoundaryProcess& a, const G4ParticleDefinition & arg0)->G4bool { return a.IsApplicable(arg0); }, jlcxx::arg("this"), jlcxx::arg("aParticleType"));
    t.method("IsApplicable", [](G4OpBoundaryProcess* a, const G4ParticleDefinition & arg0)->G4bool { return a->IsApplicable(arg0); }, jlcxx::arg("this"), jlcxx::arg("aParticleType"));

    DEBUG_MSG("Adding wrapper for G4double G4OpBoundaryProcess::GetMeanFreePath(const G4Track &, G4double, G4ForceCondition *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4OpBoundaryProcess::GetMeanFreePath(const G4Track &, G4double, G4ForceCondition *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:130:20
    t.method("GetMeanFreePath", [](G4OpBoundaryProcess& a, const G4Track & arg0, G4double arg1, G4ForceCondition * arg2)->G4double { return a.GetMeanFreePath(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("arg0"), jlcxx::arg("arg1"), jlcxx::arg("condition"));
    t.method("GetMeanFreePath", [](G4OpBoundaryProcess* a, const G4Track & arg0, G4double arg1, G4ForceCondition * arg2)->G4double { return a->GetMeanFreePath(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("arg0"), jlcxx::arg("arg1"), jlcxx::arg("condition"));

    DEBUG_MSG("Adding wrapper for G4VParticleChange * G4OpBoundaryProcess::PostStepDoIt(const G4Track &, const G4Step &) (" __HERE__ ")");
    // signature to use in the veto list: G4VParticleChange * G4OpBoundaryProcess::PostStepDoIt(const G4Track &, const G4Step &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:136:22
    t.method("PostStepDoIt", [](G4OpBoundaryProcess& a, const G4Track & arg0, const G4Step & arg1)->G4VParticleChange * { return a.PostStepDoIt(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("aTrack"), jlcxx::arg("aStep"));
    t.method("PostStepDoIt", [](G4OpBoundaryProcess* a, const G4Track & arg0, const G4Step & arg1)->G4VParticleChange * { return a->PostStepDoIt(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("aTrack"), jlcxx::arg("aStep"));

    DEBUG_MSG("Adding wrapper for G4OpBoundaryProcessStatus G4OpBoundaryProcess::GetStatus() (" __HERE__ ")");
    // signature to use in the veto list: G4OpBoundaryProcessStatus G4OpBoundaryProcess::GetStatus()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:140:37
    t.method("GetStatus", [](G4OpBoundaryProcess const& a)->G4OpBoundaryProcessStatus { return a.GetStatus(); }, jlcxx::arg("this"));
    t.method("GetStatus", [](G4OpBoundaryProcess const* a)->G4OpBoundaryProcessStatus { return a->GetStatus(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::SetInvokeSD(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4OpBoundaryProcess::SetInvokeSD(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:143:16
    t.method("SetInvokeSD", [](G4OpBoundaryProcess& a, G4bool arg0)->void { a.SetInvokeSD(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));
    t.method("SetInvokeSD", [](G4OpBoundaryProcess* a, G4bool arg0)->void { a->SetInvokeSD(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));

    DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::PreparePhysicsTable(const G4ParticleDefinition &) (" __HERE__ ")");
    // signature to use in the veto list: void G4OpBoundaryProcess::PreparePhysicsTable(const G4ParticleDefinition &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:146:16
    t.method("PreparePhysicsTable", [](G4OpBoundaryProcess& a, const G4ParticleDefinition & arg0)->void { a.PreparePhysicsTable(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));
    t.method("PreparePhysicsTable", [](G4OpBoundaryProcess* a, const G4ParticleDefinition & arg0)->void { a->PreparePhysicsTable(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));

    DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::Initialise() (" __HERE__ ")");
    // signature to use in the veto list: void G4OpBoundaryProcess::Initialise()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:148:16
    t.method("Initialise", [](G4OpBoundaryProcess& a)->void { a.Initialise(); }, jlcxx::arg("this"));
    t.method("Initialise", [](G4OpBoundaryProcess* a)->void { a->Initialise(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4OpBoundaryProcess::SetVerboseLevel(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4OpBoundaryProcess::SetVerboseLevel(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4OpBoundaryProcess.hh:150:8
    t.method("SetVerboseLevel", [](G4OpBoundaryProcess& a, G4int arg0)->void { a.SetVerboseLevel(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));
    t.method("SetVerboseLevel", [](G4OpBoundaryProcess* a, G4int arg0)->void { a->SetVerboseLevel(arg0); }, jlcxx::arg("this"), jlcxx::arg("arg0"));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4OpBoundaryProcess>> type_;
};
std::shared_ptr<Wrapper> newJlG4OpBoundaryProcess(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4OpBoundaryProcess(module));
}
