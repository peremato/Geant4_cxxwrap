// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4PrimaryVertex> : std::false_type { };
  template<> struct DefaultConstructible<G4PrimaryVertex> : std::false_type { };
}

// Class generating the wrapper for type G4PrimaryVertex
// signature to use in the veto file: G4PrimaryVertex
struct JlG4PrimaryVertex: public Wrapper {

  JlG4PrimaryVertex(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4PrimaryVertex (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:50:7
    jlcxx::TypeWrapper<G4PrimaryVertex>  t = jlModule.add_type<G4PrimaryVertex>("G4PrimaryVertex");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4PrimaryVertex>>(new jlcxx::TypeWrapper<G4PrimaryVertex>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes    );


    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::G4PrimaryVertex(G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:55:5
    t.constructor<G4double, G4double, G4double, G4double>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("x0"), jlcxx::arg("y0"), jlcxx::arg("z0"), jlcxx::arg("t0")    );


    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::G4PrimaryVertex(G4ThreeVector, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:56:5
    t.constructor<G4ThreeVector, G4double>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("xyz0"), jlcxx::arg("t0")    );


    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::G4PrimaryVertex(const G4PrimaryVertex &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:62:5
    t.constructor<const G4PrimaryVertex &>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("right")    );

    DEBUG_MSG("Adding wrapper for G4PrimaryVertex & G4PrimaryVertex::operator=(const G4PrimaryVertex &) (" __HERE__ ")");
    // signature to use in the veto list: G4PrimaryVertex & G4PrimaryVertex::operator=(const G4PrimaryVertex &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:63:22
    t.method("assign", [](G4PrimaryVertex& a, const G4PrimaryVertex & arg0)->G4PrimaryVertex & { return a.operator=(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));
    t.method("assign", [](G4PrimaryVertex* a, const G4PrimaryVertex & arg0)->G4PrimaryVertex & { return a->operator=(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4PrimaryVertex::operator==(const G4PrimaryVertex &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4PrimaryVertex::operator==(const G4PrimaryVertex &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:66:12
    t.method("==", [](G4PrimaryVertex const& a, const G4PrimaryVertex & arg0)->G4bool { return a.operator==(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));
    t.method("==", [](G4PrimaryVertex const* a, const G4PrimaryVertex & arg0)->G4bool { return a->operator==(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));

    DEBUG_MSG("Adding wrapper for G4bool G4PrimaryVertex::operator!=(const G4PrimaryVertex &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4PrimaryVertex::operator!=(const G4PrimaryVertex &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:67:12
    t.method("!=", [](G4PrimaryVertex const& a, const G4PrimaryVertex & arg0)->G4bool { return a.operator!=(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));
    t.method("!=", [](G4PrimaryVertex const* a, const G4PrimaryVertex & arg0)->G4bool { return a->operator!=(arg0); }, jlcxx::arg("this"), jlcxx::arg("right"));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4PrimaryVertex::GetPosition() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4PrimaryVertex::GetPosition()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:75:26
    t.method("GetPosition", [](G4PrimaryVertex const& a)->G4ThreeVector { return a.GetPosition(); }, jlcxx::arg("this"));
    t.method("GetPosition", [](G4PrimaryVertex const* a)->G4ThreeVector { return a->GetPosition(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetPosition(G4double, G4double, G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryVertex::SetPosition(G4double, G4double, G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:76:17
    t.method("SetPosition", [](G4PrimaryVertex& a, G4double arg0, G4double arg1, G4double arg2)->void { a.SetPosition(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("x0"), jlcxx::arg("y0"), jlcxx::arg("z0"));
    t.method("SetPosition", [](G4PrimaryVertex* a, G4double arg0, G4double arg1, G4double arg2)->void { a->SetPosition(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("x0"), jlcxx::arg("y0"), jlcxx::arg("z0"));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetX0() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryVertex::GetX0()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:77:21
    t.method("GetX0", [](G4PrimaryVertex const& a)->G4double { return a.GetX0(); }, jlcxx::arg("this"));
    t.method("GetX0", [](G4PrimaryVertex const* a)->G4double { return a->GetX0(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetY0() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryVertex::GetY0()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:78:21
    t.method("GetY0", [](G4PrimaryVertex const& a)->G4double { return a.GetY0(); }, jlcxx::arg("this"));
    t.method("GetY0", [](G4PrimaryVertex const* a)->G4double { return a->GetY0(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetZ0() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryVertex::GetZ0()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:79:21
    t.method("GetZ0", [](G4PrimaryVertex const& a)->G4double { return a.GetZ0(); }, jlcxx::arg("this"));
    t.method("GetZ0", [](G4PrimaryVertex const* a)->G4double { return a->GetZ0(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetT0() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryVertex::GetT0()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:80:21
    t.method("GetT0", [](G4PrimaryVertex const& a)->G4double { return a.GetT0(); }, jlcxx::arg("this"));
    t.method("GetT0", [](G4PrimaryVertex const* a)->G4double { return a->GetT0(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetT0(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryVertex::SetT0(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:81:17
    t.method("SetT0", [](G4PrimaryVertex& a, G4double arg0)->void { a.SetT0(arg0); }, jlcxx::arg("this"), jlcxx::arg("t0"));
    t.method("SetT0", [](G4PrimaryVertex* a, G4double arg0)->void { a->SetT0(arg0); }, jlcxx::arg("this"), jlcxx::arg("t0"));

    DEBUG_MSG("Adding wrapper for G4int G4PrimaryVertex::GetNumberOfParticle() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4PrimaryVertex::GetNumberOfParticle()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:82:18
    t.method("GetNumberOfParticle", [](G4PrimaryVertex const& a)->G4int { return a.GetNumberOfParticle(); }, jlcxx::arg("this"));
    t.method("GetNumberOfParticle", [](G4PrimaryVertex const* a)->G4int { return a->GetNumberOfParticle(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetPrimary(G4PrimaryParticle *) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryVertex::SetPrimary(G4PrimaryParticle *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:83:17
    t.method("SetPrimary", [](G4PrimaryVertex& a, G4PrimaryParticle * arg0)->void { a.SetPrimary(arg0); }, jlcxx::arg("this"), jlcxx::arg("pp"));
    t.method("SetPrimary", [](G4PrimaryVertex* a, G4PrimaryParticle * arg0)->void { a->SetPrimary(arg0); }, jlcxx::arg("this"), jlcxx::arg("pp"));

    DEBUG_MSG("Adding wrapper for G4PrimaryParticle * G4PrimaryVertex::GetPrimary(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4PrimaryParticle * G4PrimaryVertex::GetPrimary(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:84:24
    t.method("GetPrimary", [](G4PrimaryVertex const& a)->G4PrimaryParticle * { return a.GetPrimary(); }, jlcxx::arg("this"));
    t.method("GetPrimary", [](G4PrimaryVertex const& a, G4int arg0)->G4PrimaryParticle * { return a.GetPrimary(arg0); }, jlcxx::arg("this"), jlcxx::arg("i"));
    t.method("GetPrimary", [](G4PrimaryVertex const* a)->G4PrimaryParticle * { return a->GetPrimary(); }, jlcxx::arg("this"));
    t.method("GetPrimary", [](G4PrimaryVertex const* a, G4int arg0)->G4PrimaryParticle * { return a->GetPrimary(arg0); }, jlcxx::arg("this"), jlcxx::arg("i"));

    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetNext(G4PrimaryVertex *) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryVertex::SetNext(G4PrimaryVertex *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:85:17
    t.method("SetNext", [](G4PrimaryVertex& a, G4PrimaryVertex * arg0)->void { a.SetNext(arg0); }, jlcxx::arg("this"), jlcxx::arg("nv"));
    t.method("SetNext", [](G4PrimaryVertex* a, G4PrimaryVertex * arg0)->void { a->SetNext(arg0); }, jlcxx::arg("this"), jlcxx::arg("nv"));

    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::ClearNext() (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryVertex::ClearNext()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:86:17
    t.method("ClearNext", [](G4PrimaryVertex& a)->void { a.ClearNext(); }, jlcxx::arg("this"));
    t.method("ClearNext", [](G4PrimaryVertex* a)->void { a->ClearNext(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4PrimaryVertex * G4PrimaryVertex::GetNext() (" __HERE__ ")");
    // signature to use in the veto list: G4PrimaryVertex * G4PrimaryVertex::GetNext()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:87:29
    t.method("GetNext", [](G4PrimaryVertex const& a)->G4PrimaryVertex * { return a.GetNext(); }, jlcxx::arg("this"));
    t.method("GetNext", [](G4PrimaryVertex const* a)->G4PrimaryVertex * { return a->GetNext(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryVertex::GetWeight() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryVertex::GetWeight()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:88:21
    t.method("GetWeight", [](G4PrimaryVertex const& a)->G4double { return a.GetWeight(); }, jlcxx::arg("this"));
    t.method("GetWeight", [](G4PrimaryVertex const* a)->G4double { return a->GetWeight(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::SetWeight(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryVertex::SetWeight(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:89:17
    t.method("SetWeight", [](G4PrimaryVertex& a, G4double arg0)->void { a.SetWeight(arg0); }, jlcxx::arg("this"), jlcxx::arg("w"));
    t.method("SetWeight", [](G4PrimaryVertex* a, G4double arg0)->void { a->SetWeight(arg0); }, jlcxx::arg("this"), jlcxx::arg("w"));

    DEBUG_MSG("Adding wrapper for void G4PrimaryVertex::Print() (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryVertex::Print()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PrimaryVertex.hh:93:10
    t.method("Print", [](G4PrimaryVertex const& a)->void { a.Print(); }, jlcxx::arg("this"));
    t.method("Print", [](G4PrimaryVertex const* a)->void { a->Print(); }, jlcxx::arg("this"));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4PrimaryVertex>> type_;
};
std::shared_ptr<Wrapper> newJlG4PrimaryVertex(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4PrimaryVertex(module));
}
