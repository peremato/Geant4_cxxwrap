// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4DecayPhysics> : std::false_type { };
  template<> struct DefaultConstructible<G4DecayPhysics> : std::false_type { };
template<> struct SuperType<G4DecayPhysics> { typedef G4VPhysicsConstructor type; };
}

// Class generating the wrapper for type G4DecayPhysics
// signature to use in the veto file: G4DecayPhysics
struct JlG4DecayPhysics: public Wrapper {

  JlG4DecayPhysics(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4DecayPhysics (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4DecayPhysics.hh:48:7
    jlcxx::TypeWrapper<G4DecayPhysics>  t = jlModule.add_type<G4DecayPhysics>("G4DecayPhysics",
      jlcxx::julia_base_type<G4VPhysicsConstructor>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4DecayPhysics>>(new jlcxx::TypeWrapper<G4DecayPhysics>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes    );


    DEBUG_MSG("Adding wrapper for void G4DecayPhysics::G4DecayPhysics(G4int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4DecayPhysics.hh:51:3
    t.constructor<G4int>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("ver")    );


    DEBUG_MSG("Adding wrapper for void G4DecayPhysics::G4DecayPhysics(const G4String &, G4int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4DecayPhysics.hh:52:3
    t.constructor<const G4String &>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("name")    );
    t.constructor<const G4String &, G4int>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("name"), jlcxx::arg("ver")    );

    DEBUG_MSG("Adding wrapper for void G4DecayPhysics::ConstructParticle() (" __HERE__ ")");
    // signature to use in the veto list: void G4DecayPhysics::ConstructParticle()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4DecayPhysics.hh:57:8
    t.method("ConstructParticle", [](G4DecayPhysics& a)->void { a.ConstructParticle(); }, jlcxx::arg("this"));
    t.method("ConstructParticle", [](G4DecayPhysics* a)->void { a->ConstructParticle(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4DecayPhysics::ConstructProcess() (" __HERE__ ")");
    // signature to use in the veto list: void G4DecayPhysics::ConstructProcess()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4DecayPhysics.hh:62:8
    t.method("ConstructProcess", [](G4DecayPhysics& a)->void { a.ConstructProcess(); }, jlcxx::arg("this"));
    t.method("ConstructProcess", [](G4DecayPhysics* a)->void { a->ConstructProcess(); }, jlcxx::arg("this"));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4DecayPhysics>> type_;
};
std::shared_ptr<Wrapper> newJlG4DecayPhysics(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4DecayPhysics(module));
}
