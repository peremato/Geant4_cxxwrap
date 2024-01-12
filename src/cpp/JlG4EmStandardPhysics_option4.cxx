// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4EmStandardPhysics_option4> : std::false_type { };
  template<> struct DefaultConstructible<G4EmStandardPhysics_option4> : std::false_type { };
template<> struct SuperType<G4EmStandardPhysics_option4> { typedef G4VPhysicsConstructor type; };
}

// Class generating the wrapper for type G4EmStandardPhysics_option4
// signature to use in the veto file: G4EmStandardPhysics_option4
struct JlG4EmStandardPhysics_option4: public Wrapper {

  JlG4EmStandardPhysics_option4(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4EmStandardPhysics_option4 (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4EmStandardPhysics_option4.hh:52:7
    jlcxx::TypeWrapper<G4EmStandardPhysics_option4>  t = jlModule.add_type<G4EmStandardPhysics_option4>("G4EmStandardPhysics_option4",
      jlcxx::julia_base_type<G4VPhysicsConstructor>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4EmStandardPhysics_option4>>(new jlcxx::TypeWrapper<G4EmStandardPhysics_option4>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4EmStandardPhysics_option4::G4EmStandardPhysics_option4(G4int, const G4String &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4EmStandardPhysics_option4.hh:56:12
    t.constructor<G4int>(/*finalize=*/true);
    t.constructor<G4int, const G4String &>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for void G4EmStandardPhysics_option4::ConstructParticle() (" __HERE__ ")");
    // signature to use in the veto list: void G4EmStandardPhysics_option4::ConstructParticle()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4EmStandardPhysics_option4.hh:60:8
    t.method("ConstructParticle", static_cast<void (G4EmStandardPhysics_option4::*)() >(&G4EmStandardPhysics_option4::ConstructParticle));

    DEBUG_MSG("Adding wrapper for void G4EmStandardPhysics_option4::ConstructProcess() (" __HERE__ ")");
    // signature to use in the veto list: void G4EmStandardPhysics_option4::ConstructProcess()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4EmStandardPhysics_option4.hh:61:8
    t.method("ConstructProcess", static_cast<void (G4EmStandardPhysics_option4::*)() >(&G4EmStandardPhysics_option4::ConstructProcess));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4EmStandardPhysics_option4>> type_;
};
std::shared_ptr<Wrapper> newJlG4EmStandardPhysics_option4(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4EmStandardPhysics_option4(module));
}
