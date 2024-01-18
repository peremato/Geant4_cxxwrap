// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4OpticalPhysics> : std::false_type { };
  template<> struct DefaultConstructible<G4OpticalPhysics> : std::false_type { };
template<> struct SuperType<G4OpticalPhysics> { typedef G4VPhysicsConstructor type; };
}

// Class generating the wrapper for type G4OpticalPhysics
// signature to use in the veto file: G4OpticalPhysics
struct JlG4OpticalPhysics: public Wrapper {

  JlG4OpticalPhysics(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4OpticalPhysics (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4OpticalPhysics.hh:50:7
    jlcxx::TypeWrapper<G4OpticalPhysics>  t = jlModule.add_type<G4OpticalPhysics>("G4OpticalPhysics",
      jlcxx::julia_base_type<G4VPhysicsConstructor>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4OpticalPhysics>>(new jlcxx::TypeWrapper<G4OpticalPhysics>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4OpticalPhysics::G4OpticalPhysics(G4int, const G4String &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4OpticalPhysics.hh:53:3
    t.constructor<G4int>(/*finalize=*/true);
    t.constructor<G4int, const G4String &>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for void G4OpticalPhysics::PrintStatistics() (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalPhysics::PrintStatistics()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4OpticalPhysics.hh:55:8
    t.method("PrintStatistics", static_cast<void (G4OpticalPhysics::*)()  const>(&G4OpticalPhysics::PrintStatistics));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4OpticalPhysics>> type_;
};
std::shared_ptr<Wrapper> newJlG4OpticalPhysics(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4OpticalPhysics(module));
}