// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4IonisParamElm> : std::false_type { };
  template<> struct DefaultConstructible<G4IonisParamElm> : std::false_type { };
}

// Class generating the wrapper for type G4IonisParamElm
// signature to use in the veto file: G4IonisParamElm
struct JlG4IonisParamElm: public Wrapper {

  JlG4IonisParamElm(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4IonisParamElm (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4IonisParamElm.hh:40:7
    jlcxx::TypeWrapper<G4IonisParamElm>  t = jlModule.add_type<G4IonisParamElm>("G4IonisParamElm");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4IonisParamElm>>(new jlcxx::TypeWrapper<G4IonisParamElm>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4IonisParamElm>> type_;
};
std::shared_ptr<Wrapper> newJlG4IonisParamElm(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4IonisParamElm(module));
}
