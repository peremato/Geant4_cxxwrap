// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Field> : std::false_type { };
  template<> struct DefaultConstructible<G4Field> : std::false_type { };
}

// Class generating the wrapper for type G4Field
// signature to use in the veto file: G4Field
struct JlG4Field: public Wrapper {

  JlG4Field(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Field (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4Field.hh:55:7
    jlcxx::TypeWrapper<G4Field>  t = jlModule.add_type<G4Field>("G4Field");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Field>>(new jlcxx::TypeWrapper<G4Field>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Field>> type_;
};
std::shared_ptr<Wrapper> newJlG4Field(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Field(module));
}
