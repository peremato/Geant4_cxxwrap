// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4JLMagField> : std::false_type { };
  template<> struct DefaultConstructible<G4JLMagField> : std::false_type { };
template<> struct SuperType<G4JLMagField> { typedef G4MagneticField type; };
}

// Class generating the wrapper for type G4JLMagField
// signature to use in the veto file: G4JLMagField
struct JlG4JLMagField: public Wrapper {

  JlG4JLMagField(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4JLMagField (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:153:7
    jlcxx::TypeWrapper<G4JLMagField>  t = jlModule.add_type<G4JLMagField>("G4JLMagField",
      jlcxx::julia_base_type<G4MagneticField>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4JLMagField>>(new jlcxx::TypeWrapper<G4JLMagField>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4JLMagField::G4JLMagField(getfield_f, void *) (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:155:3
    t.constructor<getfield_f, void *>(/*finalize=*/true);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4JLMagField>> type_;
};
std::shared_ptr<Wrapper> newJlG4JLMagField(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4JLMagField(module));
}
