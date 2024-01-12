// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4GDMLWriteStructure> : std::false_type { };
  template<> struct DefaultConstructible<G4GDMLWriteStructure> : std::false_type { };
}

// Class generating the wrapper for type G4GDMLWriteStructure
// signature to use in the veto file: G4GDMLWriteStructure
struct JlG4GDMLWriteStructure: public Wrapper {

  JlG4GDMLWriteStructure(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4GDMLWriteStructure (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4GDMLWriteStructure.hh:51:7
    jlcxx::TypeWrapper<G4GDMLWriteStructure>  t = jlModule.add_type<G4GDMLWriteStructure>("G4GDMLWriteStructure");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4GDMLWriteStructure>>(new jlcxx::TypeWrapper<G4GDMLWriteStructure>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4GDMLWriteStructure>> type_;
};
std::shared_ptr<Wrapper> newJlG4GDMLWriteStructure(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4GDMLWriteStructure(module));
}
