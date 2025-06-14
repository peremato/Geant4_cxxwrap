// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4SPSRandomGenerator> : std::false_type { };
  template<> struct DefaultConstructible<G4SPSRandomGenerator> : std::false_type { };
}

// Class generating the wrapper for type G4SPSRandomGenerator
// signature to use in the veto file: G4SPSRandomGenerator
struct JlG4SPSRandomGenerator: public Wrapper {

  JlG4SPSRandomGenerator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4SPSRandomGenerator (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4SPSRandomGenerator.hh:60:7
    jlcxx::TypeWrapper<G4SPSRandomGenerator>  t = jlModule.add_type<G4SPSRandomGenerator>("G4SPSRandomGenerator");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4SPSRandomGenerator>>(new jlcxx::TypeWrapper<G4SPSRandomGenerator>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes    );
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4SPSRandomGenerator>> type_;
};
std::shared_ptr<Wrapper> newJlG4SPSRandomGenerator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4SPSRandomGenerator(module));
}
