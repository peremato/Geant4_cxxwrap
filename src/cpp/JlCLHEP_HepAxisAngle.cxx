// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<CLHEP::HepAxisAngle> : std::false_type { };
  template<> struct DefaultConstructible<CLHEP::HepAxisAngle> : std::false_type { };
}

// Class generating the wrapper for type CLHEP::HepAxisAngle
// signature to use in the veto file: CLHEP::HepAxisAngle
struct JlCLHEP_HepAxisAngle: public Wrapper {

  JlCLHEP_HepAxisAngle(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type CLHEP::HepAxisAngle (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/AxisAngle.h:36:7
    jlcxx::TypeWrapper<CLHEP::HepAxisAngle>  t = jlModule.add_type<CLHEP::HepAxisAngle>("CLHEP!HepAxisAngle");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<CLHEP::HepAxisAngle>>(new jlcxx::TypeWrapper<CLHEP::HepAxisAngle>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes    );
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<CLHEP::HepAxisAngle>> type_;
};
std::shared_ptr<Wrapper> newJlCLHEP_HepAxisAngle(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlCLHEP_HepAxisAngle(module));
}
