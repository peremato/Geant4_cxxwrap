// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<CLHEP::HepRotation::HepRotation_row> : std::false_type { };
  template<> struct DefaultConstructible<CLHEP::HepRotation::HepRotation_row> : std::false_type { };
}

// Class generating the wrapper for type CLHEP::HepRotation::HepRotation_row
// signature to use in the veto file: CLHEP::HepRotation::HepRotation_row
struct JlCLHEP_HepRotation_HepRotation_row: public Wrapper {

  JlCLHEP_HepRotation_HepRotation_row(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type CLHEP::HepRotation::HepRotation_row (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Vector/Rotation.h:147:9
    jlcxx::TypeWrapper<CLHEP::HepRotation::HepRotation_row>  t = jlModule.add_type<CLHEP::HepRotation::HepRotation_row>("CLHEP!HepRotation!HepRotation_row");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<CLHEP::HepRotation::HepRotation_row>>(new jlcxx::TypeWrapper<CLHEP::HepRotation::HepRotation_row>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void CLHEP::HepRotation::HepRotation_row::HepRotation_row(const CLHEP::HepRotation &, int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Vector/Rotation.h:149:12
    t.constructor<const CLHEP::HepRotation &, int>(/*finalize=*/jlcxx::finalize_policy::yes);
    module_.set_override_module(jl_base_module);


    DEBUG_MSG("Adding getindex method to wrap double CLHEP::HepRotation::HepRotation_row::operator[](int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Vector/Rotation.h:150:19
    t.method("getindex",
      [](CLHEP::HepRotation::HepRotation_row& a, int i){
      return a[i];
    });

    module_.unset_override_module();
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<CLHEP::HepRotation::HepRotation_row>> type_;
};
std::shared_ptr<Wrapper> newJlCLHEP_HepRotation_HepRotation_row(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlCLHEP_HepRotation_HepRotation_row(module));
}
