// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<CLHEP::HepEulerAngles> : std::false_type { };
  template<> struct DefaultConstructible<CLHEP::HepEulerAngles> : std::false_type { };
}

// Class generating the wrapper for type CLHEP::HepEulerAngles
// signature to use in the veto file: CLHEP::HepEulerAngles
struct JlCLHEP_HepEulerAngles: public Wrapper {

  JlCLHEP_HepEulerAngles(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type CLHEP::HepEulerAngles (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:38:7
    jlcxx::TypeWrapper<CLHEP::HepEulerAngles>  t = jlModule.add_type<CLHEP::HepEulerAngles>("CLHEP!HepEulerAngles");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<CLHEP::HepEulerAngles>>(new jlcxx::TypeWrapper<CLHEP::HepEulerAngles>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void CLHEP::HepEulerAngles::HepEulerAngles(double, double, double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:48:10
    t.constructor<double, double, double>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::getPhi() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::getPhi()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:56:19
    t.method("getPhi", static_cast<double (CLHEP::HepEulerAngles::*)()  const>(&CLHEP::HepEulerAngles::getPhi));

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::phi() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::phi()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:57:19
    t.method("phi", static_cast<double (CLHEP::HepEulerAngles::*)()  const>(&CLHEP::HepEulerAngles::phi));


    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::getTheta() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::getTheta()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:60:19
    t.method("getTheta", static_cast<double (CLHEP::HepEulerAngles::*)()  const>(&CLHEP::HepEulerAngles::getTheta));

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::theta() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::theta()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:61:19
    t.method("theta", static_cast<double (CLHEP::HepEulerAngles::*)()  const>(&CLHEP::HepEulerAngles::theta));


    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::getPsi() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::getPsi()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:64:19
    t.method("getPsi", static_cast<double (CLHEP::HepEulerAngles::*)()  const>(&CLHEP::HepEulerAngles::getPsi));

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::psi() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::psi()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:65:19
    t.method("psi", static_cast<double (CLHEP::HepEulerAngles::*)()  const>(&CLHEP::HepEulerAngles::psi));



    module_.set_override_module(jl_base_module);







    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::getTolerance() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::getTolerance()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:83:24
    module_.method("CLHEP!HepEulerAngles!getTolerance", static_cast<double (*)() >(&CLHEP::HepEulerAngles::getTolerance));

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::setTolerance(double) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::setTolerance(double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:84:24
    module_.method("CLHEP!HepEulerAngles!setTolerance", static_cast<double (*)(double) >(&CLHEP::HepEulerAngles::setTolerance));


  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<CLHEP::HepEulerAngles>> type_;
};
std::shared_ptr<Wrapper> newJlCLHEP_HepEulerAngles(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlCLHEP_HepEulerAngles(module));
}
