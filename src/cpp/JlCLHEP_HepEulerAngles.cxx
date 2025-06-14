// this file was auto-generated by wrapit v1.6.0
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
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes    );


    DEBUG_MSG("Adding wrapper for void CLHEP::HepEulerAngles::HepEulerAngles(double, double, double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:48:10
    t.constructor<double, double, double>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("phi"), jlcxx::arg("theta"), jlcxx::arg("psi")    );

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::getPhi() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::getPhi()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:56:19
    t.method("getPhi", [](CLHEP::HepEulerAngles const& a)->double { return a.getPhi(); }, jlcxx::arg("this"));
    t.method("getPhi", [](CLHEP::HepEulerAngles const* a)->double { return a->getPhi(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::phi() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::phi()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:57:19
    t.method("phi", [](CLHEP::HepEulerAngles const& a)->double { return a.phi(); }, jlcxx::arg("this"));
    t.method("phi", [](CLHEP::HepEulerAngles const* a)->double { return a->phi(); }, jlcxx::arg("this"));


    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::getTheta() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::getTheta()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:60:19
    t.method("getTheta", [](CLHEP::HepEulerAngles const& a)->double { return a.getTheta(); }, jlcxx::arg("this"));
    t.method("getTheta", [](CLHEP::HepEulerAngles const* a)->double { return a->getTheta(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::theta() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::theta()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:61:19
    t.method("theta", [](CLHEP::HepEulerAngles const& a)->double { return a.theta(); }, jlcxx::arg("this"));
    t.method("theta", [](CLHEP::HepEulerAngles const* a)->double { return a->theta(); }, jlcxx::arg("this"));


    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::getPsi() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::getPsi()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:64:19
    t.method("getPsi", [](CLHEP::HepEulerAngles const& a)->double { return a.getPsi(); }, jlcxx::arg("this"));
    t.method("getPsi", [](CLHEP::HepEulerAngles const* a)->double { return a->getPsi(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::psi() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::psi()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:65:19
    t.method("psi", [](CLHEP::HepEulerAngles const& a)->double { return a.psi(); }, jlcxx::arg("this"));
    t.method("psi", [](CLHEP::HepEulerAngles const* a)->double { return a->psi(); }, jlcxx::arg("this"));



    module_.set_override_module(jl_base_module);







    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::getTolerance() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::getTolerance()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:83:24
    module_.method("CLHEP!HepEulerAngles!getTolerance", []()->double { return CLHEP::HepEulerAngles::getTolerance(); });

    DEBUG_MSG("Adding wrapper for double CLHEP::HepEulerAngles::setTolerance(double) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepEulerAngles::setTolerance(double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Vector/EulerAngles.h:84:24
    module_.method("CLHEP!HepEulerAngles!setTolerance", [](double arg0)->double { return CLHEP::HepEulerAngles::setTolerance(arg0); }, jlcxx::arg("tol"));


  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<CLHEP::HepEulerAngles>> type_;
};
std::shared_ptr<Wrapper> newJlCLHEP_HepEulerAngles(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlCLHEP_HepEulerAngles(module));
}
