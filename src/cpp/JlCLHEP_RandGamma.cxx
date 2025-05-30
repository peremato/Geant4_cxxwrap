// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<CLHEP::RandGamma> : std::false_type { };
  template<> struct DefaultConstructible<CLHEP::RandGamma> : std::false_type { };
template<> struct SuperType<CLHEP::RandGamma> { typedef CLHEP::HepRandom type; };
}

// Class generating the wrapper for type CLHEP::RandGamma
// signature to use in the veto file: CLHEP::RandGamma
struct JlCLHEP_RandGamma: public Wrapper {

  JlCLHEP_RandGamma(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type CLHEP::RandGamma (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:36:7
    jlcxx::TypeWrapper<CLHEP::RandGamma>  t = jlModule.add_type<CLHEP::RandGamma>("CLHEP!RandGamma",
      jlcxx::julia_base_type<CLHEP::HepRandom>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<CLHEP::RandGamma>>(new jlcxx::TypeWrapper<CLHEP::RandGamma>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void CLHEP::RandGamma::RandGamma(CLHEP::HepRandomEngine &, double, double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:40:10
    t.constructor<CLHEP::HepRandomEngine &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<CLHEP::HepRandomEngine &, double>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<CLHEP::HepRandomEngine &, double, double>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void CLHEP::RandGamma::RandGamma(CLHEP::HepRandomEngine *, double, double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:42:10
    t.constructor<CLHEP::HepRandomEngine *>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<CLHEP::HepRandomEngine *, double>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<CLHEP::HepRandomEngine *, double, double>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGamma::shoot() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGamma::shoot()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:58:24
    module_.method("CLHEP!RandGamma!shoot", static_cast<double (*)() >(&CLHEP::RandGamma::shoot));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGamma::shoot(double, double) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGamma::shoot(double, double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:60:17
    module_.method("CLHEP!RandGamma!shoot", static_cast<double (*)(double, double) >(&CLHEP::RandGamma::shoot));

    DEBUG_MSG("Adding wrapper for void CLHEP::RandGamma::shootArray(const int, double *, double, double) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandGamma::shootArray(const int, double *, double, double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:62:15
    module_.method("CLHEP!RandGamma!shootArray", static_cast<void (*)(const int, double *, double, double) >(&CLHEP::RandGamma::shootArray));
    module_.method("CLHEP!RandGamma!shootArray", [](const int arg0, double * arg1)->void { CLHEP::RandGamma::shootArray(arg0, arg1); });
    module_.method("CLHEP!RandGamma!shootArray", [](const int arg0, double * arg1, double arg2)->void { CLHEP::RandGamma::shootArray(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGamma::shoot(CLHEP::HepRandomEngine *) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGamma::shoot(CLHEP::HepRandomEngine *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:68:24
    module_.method("CLHEP!RandGamma!shoot", static_cast<double (*)(CLHEP::HepRandomEngine *) >(&CLHEP::RandGamma::shoot));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGamma::shoot(CLHEP::HepRandomEngine *, double, double) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGamma::shoot(CLHEP::HepRandomEngine *, double, double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:70:17
    module_.method("CLHEP!RandGamma!shoot", static_cast<double (*)(CLHEP::HepRandomEngine *, double, double) >(&CLHEP::RandGamma::shoot));

    DEBUG_MSG("Adding wrapper for void CLHEP::RandGamma::shootArray(CLHEP::HepRandomEngine *, const int, double *, double, double) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandGamma::shootArray(CLHEP::HepRandomEngine *, const int, double *, double, double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:73:15
    module_.method("CLHEP!RandGamma!shootArray", static_cast<void (*)(CLHEP::HepRandomEngine *, const int, double *, double, double) >(&CLHEP::RandGamma::shootArray));
    module_.method("CLHEP!RandGamma!shootArray", [](CLHEP::HepRandomEngine * arg0, const int arg1, double * arg2)->void { CLHEP::RandGamma::shootArray(arg0, arg1, arg2); });
    module_.method("CLHEP!RandGamma!shootArray", [](CLHEP::HepRandomEngine * arg0, const int arg1, double * arg2, double arg3)->void { CLHEP::RandGamma::shootArray(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGamma::fire() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGamma::fire()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:80:17
    t.method("fire", static_cast<double (CLHEP::RandGamma::*)() >(&CLHEP::RandGamma::fire));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGamma::fire(double, double) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGamma::fire(double, double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:82:10
    t.method("fire", static_cast<double (CLHEP::RandGamma::*)(double, double) >(&CLHEP::RandGamma::fire));

    DEBUG_MSG("Adding wrapper for void CLHEP::RandGamma::fireArray(const int, double *) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandGamma::fireArray(const int, double *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:84:8
    t.method("fireArray", static_cast<void (CLHEP::RandGamma::*)(const int, double *) >(&CLHEP::RandGamma::fireArray));

    DEBUG_MSG("Adding wrapper for void CLHEP::RandGamma::fireArray(const int, double *, double, double) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandGamma::fireArray(const int, double *, double, double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:85:8
    t.method("fireArray", static_cast<void (CLHEP::RandGamma::*)(const int, double *, double, double) >(&CLHEP::RandGamma::fireArray));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGamma::operator()() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGamma::operator()()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:87:17
    t.method("paren", static_cast<double (CLHEP::RandGamma::*)() >(&CLHEP::RandGamma::operator()));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGamma::operator()(double, double) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGamma::operator()(double, double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:88:17
    t.method("paren", static_cast<double (CLHEP::RandGamma::*)(double, double) >(&CLHEP::RandGamma::operator()));

    DEBUG_MSG("Adding wrapper for std::string CLHEP::RandGamma::name() (" __HERE__ ")");
    // signature to use in the veto list: std::string CLHEP::RandGamma::name()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:95:15
    t.method("name", static_cast<std::string (CLHEP::RandGamma::*)()  const>(&CLHEP::RandGamma::name));

    DEBUG_MSG("Adding wrapper for CLHEP::HepRandomEngine & CLHEP::RandGamma::engine() (" __HERE__ ")");
    // signature to use in the veto list: CLHEP::HepRandomEngine & CLHEP::RandGamma::engine()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:96:21
    t.method("engine", static_cast<CLHEP::HepRandomEngine & (CLHEP::RandGamma::*)() >(&CLHEP::RandGamma::engine));

    DEBUG_MSG("Adding wrapper for std::string CLHEP::RandGamma::distributionName() (" __HERE__ ")");
    // signature to use in the veto list: std::string CLHEP::RandGamma::distributionName()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/CLHEP/Random/RandGamma.h:98:22
    module_.method("CLHEP!RandGamma!distributionName", static_cast<std::string (*)() >(&CLHEP::RandGamma::distributionName));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<CLHEP::RandGamma>> type_;
};
std::shared_ptr<Wrapper> newJlCLHEP_RandGamma(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlCLHEP_RandGamma(module));
}
