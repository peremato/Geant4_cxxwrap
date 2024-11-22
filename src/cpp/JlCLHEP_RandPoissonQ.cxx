// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<CLHEP::RandPoissonQ> : std::false_type { };
  template<> struct DefaultConstructible<CLHEP::RandPoissonQ> : std::false_type { };
}

// Class generating the wrapper for type CLHEP::RandPoissonQ
// signature to use in the veto file: CLHEP::RandPoissonQ
struct JlCLHEP_RandPoissonQ: public Wrapper {

  JlCLHEP_RandPoissonQ(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type CLHEP::RandPoissonQ (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:31:7
    jlcxx::TypeWrapper<CLHEP::RandPoissonQ>  t = jlModule.add_type<CLHEP::RandPoissonQ>("CLHEP!RandPoissonQ");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<CLHEP::RandPoissonQ>>(new jlcxx::TypeWrapper<CLHEP::RandPoissonQ>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void CLHEP::RandPoissonQ::RandPoissonQ(CLHEP::HepRandomEngine &, double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:35:10
    t.constructor<CLHEP::HepRandomEngine &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<CLHEP::HepRandomEngine &, double>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void CLHEP::RandPoissonQ::RandPoissonQ(CLHEP::HepRandomEngine *, double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:36:10
    t.constructor<CLHEP::HepRandomEngine *>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<CLHEP::HepRandomEngine *, double>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for long CLHEP::RandPoissonQ::shoot(double) (" __HERE__ ")");
    // signature to use in the veto list: long CLHEP::RandPoissonQ::shoot(double)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:74:16
    module_.method("CLHEP!RandPoissonQ!shoot", static_cast<long (*)(double) >(&CLHEP::RandPoissonQ::shoot));
    module_.method("CLHEP!RandPoissonQ!shoot", []()->long { return CLHEP::RandPoissonQ::shoot(); });

    DEBUG_MSG("Adding wrapper for void CLHEP::RandPoissonQ::shootArray(const int, long *, double) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandPoissonQ::shootArray(const int, long *, double)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:76:16
    module_.method("CLHEP!RandPoissonQ!shootArray", static_cast<void (*)(const int, long *, double) >(&CLHEP::RandPoissonQ::shootArray));
    module_.method("CLHEP!RandPoissonQ!shootArray", [](const int arg0, long * arg1)->void { CLHEP::RandPoissonQ::shootArray(arg0, arg1); });

    DEBUG_MSG("Adding wrapper for long CLHEP::RandPoissonQ::shoot(CLHEP::HepRandomEngine *, double) (" __HERE__ ")");
    // signature to use in the veto list: long CLHEP::RandPoissonQ::shoot(CLHEP::HepRandomEngine *, double)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:81:16
    module_.method("CLHEP!RandPoissonQ!shoot", static_cast<long (*)(CLHEP::HepRandomEngine *, double) >(&CLHEP::RandPoissonQ::shoot));
    module_.method("CLHEP!RandPoissonQ!shoot", [](CLHEP::HepRandomEngine * arg0)->long { return CLHEP::RandPoissonQ::shoot(arg0); });

    DEBUG_MSG("Adding wrapper for long CLHEP::RandPoissonQ::fire() (" __HERE__ ")");
    // signature to use in the veto list: long CLHEP::RandPoissonQ::fire()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:89:9
    t.method("fire", static_cast<long (CLHEP::RandPoissonQ::*)() >(&CLHEP::RandPoissonQ::fire));

    DEBUG_MSG("Adding wrapper for long CLHEP::RandPoissonQ::fire(double) (" __HERE__ ")");
    // signature to use in the veto list: long CLHEP::RandPoissonQ::fire(double)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:90:9
    t.method("fire", static_cast<long (CLHEP::RandPoissonQ::*)(double) >(&CLHEP::RandPoissonQ::fire));

    DEBUG_MSG("Adding wrapper for void CLHEP::RandPoissonQ::fireArray(const int, long *) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandPoissonQ::fireArray(const int, long *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:92:8
    t.method("fireArray", static_cast<void (CLHEP::RandPoissonQ::*)(const int, long *) >(&CLHEP::RandPoissonQ::fireArray));

    DEBUG_MSG("Adding wrapper for void CLHEP::RandPoissonQ::fireArray(const int, long *, double) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandPoissonQ::fireArray(const int, long *, double)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:93:8
    t.method("fireArray", static_cast<void (CLHEP::RandPoissonQ::*)(const int, long *, double) >(&CLHEP::RandPoissonQ::fireArray));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandPoissonQ::operator()() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandPoissonQ::operator()()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:95:10
    t.method("paren", static_cast<double (CLHEP::RandPoissonQ::*)() >(&CLHEP::RandPoissonQ::operator()));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandPoissonQ::operator()(double) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandPoissonQ::operator()(double)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:96:10
    t.method("paren", static_cast<double (CLHEP::RandPoissonQ::*)(double) >(&CLHEP::RandPoissonQ::operator()));

    DEBUG_MSG("Adding wrapper for std::string CLHEP::RandPoissonQ::name() (" __HERE__ ")");
    // signature to use in the veto list: std::string CLHEP::RandPoissonQ::name()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:98:15
    t.method("name", static_cast<std::string (CLHEP::RandPoissonQ::*)()  const>(&CLHEP::RandPoissonQ::name));

    DEBUG_MSG("Adding wrapper for CLHEP::HepRandomEngine & CLHEP::RandPoissonQ::engine() (" __HERE__ ")");
    // signature to use in the veto list: CLHEP::HepRandomEngine & CLHEP::RandPoissonQ::engine()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:99:21
    t.method("engine", static_cast<CLHEP::HepRandomEngine & (CLHEP::RandPoissonQ::*)() >(&CLHEP::RandPoissonQ::engine));

    DEBUG_MSG("Adding wrapper for std::string CLHEP::RandPoissonQ::distributionName() (" __HERE__ ")");
    // signature to use in the veto list: std::string CLHEP::RandPoissonQ::distributionName()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:101:22
    module_.method("CLHEP!RandPoissonQ!distributionName", static_cast<std::string (*)() >(&CLHEP::RandPoissonQ::distributionName));

    DEBUG_MSG("Adding wrapper for int CLHEP::RandPoissonQ::tableBoundary() (" __HERE__ ")");
    // signature to use in the veto list: int CLHEP::RandPoissonQ::tableBoundary()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/RandPoissonQ.h:110:21
    module_.method("CLHEP!RandPoissonQ!tableBoundary", static_cast<int (*)() >(&CLHEP::RandPoissonQ::tableBoundary));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<CLHEP::RandPoissonQ>> type_;
};
std::shared_ptr<Wrapper> newJlCLHEP_RandPoissonQ(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlCLHEP_RandPoissonQ(module));
}
