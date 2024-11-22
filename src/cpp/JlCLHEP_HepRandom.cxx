// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<CLHEP::HepRandom> : std::false_type { };
  template<> struct DefaultConstructible<CLHEP::HepRandom> : std::false_type { };
}

// Class generating the wrapper for type CLHEP::HepRandom
// signature to use in the veto file: CLHEP::HepRandom
struct JlCLHEP_HepRandom: public Wrapper {

  JlCLHEP_HepRandom(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type CLHEP::HepRandom (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:48:7
    jlcxx::TypeWrapper<CLHEP::HepRandom>  t = jlModule.add_type<CLHEP::HepRandom>("CLHEP!HepRandom");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<CLHEP::HepRandom>>(new jlcxx::TypeWrapper<CLHEP::HepRandom>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::HepRandom(long) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:53:3
    t.constructor<long>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::HepRandom(CLHEP::HepRandomEngine &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:57:3
    t.constructor<CLHEP::HepRandomEngine &>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::HepRandom(CLHEP::HepRandomEngine *) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:58:3
    t.constructor<CLHEP::HepRandomEngine *>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for double CLHEP::HepRandom::flat() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepRandom::flat()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:68:10
    t.method("flat", static_cast<double (CLHEP::HepRandom::*)() >(&CLHEP::HepRandom::flat));

    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::flatArray(const int, double *) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::HepRandom::flatArray(const int, double *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:71:8
    t.method("flatArray", static_cast<void (CLHEP::HepRandom::*)(const int, double *) >(&CLHEP::HepRandom::flatArray));

    DEBUG_MSG("Adding wrapper for double CLHEP::HepRandom::flat(CLHEP::HepRandomEngine *) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepRandom::flat(CLHEP::HepRandomEngine *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:74:17
    t.method("flat", static_cast<double (CLHEP::HepRandom::*)(CLHEP::HepRandomEngine *) >(&CLHEP::HepRandom::flat));

    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::flatArray(CLHEP::HepRandomEngine *, const int, double *) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::HepRandom::flatArray(CLHEP::HepRandomEngine *, const int, double *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:77:15
    t.method("flatArray", static_cast<void (CLHEP::HepRandom::*)(CLHEP::HepRandomEngine *, const int, double *) >(&CLHEP::HepRandom::flatArray));

    DEBUG_MSG("Adding wrapper for double CLHEP::HepRandom::operator()() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::HepRandom::operator()()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:82:18
    t.method("paren", static_cast<double (CLHEP::HepRandom::*)() >(&CLHEP::HepRandom::operator()));

    DEBUG_MSG("Adding wrapper for std::string CLHEP::HepRandom::name() (" __HERE__ ")");
    // signature to use in the veto list: std::string CLHEP::HepRandom::name()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:85:23
    t.method("name", static_cast<std::string (CLHEP::HepRandom::*)()  const>(&CLHEP::HepRandom::name));

    DEBUG_MSG("Adding wrapper for CLHEP::HepRandomEngine & CLHEP::HepRandom::engine() (" __HERE__ ")");
    // signature to use in the veto list: CLHEP::HepRandomEngine & CLHEP::HepRandom::engine()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:86:29
    t.method("engine", static_cast<CLHEP::HepRandomEngine & (CLHEP::HepRandom::*)() >(&CLHEP::HepRandom::engine));

    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::setTheSeed(long, int) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::HepRandom::setTheSeed(long, int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:97:15
    module_.method("CLHEP!HepRandom!setTheSeed", static_cast<void (*)(long, int) >(&CLHEP::HepRandom::setTheSeed));
    module_.method("CLHEP!HepRandom!setTheSeed", [](long arg0)->void { CLHEP::HepRandom::setTheSeed(arg0); });

    DEBUG_MSG("Adding wrapper for long CLHEP::HepRandom::getTheSeed() (" __HERE__ ")");
    // signature to use in the veto list: long CLHEP::HepRandom::getTheSeed()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:100:15
    module_.method("CLHEP!HepRandom!getTheSeed", static_cast<long (*)() >(&CLHEP::HepRandom::getTheSeed));

    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::setTheSeeds(const long *, int) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::HepRandom::setTheSeeds(const long *, int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:103:15
    module_.method("CLHEP!HepRandom!setTheSeeds", static_cast<void (*)(const long *, int) >(&CLHEP::HepRandom::setTheSeeds));
    module_.method("CLHEP!HepRandom!setTheSeeds", [](const long * arg0)->void { CLHEP::HepRandom::setTheSeeds(arg0); });

    DEBUG_MSG("Adding wrapper for const long * CLHEP::HepRandom::getTheSeeds() (" __HERE__ ")");
    // signature to use in the veto list: const long * CLHEP::HepRandom::getTheSeeds()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:106:22
    module_.method("CLHEP!HepRandom!getTheSeeds", static_cast<const long * (*)() >(&CLHEP::HepRandom::getTheSeeds));

    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::getTheTableSeeds(long *, int) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::HepRandom::getTheTableSeeds(long *, int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:109:15
    module_.method("CLHEP!HepRandom!getTheTableSeeds", static_cast<void (*)(long *, int) >(&CLHEP::HepRandom::getTheTableSeeds));

    DEBUG_MSG("Adding wrapper for CLHEP::HepRandom * CLHEP::HepRandom::getTheGenerator() (" __HERE__ ")");
    // signature to use in the veto list: CLHEP::HepRandom * CLHEP::HepRandom::getTheGenerator()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:112:22
    module_.method("CLHEP!HepRandom!getTheGenerator", static_cast<CLHEP::HepRandom * (*)() >(&CLHEP::HepRandom::getTheGenerator));

    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::setTheEngine(CLHEP::HepRandomEngine *) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::HepRandom::setTheEngine(CLHEP::HepRandomEngine *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:115:15
    module_.method("CLHEP!HepRandom!setTheEngine", static_cast<void (*)(CLHEP::HepRandomEngine *) >(&CLHEP::HepRandom::setTheEngine));

    DEBUG_MSG("Adding wrapper for CLHEP::HepRandomEngine * CLHEP::HepRandom::getTheEngine() (" __HERE__ ")");
    // signature to use in the veto list: CLHEP::HepRandomEngine * CLHEP::HepRandom::getTheEngine()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:118:28
    module_.method("CLHEP!HepRandom!getTheEngine", static_cast<CLHEP::HepRandomEngine * (*)() >(&CLHEP::HepRandom::getTheEngine));

    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::saveEngineStatus(const char []) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::HepRandom::saveEngineStatus(const char [])
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:121:15
    module_.method("CLHEP!HepRandom!saveEngineStatus", static_cast<void (*)(const char []) >(&CLHEP::HepRandom::saveEngineStatus));
    module_.method("CLHEP!HepRandom!saveEngineStatus", []()->void { CLHEP::HepRandom::saveEngineStatus(); });

    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::restoreEngineStatus(const char []) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::HepRandom::restoreEngineStatus(const char [])
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:124:15
    module_.method("CLHEP!HepRandom!restoreEngineStatus", static_cast<void (*)(const char []) >(&CLHEP::HepRandom::restoreEngineStatus));
    module_.method("CLHEP!HepRandom!restoreEngineStatus", []()->void { CLHEP::HepRandom::restoreEngineStatus(); });

    DEBUG_MSG("Adding wrapper for void CLHEP::HepRandom::showEngineStatus() (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::HepRandom::showEngineStatus()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:145:15
    module_.method("CLHEP!HepRandom!showEngineStatus", static_cast<void (*)() >(&CLHEP::HepRandom::showEngineStatus));

    DEBUG_MSG("Adding wrapper for int CLHEP::HepRandom::createInstance() (" __HERE__ ")");
    // signature to use in the veto list: int CLHEP::HepRandom::createInstance()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:148:14
    module_.method("CLHEP!HepRandom!createInstance", static_cast<int (*)() >(&CLHEP::HepRandom::createInstance));

    DEBUG_MSG("Adding wrapper for std::string CLHEP::HepRandom::distributionName() (" __HERE__ ")");
    // signature to use in the veto list: std::string CLHEP::HepRandom::distributionName()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/CLHEP/Random/Random.h:151:22
    module_.method("CLHEP!HepRandom!distributionName", static_cast<std::string (*)() >(&CLHEP::HepRandom::distributionName));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<CLHEP::HepRandom>> type_;
};
std::shared_ptr<Wrapper> newJlCLHEP_HepRandom(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlCLHEP_HepRandom(module));
}
