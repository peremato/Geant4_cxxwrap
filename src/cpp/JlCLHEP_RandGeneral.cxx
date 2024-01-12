// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<CLHEP::RandGeneral> : std::false_type { };
  template<> struct DefaultConstructible<CLHEP::RandGeneral> : std::false_type { };
template<> struct SuperType<CLHEP::RandGeneral> { typedef CLHEP::HepRandom type; };
}

// Class generating the wrapper for type CLHEP::RandGeneral
// signature to use in the veto file: CLHEP::RandGeneral
struct JlCLHEP_RandGeneral: public Wrapper {

  JlCLHEP_RandGeneral(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type CLHEP::RandGeneral (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:38:7
    jlcxx::TypeWrapper<CLHEP::RandGeneral>  t = jlModule.add_type<CLHEP::RandGeneral>("CLHEP!RandGeneral",
      jlcxx::julia_base_type<CLHEP::HepRandom>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<CLHEP::RandGeneral>>(new jlcxx::TypeWrapper<CLHEP::RandGeneral>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void CLHEP::RandGeneral::RandGeneral(const double *, int, int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:42:3
    t.constructor<const double *, int>(/*finalize=*/true);
    t.constructor<const double *, int, int>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void CLHEP::RandGeneral::RandGeneral(CLHEP::HepRandomEngine &, const double *, int, int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:45:3
    t.constructor<CLHEP::HepRandomEngine &, const double *, int>(/*finalize=*/true);
    t.constructor<CLHEP::HepRandomEngine &, const double *, int, int>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void CLHEP::RandGeneral::RandGeneral(CLHEP::HepRandomEngine *, const double *, int, int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:49:3
    t.constructor<CLHEP::HepRandomEngine *, const double *, int>(/*finalize=*/true);
    t.constructor<CLHEP::HepRandomEngine *, const double *, int, int>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGeneral::shoot() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGeneral::shoot()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:95:17
    t.method("shoot", static_cast<double (CLHEP::RandGeneral::*)() >(&CLHEP::RandGeneral::shoot));

    DEBUG_MSG("Adding wrapper for void CLHEP::RandGeneral::shootArray(const int, double *) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandGeneral::shootArray(const int, double *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:97:15
    t.method("shootArray", static_cast<void (CLHEP::RandGeneral::*)(const int, double *) >(&CLHEP::RandGeneral::shootArray));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGeneral::shoot(CLHEP::HepRandomEngine *) (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGeneral::shoot(CLHEP::HepRandomEngine *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:102:10
    t.method("shoot", static_cast<double (CLHEP::RandGeneral::*)(CLHEP::HepRandomEngine *) >(&CLHEP::RandGeneral::shoot));

    DEBUG_MSG("Adding wrapper for void CLHEP::RandGeneral::shootArray(CLHEP::HepRandomEngine *, const int, double *) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandGeneral::shootArray(CLHEP::HepRandomEngine *, const int, double *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:104:8
    t.method("shootArray", static_cast<void (CLHEP::RandGeneral::*)(CLHEP::HepRandomEngine *, const int, double *) >(&CLHEP::RandGeneral::shootArray));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGeneral::fire() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGeneral::fire()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:110:10
    t.method("fire", static_cast<double (CLHEP::RandGeneral::*)() >(&CLHEP::RandGeneral::fire));

    DEBUG_MSG("Adding wrapper for void CLHEP::RandGeneral::fireArray(const int, double *) (" __HERE__ ")");
    // signature to use in the veto list: void CLHEP::RandGeneral::fireArray(const int, double *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:112:8
    t.method("fireArray", static_cast<void (CLHEP::RandGeneral::*)(const int, double *) >(&CLHEP::RandGeneral::fireArray));

    DEBUG_MSG("Adding wrapper for double CLHEP::RandGeneral::operator()() (" __HERE__ ")");
    // signature to use in the veto list: double CLHEP::RandGeneral::operator()()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:114:10
    t.method("paren", static_cast<double (CLHEP::RandGeneral::*)() >(&CLHEP::RandGeneral::operator()));

    DEBUG_MSG("Adding wrapper for std::string CLHEP::RandGeneral::name() (" __HERE__ ")");
    // signature to use in the veto list: std::string CLHEP::RandGeneral::name()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:121:15
    t.method("name", static_cast<std::string (CLHEP::RandGeneral::*)()  const>(&CLHEP::RandGeneral::name));

    DEBUG_MSG("Adding wrapper for CLHEP::HepRandomEngine & CLHEP::RandGeneral::engine() (" __HERE__ ")");
    // signature to use in the veto list: CLHEP::HepRandomEngine & CLHEP::RandGeneral::engine()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:122:21
    t.method("engine", static_cast<CLHEP::HepRandomEngine & (CLHEP::RandGeneral::*)() >(&CLHEP::RandGeneral::engine));

    DEBUG_MSG("Adding wrapper for std::string CLHEP::RandGeneral::distributionName() (" __HERE__ ")");
    // signature to use in the veto list: std::string CLHEP::RandGeneral::distributionName()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Random/RandGeneral.h:124:22
    t.method("CLHEP!RandGeneral!distributionName", static_cast<std::string (*)() >(&CLHEP::RandGeneral::distributionName));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<CLHEP::RandGeneral>> type_;
};
std::shared_ptr<Wrapper> newJlCLHEP_RandGeneral(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlCLHEP_RandGeneral(module));
}
