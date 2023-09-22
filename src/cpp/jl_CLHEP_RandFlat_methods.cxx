// this file was auto-generated by wrapit v0.1.0-54-g4322429
#include <type_traits>
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

#include "cpp/jlGeant4.h"


#ifdef VERBOSE_IMPORT
#  define DEBUG_MSG(a) std::cerr << a << "\n"
#else
#  define DEBUG_MSG(a)
#endif
#define __HERE__  __FILE__ ":" QUOTE2(__LINE__)
#define QUOTE(arg) #arg
#define QUOTE2(arg) QUOTE(arg)
void add_methods_for_CLHEP_RandFlat(jlcxx::Module& types, jlcxx::TypeWrapper<CLHEP::RandFlat>& t128) {


  /**********************************************************************/
  /* Wrappers for the methods of class CLHEP::RandFlat
   */


  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::RandFlat(CLHEP::HepRandomEngine &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:45:10
  t128.constructor<CLHEP::HepRandomEngine &>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::RandFlat(CLHEP::HepRandomEngine &, double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:46:10
  t128.constructor<CLHEP::HepRandomEngine &, double>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::RandFlat(CLHEP::HepRandomEngine &, double, double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:47:10
  t128.constructor<CLHEP::HepRandomEngine &, double, double>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::RandFlat(CLHEP::HepRandomEngine *) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:48:10
  t128.constructor<CLHEP::HepRandomEngine *>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::RandFlat(CLHEP::HepRandomEngine *, double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:49:10
  t128.constructor<CLHEP::HepRandomEngine *, double>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::RandFlat(CLHEP::HepRandomEngine *, double, double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:50:10
  t128.constructor<CLHEP::HepRandomEngine *, double, double>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::shoot() (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::shoot()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:65:18
  t128.method("CLHEP!RandFlat!shoot", static_cast<double (*)() >(&CLHEP::RandFlat::shoot));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::shoot(double) (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::shoot(double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:67:25
  t128.method("CLHEP!RandFlat!shoot", static_cast<double (*)(double) >(&CLHEP::RandFlat::shoot));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::shoot(double, double) (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::shoot(double, double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:69:25
  t128.method("CLHEP!RandFlat!shoot", static_cast<double (*)(double, double) >(&CLHEP::RandFlat::shoot));

  DEBUG_MSG("Adding wrapper for long CLHEP::RandFlat::shootInt(long) (" __HERE__ ")");
  // signature to use in the veto list: long CLHEP::RandFlat::shootInt(long)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:71:23
  t128.method("CLHEP!RandFlat!shootInt", static_cast<long (*)(long) >(&CLHEP::RandFlat::shootInt));

  DEBUG_MSG("Adding wrapper for long CLHEP::RandFlat::shootInt(long, long) (" __HERE__ ")");
  // signature to use in the veto list: long CLHEP::RandFlat::shootInt(long, long)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:73:23
  t128.method("CLHEP!RandFlat!shootInt", static_cast<long (*)(long, long) >(&CLHEP::RandFlat::shootInt));

  DEBUG_MSG("Adding wrapper for int CLHEP::RandFlat::shootBit() (" __HERE__ ")");
  // signature to use in the veto list: int CLHEP::RandFlat::shootBit()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:75:22
  t128.method("CLHEP!RandFlat!shootBit", static_cast<int (*)() >(&CLHEP::RandFlat::shootBit));

  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::shootArray(const int, double *) (" __HERE__ ")");
  // signature to use in the veto list: void CLHEP::RandFlat::shootArray(const int, double *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:77:16
  t128.method("CLHEP!RandFlat!shootArray", static_cast<void (*)(const int, double *) >(&CLHEP::RandFlat::shootArray));

  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::shootArray(const int, double *, double, double) (" __HERE__ ")");
  // signature to use in the veto list: void CLHEP::RandFlat::shootArray(const int, double *, double, double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:79:16
  t128.method("CLHEP!RandFlat!shootArray", static_cast<void (*)(const int, double *, double, double) >(&CLHEP::RandFlat::shootArray));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::shoot(CLHEP::HepRandomEngine *) (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::shoot(CLHEP::HepRandomEngine *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:85:25
  t128.method("CLHEP!RandFlat!shoot", static_cast<double (*)(CLHEP::HepRandomEngine *) >(&CLHEP::RandFlat::shoot));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::shoot(CLHEP::HepRandomEngine *, double) (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::shoot(CLHEP::HepRandomEngine *, double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:87:25
  t128.method("CLHEP!RandFlat!shoot", static_cast<double (*)(CLHEP::HepRandomEngine *, double) >(&CLHEP::RandFlat::shoot));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::shoot(CLHEP::HepRandomEngine *, double, double) (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::shoot(CLHEP::HepRandomEngine *, double, double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:89:25
  t128.method("CLHEP!RandFlat!shoot", static_cast<double (*)(CLHEP::HepRandomEngine *, double, double) >(&CLHEP::RandFlat::shoot));

  DEBUG_MSG("Adding wrapper for long CLHEP::RandFlat::shootInt(CLHEP::HepRandomEngine *, long) (" __HERE__ ")");
  // signature to use in the veto list: long CLHEP::RandFlat::shootInt(CLHEP::HepRandomEngine *, long)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:91:23
  t128.method("CLHEP!RandFlat!shootInt", static_cast<long (*)(CLHEP::HepRandomEngine *, long) >(&CLHEP::RandFlat::shootInt));

  DEBUG_MSG("Adding wrapper for long CLHEP::RandFlat::shootInt(CLHEP::HepRandomEngine *, long, long) (" __HERE__ ")");
  // signature to use in the veto list: long CLHEP::RandFlat::shootInt(CLHEP::HepRandomEngine *, long, long)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:93:23
  t128.method("CLHEP!RandFlat!shootInt", static_cast<long (*)(CLHEP::HepRandomEngine *, long, long) >(&CLHEP::RandFlat::shootInt));

  DEBUG_MSG("Adding wrapper for int CLHEP::RandFlat::shootBit(CLHEP::HepRandomEngine *) (" __HERE__ ")");
  // signature to use in the veto list: int CLHEP::RandFlat::shootBit(CLHEP::HepRandomEngine *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:95:22
  t128.method("CLHEP!RandFlat!shootBit", static_cast<int (*)(CLHEP::HepRandomEngine *) >(&CLHEP::RandFlat::shootBit));

  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::shootArray(CLHEP::HepRandomEngine *, const int, double *) (" __HERE__ ")");
  // signature to use in the veto list: void CLHEP::RandFlat::shootArray(CLHEP::HepRandomEngine *, const int, double *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:97:23
  t128.method("CLHEP!RandFlat!shootArray", static_cast<void (*)(CLHEP::HepRandomEngine *, const int, double *) >(&CLHEP::RandFlat::shootArray));

  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::shootArray(CLHEP::HepRandomEngine *, const int, double *, double, double) (" __HERE__ ")");
  // signature to use in the veto list: void CLHEP::RandFlat::shootArray(CLHEP::HepRandomEngine *, const int, double *, double, double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:100:16
  t128.method("CLHEP!RandFlat!shootArray", static_cast<void (*)(CLHEP::HepRandomEngine *, const int, double *, double, double) >(&CLHEP::RandFlat::shootArray));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::fire() (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::fire()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:107:17
  t128.method("fire", static_cast<double (CLHEP::RandFlat::*)() >(&CLHEP::RandFlat::fire));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::fire(double) (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::fire(double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:109:17
  t128.method("fire", static_cast<double (CLHEP::RandFlat::*)(double) >(&CLHEP::RandFlat::fire));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::fire(double, double) (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::fire(double, double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:111:17
  t128.method("fire", static_cast<double (CLHEP::RandFlat::*)(double, double) >(&CLHEP::RandFlat::fire));

  DEBUG_MSG("Adding wrapper for long CLHEP::RandFlat::fireInt(long) (" __HERE__ ")");
  // signature to use in the veto list: long CLHEP::RandFlat::fireInt(long)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:113:15
  t128.method("fireInt", static_cast<long (CLHEP::RandFlat::*)(long) >(&CLHEP::RandFlat::fireInt));

  DEBUG_MSG("Adding wrapper for long CLHEP::RandFlat::fireInt(long, long) (" __HERE__ ")");
  // signature to use in the veto list: long CLHEP::RandFlat::fireInt(long, long)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:115:15
  t128.method("fireInt", static_cast<long (CLHEP::RandFlat::*)(long, long) >(&CLHEP::RandFlat::fireInt));

  DEBUG_MSG("Adding wrapper for int CLHEP::RandFlat::fireBit() (" __HERE__ ")");
  // signature to use in the veto list: int CLHEP::RandFlat::fireBit()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:117:14
  t128.method("fireBit", static_cast<int (CLHEP::RandFlat::*)() >(&CLHEP::RandFlat::fireBit));

  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::fireArray(const int, double *) (" __HERE__ ")");
  // signature to use in the veto list: void CLHEP::RandFlat::fireArray(const int, double *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:119:8
  t128.method("fireArray", static_cast<void (CLHEP::RandFlat::*)(const int, double *) >(&CLHEP::RandFlat::fireArray));

  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::fireArray(const int, double *, double, double) (" __HERE__ ")");
  // signature to use in the veto list: void CLHEP::RandFlat::fireArray(const int, double *, double, double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:121:8
  t128.method("fireArray", static_cast<void (CLHEP::RandFlat::*)(const int, double *, double, double) >(&CLHEP::RandFlat::fireArray));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::operator()() (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::operator()()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:124:10
  t128.method("paren", static_cast<double (CLHEP::RandFlat::*)() >(&CLHEP::RandFlat::operator()));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::operator()(double) (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::operator()(double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:125:10
  t128.method("paren", static_cast<double (CLHEP::RandFlat::*)(double) >(&CLHEP::RandFlat::operator()));

  DEBUG_MSG("Adding wrapper for double CLHEP::RandFlat::operator()(double, double) (" __HERE__ ")");
  // signature to use in the veto list: double CLHEP::RandFlat::operator()(double, double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:126:10
  t128.method("paren", static_cast<double (CLHEP::RandFlat::*)(double, double) >(&CLHEP::RandFlat::operator()));

  DEBUG_MSG("Adding wrapper for CLHEP::HepRandomEngine & CLHEP::RandFlat::engine() (" __HERE__ ")");
  // signature to use in the veto list: CLHEP::HepRandomEngine & CLHEP::RandFlat::engine()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:134:21
  t128.method("engine", static_cast<CLHEP::HepRandomEngine & (CLHEP::RandFlat::*)() >(&CLHEP::RandFlat::engine));

  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::saveEngineStatus(const char []) (" __HERE__ ")");
  // signature to use in the veto list: void CLHEP::RandFlat::saveEngineStatus(const char [])
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:144:15
  t128.method("CLHEP!RandFlat!saveEngineStatus", static_cast<void (*)(const char []) >(&CLHEP::RandFlat::saveEngineStatus));
  t128.method("CLHEP!RandFlat!saveEngineStatus", []()->void{ CLHEP::RandFlat::saveEngineStatus(); });

  DEBUG_MSG("Adding wrapper for void CLHEP::RandFlat::restoreEngineStatus(const char []) (" __HERE__ ")");
  // signature to use in the veto list: void CLHEP::RandFlat::restoreEngineStatus(const char [])
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/CLHEP/Random/RandFlat.h:147:15
  t128.method("CLHEP!RandFlat!restoreEngineStatus", static_cast<void (*)(const char []) >(&CLHEP::RandFlat::restoreEngineStatus));
  t128.method("CLHEP!RandFlat!restoreEngineStatus", []()->void{ CLHEP::RandFlat::restoreEngineStatus(); });

  /* End of CLHEP::RandFlat class method wrappers
   **********************************************************************/

}