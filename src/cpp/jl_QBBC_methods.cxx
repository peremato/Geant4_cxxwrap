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
void add_methods_for_QBBC(jlcxx::Module& types, jlcxx::TypeWrapper<QBBC>& t190) {


  /**********************************************************************/
  /* Wrappers for the methods of class QBBC
   */


  DEBUG_MSG("Adding wrapper for void QBBC::QBBC(G4int, const G4String &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/QBBC.hh:47:12
  t190.constructor<G4int>(/*finalize=*/true);
  t190.constructor<G4int, const G4String &>(/*finalize=*/true);

  /* End of QBBC class method wrappers
   **********************************************************************/

}
