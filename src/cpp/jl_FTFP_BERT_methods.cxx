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
void add_methods_for_FTFP_BERT(jlcxx::Module& types, jlcxx::TypeWrapper<FTFP_BERT>& t175) {


  /**********************************************************************/
  /* Wrappers for the methods of class FTFP_BERT
   */


  DEBUG_MSG("Adding wrapper for void FTFP_BERT::FTFP_BERT(G4int) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/FTFP_BERT.hh:49:3
  t175.constructor<G4int>(/*finalize=*/true);

  /* End of FTFP_BERT class method wrappers
   **********************************************************************/

}
