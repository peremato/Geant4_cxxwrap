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
void add_methods_for_G4MTRunManagerKernel(jlcxx::Module& types, jlcxx::TypeWrapper<G4MTRunManagerKernel>& t140) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4MTRunManagerKernel
   */

  DEBUG_MSG("Adding wrapper for void G4MTRunManagerKernel::SetUpDecayChannels() (" __HERE__ ")");
  // signature to use in the veto list: void G4MTRunManagerKernel::SetUpDecayChannels()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MTRunManagerKernel.hh:79:10
  t140.method("SetUpDecayChannels", static_cast<void (G4MTRunManagerKernel::*)() >(&G4MTRunManagerKernel::SetUpDecayChannels));

  DEBUG_MSG("Adding wrapper for void G4MTRunManagerKernel::BroadcastAbortRun(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4MTRunManagerKernel::BroadcastAbortRun(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MTRunManagerKernel.hh:84:10
  t140.method("BroadcastAbortRun", static_cast<void (G4MTRunManagerKernel::*)(G4bool) >(&G4MTRunManagerKernel::BroadcastAbortRun));

  /* End of G4MTRunManagerKernel class method wrappers
   **********************************************************************/

}