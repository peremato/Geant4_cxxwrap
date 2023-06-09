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
void add_methods_for_G4UserSteppingAction(jlcxx::Module& types, jlcxx::TypeWrapper<G4UserSteppingAction>& t92) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4UserSteppingAction
   */

  DEBUG_MSG("Adding wrapper for void G4UserSteppingAction::SetSteppingManagerPointer(G4SteppingManager *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UserSteppingAction::SetSteppingManagerPointer(G4SteppingManager *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UserSteppingAction.hh:55:18
  t92.method("SetSteppingManagerPointer", static_cast<void (G4UserSteppingAction::*)(G4SteppingManager *) >(&G4UserSteppingAction::SetSteppingManagerPointer));

  DEBUG_MSG("Adding wrapper for void G4UserSteppingAction::UserSteppingAction(const G4Step *) (" __HERE__ ")");
  // signature to use in the veto list: void G4UserSteppingAction::UserSteppingAction(const G4Step *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4UserSteppingAction.hh:56:18
  t92.method("UserSteppingAction", static_cast<void (G4UserSteppingAction::*)(const G4Step *) >(&G4UserSteppingAction::UserSteppingAction));

  /* End of G4UserSteppingAction class method wrappers
   **********************************************************************/

}
