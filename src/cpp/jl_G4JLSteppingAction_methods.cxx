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
void add_methods_for_G4JLSteppingAction(jlcxx::Module& types, jlcxx::TypeWrapper<G4JLSteppingAction>& t111) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4JLSteppingAction
   */


  DEBUG_MSG("Adding wrapper for void G4JLSteppingAction::G4JLSteppingAction(stepaction_f) (" __HERE__ ")");
  // defined in ./cpp/Geant4Wrap.h:135:3
  t111.constructor<stepaction_f>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for void G4JLSteppingAction::UserSteppingAction(const G4Step *) (" __HERE__ ")");
  // signature to use in the veto list: void G4JLSteppingAction::UserSteppingAction(const G4Step *)
  // defined in ./cpp/Geant4Wrap.h:137:16
  t111.method("UserSteppingAction", static_cast<void (G4JLSteppingAction::*)(const G4Step *) >(&G4JLSteppingAction::UserSteppingAction));

  /* End of G4JLSteppingAction class method wrappers
   **********************************************************************/

}