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
void add_methods_for_G4VUserPrimaryGeneratorAction(jlcxx::Module& types, jlcxx::TypeWrapper<G4VUserPrimaryGeneratorAction>& t93) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4VUserPrimaryGeneratorAction
   */

  DEBUG_MSG("Adding wrapper for void G4VUserPrimaryGeneratorAction::GeneratePrimaries(G4Event *) (" __HERE__ ")");
  // signature to use in the veto list: void G4VUserPrimaryGeneratorAction::GeneratePrimaries(G4Event *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4VUserPrimaryGeneratorAction.hh:54:18
  t93.method("GeneratePrimaries", static_cast<void (G4VUserPrimaryGeneratorAction::*)(G4Event *) >(&G4VUserPrimaryGeneratorAction::GeneratePrimaries));

  /* End of G4VUserPrimaryGeneratorAction class method wrappers
   **********************************************************************/

}
