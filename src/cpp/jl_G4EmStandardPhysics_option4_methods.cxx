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
void add_methods_for_G4EmStandardPhysics_option4(jlcxx::Module& types, jlcxx::TypeWrapper<G4EmStandardPhysics_option4>& t215) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4EmStandardPhysics_option4
   */


  DEBUG_MSG("Adding wrapper for void G4EmStandardPhysics_option4::G4EmStandardPhysics_option4(G4int, const G4String &) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EmStandardPhysics_option4.hh:56:12
  t215.constructor<G4int>(/*finalize=*/true);
  t215.constructor<G4int, const G4String &>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for void G4EmStandardPhysics_option4::ConstructParticle() (" __HERE__ ")");
  // signature to use in the veto list: void G4EmStandardPhysics_option4::ConstructParticle()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EmStandardPhysics_option4.hh:60:8
  t215.method("ConstructParticle", static_cast<void (G4EmStandardPhysics_option4::*)() >(&G4EmStandardPhysics_option4::ConstructParticle));

  DEBUG_MSG("Adding wrapper for void G4EmStandardPhysics_option4::ConstructProcess() (" __HERE__ ")");
  // signature to use in the veto list: void G4EmStandardPhysics_option4::ConstructProcess()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4EmStandardPhysics_option4.hh:61:8
  t215.method("ConstructProcess", static_cast<void (G4EmStandardPhysics_option4::*)() >(&G4EmStandardPhysics_option4::ConstructProcess));

  /* End of G4EmStandardPhysics_option4 class method wrappers
   **********************************************************************/

}
