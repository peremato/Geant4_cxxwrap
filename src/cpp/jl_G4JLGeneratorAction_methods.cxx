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
void add_methods_for_G4JLGeneratorAction(jlcxx::Module& types, jlcxx::TypeWrapper<G4JLGeneratorAction>& t110) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4JLGeneratorAction
   */


  DEBUG_MSG("Adding wrapper for void G4JLGeneratorAction::G4JLGeneratorAction(generate_f, void *) (" __HERE__ ")");
  // defined in ./cpp/Geant4Wrap.h:138:3
  t110.constructor<generate_f, void *>(/*finalize=*/false);

  DEBUG_MSG("Adding wrapper for void G4JLGeneratorAction::GeneratePrimaries(G4Event *) (" __HERE__ ")");
  // signature to use in the veto list: void G4JLGeneratorAction::GeneratePrimaries(G4Event *)
  // defined in ./cpp/Geant4Wrap.h:140:8
  t110.method("GeneratePrimaries", static_cast<void (G4JLGeneratorAction::*)(G4Event *) >(&G4JLGeneratorAction::GeneratePrimaries));

  /* End of G4JLGeneratorAction class method wrappers
   **********************************************************************/

}
