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
void add_methods_for_G4JLRunAction(jlcxx::Module& types, jlcxx::TypeWrapper<G4JLRunAction>& t116) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4JLRunAction
   */


  DEBUG_MSG("Adding wrapper for void G4JLRunAction::G4JLRunAction(runaction_f, void *, runaction_f, void *) (" __HERE__ ")");
  // defined in ./cpp/Geant4Wrap.h:193:5
  t116.constructor<runaction_f>(/*finalize=*/false);
  t116.constructor<runaction_f, void *>(/*finalize=*/false);
  t116.constructor<runaction_f, void *, runaction_f>(/*finalize=*/false);
  t116.constructor<runaction_f, void *, runaction_f, void *>(/*finalize=*/false);

  DEBUG_MSG("Adding wrapper for void G4JLRunAction::BeginOfRunAction(const G4Run *) (" __HERE__ ")");
  // signature to use in the veto list: void G4JLRunAction::BeginOfRunAction(const G4Run *)
  // defined in ./cpp/Geant4Wrap.h:197:18
  t116.method("BeginOfRunAction", static_cast<void (G4JLRunAction::*)(const G4Run *) >(&G4JLRunAction::BeginOfRunAction));

  DEBUG_MSG("Adding wrapper for void G4JLRunAction::EndOfRunAction(const G4Run *) (" __HERE__ ")");
  // signature to use in the veto list: void G4JLRunAction::EndOfRunAction(const G4Run *)
  // defined in ./cpp/Geant4Wrap.h:198:20
  t116.method("EndOfRunAction", static_cast<void (G4JLRunAction::*)(const G4Run *) >(&G4JLRunAction::EndOfRunAction));

  /* End of G4JLRunAction class method wrappers
   **********************************************************************/

}
