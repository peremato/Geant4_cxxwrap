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
void add_methods_for_G4JLWorkerInitialization(jlcxx::Module& types, jlcxx::TypeWrapper<G4JLWorkerInitialization>& t109) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4JLWorkerInitialization
   */

  DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerInitialize() (" __HERE__ ")");
  // signature to use in the veto list: void G4JLWorkerInitialization::WorkerInitialize()
  // defined in ./cpp/Geant4Wrap.h:128:18
  t109.method("WorkerInitialize", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerInitialize));

  DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerStart() (" __HERE__ ")");
  // signature to use in the veto list: void G4JLWorkerInitialization::WorkerStart()
  // defined in ./cpp/Geant4Wrap.h:129:18
  t109.method("WorkerStart", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerStart));

  DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerRunStart() (" __HERE__ ")");
  // signature to use in the veto list: void G4JLWorkerInitialization::WorkerRunStart()
  // defined in ./cpp/Geant4Wrap.h:130:18
  t109.method("WorkerRunStart", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerRunStart));

  DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerRunEnd() (" __HERE__ ")");
  // signature to use in the veto list: void G4JLWorkerInitialization::WorkerRunEnd()
  // defined in ./cpp/Geant4Wrap.h:131:18
  t109.method("WorkerRunEnd", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerRunEnd));

  DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerStop() (" __HERE__ ")");
  // signature to use in the veto list: void G4JLWorkerInitialization::WorkerStop()
  // defined in ./cpp/Geant4Wrap.h:132:18
  t109.method("WorkerStop", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerStop));

  /* End of G4JLWorkerInitialization class method wrappers
   **********************************************************************/

}
