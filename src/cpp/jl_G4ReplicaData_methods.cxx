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
void add_methods_for_G4ReplicaData(jlcxx::Module& types, jlcxx::TypeWrapper<G4ReplicaData>& t152) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4ReplicaData
   */

  DEBUG_MSG("Adding wrapper for void G4ReplicaData::initialize() (" __HERE__ ")");
  // signature to use in the veto list: void G4ReplicaData::initialize()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:82:8
  t152.method("initialize", static_cast<void (G4ReplicaData::*)() >(&G4ReplicaData::initialize));

  DEBUG_MSG("Adding fcopyNo methods  to provide read access to the field fcopyNo (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:84:9
  // signature to use in the veto list: G4ReplicaData::fcopyNo
  t152.method("fcopyNo", [](const G4ReplicaData& a) -> G4int { return a.fcopyNo; });
  t152.method("fcopyNo", [](G4ReplicaData& a) -> G4int { return a.fcopyNo; });
  t152.method("fcopyNo", [](const G4ReplicaData* a) -> G4int { return a->fcopyNo; });
  t152.method("fcopyNo", [](G4ReplicaData* a) -> G4int { return a->fcopyNo; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4PVReplica.hh:84:9
  // signature to use in the veto list: G4ReplicaData::fcopyNo
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fcopyNo! methods to provide write access to the field fcopyNo (" __HERE__ ")");
  t152.method("fcopyNo!", [](G4ReplicaData& a, G4int val) -> G4int { return a.fcopyNo = val; });

  DEBUG_MSG("Adding fcopyNo! methods to provide write access to the field fcopyNo (" __HERE__ ")");
  t152.method("fcopyNo!", [](G4ReplicaData* a, G4int val) -> G4int { return a->fcopyNo = val; });

  /* End of G4ReplicaData class method wrappers
   **********************************************************************/

}
