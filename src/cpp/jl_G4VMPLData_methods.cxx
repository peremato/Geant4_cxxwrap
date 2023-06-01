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
void add_methods_for_G4VMPLData(jlcxx::Module& types, jlcxx::TypeWrapper<G4VMPLData>& t165) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4VMPLData
   */

  DEBUG_MSG("Adding wrapper for void G4VMPLData::initialize() (" __HERE__ ")");
  // signature to use in the veto list: void G4VMPLData::initialize()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VModularPhysicsList.hh:61:10
  t165.method("initialize", static_cast<void (G4VMPLData::*)() >(&G4VMPLData::initialize));

  DEBUG_MSG("Adding physicsVector methods  to provide read access to the field physicsVector (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VModularPhysicsList.hh:64:28
  // signature to use in the veto list: G4VMPLData::physicsVector
  t165.method("physicsVector", [](const G4VMPLData& a) -> G4VMPLData::G4PhysConstVectorData * { return a.physicsVector; });
  t165.method("physicsVector", [](G4VMPLData& a) -> G4VMPLData::G4PhysConstVectorData * { return a.physicsVector; });
  t165.method("physicsVector", [](const G4VMPLData* a) -> G4VMPLData::G4PhysConstVectorData * { return a->physicsVector; });
  t165.method("physicsVector", [](G4VMPLData* a) -> G4VMPLData::G4PhysConstVectorData * { return a->physicsVector; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VModularPhysicsList.hh:64:28
  // signature to use in the veto list: G4VMPLData::physicsVector
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding physicsVector! methods to provide write access to the field physicsVector (" __HERE__ ")");
  t165.method("physicsVector!", [](G4VMPLData& a, G4VMPLData::G4PhysConstVectorData * val) -> G4VMPLData::G4PhysConstVectorData * { return a.physicsVector = val; });

  DEBUG_MSG("Adding physicsVector! methods to provide write access to the field physicsVector (" __HERE__ ")");
  t165.method("physicsVector!", [](G4VMPLData* a, G4VMPLData::G4PhysConstVectorData * val) -> G4VMPLData::G4PhysConstVectorData * { return a->physicsVector = val; });

  /* End of G4VMPLData class method wrappers
   **********************************************************************/

}
