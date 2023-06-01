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
void add_methods_for_G4PVData(jlcxx::Module& types, jlcxx::TypeWrapper<G4PVData>& t23) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4PVData
   */

  DEBUG_MSG("Adding wrapper for void G4PVData::initialize() (" __HERE__ ")");
  // signature to use in the veto list: void G4PVData::initialize()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VPhysicalVolume.hh:65:10
  t23.method("initialize", static_cast<void (G4PVData::*)() >(&G4PVData::initialize));

  DEBUG_MSG("Adding frot methods  to provide read access to the field frot (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VPhysicalVolume.hh:71:23
  // signature to use in the veto list: G4PVData::frot
  t23.method("frot", [](const G4PVData& a) -> G4RotationMatrix * { return a.frot; });
  t23.method("frot", [](G4PVData& a) -> G4RotationMatrix * { return a.frot; });
  t23.method("frot", [](const G4PVData* a) -> G4RotationMatrix * { return a->frot; });
  t23.method("frot", [](G4PVData* a) -> G4RotationMatrix * { return a->frot; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VPhysicalVolume.hh:71:23
  // signature to use in the veto list: G4PVData::frot
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding frot! methods to provide write access to the field frot (" __HERE__ ")");
  t23.method("frot!", [](G4PVData& a, G4RotationMatrix * val) -> G4RotationMatrix * { return a.frot = val; });

  DEBUG_MSG("Adding frot! methods to provide write access to the field frot (" __HERE__ ")");
  t23.method("frot!", [](G4PVData* a, G4RotationMatrix * val) -> G4RotationMatrix * { return a->frot = val; });

  DEBUG_MSG("Adding tx methods  to provide read access to the field tx (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VPhysicalVolume.hh:72:14
  // signature to use in the veto list: G4PVData::tx
  t23.method("tx", [](const G4PVData& a) -> G4double { return a.tx; });
  t23.method("tx", [](G4PVData& a) -> G4double { return a.tx; });
  t23.method("tx", [](const G4PVData* a) -> G4double { return a->tx; });
  t23.method("tx", [](G4PVData* a) -> G4double { return a->tx; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VPhysicalVolume.hh:72:14
  // signature to use in the veto list: G4PVData::tx
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding tx! methods to provide write access to the field tx (" __HERE__ ")");
  t23.method("tx!", [](G4PVData& a, G4double val) -> G4double { return a.tx = val; });

  DEBUG_MSG("Adding tx! methods to provide write access to the field tx (" __HERE__ ")");
  t23.method("tx!", [](G4PVData* a, G4double val) -> G4double { return a->tx = val; });

  DEBUG_MSG("Adding ty methods  to provide read access to the field ty (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VPhysicalVolume.hh:72:23
  // signature to use in the veto list: G4PVData::ty
  t23.method("ty", [](const G4PVData& a) -> G4double { return a.ty; });
  t23.method("ty", [](G4PVData& a) -> G4double { return a.ty; });
  t23.method("ty", [](const G4PVData* a) -> G4double { return a->ty; });
  t23.method("ty", [](G4PVData* a) -> G4double { return a->ty; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VPhysicalVolume.hh:72:23
  // signature to use in the veto list: G4PVData::ty
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding ty! methods to provide write access to the field ty (" __HERE__ ")");
  t23.method("ty!", [](G4PVData& a, G4double val) -> G4double { return a.ty = val; });

  DEBUG_MSG("Adding ty! methods to provide write access to the field ty (" __HERE__ ")");
  t23.method("ty!", [](G4PVData* a, G4double val) -> G4double { return a->ty = val; });

  DEBUG_MSG("Adding tz methods  to provide read access to the field tz (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VPhysicalVolume.hh:72:32
  // signature to use in the veto list: G4PVData::tz
  t23.method("tz", [](const G4PVData& a) -> G4double { return a.tz; });
  t23.method("tz", [](G4PVData& a) -> G4double { return a.tz; });
  t23.method("tz", [](const G4PVData* a) -> G4double { return a->tz; });
  t23.method("tz", [](G4PVData* a) -> G4double { return a->tz; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4VPhysicalVolume.hh:72:32
  // signature to use in the veto list: G4PVData::tz
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding tz! methods to provide write access to the field tz (" __HERE__ ")");
  t23.method("tz!", [](G4PVData& a, G4double val) -> G4double { return a.tz = val; });

  DEBUG_MSG("Adding tz! methods to provide write access to the field tz (" __HERE__ ")");
  t23.method("tz!", [](G4PVData* a, G4double val) -> G4double { return a->tz = val; });

  /* End of G4PVData class method wrappers
   **********************************************************************/

}
