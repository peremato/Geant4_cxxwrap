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
void add_methods_for_G4LVData(jlcxx::Module& types, jlcxx::TypeWrapper<G4LVData>& t41) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4LVData
   */

  DEBUG_MSG("Adding wrapper for void G4LVData::initialize() (" __HERE__ ")");
  // signature to use in the veto list: void G4LVData::initialize()
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:136:10
  t41.method("initialize", static_cast<void (G4LVData::*)() >(&G4LVData::initialize));

  DEBUG_MSG("Adding fSolid methods  to provide read access to the field fSolid (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:148:15
  // signature to use in the veto list: G4LVData::fSolid
  t41.method("fSolid", [](const G4LVData& a) -> G4VSolid * { return a.fSolid; });
  t41.method("fSolid", [](G4LVData& a) -> G4VSolid * { return a.fSolid; });
  t41.method("fSolid", [](const G4LVData* a) -> G4VSolid * { return a->fSolid; });
  t41.method("fSolid", [](G4LVData* a) -> G4VSolid * { return a->fSolid; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:148:15
  // signature to use in the veto list: G4LVData::fSolid
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fSolid! methods to provide write access to the field fSolid (" __HERE__ ")");
  t41.method("fSolid!", [](G4LVData& a, G4VSolid * val) -> G4VSolid * { return a.fSolid = val; });

  DEBUG_MSG("Adding fSolid! methods to provide write access to the field fSolid (" __HERE__ ")");
  t41.method("fSolid!", [](G4LVData* a, G4VSolid * val) -> G4VSolid * { return a->fSolid = val; });

  DEBUG_MSG("Adding fSensitiveDetector methods  to provide read access to the field fSensitiveDetector (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:150:27
  // signature to use in the veto list: G4LVData::fSensitiveDetector
  t41.method("fSensitiveDetector", [](const G4LVData& a) -> G4VSensitiveDetector * { return a.fSensitiveDetector; });
  t41.method("fSensitiveDetector", [](G4LVData& a) -> G4VSensitiveDetector * { return a.fSensitiveDetector; });
  t41.method("fSensitiveDetector", [](const G4LVData* a) -> G4VSensitiveDetector * { return a->fSensitiveDetector; });
  t41.method("fSensitiveDetector", [](G4LVData* a) -> G4VSensitiveDetector * { return a->fSensitiveDetector; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:150:27
  // signature to use in the veto list: G4LVData::fSensitiveDetector
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fSensitiveDetector! methods to provide write access to the field fSensitiveDetector (" __HERE__ ")");
  t41.method("fSensitiveDetector!", [](G4LVData& a, G4VSensitiveDetector * val) -> G4VSensitiveDetector * { return a.fSensitiveDetector = val; });

  DEBUG_MSG("Adding fSensitiveDetector! methods to provide write access to the field fSensitiveDetector (" __HERE__ ")");
  t41.method("fSensitiveDetector!", [](G4LVData* a, G4VSensitiveDetector * val) -> G4VSensitiveDetector * { return a->fSensitiveDetector = val; });

  DEBUG_MSG("Adding fFieldManager methods  to provide read access to the field fFieldManager (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:152:21
  // signature to use in the veto list: G4LVData::fFieldManager
  t41.method("fFieldManager", [](const G4LVData& a) -> G4FieldManager * { return a.fFieldManager; });
  t41.method("fFieldManager", [](G4LVData& a) -> G4FieldManager * { return a.fFieldManager; });
  t41.method("fFieldManager", [](const G4LVData* a) -> G4FieldManager * { return a->fFieldManager; });
  t41.method("fFieldManager", [](G4LVData* a) -> G4FieldManager * { return a->fFieldManager; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:152:21
  // signature to use in the veto list: G4LVData::fFieldManager
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fFieldManager! methods to provide write access to the field fFieldManager (" __HERE__ ")");
  t41.method("fFieldManager!", [](G4LVData& a, G4FieldManager * val) -> G4FieldManager * { return a.fFieldManager = val; });

  DEBUG_MSG("Adding fFieldManager! methods to provide write access to the field fFieldManager (" __HERE__ ")");
  t41.method("fFieldManager!", [](G4LVData* a, G4FieldManager * val) -> G4FieldManager * { return a->fFieldManager = val; });

  DEBUG_MSG("Adding fMaterial methods  to provide read access to the field fMaterial (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:154:17
  // signature to use in the veto list: G4LVData::fMaterial
  t41.method("fMaterial", [](const G4LVData& a) -> G4Material * { return a.fMaterial; });
  t41.method("fMaterial", [](G4LVData& a) -> G4Material * { return a.fMaterial; });
  t41.method("fMaterial", [](const G4LVData* a) -> G4Material * { return a->fMaterial; });
  t41.method("fMaterial", [](G4LVData* a) -> G4Material * { return a->fMaterial; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:154:17
  // signature to use in the veto list: G4LVData::fMaterial
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fMaterial! methods to provide write access to the field fMaterial (" __HERE__ ")");
  t41.method("fMaterial!", [](G4LVData& a, G4Material * val) -> G4Material * { return a.fMaterial = val; });

  DEBUG_MSG("Adding fMaterial! methods to provide write access to the field fMaterial (" __HERE__ ")");
  t41.method("fMaterial!", [](G4LVData* a, G4Material * val) -> G4Material * { return a->fMaterial = val; });

  DEBUG_MSG("Adding fMass methods  to provide read access to the field fMass (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:156:14
  // signature to use in the veto list: G4LVData::fMass
  t41.method("fMass", [](const G4LVData& a) -> G4double { return a.fMass; });
  t41.method("fMass", [](G4LVData& a) -> G4double { return a.fMass; });
  t41.method("fMass", [](const G4LVData* a) -> G4double { return a->fMass; });
  t41.method("fMass", [](G4LVData* a) -> G4double { return a->fMass; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:156:14
  // signature to use in the veto list: G4LVData::fMass
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fMass! methods to provide write access to the field fMass (" __HERE__ ")");
  t41.method("fMass!", [](G4LVData& a, G4double val) -> G4double { return a.fMass = val; });

  DEBUG_MSG("Adding fMass! methods to provide write access to the field fMass (" __HERE__ ")");
  t41.method("fMass!", [](G4LVData* a, G4double val) -> G4double { return a->fMass = val; });

  DEBUG_MSG("Adding fCutsCouple methods  to provide read access to the field fCutsCouple (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:158:27
  // signature to use in the veto list: G4LVData::fCutsCouple
  t41.method("fCutsCouple", [](const G4LVData& a) -> G4MaterialCutsCouple * { return a.fCutsCouple; });
  t41.method("fCutsCouple", [](G4LVData& a) -> G4MaterialCutsCouple * { return a.fCutsCouple; });
  t41.method("fCutsCouple", [](const G4LVData* a) -> G4MaterialCutsCouple * { return a->fCutsCouple; });
  t41.method("fCutsCouple", [](G4LVData* a) -> G4MaterialCutsCouple * { return a->fCutsCouple; });
  // defined in /Users/mato/.julia/artifacts/3ddffb81697f6dd4742c75d8b8d14865fe8a388c/include/Geant4/G4LogicalVolume.hh:158:27
  // signature to use in the veto list: G4LVData::fCutsCouple
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fCutsCouple! methods to provide write access to the field fCutsCouple (" __HERE__ ")");
  t41.method("fCutsCouple!", [](G4LVData& a, G4MaterialCutsCouple * val) -> G4MaterialCutsCouple * { return a.fCutsCouple = val; });

  DEBUG_MSG("Adding fCutsCouple! methods to provide write access to the field fCutsCouple (" __HERE__ ")");
  t41.method("fCutsCouple!", [](G4LVData* a, G4MaterialCutsCouple * val) -> G4MaterialCutsCouple * { return a->fCutsCouple = val; });

  /* End of G4LVData class method wrappers
   **********************************************************************/

}
