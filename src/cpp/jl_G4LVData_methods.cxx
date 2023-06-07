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
void add_methods_for_G4LVData(jlcxx::Module& types, jlcxx::TypeWrapper<G4LVData>& t45) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4LVData
   */

  DEBUG_MSG("Adding wrapper for void G4LVData::initialize() (" __HERE__ ")");
  // signature to use in the veto list: void G4LVData::initialize()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:136:10
  t45.method("initialize", static_cast<void (G4LVData::*)() >(&G4LVData::initialize));

  DEBUG_MSG("Adding fSolid methods  to provide read access to the field fSolid (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:148:15
  // signature to use in the veto list: G4LVData::fSolid
  t45.method("fSolid", [](const G4LVData& a) -> G4VSolid * { return a.fSolid; });
  t45.method("fSolid", [](G4LVData& a) -> G4VSolid * { return a.fSolid; });
  t45.method("fSolid", [](const G4LVData* a) -> G4VSolid * { return a->fSolid; });
  t45.method("fSolid", [](G4LVData* a) -> G4VSolid * { return a->fSolid; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:148:15
  // signature to use in the veto list: G4LVData::fSolid
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fSolid! methods to provide write access to the field fSolid (" __HERE__ ")");
  t45.method("fSolid!", [](G4LVData& a, G4VSolid * val) -> G4VSolid * { return a.fSolid = val; });

  DEBUG_MSG("Adding fSolid! methods to provide write access to the field fSolid (" __HERE__ ")");
  t45.method("fSolid!", [](G4LVData* a, G4VSolid * val) -> G4VSolid * { return a->fSolid = val; });

  DEBUG_MSG("Adding fSensitiveDetector methods  to provide read access to the field fSensitiveDetector (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:150:27
  // signature to use in the veto list: G4LVData::fSensitiveDetector
  t45.method("fSensitiveDetector", [](const G4LVData& a) -> G4VSensitiveDetector * { return a.fSensitiveDetector; });
  t45.method("fSensitiveDetector", [](G4LVData& a) -> G4VSensitiveDetector * { return a.fSensitiveDetector; });
  t45.method("fSensitiveDetector", [](const G4LVData* a) -> G4VSensitiveDetector * { return a->fSensitiveDetector; });
  t45.method("fSensitiveDetector", [](G4LVData* a) -> G4VSensitiveDetector * { return a->fSensitiveDetector; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:150:27
  // signature to use in the veto list: G4LVData::fSensitiveDetector
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fSensitiveDetector! methods to provide write access to the field fSensitiveDetector (" __HERE__ ")");
  t45.method("fSensitiveDetector!", [](G4LVData& a, G4VSensitiveDetector * val) -> G4VSensitiveDetector * { return a.fSensitiveDetector = val; });

  DEBUG_MSG("Adding fSensitiveDetector! methods to provide write access to the field fSensitiveDetector (" __HERE__ ")");
  t45.method("fSensitiveDetector!", [](G4LVData* a, G4VSensitiveDetector * val) -> G4VSensitiveDetector * { return a->fSensitiveDetector = val; });

  DEBUG_MSG("Adding fFieldManager methods  to provide read access to the field fFieldManager (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:152:21
  // signature to use in the veto list: G4LVData::fFieldManager
  t45.method("fFieldManager", [](const G4LVData& a) -> G4FieldManager * { return a.fFieldManager; });
  t45.method("fFieldManager", [](G4LVData& a) -> G4FieldManager * { return a.fFieldManager; });
  t45.method("fFieldManager", [](const G4LVData* a) -> G4FieldManager * { return a->fFieldManager; });
  t45.method("fFieldManager", [](G4LVData* a) -> G4FieldManager * { return a->fFieldManager; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:152:21
  // signature to use in the veto list: G4LVData::fFieldManager
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fFieldManager! methods to provide write access to the field fFieldManager (" __HERE__ ")");
  t45.method("fFieldManager!", [](G4LVData& a, G4FieldManager * val) -> G4FieldManager * { return a.fFieldManager = val; });

  DEBUG_MSG("Adding fFieldManager! methods to provide write access to the field fFieldManager (" __HERE__ ")");
  t45.method("fFieldManager!", [](G4LVData* a, G4FieldManager * val) -> G4FieldManager * { return a->fFieldManager = val; });

  DEBUG_MSG("Adding fMaterial methods  to provide read access to the field fMaterial (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:154:17
  // signature to use in the veto list: G4LVData::fMaterial
  t45.method("fMaterial", [](const G4LVData& a) -> G4Material * { return a.fMaterial; });
  t45.method("fMaterial", [](G4LVData& a) -> G4Material * { return a.fMaterial; });
  t45.method("fMaterial", [](const G4LVData* a) -> G4Material * { return a->fMaterial; });
  t45.method("fMaterial", [](G4LVData* a) -> G4Material * { return a->fMaterial; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:154:17
  // signature to use in the veto list: G4LVData::fMaterial
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fMaterial! methods to provide write access to the field fMaterial (" __HERE__ ")");
  t45.method("fMaterial!", [](G4LVData& a, G4Material * val) -> G4Material * { return a.fMaterial = val; });

  DEBUG_MSG("Adding fMaterial! methods to provide write access to the field fMaterial (" __HERE__ ")");
  t45.method("fMaterial!", [](G4LVData* a, G4Material * val) -> G4Material * { return a->fMaterial = val; });

  DEBUG_MSG("Adding fMass methods  to provide read access to the field fMass (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:156:14
  // signature to use in the veto list: G4LVData::fMass
  t45.method("fMass", [](const G4LVData& a) -> G4double { return a.fMass; });
  t45.method("fMass", [](G4LVData& a) -> G4double { return a.fMass; });
  t45.method("fMass", [](const G4LVData* a) -> G4double { return a->fMass; });
  t45.method("fMass", [](G4LVData* a) -> G4double { return a->fMass; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:156:14
  // signature to use in the veto list: G4LVData::fMass
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fMass! methods to provide write access to the field fMass (" __HERE__ ")");
  t45.method("fMass!", [](G4LVData& a, G4double val) -> G4double { return a.fMass = val; });

  DEBUG_MSG("Adding fMass! methods to provide write access to the field fMass (" __HERE__ ")");
  t45.method("fMass!", [](G4LVData* a, G4double val) -> G4double { return a->fMass = val; });

  DEBUG_MSG("Adding fCutsCouple methods  to provide read access to the field fCutsCouple (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:158:27
  // signature to use in the veto list: G4LVData::fCutsCouple
  t45.method("fCutsCouple", [](const G4LVData& a) -> G4MaterialCutsCouple * { return a.fCutsCouple; });
  t45.method("fCutsCouple", [](G4LVData& a) -> G4MaterialCutsCouple * { return a.fCutsCouple; });
  t45.method("fCutsCouple", [](const G4LVData* a) -> G4MaterialCutsCouple * { return a->fCutsCouple; });
  t45.method("fCutsCouple", [](G4LVData* a) -> G4MaterialCutsCouple * { return a->fCutsCouple; });
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4LogicalVolume.hh:158:27
  // signature to use in the veto list: G4LVData::fCutsCouple
  // with ! suffix to veto the setter only

  DEBUG_MSG("Adding fCutsCouple! methods to provide write access to the field fCutsCouple (" __HERE__ ")");
  t45.method("fCutsCouple!", [](G4LVData& a, G4MaterialCutsCouple * val) -> G4MaterialCutsCouple * { return a.fCutsCouple = val; });

  DEBUG_MSG("Adding fCutsCouple! methods to provide write access to the field fCutsCouple (" __HERE__ ")");
  t45.method("fCutsCouple!", [](G4LVData* a, G4MaterialCutsCouple * val) -> G4MaterialCutsCouple * { return a->fCutsCouple = val; });

  /* End of G4LVData class method wrappers
   **********************************************************************/

}
