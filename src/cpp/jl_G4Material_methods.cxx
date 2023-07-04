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
void add_methods_for_G4Material(jlcxx::Module& types, jlcxx::TypeWrapper<G4Material>& t41) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4Material
   */


  DEBUG_MSG("Adding wrapper for void G4Material::G4Material(const G4String &, G4double, G4double, G4double, G4State, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:122:3
  t41.constructor<const G4String &, G4double, G4double, G4double>(/*finalize=*/false);
  t41.constructor<const G4String &, G4double, G4double, G4double, G4State>(/*finalize=*/false);
  t41.constructor<const G4String &, G4double, G4double, G4double, G4State, G4double>(/*finalize=*/false);
  t41.constructor<const G4String &, G4double, G4double, G4double, G4State, G4double, G4double>(/*finalize=*/false);


  DEBUG_MSG("Adding wrapper for void G4Material::G4Material(const G4String &, G4double, G4int, G4State, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:134:3
  t41.constructor<const G4String &, G4double, G4int>(/*finalize=*/false);
  t41.constructor<const G4String &, G4double, G4int, G4State>(/*finalize=*/false);
  t41.constructor<const G4String &, G4double, G4int, G4State, G4double>(/*finalize=*/false);
  t41.constructor<const G4String &, G4double, G4int, G4State, G4double, G4double>(/*finalize=*/false);


  DEBUG_MSG("Adding wrapper for void G4Material::G4Material(const G4String &, G4double, const G4Material *, G4State, G4double, G4double) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:144:3
  t41.constructor<const G4String &, G4double, const G4Material *>(/*finalize=*/false);
  t41.constructor<const G4String &, G4double, const G4Material *, G4State>(/*finalize=*/false);
  t41.constructor<const G4String &, G4double, const G4Material *, G4State, G4double>(/*finalize=*/false);
  t41.constructor<const G4String &, G4double, const G4Material *, G4State, G4double, G4double>(/*finalize=*/false);

  DEBUG_MSG("Adding wrapper for void G4Material::AddElementByNumberOfAtoms(const G4Element *, G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::AddElementByNumberOfAtoms(const G4Element *, G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:154:8
  t41.method("AddElementByNumberOfAtoms", static_cast<void (G4Material::*)(const G4Element *, G4int) >(&G4Material::AddElementByNumberOfAtoms));

  DEBUG_MSG("Adding wrapper for void G4Material::AddElement(G4Element *, G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::AddElement(G4Element *, G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:156:8
  t41.method("AddElement", static_cast<void (G4Material::*)(G4Element *, G4int) >(&G4Material::AddElement));

  DEBUG_MSG("Adding wrapper for void G4Material::AddElementByMassFraction(const G4Element *, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::AddElementByMassFraction(const G4Element *, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:162:8
  t41.method("AddElementByMassFraction", static_cast<void (G4Material::*)(const G4Element *, G4double) >(&G4Material::AddElementByMassFraction));

  DEBUG_MSG("Adding wrapper for void G4Material::AddElement(G4Element *, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::AddElement(G4Element *, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:163:15
  t41.method("AddElement", static_cast<void (G4Material::*)(G4Element *, G4double) >(&G4Material::AddElement));

  DEBUG_MSG("Adding wrapper for void G4Material::AddMaterial(G4Material *, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::AddMaterial(G4Material *, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:166:8
  t41.method("AddMaterial", static_cast<void (G4Material::*)(G4Material *, G4double) >(&G4Material::AddMaterial));

  DEBUG_MSG("Adding wrapper for const G4String & G4Material::GetName() (" __HERE__ ")");
  // signature to use in the veto list: const G4String & G4Material::GetName()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:172:26
  t41.method("GetName", static_cast<const G4String & (G4Material::*)()  const>(&G4Material::GetName));

  DEBUG_MSG("Adding wrapper for const G4String & G4Material::GetChemicalFormula() (" __HERE__ ")");
  // signature to use in the veto list: const G4String & G4Material::GetChemicalFormula()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:173:26
  t41.method("GetChemicalFormula", static_cast<const G4String & (G4Material::*)()  const>(&G4Material::GetChemicalFormula));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetFreeElectronDensity() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetFreeElectronDensity()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:174:19
  t41.method("GetFreeElectronDensity", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetFreeElectronDensity));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetDensity() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetDensity()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:175:19
  t41.method("GetDensity", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetDensity));

  DEBUG_MSG("Adding wrapper for G4State G4Material::GetState() (" __HERE__ ")");
  // signature to use in the veto list: G4State G4Material::GetState()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:176:19
  t41.method("GetState", static_cast<G4State (G4Material::*)()  const>(&G4Material::GetState));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetTemperature() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetTemperature()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:177:19
  t41.method("GetTemperature", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetTemperature));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetPressure() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetPressure()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:178:19
  t41.method("GetPressure", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetPressure));

  DEBUG_MSG("Adding wrapper for size_t G4Material::GetNumberOfElements() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4Material::GetNumberOfElements()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:181:17
  t41.method("GetNumberOfElements", static_cast<size_t (G4Material::*)()  const>(&G4Material::GetNumberOfElements));

  DEBUG_MSG("Adding wrapper for const G4ElementVector * G4Material::GetElementVector() (" __HERE__ ")");
  // signature to use in the veto list: const G4ElementVector * G4Material::GetElementVector()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:185:20
  t41.method("GetElementVector", static_cast<const G4ElementVector * (G4Material::*)()  const>(&G4Material::GetElementVector));

  DEBUG_MSG("Adding wrapper for const G4double * G4Material::GetFractionVector() (" __HERE__ ")");
  // signature to use in the veto list: const G4double * G4Material::GetFractionVector()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:189:13
  t41.method("GetFractionVector", static_cast<const G4double * (G4Material::*)()  const>(&G4Material::GetFractionVector));

  DEBUG_MSG("Adding wrapper for const G4int * G4Material::GetAtomsVector() (" __HERE__ ")");
  // signature to use in the veto list: const G4int * G4Material::GetAtomsVector()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:193:13
  t41.method("GetAtomsVector", static_cast<const G4int * (G4Material::*)()  const>(&G4Material::GetAtomsVector));

  DEBUG_MSG("Adding wrapper for const G4Element * G4Material::GetElement(G4int) (" __HERE__ ")");
  // signature to use in the veto list: const G4Element * G4Material::GetElement(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:197:14
  t41.method("GetElement", static_cast<const G4Element * (G4Material::*)(G4int)  const>(&G4Material::GetElement));

  DEBUG_MSG("Adding wrapper for const G4double * G4Material::GetVecNbOfAtomsPerVolume() (" __HERE__ ")");
  // signature to use in the veto list: const G4double * G4Material::GetVecNbOfAtomsPerVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:201:13
  t41.method("GetVecNbOfAtomsPerVolume", static_cast<const G4double * (G4Material::*)()  const>(&G4Material::GetVecNbOfAtomsPerVolume));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetTotNbOfAtomsPerVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetTotNbOfAtomsPerVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:204:13
  t41.method("GetTotNbOfAtomsPerVolume", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetTotNbOfAtomsPerVolume));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetTotNbOfElectPerVolume() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetTotNbOfElectPerVolume()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:207:13
  t41.method("GetTotNbOfElectPerVolume", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetTotNbOfElectPerVolume));

  DEBUG_MSG("Adding wrapper for const G4double * G4Material::GetAtomicNumDensityVector() (" __HERE__ ")");
  // signature to use in the veto list: const G4double * G4Material::GetAtomicNumDensityVector()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:211:13
  t41.method("GetAtomicNumDensityVector", static_cast<const G4double * (G4Material::*)()  const>(&G4Material::GetAtomicNumDensityVector));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetElectronDensity() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetElectronDensity()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:212:20
  t41.method("GetElectronDensity", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetElectronDensity));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetRadlen() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetRadlen()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:215:20
  t41.method("GetRadlen", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetRadlen));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetNuclearInterLength() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetNuclearInterLength()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:218:19
  t41.method("GetNuclearInterLength", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetNuclearInterLength));

  DEBUG_MSG("Adding wrapper for G4IonisParamMat * G4Material::GetIonisation() (" __HERE__ ")");
  // signature to use in the veto list: G4IonisParamMat * G4Material::GetIonisation()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:221:27
  t41.method("GetIonisation", static_cast<G4IonisParamMat * (G4Material::*)()  const>(&G4Material::GetIonisation));

  DEBUG_MSG("Adding wrapper for G4SandiaTable * G4Material::GetSandiaTable() (" __HERE__ ")");
  // signature to use in the veto list: G4SandiaTable * G4Material::GetSandiaTable()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:224:25
  t41.method("GetSandiaTable", static_cast<G4SandiaTable * (G4Material::*)()  const>(&G4Material::GetSandiaTable));

  DEBUG_MSG("Adding wrapper for const G4Material * G4Material::GetBaseMaterial() (" __HERE__ ")");
  // signature to use in the veto list: const G4Material * G4Material::GetBaseMaterial()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:228:21
  t41.method("GetBaseMaterial", static_cast<const G4Material * (G4Material::*)()  const>(&G4Material::GetBaseMaterial));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetMassOfMolecule() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetMassOfMolecule()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:236:19
  t41.method("GetMassOfMolecule", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetMassOfMolecule));

  DEBUG_MSG("Adding wrapper for void G4Material::SetChemicalFormula(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::SetChemicalFormula(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:238:8
  t41.method("SetChemicalFormula", static_cast<void (G4Material::*)(const G4String &) >(&G4Material::SetChemicalFormula));

  DEBUG_MSG("Adding wrapper for void G4Material::SetFreeElectronDensity(G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::SetFreeElectronDensity(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:240:8
  t41.method("SetFreeElectronDensity", static_cast<void (G4Material::*)(G4double) >(&G4Material::SetFreeElectronDensity));

  DEBUG_MSG("Adding wrapper for void G4Material::ComputeDensityEffectOnFly(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::ComputeDensityEffectOnFly(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:242:8
  t41.method("ComputeDensityEffectOnFly", static_cast<void (G4Material::*)(G4bool) >(&G4Material::ComputeDensityEffectOnFly));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetZ() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetZ()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:245:12
  t41.method("GetZ", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetZ));

  DEBUG_MSG("Adding wrapper for G4double G4Material::GetA() (" __HERE__ ")");
  // signature to use in the veto list: G4double G4Material::GetA()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:246:12
  t41.method("GetA", static_cast<G4double (G4Material::*)()  const>(&G4Material::GetA));

  DEBUG_MSG("Adding wrapper for void G4Material::SetMaterialPropertiesTable(G4MaterialPropertiesTable *) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::SetMaterialPropertiesTable(G4MaterialPropertiesTable *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:249:8
  t41.method("SetMaterialPropertiesTable", static_cast<void (G4Material::*)(G4MaterialPropertiesTable *) >(&G4Material::SetMaterialPropertiesTable));

  DEBUG_MSG("Adding wrapper for G4MaterialPropertiesTable * G4Material::GetMaterialPropertiesTable() (" __HERE__ ")");
  // signature to use in the veto list: G4MaterialPropertiesTable * G4Material::GetMaterialPropertiesTable()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:251:37
  t41.method("GetMaterialPropertiesTable", static_cast<G4MaterialPropertiesTable * (G4Material::*)()  const>(&G4Material::GetMaterialPropertiesTable));

  DEBUG_MSG("Adding wrapper for size_t G4Material::GetIndex() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4Material::GetIndex()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:255:17
  t41.method("GetIndex", static_cast<size_t (G4Material::*)()  const>(&G4Material::GetIndex));

  DEBUG_MSG("Adding wrapper for G4MaterialTable * G4Material::GetMaterialTable() (" __HERE__ ")");
  // signature to use in the veto list: G4MaterialTable * G4Material::GetMaterialTable()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:259:27
  t41.method("G4Material!GetMaterialTable", static_cast<G4MaterialTable * (*)() >(&G4Material::GetMaterialTable));

  DEBUG_MSG("Adding wrapper for size_t G4Material::GetNumberOfMaterials() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4Material::GetNumberOfMaterials()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:261:17
  t41.method("G4Material!GetNumberOfMaterials", static_cast<size_t (*)() >(&G4Material::GetNumberOfMaterials));

  DEBUG_MSG("Adding wrapper for G4Material * G4Material::GetMaterial(const G4String &, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4Material::GetMaterial(const G4String &, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:264:22
  t41.method("G4Material!GetMaterial", static_cast<G4Material * (*)(const G4String &, G4bool) >(&G4Material::GetMaterial));
  t41.method("G4Material!GetMaterial", [](const G4String & arg0)->G4Material *{ return G4Material::GetMaterial(arg0); });

  DEBUG_MSG("Adding wrapper for G4Material * G4Material::GetMaterial(G4double, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4Material::GetMaterial(G4double, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:267:22
  t41.method("G4Material!GetMaterial", static_cast<G4Material * (*)(G4double, G4double, G4double) >(&G4Material::GetMaterial));

  DEBUG_MSG("Adding wrapper for G4Material * G4Material::GetMaterial(size_t, G4double) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4Material::GetMaterial(size_t, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:270:22
  t41.method("G4Material!GetMaterial", static_cast<G4Material * (*)(size_t, G4double) >(&G4Material::GetMaterial));

  DEBUG_MSG("Adding wrapper for void G4Material::SetName(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4Material::SetName(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:284:15
  t41.method("SetName", static_cast<void (G4Material::*)(const G4String &) >(&G4Material::SetName));

  DEBUG_MSG("Adding wrapper for G4bool G4Material::IsExtended() (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4Material::IsExtended()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4Material.hh:286:18
  t41.method("IsExtended", static_cast<G4bool (G4Material::*)()  const>(&G4Material::IsExtended));

  /* End of G4Material class method wrappers
   **********************************************************************/

}