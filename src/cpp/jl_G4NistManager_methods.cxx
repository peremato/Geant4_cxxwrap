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
void add_methods_for_G4NistManager(jlcxx::Module& types, jlcxx::TypeWrapper<G4NistManager>& t150) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4NistManager
   */

  DEBUG_MSG("Adding wrapper for G4NistManager * G4NistManager::Instance() (" __HERE__ ")");
  // signature to use in the veto list: G4NistManager * G4NistManager::Instance()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:87:25
  t150.method("G4NistManager!Instance", static_cast<G4NistManager * (*)() >(&G4NistManager::Instance));

  DEBUG_MSG("Adding wrapper for G4Element * G4NistManager::GetElement(size_t) (" __HERE__ ")");
  // signature to use in the veto list: G4Element * G4NistManager::GetElement(size_t)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:92:21
  t150.method("GetElement", static_cast<G4Element * (G4NistManager::*)(size_t)  const>(&G4NistManager::GetElement));

  DEBUG_MSG("Adding wrapper for G4Element * G4NistManager::FindElement(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4Element * G4NistManager::FindElement(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:96:21
  t150.method("FindElement", static_cast<G4Element * (G4NistManager::*)(G4int)  const>(&G4NistManager::FindElement));

  DEBUG_MSG("Adding wrapper for G4Element * G4NistManager::FindOrBuildElement(G4int, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4Element * G4NistManager::FindOrBuildElement(G4int, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:97:21
  t150.method("FindOrBuildElement", static_cast<G4Element * (G4NistManager::*)(G4int, G4bool) >(&G4NistManager::FindOrBuildElement));
  t150.method("FindOrBuildElement", [](G4NistManager& a, G4int arg0)->G4Element *{ return a.FindOrBuildElement(arg0); });
  t150.method("FindOrBuildElement", [](G4NistManager* a, G4int arg0)->G4Element *{ return a->FindOrBuildElement(arg0); });

  DEBUG_MSG("Adding wrapper for G4Element * G4NistManager::FindOrBuildElement(const G4String &, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4Element * G4NistManager::FindOrBuildElement(const G4String &, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:101:21
  t150.method("FindOrBuildElement", static_cast<G4Element * (G4NistManager::*)(const G4String &, G4bool) >(&G4NistManager::FindOrBuildElement));
  t150.method("FindOrBuildElement", [](G4NistManager& a, const G4String & arg0)->G4Element *{ return a.FindOrBuildElement(arg0); });
  t150.method("FindOrBuildElement", [](G4NistManager* a, const G4String & arg0)->G4Element *{ return a->FindOrBuildElement(arg0); });

  DEBUG_MSG("Adding wrapper for size_t G4NistManager::GetNumberOfElements() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4NistManager::GetNumberOfElements()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:106:17
  t150.method("GetNumberOfElements", static_cast<size_t (G4NistManager::*)()  const>(&G4NistManager::GetNumberOfElements));

  DEBUG_MSG("Adding wrapper for G4int G4NistManager::GetZ(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4NistManager::GetZ(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:110:16
  t150.method("GetZ", static_cast<G4int (G4NistManager::*)(const G4String &)  const>(&G4NistManager::GetZ));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetAtomicMassAmu(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetAtomicMassAmu(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:115:19
  t150.method("GetAtomicMassAmu", static_cast<G4double (G4NistManager::*)(const G4String &)  const>(&G4NistManager::GetAtomicMassAmu));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetAtomicMassAmu(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetAtomicMassAmu(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:120:19
  t150.method("GetAtomicMassAmu", static_cast<G4double (G4NistManager::*)(G4int)  const>(&G4NistManager::GetAtomicMassAmu));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetIsotopeMass(G4int, G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetIsotopeMass(G4int, G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:124:19
  t150.method("GetIsotopeMass", static_cast<G4double (G4NistManager::*)(G4int, G4int)  const>(&G4NistManager::GetIsotopeMass));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetAtomicMass(G4int, G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetAtomicMass(G4int, G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:129:19
  t150.method("GetAtomicMass", static_cast<G4double (G4NistManager::*)(G4int, G4int)  const>(&G4NistManager::GetAtomicMass));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetTotalElectronBindingEnergy(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetTotalElectronBindingEnergy(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:133:19
  t150.method("GetTotalElectronBindingEnergy", static_cast<G4double (G4NistManager::*)(G4int)  const>(&G4NistManager::GetTotalElectronBindingEnergy));

  DEBUG_MSG("Adding wrapper for G4int G4NistManager::GetNistFirstIsotopeN(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4NistManager::GetNistFirstIsotopeN(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:137:16
  t150.method("GetNistFirstIsotopeN", static_cast<G4int (G4NistManager::*)(G4int)  const>(&G4NistManager::GetNistFirstIsotopeN));

  DEBUG_MSG("Adding wrapper for G4int G4NistManager::GetNumberOfNistIsotopes(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4NistManager::GetNumberOfNistIsotopes(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:141:16
  t150.method("GetNumberOfNistIsotopes", static_cast<G4int (G4NistManager::*)(G4int)  const>(&G4NistManager::GetNumberOfNistIsotopes));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetIsotopeAbundance(G4int, G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetIsotopeAbundance(G4int, G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:145:19
  t150.method("GetIsotopeAbundance", static_cast<G4double (G4NistManager::*)(G4int, G4int)  const>(&G4NistManager::GetIsotopeAbundance));

  DEBUG_MSG("Adding wrapper for void G4NistManager::PrintElement(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4NistManager::PrintElement(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:149:15
  t150.method("PrintElement", static_cast<void (G4NistManager::*)(G4int)  const>(&G4NistManager::PrintElement));

  DEBUG_MSG("Adding wrapper for void G4NistManager::PrintElement(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4NistManager::PrintElement(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:153:8
  t150.method("PrintElement", static_cast<void (G4NistManager::*)(const G4String &)  const>(&G4NistManager::PrintElement));

  DEBUG_MSG("Adding wrapper for void G4NistManager::PrintG4Element(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4NistManager::PrintG4Element(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:157:8
  t150.method("PrintG4Element", static_cast<void (G4NistManager::*)(const G4String &)  const>(&G4NistManager::PrintG4Element));

  DEBUG_MSG("Adding wrapper for const std::vector<G4String> & G4NistManager::GetNistElementNames() (" __HERE__ ")");
  // signature to use in the veto list: const std::vector<G4String> & G4NistManager::GetNistElementNames()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:161:39
  t150.method("GetNistElementNames", static_cast<const std::vector<G4String> & (G4NistManager::*)()  const>(&G4NistManager::GetNistElementNames));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetMeanIonisationEnergy(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetMeanIonisationEnergy(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:165:19
  t150.method("GetMeanIonisationEnergy", static_cast<G4double (G4NistManager::*)(G4int)  const>(&G4NistManager::GetMeanIonisationEnergy));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetNominalDensity(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetNominalDensity(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:170:19
  t150.method("GetNominalDensity", static_cast<G4double (G4NistManager::*)(G4int)  const>(&G4NistManager::GetNominalDensity));

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::GetMaterial(size_t) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::GetMaterial(size_t)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:174:22
  t150.method("GetMaterial", static_cast<G4Material * (G4NistManager::*)(size_t)  const>(&G4NistManager::GetMaterial));

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::FindMaterial(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::FindMaterial(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:178:22
  t150.method("FindMaterial", static_cast<G4Material * (G4NistManager::*)(const G4String &)  const>(&G4NistManager::FindMaterial));

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::FindOrBuildMaterial(const G4String &, G4bool, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::FindOrBuildMaterial(const G4String &, G4bool, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:179:22
  t150.method("FindOrBuildMaterial", static_cast<G4Material * (G4NistManager::*)(const G4String &, G4bool, G4bool) >(&G4NistManager::FindOrBuildMaterial));
  t150.method("FindOrBuildMaterial", [](G4NistManager& a, const G4String & arg0)->G4Material *{ return a.FindOrBuildMaterial(arg0); });
  t150.method("FindOrBuildMaterial", [](G4NistManager& a, const G4String & arg0, G4bool arg1)->G4Material *{ return a.FindOrBuildMaterial(arg0, arg1); });
  t150.method("FindOrBuildMaterial", [](G4NistManager* a, const G4String & arg0)->G4Material *{ return a->FindOrBuildMaterial(arg0); });
  t150.method("FindOrBuildMaterial", [](G4NistManager* a, const G4String & arg0, G4bool arg1)->G4Material *{ return a->FindOrBuildMaterial(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::FindSimpleMaterial(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::FindSimpleMaterial(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:185:22
  t150.method("FindSimpleMaterial", static_cast<G4Material * (G4NistManager::*)(G4int)  const>(&G4NistManager::FindSimpleMaterial));

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::FindOrBuildSimpleMaterial(G4int, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::FindOrBuildSimpleMaterial(G4int, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:186:22
  t150.method("FindOrBuildSimpleMaterial", static_cast<G4Material * (G4NistManager::*)(G4int, G4bool) >(&G4NistManager::FindOrBuildSimpleMaterial));
  t150.method("FindOrBuildSimpleMaterial", [](G4NistManager& a, G4int arg0)->G4Material *{ return a.FindOrBuildSimpleMaterial(arg0); });
  t150.method("FindOrBuildSimpleMaterial", [](G4NistManager* a, G4int arg0)->G4Material *{ return a->FindOrBuildSimpleMaterial(arg0); });

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::BuildMaterialWithNewDensity(const G4String &, const G4String &, G4double, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::BuildMaterialWithNewDensity(const G4String &, const G4String &, G4double, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:192:15
  t150.method("BuildMaterialWithNewDensity", static_cast<G4Material * (G4NistManager::*)(const G4String &, const G4String &, G4double, G4double, G4double) >(&G4NistManager::BuildMaterialWithNewDensity));
  t150.method("BuildMaterialWithNewDensity", [](G4NistManager& a, const G4String & arg0, const G4String & arg1)->G4Material *{ return a.BuildMaterialWithNewDensity(arg0, arg1); });
  t150.method("BuildMaterialWithNewDensity", [](G4NistManager& a, const G4String & arg0, const G4String & arg1, G4double arg2)->G4Material *{ return a.BuildMaterialWithNewDensity(arg0, arg1, arg2); });
  t150.method("BuildMaterialWithNewDensity", [](G4NistManager& a, const G4String & arg0, const G4String & arg1, G4double arg2, G4double arg3)->G4Material *{ return a.BuildMaterialWithNewDensity(arg0, arg1, arg2, arg3); });
  t150.method("BuildMaterialWithNewDensity", [](G4NistManager* a, const G4String & arg0, const G4String & arg1)->G4Material *{ return a->BuildMaterialWithNewDensity(arg0, arg1); });
  t150.method("BuildMaterialWithNewDensity", [](G4NistManager* a, const G4String & arg0, const G4String & arg1, G4double arg2)->G4Material *{ return a->BuildMaterialWithNewDensity(arg0, arg1, arg2); });
  t150.method("BuildMaterialWithNewDensity", [](G4NistManager* a, const G4String & arg0, const G4String & arg1, G4double arg2, G4double arg3)->G4Material *{ return a->BuildMaterialWithNewDensity(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::ConstructNewMaterial(const G4String &, const std::vector<G4String> &, const std::vector<G4int> &, G4double, G4bool, G4State, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::ConstructNewMaterial(const G4String &, const std::vector<G4String> &, const std::vector<G4int> &, G4double, G4bool, G4State, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:201:22
  t150.method("ConstructNewMaterial", static_cast<G4Material * (G4NistManager::*)(const G4String &, const std::vector<G4String> &, const std::vector<G4int> &, G4double, G4bool, G4State, G4double, G4double) >(&G4NistManager::ConstructNewMaterial));
  t150.method("ConstructNewMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4double arg3)->G4Material *{ return a.ConstructNewMaterial(arg0, arg1, arg2, arg3); });
  t150.method("ConstructNewMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4double arg3, G4bool arg4)->G4Material *{ return a.ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4); });
  t150.method("ConstructNewMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4double arg3, G4bool arg4, G4State arg5)->G4Material *{ return a.ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4, arg5); });
  t150.method("ConstructNewMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4double arg3, G4bool arg4, G4State arg5, G4double arg6)->G4Material *{ return a.ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4, arg5, arg6); });
  t150.method("ConstructNewMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4double arg3)->G4Material *{ return a->ConstructNewMaterial(arg0, arg1, arg2, arg3); });
  t150.method("ConstructNewMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4double arg3, G4bool arg4)->G4Material *{ return a->ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4); });
  t150.method("ConstructNewMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4double arg3, G4bool arg4, G4State arg5)->G4Material *{ return a->ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4, arg5); });
  t150.method("ConstructNewMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4double arg3, G4bool arg4, G4State arg5, G4double arg6)->G4Material *{ return a->ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4, arg5, arg6); });

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::ConstructNewMaterial(const G4String &, const std::vector<G4String> &, const std::vector<G4double> &, G4double, G4bool, G4State, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::ConstructNewMaterial(const G4String &, const std::vector<G4String> &, const std::vector<G4double> &, G4double, G4bool, G4State, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:214:22
  t150.method("ConstructNewMaterial", static_cast<G4Material * (G4NistManager::*)(const G4String &, const std::vector<G4String> &, const std::vector<G4double> &, G4double, G4bool, G4State, G4double, G4double) >(&G4NistManager::ConstructNewMaterial));
  t150.method("ConstructNewMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4double> & arg2, G4double arg3)->G4Material *{ return a.ConstructNewMaterial(arg0, arg1, arg2, arg3); });
  t150.method("ConstructNewMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4double> & arg2, G4double arg3, G4bool arg4)->G4Material *{ return a.ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4); });
  t150.method("ConstructNewMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4double> & arg2, G4double arg3, G4bool arg4, G4State arg5)->G4Material *{ return a.ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4, arg5); });
  t150.method("ConstructNewMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4double> & arg2, G4double arg3, G4bool arg4, G4State arg5, G4double arg6)->G4Material *{ return a.ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4, arg5, arg6); });
  t150.method("ConstructNewMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4double> & arg2, G4double arg3)->G4Material *{ return a->ConstructNewMaterial(arg0, arg1, arg2, arg3); });
  t150.method("ConstructNewMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4double> & arg2, G4double arg3, G4bool arg4)->G4Material *{ return a->ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4); });
  t150.method("ConstructNewMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4double> & arg2, G4double arg3, G4bool arg4, G4State arg5)->G4Material *{ return a->ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4, arg5); });
  t150.method("ConstructNewMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4double> & arg2, G4double arg3, G4bool arg4, G4State arg5, G4double arg6)->G4Material *{ return a->ConstructNewMaterial(arg0, arg1, arg2, arg3, arg4, arg5, arg6); });

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::ConstructNewGasMaterial(const G4String &, const G4String &, G4double, G4double, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::ConstructNewGasMaterial(const G4String &, const G4String &, G4double, G4double, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:226:22
  t150.method("ConstructNewGasMaterial", static_cast<G4Material * (G4NistManager::*)(const G4String &, const G4String &, G4double, G4double, G4bool) >(&G4NistManager::ConstructNewGasMaterial));
  t150.method("ConstructNewGasMaterial", [](G4NistManager& a, const G4String & arg0, const G4String & arg1, G4double arg2, G4double arg3)->G4Material *{ return a.ConstructNewGasMaterial(arg0, arg1, arg2, arg3); });
  t150.method("ConstructNewGasMaterial", [](G4NistManager* a, const G4String & arg0, const G4String & arg1, G4double arg2, G4double arg3)->G4Material *{ return a->ConstructNewGasMaterial(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for G4Material * G4NistManager::ConstructNewIdealGasMaterial(const G4String &, const std::vector<G4String> &, const std::vector<G4int> &, G4bool, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: G4Material * G4NistManager::ConstructNewIdealGasMaterial(const G4String &, const std::vector<G4String> &, const std::vector<G4int> &, G4bool, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:234:22
  t150.method("ConstructNewIdealGasMaterial", static_cast<G4Material * (G4NistManager::*)(const G4String &, const std::vector<G4String> &, const std::vector<G4int> &, G4bool, G4double, G4double) >(&G4NistManager::ConstructNewIdealGasMaterial));
  t150.method("ConstructNewIdealGasMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2)->G4Material *{ return a.ConstructNewIdealGasMaterial(arg0, arg1, arg2); });
  t150.method("ConstructNewIdealGasMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4bool arg3)->G4Material *{ return a.ConstructNewIdealGasMaterial(arg0, arg1, arg2, arg3); });
  t150.method("ConstructNewIdealGasMaterial", [](G4NistManager& a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4bool arg3, G4double arg4)->G4Material *{ return a.ConstructNewIdealGasMaterial(arg0, arg1, arg2, arg3, arg4); });
  t150.method("ConstructNewIdealGasMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2)->G4Material *{ return a->ConstructNewIdealGasMaterial(arg0, arg1, arg2); });
  t150.method("ConstructNewIdealGasMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4bool arg3)->G4Material *{ return a->ConstructNewIdealGasMaterial(arg0, arg1, arg2, arg3); });
  t150.method("ConstructNewIdealGasMaterial", [](G4NistManager* a, const G4String & arg0, const std::vector<G4String> & arg1, const std::vector<G4int> & arg2, G4bool arg3, G4double arg4)->G4Material *{ return a->ConstructNewIdealGasMaterial(arg0, arg1, arg2, arg3, arg4); });

  DEBUG_MSG("Adding wrapper for void G4NistManager::SetDensityEffectCalculatorFlag(const G4String &, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4NistManager::SetDensityEffectCalculatorFlag(const G4String &, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:244:8
  t150.method("SetDensityEffectCalculatorFlag", static_cast<void (G4NistManager::*)(const G4String &, G4bool) >(&G4NistManager::SetDensityEffectCalculatorFlag));

  DEBUG_MSG("Adding wrapper for void G4NistManager::SetDensityEffectCalculatorFlag(G4Material *, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4NistManager::SetDensityEffectCalculatorFlag(G4Material *, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:248:8
  t150.method("SetDensityEffectCalculatorFlag", static_cast<void (G4NistManager::*)(G4Material *, G4bool) >(&G4NistManager::SetDensityEffectCalculatorFlag));

  DEBUG_MSG("Adding wrapper for size_t G4NistManager::GetNumberOfMaterials() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4NistManager::GetNumberOfMaterials()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:252:17
  t150.method("GetNumberOfMaterials", static_cast<size_t (G4NistManager::*)()  const>(&G4NistManager::GetNumberOfMaterials));

  DEBUG_MSG("Adding wrapper for G4int G4NistManager::GetVerbose() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4NistManager::GetVerbose()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:254:16
  t150.method("GetVerbose", static_cast<G4int (G4NistManager::*)()  const>(&G4NistManager::GetVerbose));

  DEBUG_MSG("Adding wrapper for void G4NistManager::SetVerbose(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4NistManager::SetVerbose(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:256:8
  t150.method("SetVerbose", static_cast<void (G4NistManager::*)(G4int) >(&G4NistManager::SetVerbose));

  DEBUG_MSG("Adding wrapper for void G4NistManager::PrintG4Material(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4NistManager::PrintG4Material(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:260:8
  t150.method("PrintG4Material", static_cast<void (G4NistManager::*)(const G4String &)  const>(&G4NistManager::PrintG4Material));

  DEBUG_MSG("Adding wrapper for void G4NistManager::ListMaterials(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4NistManager::ListMaterials(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:269:15
  t150.method("ListMaterials", static_cast<void (G4NistManager::*)(const G4String &)  const>(&G4NistManager::ListMaterials));

  DEBUG_MSG("Adding wrapper for const std::vector<G4String> & G4NistManager::GetNistMaterialNames() (" __HERE__ ")");
  // signature to use in the veto list: const std::vector<G4String> & G4NistManager::GetNistMaterialNames()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:273:39
  t150.method("GetNistMaterialNames", static_cast<const std::vector<G4String> & (G4NistManager::*)()  const>(&G4NistManager::GetNistMaterialNames));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetZ13(G4double) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetZ13(G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:277:19
  t150.method("GetZ13", static_cast<G4double (G4NistManager::*)(G4double)  const>(&G4NistManager::GetZ13));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetZ13(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetZ13(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:278:19
  t150.method("GetZ13", static_cast<G4double (G4NistManager::*)(G4int)  const>(&G4NistManager::GetZ13));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetA27(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetA27(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:282:19
  t150.method("GetA27", static_cast<G4double (G4NistManager::*)(G4int)  const>(&G4NistManager::GetA27));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetLOGZ(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetLOGZ(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:286:19
  t150.method("GetLOGZ", static_cast<G4double (G4NistManager::*)(G4int)  const>(&G4NistManager::GetLOGZ));

  DEBUG_MSG("Adding wrapper for G4double G4NistManager::GetLOGAMU(G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4NistManager::GetLOGAMU(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:287:19
  t150.method("GetLOGAMU", static_cast<G4double (G4NistManager::*)(G4int)  const>(&G4NistManager::GetLOGAMU));

  DEBUG_MSG("Adding wrapper for G4ICRU90StoppingData * G4NistManager::GetICRU90StoppingData() (" __HERE__ ")");
  // signature to use in the veto list: G4ICRU90StoppingData * G4NistManager::GetICRU90StoppingData()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4NistManager.hh:289:25
  t150.method("GetICRU90StoppingData", static_cast<G4ICRU90StoppingData * (G4NistManager::*)() >(&G4NistManager::GetICRU90StoppingData));

  /* End of G4NistManager class method wrappers
   **********************************************************************/

}
