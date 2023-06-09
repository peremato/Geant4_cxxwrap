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
void add_methods_for_G4GDMLParser(jlcxx::Module& types, jlcxx::TypeWrapper<G4GDMLParser>& t178) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4GDMLParser
   */


  DEBUG_MSG("Adding wrapper for void G4GDMLParser::G4GDMLParser(G4GDMLReadStructure *) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:56:5
  t178.constructor<G4GDMLReadStructure *>(/*finalize=*/true);


  DEBUG_MSG("Adding wrapper for void G4GDMLParser::G4GDMLParser(G4GDMLReadStructure *, G4GDMLWriteStructure *) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:57:5
  t178.constructor<G4GDMLReadStructure *, G4GDMLWriteStructure *>(/*finalize=*/true);

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::Read(const G4String &, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::Read(const G4String &, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:62:17
  t178.method("Read", static_cast<void (G4GDMLParser::*)(const G4String &, G4bool) >(&G4GDMLParser::Read));
  t178.method("Read", [](G4GDMLParser& a, const G4String & arg0)->void{ a.Read(arg0); });
  t178.method("Read", [](G4GDMLParser* a, const G4String & arg0)->void{ a->Read(arg0); });

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::ReadModule(const G4String &, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::ReadModule(const G4String &, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:67:17
  t178.method("ReadModule", static_cast<void (G4GDMLParser::*)(const G4String &, G4bool) >(&G4GDMLParser::ReadModule));
  t178.method("ReadModule", [](G4GDMLParser& a, const G4String & arg0)->void{ a.ReadModule(arg0); });
  t178.method("ReadModule", [](G4GDMLParser* a, const G4String & arg0)->void{ a->ReadModule(arg0); });

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::Write(const G4String &, const G4VPhysicalVolume *, G4bool, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::Write(const G4String &, const G4VPhysicalVolume *, G4bool, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:72:17
  t178.method("Write", static_cast<void (G4GDMLParser::*)(const G4String &, const G4VPhysicalVolume *, G4bool, const G4String &) >(&G4GDMLParser::Write));
  t178.method("Write", [](G4GDMLParser& a, const G4String & arg0)->void{ a.Write(arg0); });
  t178.method("Write", [](G4GDMLParser& a, const G4String & arg0, const G4VPhysicalVolume * arg1)->void{ a.Write(arg0, arg1); });
  t178.method("Write", [](G4GDMLParser& a, const G4String & arg0, const G4VPhysicalVolume * arg1, G4bool arg2)->void{ a.Write(arg0, arg1, arg2); });
  t178.method("Write", [](G4GDMLParser* a, const G4String & arg0)->void{ a->Write(arg0); });
  t178.method("Write", [](G4GDMLParser* a, const G4String & arg0, const G4VPhysicalVolume * arg1)->void{ a->Write(arg0, arg1); });
  t178.method("Write", [](G4GDMLParser* a, const G4String & arg0, const G4VPhysicalVolume * arg1, G4bool arg2)->void{ a->Write(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::Write(const G4String &, const G4LogicalVolume *, G4bool, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::Write(const G4String &, const G4LogicalVolume *, G4bool, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:83:17
  t178.method("Write", static_cast<void (G4GDMLParser::*)(const G4String &, const G4LogicalVolume *, G4bool, const G4String &) >(&G4GDMLParser::Write));
  t178.method("Write", [](G4GDMLParser& a, const G4String & arg0, const G4LogicalVolume * arg1)->void{ a.Write(arg0, arg1); });
  t178.method("Write", [](G4GDMLParser& a, const G4String & arg0, const G4LogicalVolume * arg1, G4bool arg2)->void{ a.Write(arg0, arg1, arg2); });
  t178.method("Write", [](G4GDMLParser* a, const G4String & arg0, const G4LogicalVolume * arg1)->void{ a->Write(arg0, arg1); });
  t178.method("Write", [](G4GDMLParser* a, const G4String & arg0, const G4LogicalVolume * arg1, G4bool arg2)->void{ a->Write(arg0, arg1, arg2); });

  DEBUG_MSG("Adding wrapper for G4LogicalVolume * G4GDMLParser::ParseST(const G4String &, G4Material *, G4Material *) (" __HERE__ ")");
  // signature to use in the veto list: G4LogicalVolume * G4GDMLParser::ParseST(const G4String &, G4Material *, G4Material *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:94:29
  t178.method("ParseST", static_cast<G4LogicalVolume * (G4GDMLParser::*)(const G4String &, G4Material *, G4Material *) >(&G4GDMLParser::ParseST));

  DEBUG_MSG("Adding wrapper for G4bool G4GDMLParser::IsValid(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4GDMLParser::IsValid(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:104:19
  t178.method("IsValid", static_cast<G4bool (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::IsValid));

  DEBUG_MSG("Adding wrapper for G4double G4GDMLParser::GetConstant(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4GDMLParser::GetConstant(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:105:21
  t178.method("GetConstant", static_cast<G4double (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetConstant));

  DEBUG_MSG("Adding wrapper for G4double G4GDMLParser::GetVariable(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4GDMLParser::GetVariable(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:106:21
  t178.method("GetVariable", static_cast<G4double (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetVariable));

  DEBUG_MSG("Adding wrapper for G4double G4GDMLParser::GetQuantity(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4GDMLParser::GetQuantity(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:107:21
  t178.method("GetQuantity", static_cast<G4double (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetQuantity));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4GDMLParser::GetPosition(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4GDMLParser::GetPosition(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:108:26
  t178.method("GetPosition", static_cast<G4ThreeVector (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetPosition));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4GDMLParser::GetRotation(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4GDMLParser::GetRotation(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:109:26
  t178.method("GetRotation", static_cast<G4ThreeVector (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetRotation));

  DEBUG_MSG("Adding wrapper for G4ThreeVector G4GDMLParser::GetScale(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4ThreeVector G4GDMLParser::GetScale(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:110:26
  t178.method("GetScale", static_cast<G4ThreeVector (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetScale));

  DEBUG_MSG("Adding wrapper for G4GDMLMatrix G4GDMLParser::GetMatrix(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4GDMLMatrix G4GDMLParser::GetMatrix(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:111:25
  t178.method("GetMatrix", static_cast<G4GDMLMatrix (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetMatrix));

  DEBUG_MSG("Adding wrapper for G4LogicalVolume * G4GDMLParser::GetVolume(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4LogicalVolume * G4GDMLParser::GetVolume(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:112:29
  t178.method("GetVolume", static_cast<G4LogicalVolume * (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetVolume));

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4GDMLParser::GetPhysVolume(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4GDMLParser::GetPhysVolume(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:113:31
  t178.method("GetPhysVolume", static_cast<G4VPhysicalVolume * (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetPhysVolume));

  DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4GDMLParser::GetWorldVolume(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4VPhysicalVolume * G4GDMLParser::GetWorldVolume(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:115:12
  t178.method("GetWorldVolume", static_cast<G4VPhysicalVolume * (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetWorldVolume));
  t178.method("GetWorldVolume", [](G4GDMLParser const& a)->G4VPhysicalVolume *{ return a.GetWorldVolume(); });
  t178.method("GetWorldVolume", [](G4GDMLParser const* a)->G4VPhysicalVolume *{ return a->GetWorldVolume(); });

  DEBUG_MSG("Adding wrapper for G4GDMLAuxListType G4GDMLParser::GetVolumeAuxiliaryInformation(G4LogicalVolume *) (" __HERE__ ")");
  // signature to use in the veto list: G4GDMLAuxListType G4GDMLParser::GetVolumeAuxiliaryInformation(G4LogicalVolume *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:117:12
  t178.method("GetVolumeAuxiliaryInformation", static_cast<G4GDMLAuxListType (G4GDMLParser::*)(G4LogicalVolume *)  const>(&G4GDMLParser::GetVolumeAuxiliaryInformation));

  DEBUG_MSG("Adding wrapper for const G4GDMLAuxListType * G4GDMLParser::GetAuxList() (" __HERE__ ")");
  // signature to use in the veto list: const G4GDMLAuxListType * G4GDMLParser::GetAuxList()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:119:37
  t178.method("GetAuxList", static_cast<const G4GDMLAuxListType * (G4GDMLParser::*)()  const>(&G4GDMLParser::GetAuxList));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::AddAuxiliary(G4GDMLAuxStructType) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::AddAuxiliary(G4GDMLAuxStructType)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:120:17
  t178.method("AddAuxiliary", static_cast<void (G4GDMLParser::*)(G4GDMLAuxStructType) >(&G4GDMLParser::AddAuxiliary));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::StripNamePointers() (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::StripNamePointers()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:121:17
  t178.method("StripNamePointers", static_cast<void (G4GDMLParser::*)()  const>(&G4GDMLParser::StripNamePointers));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetStripFlag(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::SetStripFlag(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:122:17
  t178.method("SetStripFlag", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetStripFlag));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetOverlapCheck(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::SetOverlapCheck(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:123:17
  t178.method("SetOverlapCheck", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetOverlapCheck));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetRegionExport(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::SetRegionExport(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:124:17
  t178.method("SetRegionExport", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetRegionExport));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetEnergyCutsExport(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::SetEnergyCutsExport(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:125:17
  t178.method("SetEnergyCutsExport", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetEnergyCutsExport));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetSDExport(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::SetSDExport(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:126:17
  t178.method("SetSDExport", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetSDExport));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetReverseSearch(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::SetReverseSearch(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:127:17
  t178.method("SetReverseSearch", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetReverseSearch));

  DEBUG_MSG("Adding wrapper for G4int G4GDMLParser::GetMaxExportLevel() (" __HERE__ ")");
  // signature to use in the veto list: G4int G4GDMLParser::GetMaxExportLevel()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:129:18
  t178.method("GetMaxExportLevel", static_cast<G4int (G4GDMLParser::*)()  const>(&G4GDMLParser::GetMaxExportLevel));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetMaxExportLevel(G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::SetMaxExportLevel(G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:130:17
  t178.method("SetMaxExportLevel", static_cast<void (G4GDMLParser::*)(G4int) >(&G4GDMLParser::SetMaxExportLevel));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::Clear() (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::Clear()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:132:17
  t178.method("Clear", static_cast<void (G4GDMLParser::*)() >(&G4GDMLParser::Clear));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::AddModule(const G4VPhysicalVolume *const) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::AddModule(const G4VPhysicalVolume *const)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:136:17
  t178.method("AddModule", static_cast<void (G4GDMLParser::*)(const G4VPhysicalVolume *const) >(&G4GDMLParser::AddModule));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::AddModule(const G4int) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::AddModule(const G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:137:17
  t178.method("AddModule", static_cast<void (G4GDMLParser::*)(const G4int) >(&G4GDMLParser::AddModule));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetAddPointerToName(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::SetAddPointerToName(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:138:17
  t178.method("SetAddPointerToName", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetAddPointerToName));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::AddVolumeAuxiliary(G4GDMLAuxStructType, const G4LogicalVolume *const) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::AddVolumeAuxiliary(G4GDMLAuxStructType, const G4LogicalVolume *const)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:139:17
  t178.method("AddVolumeAuxiliary", static_cast<void (G4GDMLParser::*)(G4GDMLAuxStructType, const G4LogicalVolume *const) >(&G4GDMLParser::AddVolumeAuxiliary));

  DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetOutputFileOverwrite(G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4GDMLParser::SetOutputFileOverwrite(G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4GDMLParser.hh:141:17
  t178.method("SetOutputFileOverwrite", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetOutputFileOverwrite));

  /* End of G4GDMLParser class method wrappers
   **********************************************************************/

}
