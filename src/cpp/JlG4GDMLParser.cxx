// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4GDMLParser> : std::false_type { };
  template<> struct DefaultConstructible<G4GDMLParser> : std::false_type { };
}

// Class generating the wrapper for type G4GDMLParser
// signature to use in the veto file: G4GDMLParser
struct JlG4GDMLParser: public Wrapper {

  JlG4GDMLParser(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4GDMLParser (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:57:7
    jlcxx::TypeWrapper<G4GDMLParser>  t = jlModule.add_type<G4GDMLParser>("G4GDMLParser");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4GDMLParser>>(new jlcxx::TypeWrapper<G4GDMLParser>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void G4GDMLParser::G4GDMLParser(G4GDMLReadStructure *) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:62:5
    t.constructor<G4GDMLReadStructure *>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void G4GDMLParser::G4GDMLParser(G4GDMLReadStructure *, G4GDMLWriteStructure *) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:63:5
    t.constructor<G4GDMLReadStructure *, G4GDMLWriteStructure *>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::Read(const G4String &, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::Read(const G4String &, G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:68:17
    t.method("Read", static_cast<void (G4GDMLParser::*)(const G4String &, G4bool) >(&G4GDMLParser::Read));
    t.method("Read", [](G4GDMLParser& a, const G4String & arg0)->void { a.Read(arg0); });
    t.method("Read", [](G4GDMLParser* a, const G4String & arg0)->void { a->Read(arg0); });

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::ReadModule(const G4String &, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::ReadModule(const G4String &, G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:74:17
    t.method("ReadModule", static_cast<void (G4GDMLParser::*)(const G4String &, G4bool) >(&G4GDMLParser::ReadModule));
    t.method("ReadModule", [](G4GDMLParser& a, const G4String & arg0)->void { a.ReadModule(arg0); });
    t.method("ReadModule", [](G4GDMLParser* a, const G4String & arg0)->void { a->ReadModule(arg0); });

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::Write(const G4String &, const G4VPhysicalVolume *, G4bool, const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::Write(const G4String &, const G4VPhysicalVolume *, G4bool, const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:80:17
    t.method("Write", static_cast<void (G4GDMLParser::*)(const G4String &, const G4VPhysicalVolume *, G4bool, const G4String &) >(&G4GDMLParser::Write));
    t.method("Write", [](G4GDMLParser& a, const G4String & arg0)->void { a.Write(arg0); });
    t.method("Write", [](G4GDMLParser& a, const G4String & arg0, const G4VPhysicalVolume * arg1)->void { a.Write(arg0, arg1); });
    t.method("Write", [](G4GDMLParser& a, const G4String & arg0, const G4VPhysicalVolume * arg1, G4bool arg2)->void { a.Write(arg0, arg1, arg2); });
    t.method("Write", [](G4GDMLParser* a, const G4String & arg0)->void { a->Write(arg0); });
    t.method("Write", [](G4GDMLParser* a, const G4String & arg0, const G4VPhysicalVolume * arg1)->void { a->Write(arg0, arg1); });
    t.method("Write", [](G4GDMLParser* a, const G4String & arg0, const G4VPhysicalVolume * arg1, G4bool arg2)->void { a->Write(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::Write(const G4String &, const G4LogicalVolume *, G4bool, const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::Write(const G4String &, const G4LogicalVolume *, G4bool, const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:91:17
    t.method("Write", static_cast<void (G4GDMLParser::*)(const G4String &, const G4LogicalVolume *, G4bool, const G4String &) >(&G4GDMLParser::Write));
    t.method("Write", [](G4GDMLParser& a, const G4String & arg0, const G4LogicalVolume * arg1)->void { a.Write(arg0, arg1); });
    t.method("Write", [](G4GDMLParser& a, const G4String & arg0, const G4LogicalVolume * arg1, G4bool arg2)->void { a.Write(arg0, arg1, arg2); });
    t.method("Write", [](G4GDMLParser* a, const G4String & arg0, const G4LogicalVolume * arg1)->void { a->Write(arg0, arg1); });
    t.method("Write", [](G4GDMLParser* a, const G4String & arg0, const G4LogicalVolume * arg1, G4bool arg2)->void { a->Write(arg0, arg1, arg2); });

    DEBUG_MSG("Adding wrapper for G4LogicalVolume * G4GDMLParser::ParseST(const G4String &, G4Material *, G4Material *) (" __HERE__ ")");
    // signature to use in the veto list: G4LogicalVolume * G4GDMLParser::ParseST(const G4String &, G4Material *, G4Material *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:102:29
    t.method("ParseST", static_cast<G4LogicalVolume * (G4GDMLParser::*)(const G4String &, G4Material *, G4Material *) >(&G4GDMLParser::ParseST));

    DEBUG_MSG("Adding wrapper for G4bool G4GDMLParser::IsValid(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4GDMLParser::IsValid(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:112:19
    t.method("IsValid", static_cast<G4bool (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::IsValid));

    DEBUG_MSG("Adding wrapper for G4double G4GDMLParser::GetConstant(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GDMLParser::GetConstant(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:113:21
    t.method("GetConstant", static_cast<G4double (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetConstant));

    DEBUG_MSG("Adding wrapper for G4double G4GDMLParser::GetVariable(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GDMLParser::GetVariable(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:114:21
    t.method("GetVariable", static_cast<G4double (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetVariable));

    DEBUG_MSG("Adding wrapper for G4double G4GDMLParser::GetQuantity(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4GDMLParser::GetQuantity(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:115:21
    t.method("GetQuantity", static_cast<G4double (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetQuantity));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4GDMLParser::GetPosition(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4GDMLParser::GetPosition(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:116:26
    t.method("GetPosition", static_cast<G4ThreeVector (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetPosition));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4GDMLParser::GetRotation(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4GDMLParser::GetRotation(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:117:26
    t.method("GetRotation", static_cast<G4ThreeVector (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetRotation));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4GDMLParser::GetScale(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4GDMLParser::GetScale(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:118:26
    t.method("GetScale", static_cast<G4ThreeVector (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetScale));

    DEBUG_MSG("Adding wrapper for G4GDMLMatrix G4GDMLParser::GetMatrix(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4GDMLMatrix G4GDMLParser::GetMatrix(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:119:25
    t.method("GetMatrix", static_cast<G4GDMLMatrix (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetMatrix));

    DEBUG_MSG("Adding wrapper for G4LogicalVolume * G4GDMLParser::GetVolume(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4LogicalVolume * G4GDMLParser::GetVolume(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:120:29
    t.method("GetVolume", static_cast<G4LogicalVolume * (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetVolume));

    DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4GDMLParser::GetPhysVolume(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4VPhysicalVolume * G4GDMLParser::GetPhysVolume(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:121:31
    t.method("GetPhysVolume", static_cast<G4VPhysicalVolume * (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetPhysVolume));

    DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4GDMLParser::GetWorldVolume(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: G4VPhysicalVolume * G4GDMLParser::GetWorldVolume(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:123:12
    t.method("GetWorldVolume", static_cast<G4VPhysicalVolume * (G4GDMLParser::*)(const G4String &)  const>(&G4GDMLParser::GetWorldVolume));
    t.method("GetWorldVolume", [](G4GDMLParser const& a)->G4VPhysicalVolume * { return a.GetWorldVolume(); });
    t.method("GetWorldVolume", [](G4GDMLParser const* a)->G4VPhysicalVolume * { return a->GetWorldVolume(); });

    DEBUG_MSG("Adding wrapper for G4GDMLAuxListType G4GDMLParser::GetVolumeAuxiliaryInformation(G4LogicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: G4GDMLAuxListType G4GDMLParser::GetVolumeAuxiliaryInformation(G4LogicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:125:12
    t.method("GetVolumeAuxiliaryInformation", static_cast<G4GDMLAuxListType (G4GDMLParser::*)(G4LogicalVolume *)  const>(&G4GDMLParser::GetVolumeAuxiliaryInformation));

    DEBUG_MSG("Adding wrapper for const G4GDMLAuxListType * G4GDMLParser::GetAuxList() (" __HERE__ ")");
    // signature to use in the veto list: const G4GDMLAuxListType * G4GDMLParser::GetAuxList()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:127:37
    t.method("GetAuxList", static_cast<const G4GDMLAuxListType * (G4GDMLParser::*)()  const>(&G4GDMLParser::GetAuxList));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::AddAuxiliary(G4GDMLAuxStructType) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::AddAuxiliary(G4GDMLAuxStructType)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:128:17
    t.method("AddAuxiliary", static_cast<void (G4GDMLParser::*)(G4GDMLAuxStructType) >(&G4GDMLParser::AddAuxiliary));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::StripNamePointers() (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::StripNamePointers()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:129:17
    t.method("StripNamePointers", static_cast<void (G4GDMLParser::*)()  const>(&G4GDMLParser::StripNamePointers));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetStripFlag(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetStripFlag(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:130:17
    t.method("SetStripFlag", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetStripFlag));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetOverlapCheck(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetOverlapCheck(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:131:17
    t.method("SetOverlapCheck", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetOverlapCheck));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetRegionExport(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetRegionExport(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:132:17
    t.method("SetRegionExport", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetRegionExport));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetEnergyCutsExport(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetEnergyCutsExport(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:133:17
    t.method("SetEnergyCutsExport", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetEnergyCutsExport));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetSDExport(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetSDExport(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:134:17
    t.method("SetSDExport", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetSDExport));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetReverseSearch(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetReverseSearch(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:135:17
    t.method("SetReverseSearch", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetReverseSearch));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetImportSchema(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetImportSchema(const G4String &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:136:17
    t.method("SetImportSchema", static_cast<void (G4GDMLParser::*)(const G4String &) >(&G4GDMLParser::SetImportSchema));

    DEBUG_MSG("Adding wrapper for G4int G4GDMLParser::GetMaxExportLevel() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4GDMLParser::GetMaxExportLevel()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:138:18
    t.method("GetMaxExportLevel", static_cast<G4int (G4GDMLParser::*)()  const>(&G4GDMLParser::GetMaxExportLevel));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetMaxExportLevel(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetMaxExportLevel(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:139:17
    t.method("SetMaxExportLevel", static_cast<void (G4GDMLParser::*)(G4int) >(&G4GDMLParser::SetMaxExportLevel));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::Clear() (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::Clear()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:141:17
    t.method("Clear", static_cast<void (G4GDMLParser::*)() >(&G4GDMLParser::Clear));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::AddModule(const G4VPhysicalVolume *const) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::AddModule(const G4VPhysicalVolume *const)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:145:17
    t.method("AddModule", static_cast<void (G4GDMLParser::*)(const G4VPhysicalVolume *const) >(&G4GDMLParser::AddModule));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::AddModule(const G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::AddModule(const G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:146:17
    t.method("AddModule", static_cast<void (G4GDMLParser::*)(const G4int) >(&G4GDMLParser::AddModule));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetAddPointerToName(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetAddPointerToName(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:147:17
    t.method("SetAddPointerToName", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetAddPointerToName));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::AddVolumeAuxiliary(G4GDMLAuxStructType, const G4LogicalVolume *const) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::AddVolumeAuxiliary(G4GDMLAuxStructType, const G4LogicalVolume *const)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:148:17
    t.method("AddVolumeAuxiliary", static_cast<void (G4GDMLParser::*)(G4GDMLAuxStructType, const G4LogicalVolume *const) >(&G4GDMLParser::AddVolumeAuxiliary));

    DEBUG_MSG("Adding wrapper for void G4GDMLParser::SetOutputFileOverwrite(G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4GDMLParser::SetOutputFileOverwrite(G4bool)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4GDMLParser.hh:150:17
    t.method("SetOutputFileOverwrite", static_cast<void (G4GDMLParser::*)(G4bool) >(&G4GDMLParser::SetOutputFileOverwrite));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4GDMLParser>> type_;
};
std::shared_ptr<Wrapper> newJlG4GDMLParser(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4GDMLParser(module));
}
