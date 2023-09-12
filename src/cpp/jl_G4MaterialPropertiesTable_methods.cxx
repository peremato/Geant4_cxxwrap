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
void add_methods_for_G4MaterialPropertiesTable(jlcxx::Module& types, jlcxx::TypeWrapper<G4MaterialPropertiesTable>& t40) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4MaterialPropertiesTable
   */

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::AddConstProperty(const G4String &, G4double, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::AddConstProperty(const G4String &, G4double, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:67:8
  t40.method("AddConstProperty", static_cast<void (G4MaterialPropertiesTable::*)(const G4String &, G4double, G4bool) >(&G4MaterialPropertiesTable::AddConstProperty));
  t40.method("AddConstProperty", [](G4MaterialPropertiesTable& a, const G4String & arg0, G4double arg1)->void{ a.AddConstProperty(arg0, arg1); });
  t40.method("AddConstProperty", [](G4MaterialPropertiesTable* a, const G4String & arg0, G4double arg1)->void{ a->AddConstProperty(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for G4MaterialPropertyVector * G4MaterialPropertiesTable::AddProperty(const G4String &, const std::vector<G4double> &, const std::vector<G4double> &, G4bool, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: G4MaterialPropertyVector * G4MaterialPropertiesTable::AddProperty(const G4String &, const std::vector<G4double> &, const std::vector<G4double> &, G4bool, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:73:29
  t40.method("AddProperty", static_cast<G4MaterialPropertyVector * (G4MaterialPropertiesTable::*)(const G4String &, const std::vector<G4double> &, const std::vector<G4double> &, G4bool, G4bool) >(&G4MaterialPropertiesTable::AddProperty));
  t40.method("AddProperty", [](G4MaterialPropertiesTable& a, const G4String & arg0, const std::vector<G4double> & arg1, const std::vector<G4double> & arg2)->G4MaterialPropertyVector *{ return a.AddProperty(arg0, arg1, arg2); });
  t40.method("AddProperty", [](G4MaterialPropertiesTable& a, const G4String & arg0, const std::vector<G4double> & arg1, const std::vector<G4double> & arg2, G4bool arg3)->G4MaterialPropertyVector *{ return a.AddProperty(arg0, arg1, arg2, arg3); });
  t40.method("AddProperty", [](G4MaterialPropertiesTable* a, const G4String & arg0, const std::vector<G4double> & arg1, const std::vector<G4double> & arg2)->G4MaterialPropertyVector *{ return a->AddProperty(arg0, arg1, arg2); });
  t40.method("AddProperty", [](G4MaterialPropertiesTable* a, const G4String & arg0, const std::vector<G4double> & arg1, const std::vector<G4double> & arg2, G4bool arg3)->G4MaterialPropertyVector *{ return a->AddProperty(arg0, arg1, arg2, arg3); });

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::AddProperty(const G4String &, G4MaterialPropertyVector *, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::AddProperty(const G4String &, G4MaterialPropertyVector *, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:86:8
  t40.method("AddProperty", static_cast<void (G4MaterialPropertiesTable::*)(const G4String &, G4MaterialPropertyVector *, G4bool) >(&G4MaterialPropertiesTable::AddProperty));
  t40.method("AddProperty", [](G4MaterialPropertiesTable& a, const G4String & arg0, G4MaterialPropertyVector * arg1)->void{ a.AddProperty(arg0, arg1); });
  t40.method("AddProperty", [](G4MaterialPropertiesTable* a, const G4String & arg0, G4MaterialPropertyVector * arg1)->void{ a->AddProperty(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::AddProperty(const char *, G4MaterialPropertyVector *, G4bool) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::AddProperty(const char *, G4MaterialPropertyVector *, G4bool)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:88:8
  t40.method("AddProperty", static_cast<void (G4MaterialPropertiesTable::*)(const char *, G4MaterialPropertyVector *, G4bool) >(&G4MaterialPropertiesTable::AddProperty));
  t40.method("AddProperty", [](G4MaterialPropertiesTable& a, const char * arg0, G4MaterialPropertyVector * arg1)->void{ a.AddProperty(arg0, arg1); });
  t40.method("AddProperty", [](G4MaterialPropertiesTable* a, const char * arg0, G4MaterialPropertyVector * arg1)->void{ a->AddProperty(arg0, arg1); });

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::AddProperty(const G4String &, const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::AddProperty(const G4String &, const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:93:8
  t40.method("AddProperty", static_cast<void (G4MaterialPropertiesTable::*)(const G4String &, const G4String &) >(&G4MaterialPropertiesTable::AddProperty));

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::RemoveConstProperty(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::RemoveConstProperty(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:98:8
  t40.method("RemoveConstProperty", static_cast<void (G4MaterialPropertiesTable::*)(const G4String &) >(&G4MaterialPropertiesTable::RemoveConstProperty));

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::RemoveConstProperty(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::RemoveConstProperty(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:99:8
  t40.method("RemoveConstProperty", static_cast<void (G4MaterialPropertiesTable::*)(const char *) >(&G4MaterialPropertiesTable::RemoveConstProperty));

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::RemoveProperty(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::RemoveProperty(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:102:8
  t40.method("RemoveProperty", static_cast<void (G4MaterialPropertiesTable::*)(const G4String &) >(&G4MaterialPropertiesTable::RemoveProperty));

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::RemoveProperty(const char *) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::RemoveProperty(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:103:8
  t40.method("RemoveProperty", static_cast<void (G4MaterialPropertiesTable::*)(const char *) >(&G4MaterialPropertiesTable::RemoveProperty));

  DEBUG_MSG("Adding wrapper for G4double G4MaterialPropertiesTable::GetConstProperty(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4MaterialPropertiesTable::GetConstProperty(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:106:12
  t40.method("GetConstProperty", static_cast<G4double (G4MaterialPropertiesTable::*)(const G4String &)  const>(&G4MaterialPropertiesTable::GetConstProperty));

  DEBUG_MSG("Adding wrapper for G4double G4MaterialPropertiesTable::GetConstProperty(const char *) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4MaterialPropertiesTable::GetConstProperty(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:107:12
  t40.method("GetConstProperty", static_cast<G4double (G4MaterialPropertiesTable::*)(const char *)  const>(&G4MaterialPropertiesTable::GetConstProperty));

  DEBUG_MSG("Adding wrapper for G4double G4MaterialPropertiesTable::GetConstProperty(const G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4double G4MaterialPropertiesTable::GetConstProperty(const G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:108:12
  t40.method("GetConstProperty", static_cast<G4double (G4MaterialPropertiesTable::*)(const G4int)  const>(&G4MaterialPropertiesTable::GetConstProperty));

  DEBUG_MSG("Adding wrapper for G4bool G4MaterialPropertiesTable::ConstPropertyExists(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4MaterialPropertiesTable::ConstPropertyExists(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:113:10
  t40.method("ConstPropertyExists", static_cast<G4bool (G4MaterialPropertiesTable::*)(const G4String &)  const>(&G4MaterialPropertiesTable::ConstPropertyExists));

  DEBUG_MSG("Adding wrapper for G4bool G4MaterialPropertiesTable::ConstPropertyExists(const char *) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4MaterialPropertiesTable::ConstPropertyExists(const char *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:114:10
  t40.method("ConstPropertyExists", static_cast<G4bool (G4MaterialPropertiesTable::*)(const char *)  const>(&G4MaterialPropertiesTable::ConstPropertyExists));

  DEBUG_MSG("Adding wrapper for G4bool G4MaterialPropertiesTable::ConstPropertyExists(const G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4MaterialPropertiesTable::ConstPropertyExists(const G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:115:10
  t40.method("ConstPropertyExists", static_cast<G4bool (G4MaterialPropertiesTable::*)(const G4int)  const>(&G4MaterialPropertiesTable::ConstPropertyExists));

  DEBUG_MSG("Adding wrapper for G4MaterialPropertyVector * G4MaterialPropertiesTable::GetProperty(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4MaterialPropertyVector * G4MaterialPropertiesTable::GetProperty(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:122:29
  t40.method("GetProperty", static_cast<G4MaterialPropertyVector * (G4MaterialPropertiesTable::*)(const G4String &)  const>(&G4MaterialPropertiesTable::GetProperty));

  DEBUG_MSG("Adding wrapper for G4MaterialPropertyVector * G4MaterialPropertiesTable::GetProperty(const G4int) (" __HERE__ ")");
  // signature to use in the veto list: G4MaterialPropertyVector * G4MaterialPropertiesTable::GetProperty(const G4int)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:123:29
  t40.method("GetProperty", static_cast<G4MaterialPropertyVector * (G4MaterialPropertiesTable::*)(const G4int)  const>(&G4MaterialPropertiesTable::GetProperty));

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::AddEntry(const G4String &, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::AddEntry(const G4String &, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:127:8
  t40.method("AddEntry", static_cast<void (G4MaterialPropertiesTable::*)(const G4String &, G4double, G4double) >(&G4MaterialPropertiesTable::AddEntry));

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::AddEntry(const char *, G4double, G4double) (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::AddEntry(const char *, G4double, G4double)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:129:8
  t40.method("AddEntry", static_cast<void (G4MaterialPropertiesTable::*)(const char *, G4double, G4double) >(&G4MaterialPropertiesTable::AddEntry));

  DEBUG_MSG("Adding wrapper for G4int G4MaterialPropertiesTable::GetConstPropertyIndex(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4MaterialPropertiesTable::GetConstPropertyIndex(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:133:9
  t40.method("GetConstPropertyIndex", static_cast<G4int (G4MaterialPropertiesTable::*)(const G4String &)  const>(&G4MaterialPropertiesTable::GetConstPropertyIndex));

  DEBUG_MSG("Adding wrapper for G4int G4MaterialPropertiesTable::GetPropertyIndex(const G4String &) (" __HERE__ ")");
  // signature to use in the veto list: G4int G4MaterialPropertiesTable::GetPropertyIndex(const G4String &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:138:9
  t40.method("GetPropertyIndex", static_cast<G4int (G4MaterialPropertiesTable::*)(const G4String &)  const>(&G4MaterialPropertiesTable::GetPropertyIndex));

  DEBUG_MSG("Adding wrapper for void G4MaterialPropertiesTable::DumpTable() (" __HERE__ ")");
  // signature to use in the veto list: void G4MaterialPropertiesTable::DumpTable()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:143:8
  t40.method("DumpTable", static_cast<void (G4MaterialPropertiesTable::*)()  const>(&G4MaterialPropertiesTable::DumpTable));

  DEBUG_MSG("Adding wrapper for const std::vector<G4String> & G4MaterialPropertiesTable::GetMaterialPropertyNames() (" __HERE__ ")");
  // signature to use in the veto list: const std::vector<G4String> & G4MaterialPropertiesTable::GetMaterialPropertyNames()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:147:32
  t40.method("GetMaterialPropertyNames", static_cast<const std::vector<G4String> & (G4MaterialPropertiesTable::*)()  const>(&G4MaterialPropertiesTable::GetMaterialPropertyNames));

  DEBUG_MSG("Adding wrapper for const std::vector<G4String> & G4MaterialPropertiesTable::GetMaterialConstPropertyNames() (" __HERE__ ")");
  // signature to use in the veto list: const std::vector<G4String> & G4MaterialPropertiesTable::GetMaterialConstPropertyNames()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:151:32
  t40.method("GetMaterialConstPropertyNames", static_cast<const std::vector<G4String> & (G4MaterialPropertiesTable::*)()  const>(&G4MaterialPropertiesTable::GetMaterialConstPropertyNames));

  DEBUG_MSG("Adding wrapper for const std::vector<G4MaterialPropertyVector *> & G4MaterialPropertiesTable::GetProperties() (" __HERE__ ")");
  // signature to use in the veto list: const std::vector<G4MaterialPropertyVector *> & G4MaterialPropertiesTable::GetProperties()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4MaterialPropertiesTable.hh:155:49
  t40.method("GetProperties", static_cast<const std::vector<G4MaterialPropertyVector *> & (G4MaterialPropertiesTable::*)()  const>(&G4MaterialPropertiesTable::GetProperties));

  /* End of G4MaterialPropertiesTable class method wrappers
   **********************************************************************/

}
