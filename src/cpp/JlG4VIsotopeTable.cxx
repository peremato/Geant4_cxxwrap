// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VIsotopeTable> : std::false_type { };
  template<> struct DefaultConstructible<G4VIsotopeTable> : std::false_type { };
}

// Class generating the wrapper for type G4VIsotopeTable
// signature to use in the veto file: G4VIsotopeTable
struct JlG4VIsotopeTable: public Wrapper {

  JlG4VIsotopeTable(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VIsotopeTable (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VIsotopeTable.hh:44:7
    jlcxx::TypeWrapper<G4VIsotopeTable>  t = jlModule.add_type<G4VIsotopeTable>("G4VIsotopeTable");
    jlcxx::stl::apply_stl<G4VIsotopeTable*>(jlModule);
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VIsotopeTable>>(new jlcxx::TypeWrapper<G4VIsotopeTable>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;



    DEBUG_MSG("Adding wrapper for G4VIsotopeTable & G4VIsotopeTable::operator=(const G4VIsotopeTable &) (" __HERE__ ")");
    // signature to use in the veto list: G4VIsotopeTable & G4VIsotopeTable::operator=(const G4VIsotopeTable &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VIsotopeTable.hh:53:22
    t.method("assign", static_cast<G4VIsotopeTable & (G4VIsotopeTable::*)(const G4VIsotopeTable &) >(&G4VIsotopeTable::operator=));

    DEBUG_MSG("Adding wrapper for G4int G4VIsotopeTable::GetVerboseLevel() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4VIsotopeTable::GetVerboseLevel()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VIsotopeTable.hh:81:11
    t.method("GetVerboseLevel", static_cast<G4int (G4VIsotopeTable::*)()  const>(&G4VIsotopeTable::GetVerboseLevel));

    DEBUG_MSG("Adding wrapper for void G4VIsotopeTable::SetVerboseLevel(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4VIsotopeTable::SetVerboseLevel(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VIsotopeTable.hh:82:10
    t.method("SetVerboseLevel", static_cast<void (G4VIsotopeTable::*)(G4int) >(&G4VIsotopeTable::SetVerboseLevel));

    DEBUG_MSG("Adding wrapper for void G4VIsotopeTable::DumpTable(G4int, G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4VIsotopeTable::DumpTable(G4int, G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VIsotopeTable.hh:85:10
    t.method("DumpTable", static_cast<void (G4VIsotopeTable::*)(G4int, G4int) >(&G4VIsotopeTable::DumpTable));
    t.method("DumpTable", [](G4VIsotopeTable& a)->void { a.DumpTable(); });
    t.method("DumpTable", [](G4VIsotopeTable& a, G4int arg0)->void { a.DumpTable(arg0); });
    t.method("DumpTable", [](G4VIsotopeTable* a)->void { a->DumpTable(); });
    t.method("DumpTable", [](G4VIsotopeTable* a, G4int arg0)->void { a->DumpTable(arg0); });

    DEBUG_MSG("Adding wrapper for const G4String & G4VIsotopeTable::GetName() (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4VIsotopeTable::GetName()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VIsotopeTable.hh:87:21
    t.method("GetName", static_cast<const G4String & (G4VIsotopeTable::*)()  const>(&G4VIsotopeTable::GetName));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VIsotopeTable>> type_;
};
std::shared_ptr<Wrapper> newJlG4VIsotopeTable(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VIsotopeTable(module));
}