// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4BooleanSolid> : std::false_type { };
  template<> struct DefaultConstructible<G4BooleanSolid> : std::false_type { };
template<> struct SuperType<G4BooleanSolid> { typedef G4VSolid type; };
}

// Class generating the wrapper for type G4BooleanSolid
// signature to use in the veto file: G4BooleanSolid
struct JlG4BooleanSolid: public Wrapper {

  JlG4BooleanSolid(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4BooleanSolid (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:49:7
    jlcxx::TypeWrapper<G4BooleanSolid>  t = jlModule.add_type<G4BooleanSolid>("G4BooleanSolid",
      jlcxx::julia_base_type<G4VSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4BooleanSolid>>(new jlcxx::TypeWrapper<G4BooleanSolid>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;




    DEBUG_MSG("Adding wrapper for const G4VSolid * G4BooleanSolid::GetConstituentSolid(G4int) (" __HERE__ ")");
    // signature to use in the veto list: const G4VSolid * G4BooleanSolid::GetConstituentSolid(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:70:21
    t.method("GetConstituentSolid", static_cast<const G4VSolid * (G4BooleanSolid::*)(G4int)  const>(&G4BooleanSolid::GetConstituentSolid));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4BooleanSolid::GetConstituentSolid(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4BooleanSolid::GetConstituentSolid(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:71:21
    t.method("GetConstituentSolid", static_cast<G4VSolid * (G4BooleanSolid::*)(G4int) >(&G4BooleanSolid::GetConstituentSolid));

    DEBUG_MSG("Adding wrapper for G4double G4BooleanSolid::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4BooleanSolid::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:76:14
    t.method("GetCubicVolume", static_cast<G4double (G4BooleanSolid::*)() >(&G4BooleanSolid::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4BooleanSolid::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4BooleanSolid::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:77:21
    t.method("GetSurfaceArea", static_cast<G4double (G4BooleanSolid::*)() >(&G4BooleanSolid::GetSurfaceArea));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4BooleanSolid::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4BooleanSolid::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:79:20
    t.method("GetEntityType", static_cast<G4GeometryType (G4BooleanSolid::*)()  const>(&G4BooleanSolid::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4BooleanSolid::GetPolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4BooleanSolid::GetPolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:80:19
    t.method("GetPolyhedron", static_cast<G4Polyhedron * (G4BooleanSolid::*)()  const>(&G4BooleanSolid::GetPolyhedron));

    DEBUG_MSG("Adding wrapper for G4int G4BooleanSolid::GetCubVolStatistics() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4BooleanSolid::GetCubVolStatistics()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:84:18
    t.method("GetCubVolStatistics", static_cast<G4int (G4BooleanSolid::*)()  const>(&G4BooleanSolid::GetCubVolStatistics));

    DEBUG_MSG("Adding wrapper for G4double G4BooleanSolid::GetCubVolEpsilon() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4BooleanSolid::GetCubVolEpsilon()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:85:21
    t.method("GetCubVolEpsilon", static_cast<G4double (G4BooleanSolid::*)()  const>(&G4BooleanSolid::GetCubVolEpsilon));

    DEBUG_MSG("Adding wrapper for void G4BooleanSolid::SetCubVolStatistics(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4BooleanSolid::SetCubVolStatistics(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:86:10
    t.method("SetCubVolStatistics", static_cast<void (G4BooleanSolid::*)(G4int) >(&G4BooleanSolid::SetCubVolStatistics));

    DEBUG_MSG("Adding wrapper for void G4BooleanSolid::SetCubVolEpsilon(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4BooleanSolid::SetCubVolEpsilon(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:87:10
    t.method("SetCubVolEpsilon", static_cast<void (G4BooleanSolid::*)(G4double) >(&G4BooleanSolid::SetCubVolEpsilon));

    DEBUG_MSG("Adding wrapper for G4int G4BooleanSolid::GetAreaStatistics() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4BooleanSolid::GetAreaStatistics()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:89:18
    t.method("GetAreaStatistics", static_cast<G4int (G4BooleanSolid::*)()  const>(&G4BooleanSolid::GetAreaStatistics));

    DEBUG_MSG("Adding wrapper for G4double G4BooleanSolid::GetAreaAccuracy() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4BooleanSolid::GetAreaAccuracy()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:90:21
    t.method("GetAreaAccuracy", static_cast<G4double (G4BooleanSolid::*)()  const>(&G4BooleanSolid::GetAreaAccuracy));

    DEBUG_MSG("Adding wrapper for void G4BooleanSolid::SetAreaStatistics(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4BooleanSolid::SetAreaStatistics(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:91:17
    t.method("SetAreaStatistics", static_cast<void (G4BooleanSolid::*)(G4int) >(&G4BooleanSolid::SetAreaStatistics));

    DEBUG_MSG("Adding wrapper for void G4BooleanSolid::SetAreaAccuracy(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4BooleanSolid::SetAreaAccuracy(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:92:17
    t.method("SetAreaAccuracy", static_cast<void (G4BooleanSolid::*)(G4double) >(&G4BooleanSolid::SetAreaAccuracy));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4BooleanSolid::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4BooleanSolid::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:94:19
    t.method("GetPointOnSurface", static_cast<G4ThreeVector (G4BooleanSolid::*)()  const>(&G4BooleanSolid::GetPointOnSurface));

    DEBUG_MSG("Adding wrapper for G4int G4BooleanSolid::GetNumOfConstituents() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4BooleanSolid::GetNumOfConstituents()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:96:11
    t.method("GetNumOfConstituents", static_cast<G4int (G4BooleanSolid::*)()  const>(&G4BooleanSolid::GetNumOfConstituents));

    DEBUG_MSG("Adding wrapper for G4bool G4BooleanSolid::IsFaceted() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4BooleanSolid::IsFaceted()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:97:12
    t.method("IsFaceted", static_cast<G4bool (G4BooleanSolid::*)()  const>(&G4BooleanSolid::IsFaceted));


    DEBUG_MSG("Adding wrapper for G4BooleanSolid & G4BooleanSolid::operator=(const G4BooleanSolid &) (" __HERE__ ")");
    // signature to use in the veto list: G4BooleanSolid & G4BooleanSolid::operator=(const G4BooleanSolid &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:105:21
    t.method("assign", static_cast<G4BooleanSolid & (G4BooleanSolid::*)(const G4BooleanSolid &) >(&G4BooleanSolid::operator=));

    DEBUG_MSG("Adding wrapper for void G4BooleanSolid::SetExternalBooleanProcessor(G4VBooleanProcessor *) (" __HERE__ ")");
    // signature to use in the veto list: void G4BooleanSolid::SetExternalBooleanProcessor(G4VBooleanProcessor *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:108:17
    module_.method("G4BooleanSolid!SetExternalBooleanProcessor", static_cast<void (*)(G4VBooleanProcessor *) >(&G4BooleanSolid::SetExternalBooleanProcessor));

    DEBUG_MSG("Adding wrapper for G4VBooleanProcessor * G4BooleanSolid::GetExternalBooleanProcessor() (" __HERE__ ")");
    // signature to use in the veto list: G4VBooleanProcessor * G4BooleanSolid::GetExternalBooleanProcessor()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4BooleanSolid.hh:110:33
    module_.method("G4BooleanSolid!GetExternalBooleanProcessor", static_cast<G4VBooleanProcessor * (*)() >(&G4BooleanSolid::GetExternalBooleanProcessor));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4BooleanSolid>> type_;
};
std::shared_ptr<Wrapper> newJlG4BooleanSolid(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4BooleanSolid(module));
}
