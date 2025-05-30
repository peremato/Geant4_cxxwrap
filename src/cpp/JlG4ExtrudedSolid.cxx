// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4ExtrudedSolid> : std::false_type { };
  template<> struct DefaultConstructible<G4ExtrudedSolid> : std::false_type { };
}

// Class generating the wrapper for type G4ExtrudedSolid
// signature to use in the veto file: G4ExtrudedSolid
struct JlG4ExtrudedSolid: public Wrapper {

  JlG4ExtrudedSolid(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4ExtrudedSolid (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:71:7
    jlcxx::TypeWrapper<G4ExtrudedSolid>  t = jlModule.add_type<G4ExtrudedSolid>("G4ExtrudedSolid");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4ExtrudedSolid>>(new jlcxx::TypeWrapper<G4ExtrudedSolid>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4ExtrudedSolid::G4ExtrudedSolid(const G4String &, const std::vector<G4TwoVector> &, const std::vector<G4ExtrudedSolid::ZSection> &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:87:5
    t.constructor<const G4String &, const std::vector<G4TwoVector> &, const std::vector<G4ExtrudedSolid::ZSection> &>(/*finalize=*/jlcxx::finalize_policy::yes);


    DEBUG_MSG("Adding wrapper for void G4ExtrudedSolid::G4ExtrudedSolid(const G4String &, const std::vector<G4TwoVector> &, G4double, const G4TwoVector &, G4double, const G4TwoVector &, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:92:5
    t.constructor<const G4String &, const std::vector<G4TwoVector> &, G4double>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const G4String &, const std::vector<G4TwoVector> &, G4double, const G4TwoVector &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const G4String &, const std::vector<G4TwoVector> &, G4double, const G4TwoVector &, G4double>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const G4String &, const std::vector<G4TwoVector> &, G4double, const G4TwoVector &, G4double, const G4TwoVector &>(/*finalize=*/jlcxx::finalize_policy::yes);
    t.constructor<const G4String &, const std::vector<G4TwoVector> &, G4double, const G4TwoVector &, G4double, const G4TwoVector &, G4double>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for G4int G4ExtrudedSolid::GetNofVertices() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ExtrudedSolid::GetNofVertices()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:106:24
    t.method("GetNofVertices", static_cast<G4int (G4ExtrudedSolid::*)()  const>(&G4ExtrudedSolid::GetNofVertices));

    DEBUG_MSG("Adding wrapper for G4TwoVector G4ExtrudedSolid::GetVertex(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4TwoVector G4ExtrudedSolid::GetVertex(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:107:24
    t.method("GetVertex", static_cast<G4TwoVector (G4ExtrudedSolid::*)(G4int)  const>(&G4ExtrudedSolid::GetVertex));

    DEBUG_MSG("Adding wrapper for std::vector<G4TwoVector> G4ExtrudedSolid::GetPolygon() (" __HERE__ ")");
    // signature to use in the veto list: std::vector<G4TwoVector> G4ExtrudedSolid::GetPolygon()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:108:37
    t.method("GetPolygon", static_cast<std::vector<G4TwoVector> (G4ExtrudedSolid::*)()  const>(&G4ExtrudedSolid::GetPolygon));

    DEBUG_MSG("Adding wrapper for G4int G4ExtrudedSolid::GetNofZSections() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4ExtrudedSolid::GetNofZSections()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:110:24
    t.method("GetNofZSections", static_cast<G4int (G4ExtrudedSolid::*)()  const>(&G4ExtrudedSolid::GetNofZSections));

    DEBUG_MSG("Adding wrapper for G4ExtrudedSolid::ZSection G4ExtrudedSolid::GetZSection(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4ExtrudedSolid::ZSection G4ExtrudedSolid::GetZSection(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:111:24
    t.method("GetZSection", static_cast<G4ExtrudedSolid::ZSection (G4ExtrudedSolid::*)(G4int)  const>(&G4ExtrudedSolid::GetZSection));

    DEBUG_MSG("Adding wrapper for std::vector<G4ExtrudedSolid::ZSection> G4ExtrudedSolid::GetZSections() (" __HERE__ ")");
    // signature to use in the veto list: std::vector<G4ExtrudedSolid::ZSection> G4ExtrudedSolid::GetZSections()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:112:34
    t.method("GetZSections", static_cast<std::vector<G4ExtrudedSolid::ZSection> (G4ExtrudedSolid::*)()  const>(&G4ExtrudedSolid::GetZSections));

    DEBUG_MSG("Adding wrapper for EInside G4ExtrudedSolid::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4ExtrudedSolid::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:116:14
    t.method("Inside", static_cast<EInside (G4ExtrudedSolid::*)(const G4ThreeVector &)  const>(&G4ExtrudedSolid::Inside));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4ExtrudedSolid::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4ExtrudedSolid::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:117:19
    t.method("SurfaceNormal", static_cast<G4ThreeVector (G4ExtrudedSolid::*)(const G4ThreeVector &)  const>(&G4ExtrudedSolid::SurfaceNormal));

    DEBUG_MSG("Adding wrapper for G4double G4ExtrudedSolid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4ExtrudedSolid::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:118:14
    t.method("DistanceToIn", static_cast<G4double (G4ExtrudedSolid::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4ExtrudedSolid::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4ExtrudedSolid::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4ExtrudedSolid::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:120:14
    t.method("DistanceToIn", static_cast<G4double (G4ExtrudedSolid::*)(const G4ThreeVector &)  const>(&G4ExtrudedSolid::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4ExtrudedSolid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4ExtrudedSolid::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:121:14
    t.method("DistanceToOut", static_cast<G4double (G4ExtrudedSolid::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4ExtrudedSolid::DistanceToOut));
    t.method("DistanceToOut", [](G4ExtrudedSolid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4ExtrudedSolid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4ExtrudedSolid const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); });
    t.method("DistanceToOut", [](G4ExtrudedSolid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4ExtrudedSolid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4ExtrudedSolid const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for G4double G4ExtrudedSolid::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4ExtrudedSolid::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:126:14
    t.method("DistanceToOut", static_cast<G4double (G4ExtrudedSolid::*)(const G4ThreeVector &)  const>(&G4ExtrudedSolid::DistanceToOut));

    DEBUG_MSG("Adding wrapper for void G4ExtrudedSolid::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4ExtrudedSolid::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:128:10
    t.method("BoundingLimits", static_cast<void (G4ExtrudedSolid::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4ExtrudedSolid::BoundingLimits));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4ExtrudedSolid::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4ExtrudedSolid::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:133:20
    t.method("GetEntityType", static_cast<G4GeometryType (G4ExtrudedSolid::*)()  const>(&G4ExtrudedSolid::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4bool G4ExtrudedSolid::IsFaceted() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4ExtrudedSolid::IsFaceted()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:134:12
    t.method("IsFaceted", static_cast<G4bool (G4ExtrudedSolid::*)()  const>(&G4ExtrudedSolid::IsFaceted));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4ExtrudedSolid::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4ExtrudedSolid::Clone()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:135:15
    t.method("Clone", static_cast<G4VSolid * (G4ExtrudedSolid::*)()  const>(&G4ExtrudedSolid::Clone));


    DEBUG_MSG("Adding wrapper for void G4ExtrudedSolid::G4ExtrudedSolid(const G4ExtrudedSolid &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:144:5
    t.constructor<const G4ExtrudedSolid &>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for G4ExtrudedSolid & G4ExtrudedSolid::operator=(const G4ExtrudedSolid &) (" __HERE__ ")");
    // signature to use in the veto list: G4ExtrudedSolid & G4ExtrudedSolid::operator=(const G4ExtrudedSolid &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4ExtrudedSolid.hh:145:22
    t.method("assign", static_cast<G4ExtrudedSolid & (G4ExtrudedSolid::*)(const G4ExtrudedSolid &) >(&G4ExtrudedSolid::operator=));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4ExtrudedSolid>> type_;
};
std::shared_ptr<Wrapper> newJlG4ExtrudedSolid(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4ExtrudedSolid(module));
}
