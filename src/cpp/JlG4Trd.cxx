// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Trd> : std::false_type { };
  template<> struct DefaultConstructible<G4Trd> : std::false_type { };
template<> struct SuperType<G4Trd> { typedef G4CSGSolid type; };
}

// Class generating the wrapper for type G4Trd
// signature to use in the veto file: G4Trd
struct JlG4Trd: public Wrapper {

  JlG4Trd(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Trd (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:62:7
    jlcxx::TypeWrapper<G4Trd>  t = jlModule.add_type<G4Trd>("G4Trd",
      jlcxx::julia_base_type<G4CSGSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Trd>>(new jlcxx::TypeWrapper<G4Trd>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4Trd::G4Trd(const G4String &, G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:66:5
    t.constructor<const G4String &, G4double, G4double, G4double, G4double, G4double>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("pName"), jlcxx::arg("pdx1"), jlcxx::arg("pdx2"), jlcxx::arg("pdy1"), jlcxx::arg("pdy2"), jlcxx::arg("pdz")    );

    DEBUG_MSG("Adding wrapper for G4double G4Trd::GetXHalfLength1() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::GetXHalfLength1()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:79:21
    t.method("GetXHalfLength1", [](G4Trd const& a)->G4double { return a.GetXHalfLength1(); }, jlcxx::arg("this"));
    t.method("GetXHalfLength1", [](G4Trd const* a)->G4double { return a->GetXHalfLength1(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::GetXHalfLength2() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::GetXHalfLength2()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:80:21
    t.method("GetXHalfLength2", [](G4Trd const& a)->G4double { return a.GetXHalfLength2(); }, jlcxx::arg("this"));
    t.method("GetXHalfLength2", [](G4Trd const* a)->G4double { return a->GetXHalfLength2(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::GetYHalfLength1() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::GetYHalfLength1()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:81:21
    t.method("GetYHalfLength1", [](G4Trd const& a)->G4double { return a.GetYHalfLength1(); }, jlcxx::arg("this"));
    t.method("GetYHalfLength1", [](G4Trd const* a)->G4double { return a->GetYHalfLength1(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::GetYHalfLength2() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::GetYHalfLength2()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:82:21
    t.method("GetYHalfLength2", [](G4Trd const& a)->G4double { return a.GetYHalfLength2(); }, jlcxx::arg("this"));
    t.method("GetYHalfLength2", [](G4Trd const* a)->G4double { return a->GetYHalfLength2(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::GetZHalfLength() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::GetZHalfLength()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:83:21
    t.method("GetZHalfLength", [](G4Trd const& a)->G4double { return a.GetZHalfLength(); }, jlcxx::arg("this"));
    t.method("GetZHalfLength", [](G4Trd const* a)->G4double { return a->GetZHalfLength(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Trd::SetXHalfLength1(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Trd::SetXHalfLength1(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:87:17
    t.method("SetXHalfLength1", [](G4Trd& a, G4double arg0)->void { a.SetXHalfLength1(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));
    t.method("SetXHalfLength1", [](G4Trd* a, G4double arg0)->void { a->SetXHalfLength1(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));

    DEBUG_MSG("Adding wrapper for void G4Trd::SetXHalfLength2(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Trd::SetXHalfLength2(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:88:17
    t.method("SetXHalfLength2", [](G4Trd& a, G4double arg0)->void { a.SetXHalfLength2(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));
    t.method("SetXHalfLength2", [](G4Trd* a, G4double arg0)->void { a->SetXHalfLength2(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));

    DEBUG_MSG("Adding wrapper for void G4Trd::SetYHalfLength1(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Trd::SetYHalfLength1(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:89:17
    t.method("SetYHalfLength1", [](G4Trd& a, G4double arg0)->void { a.SetYHalfLength1(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));
    t.method("SetYHalfLength1", [](G4Trd* a, G4double arg0)->void { a->SetYHalfLength1(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));

    DEBUG_MSG("Adding wrapper for void G4Trd::SetYHalfLength2(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Trd::SetYHalfLength2(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:90:17
    t.method("SetYHalfLength2", [](G4Trd& a, G4double arg0)->void { a.SetYHalfLength2(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));
    t.method("SetYHalfLength2", [](G4Trd* a, G4double arg0)->void { a->SetYHalfLength2(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));

    DEBUG_MSG("Adding wrapper for void G4Trd::SetZHalfLength(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Trd::SetZHalfLength(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:91:17
    t.method("SetZHalfLength", [](G4Trd& a, G4double arg0)->void { a.SetZHalfLength(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));
    t.method("SetZHalfLength", [](G4Trd* a, G4double arg0)->void { a->SetZHalfLength(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));

    DEBUG_MSG("Adding wrapper for void G4Trd::SetAllParameters(G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Trd::SetAllParameters(G4double, G4double, G4double, G4double, G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:93:10
    t.method("SetAllParameters", [](G4Trd& a, G4double arg0, G4double arg1, G4double arg2, G4double arg3, G4double arg4)->void { a.SetAllParameters(arg0, arg1, arg2, arg3, arg4); }, jlcxx::arg("this"), jlcxx::arg("pdx1"), jlcxx::arg("pdx2"), jlcxx::arg("pdy1"), jlcxx::arg("pdy2"), jlcxx::arg("pdz"));
    t.method("SetAllParameters", [](G4Trd* a, G4double arg0, G4double arg1, G4double arg2, G4double arg3, G4double arg4)->void { a->SetAllParameters(arg0, arg1, arg2, arg3, arg4); }, jlcxx::arg("this"), jlcxx::arg("pdx1"), jlcxx::arg("pdx2"), jlcxx::arg("pdy1"), jlcxx::arg("pdy2"), jlcxx::arg("pdz"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:99:14
    t.method("GetCubicVolume", [](G4Trd& a)->G4double { return a.GetCubicVolume(); }, jlcxx::arg("this"));
    t.method("GetCubicVolume", [](G4Trd* a)->G4double { return a->GetCubicVolume(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:100:14
    t.method("GetSurfaceArea", [](G4Trd& a)->G4double { return a.GetSurfaceArea(); }, jlcxx::arg("this"));
    t.method("GetSurfaceArea", [](G4Trd* a)->G4double { return a->GetSurfaceArea(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Trd::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Trd::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:102:10
    t.method("ComputeDimensions", [](G4Trd& a, G4VPVParameterisation * arg0, const G4int arg1, const G4VPhysicalVolume * arg2)->void { a.ComputeDimensions(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("n"), jlcxx::arg("pRep"));
    t.method("ComputeDimensions", [](G4Trd* a, G4VPVParameterisation * arg0, const G4int arg1, const G4VPhysicalVolume * arg2)->void { a->ComputeDimensions(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("n"), jlcxx::arg("pRep"));

    DEBUG_MSG("Adding wrapper for void G4Trd::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Trd::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:106:10
    t.method("BoundingLimits", [](G4Trd const& a, G4ThreeVector & arg0, G4ThreeVector & arg1)->void { a.BoundingLimits(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("pMin"), jlcxx::arg("pMax"));
    t.method("BoundingLimits", [](G4Trd const* a, G4ThreeVector & arg0, G4ThreeVector & arg1)->void { a->BoundingLimits(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("pMin"), jlcxx::arg("pMax"));

    DEBUG_MSG("Adding wrapper for EInside G4Trd::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4Trd::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:113:13
    t.method("Inside", [](G4Trd const& a, const G4ThreeVector & arg0)->EInside { return a.Inside(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));
    t.method("Inside", [](G4Trd const* a, const G4ThreeVector & arg0)->EInside { return a->Inside(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Trd::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Trd::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:115:19
    t.method("SurfaceNormal", [](G4Trd const& a, const G4ThreeVector & arg0)->G4ThreeVector { return a.SurfaceNormal(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));
    t.method("SurfaceNormal", [](G4Trd const* a, const G4ThreeVector & arg0)->G4ThreeVector { return a->SurfaceNormal(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:117:14
    t.method("DistanceToIn", [](G4Trd const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToIn(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"));
    t.method("DistanceToIn", [](G4Trd const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToIn(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:120:14
    t.method("DistanceToIn", [](G4Trd const& a, const G4ThreeVector & arg0)->G4double { return a.DistanceToIn(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));
    t.method("DistanceToIn", [](G4Trd const* a, const G4ThreeVector & arg0)->G4double { return a->DistanceToIn(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:122:14
    t.method("DistanceToOut", [](G4Trd const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"));
    t.method("DistanceToOut", [](G4Trd const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"));
    t.method("DistanceToOut", [](G4Trd const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"), jlcxx::arg("validNorm"));
    t.method("DistanceToOut", [](G4Trd const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3, G4ThreeVector * arg4)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3, arg4); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"), jlcxx::arg("validNorm"), jlcxx::arg("n"));
    t.method("DistanceToOut", [](G4Trd const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"));
    t.method("DistanceToOut", [](G4Trd const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"));
    t.method("DistanceToOut", [](G4Trd const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"), jlcxx::arg("validNorm"));
    t.method("DistanceToOut", [](G4Trd const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3, G4ThreeVector * arg4)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3, arg4); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"), jlcxx::arg("validNorm"), jlcxx::arg("n"));

    DEBUG_MSG("Adding wrapper for G4double G4Trd::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Trd::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:128:14
    t.method("DistanceToOut", [](G4Trd const& a, const G4ThreeVector & arg0)->G4double { return a.DistanceToOut(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));
    t.method("DistanceToOut", [](G4Trd const* a, const G4ThreeVector & arg0)->G4double { return a->DistanceToOut(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4Trd::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4Trd::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:130:20
    t.method("GetEntityType", [](G4Trd const& a)->G4GeometryType { return a.GetEntityType(); }, jlcxx::arg("this"));
    t.method("GetEntityType", [](G4Trd const* a)->G4GeometryType { return a->GetEntityType(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Trd::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Trd::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:132:19
    t.method("GetPointOnSurface", [](G4Trd const& a)->G4ThreeVector { return a.GetPointOnSurface(); }, jlcxx::arg("this"));
    t.method("GetPointOnSurface", [](G4Trd const* a)->G4ThreeVector { return a->GetPointOnSurface(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4bool G4Trd::IsFaceted() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Trd::IsFaceted()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:134:12
    t.method("IsFaceted", [](G4Trd const& a)->G4bool { return a.IsFaceted(); }, jlcxx::arg("this"));
    t.method("IsFaceted", [](G4Trd const* a)->G4bool { return a->IsFaceted(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4Trd::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4Trd::Clone()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:136:15
    t.method("Clone", [](G4Trd const& a)->G4VSolid * { return a.Clone(); }, jlcxx::arg("this"));
    t.method("Clone", [](G4Trd const* a)->G4VSolid * { return a->Clone(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Trd::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Trd::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:143:19
    t.method("CreatePolyhedron", [](G4Trd const& a)->G4Polyhedron * { return a.CreatePolyhedron(); }, jlcxx::arg("this"));
    t.method("CreatePolyhedron", [](G4Trd const* a)->G4Polyhedron * { return a->CreatePolyhedron(); }, jlcxx::arg("this"));


    DEBUG_MSG("Adding wrapper for void G4Trd::G4Trd(const G4Trd &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:150:5
    t.constructor<const G4Trd &>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("rhs")    );

    DEBUG_MSG("Adding wrapper for G4Trd & G4Trd::operator=(const G4Trd &) (" __HERE__ ")");
    // signature to use in the veto list: G4Trd & G4Trd::operator=(const G4Trd &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Trd.hh:151:12
    t.method("assign", [](G4Trd& a, const G4Trd & arg0)->G4Trd & { return a.operator=(arg0); }, jlcxx::arg("this"), jlcxx::arg("rhs"));
    t.method("assign", [](G4Trd* a, const G4Trd & arg0)->G4Trd & { return a->operator=(arg0); }, jlcxx::arg("this"), jlcxx::arg("rhs"));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Trd>> type_;
};
std::shared_ptr<Wrapper> newJlG4Trd(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Trd(module));
}
