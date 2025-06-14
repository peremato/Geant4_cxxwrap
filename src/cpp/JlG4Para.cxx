// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Para> : std::false_type { };
  template<> struct DefaultConstructible<G4Para> : std::false_type { };
template<> struct SuperType<G4Para> { typedef G4CSGSolid type; };
}

// Class generating the wrapper for type G4Para
// signature to use in the veto file: G4Para
struct JlG4Para: public Wrapper {

  JlG4Para(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Para (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:78:7
    jlcxx::TypeWrapper<G4Para>  t = jlModule.add_type<G4Para>("G4Para",
      jlcxx::julia_base_type<G4CSGSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Para>>(new jlcxx::TypeWrapper<G4Para>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4Para::G4Para(const G4String &, G4double, G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:82:5
    t.constructor<const G4String &, G4double, G4double, G4double, G4double, G4double, G4double>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("pName"), jlcxx::arg("pDx"), jlcxx::arg("pDy"), jlcxx::arg("pDz"), jlcxx::arg("pAlpha"), jlcxx::arg("pTheta"), jlcxx::arg("pPhi")    );


    DEBUG_MSG("Adding wrapper for void G4Para::G4Para(const G4String &, const G4ThreeVector[8]) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:86:5
    t.constructor<const G4String &, const G4ThreeVector[8]>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("pName"), jlcxx::arg("pt")    );

    DEBUG_MSG("Adding wrapper for G4double G4Para::GetZHalfLength() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::GetZHalfLength()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:93:21
    t.method("GetZHalfLength", [](G4Para const& a)->G4double { return a.GetZHalfLength(); }, jlcxx::arg("this"));
    t.method("GetZHalfLength", [](G4Para const* a)->G4double { return a->GetZHalfLength(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Para::GetSymAxis() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Para::GetSymAxis()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:94:26
    t.method("GetSymAxis", [](G4Para const& a)->G4ThreeVector { return a.GetSymAxis(); }, jlcxx::arg("this"));
    t.method("GetSymAxis", [](G4Para const* a)->G4ThreeVector { return a->GetSymAxis(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::GetYHalfLength() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::GetYHalfLength()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:95:21
    t.method("GetYHalfLength", [](G4Para const& a)->G4double { return a.GetYHalfLength(); }, jlcxx::arg("this"));
    t.method("GetYHalfLength", [](G4Para const* a)->G4double { return a->GetYHalfLength(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::GetXHalfLength() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::GetXHalfLength()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:96:21
    t.method("GetXHalfLength", [](G4Para const& a)->G4double { return a.GetXHalfLength(); }, jlcxx::arg("this"));
    t.method("GetXHalfLength", [](G4Para const* a)->G4double { return a->GetXHalfLength(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::GetTanAlpha() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::GetTanAlpha()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:97:21
    t.method("GetTanAlpha", [](G4Para const& a)->G4double { return a.GetTanAlpha(); }, jlcxx::arg("this"));
    t.method("GetTanAlpha", [](G4Para const* a)->G4double { return a->GetTanAlpha(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::GetAlpha() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::GetAlpha()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:99:21
    t.method("GetAlpha", [](G4Para const& a)->G4double { return a.GetAlpha(); }, jlcxx::arg("this"));
    t.method("GetAlpha", [](G4Para const* a)->G4double { return a->GetAlpha(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::GetTheta() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::GetTheta()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:100:21
    t.method("GetTheta", [](G4Para const& a)->G4double { return a.GetTheta(); }, jlcxx::arg("this"));
    t.method("GetTheta", [](G4Para const* a)->G4double { return a->GetTheta(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::GetPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::GetPhi()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:101:21
    t.method("GetPhi", [](G4Para const& a)->G4double { return a.GetPhi(); }, jlcxx::arg("this"));
    t.method("GetPhi", [](G4Para const* a)->G4double { return a->GetPhi(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Para::SetXHalfLength(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Para::SetXHalfLength(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:106:17
    t.method("SetXHalfLength", [](G4Para& a, G4double arg0)->void { a.SetXHalfLength(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));
    t.method("SetXHalfLength", [](G4Para* a, G4double arg0)->void { a->SetXHalfLength(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));

    DEBUG_MSG("Adding wrapper for void G4Para::SetYHalfLength(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Para::SetYHalfLength(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:107:17
    t.method("SetYHalfLength", [](G4Para& a, G4double arg0)->void { a.SetYHalfLength(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));
    t.method("SetYHalfLength", [](G4Para* a, G4double arg0)->void { a->SetYHalfLength(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));

    DEBUG_MSG("Adding wrapper for void G4Para::SetZHalfLength(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Para::SetZHalfLength(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:108:17
    t.method("SetZHalfLength", [](G4Para& a, G4double arg0)->void { a.SetZHalfLength(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));
    t.method("SetZHalfLength", [](G4Para* a, G4double arg0)->void { a->SetZHalfLength(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));

    DEBUG_MSG("Adding wrapper for void G4Para::SetAlpha(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Para::SetAlpha(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:109:17
    t.method("SetAlpha", [](G4Para& a, G4double arg0)->void { a.SetAlpha(arg0); }, jlcxx::arg("this"), jlcxx::arg("alpha"));
    t.method("SetAlpha", [](G4Para* a, G4double arg0)->void { a->SetAlpha(arg0); }, jlcxx::arg("this"), jlcxx::arg("alpha"));

    DEBUG_MSG("Adding wrapper for void G4Para::SetTanAlpha(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Para::SetTanAlpha(G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:110:17
    t.method("SetTanAlpha", [](G4Para& a, G4double arg0)->void { a.SetTanAlpha(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));
    t.method("SetTanAlpha", [](G4Para* a, G4double arg0)->void { a->SetTanAlpha(arg0); }, jlcxx::arg("this"), jlcxx::arg("val"));

    DEBUG_MSG("Adding wrapper for void G4Para::SetThetaAndPhi(G4double, G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Para::SetThetaAndPhi(G4double, G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:111:17
    t.method("SetThetaAndPhi", [](G4Para& a, G4double arg0, G4double arg1)->void { a.SetThetaAndPhi(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("pTheta"), jlcxx::arg("pPhi"));
    t.method("SetThetaAndPhi", [](G4Para* a, G4double arg0, G4double arg1)->void { a->SetThetaAndPhi(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("pTheta"), jlcxx::arg("pPhi"));

    DEBUG_MSG("Adding wrapper for void G4Para::SetAllParameters(G4double, G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Para::SetAllParameters(G4double, G4double, G4double, G4double, G4double, G4double)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:113:10
    t.method("SetAllParameters", [](G4Para& a, G4double arg0, G4double arg1, G4double arg2, G4double arg3, G4double arg4, G4double arg5)->void { a.SetAllParameters(arg0, arg1, arg2, arg3, arg4, arg5); }, jlcxx::arg("this"), jlcxx::arg("pDx"), jlcxx::arg("pDy"), jlcxx::arg("pDz"), jlcxx::arg("pAlpha"), jlcxx::arg("pTheta"), jlcxx::arg("pPhi"));
    t.method("SetAllParameters", [](G4Para* a, G4double arg0, G4double arg1, G4double arg2, G4double arg3, G4double arg4, G4double arg5)->void { a->SetAllParameters(arg0, arg1, arg2, arg3, arg4, arg5); }, jlcxx::arg("this"), jlcxx::arg("pDx"), jlcxx::arg("pDy"), jlcxx::arg("pDz"), jlcxx::arg("pAlpha"), jlcxx::arg("pTheta"), jlcxx::arg("pPhi"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:118:14
    t.method("GetCubicVolume", [](G4Para& a)->G4double { return a.GetCubicVolume(); }, jlcxx::arg("this"));
    t.method("GetCubicVolume", [](G4Para* a)->G4double { return a->GetCubicVolume(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:119:14
    t.method("GetSurfaceArea", [](G4Para& a)->G4double { return a.GetSurfaceArea(); }, jlcxx::arg("this"));
    t.method("GetSurfaceArea", [](G4Para* a)->G4double { return a->GetSurfaceArea(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4Para::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Para::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:121:10
    t.method("ComputeDimensions", [](G4Para& a, G4VPVParameterisation * arg0, const G4int arg1, const G4VPhysicalVolume * arg2)->void { a.ComputeDimensions(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("n"), jlcxx::arg("pRep"));
    t.method("ComputeDimensions", [](G4Para* a, G4VPVParameterisation * arg0, const G4int arg1, const G4VPhysicalVolume * arg2)->void { a->ComputeDimensions(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("n"), jlcxx::arg("pRep"));

    DEBUG_MSG("Adding wrapper for void G4Para::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Para::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:125:10
    t.method("BoundingLimits", [](G4Para const& a, G4ThreeVector & arg0, G4ThreeVector & arg1)->void { a.BoundingLimits(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("pMin"), jlcxx::arg("pMax"));
    t.method("BoundingLimits", [](G4Para const* a, G4ThreeVector & arg0, G4ThreeVector & arg1)->void { a->BoundingLimits(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("pMin"), jlcxx::arg("pMax"));

    DEBUG_MSG("Adding wrapper for EInside G4Para::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4Para::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:132:13
    t.method("Inside", [](G4Para const& a, const G4ThreeVector & arg0)->EInside { return a.Inside(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));
    t.method("Inside", [](G4Para const* a, const G4ThreeVector & arg0)->EInside { return a->Inside(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Para::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Para::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:134:19
    t.method("SurfaceNormal", [](G4Para const& a, const G4ThreeVector & arg0)->G4ThreeVector { return a.SurfaceNormal(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));
    t.method("SurfaceNormal", [](G4Para const* a, const G4ThreeVector & arg0)->G4ThreeVector { return a->SurfaceNormal(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:136:14
    t.method("DistanceToIn", [](G4Para const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToIn(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"));
    t.method("DistanceToIn", [](G4Para const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToIn(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:138:14
    t.method("DistanceToIn", [](G4Para const& a, const G4ThreeVector & arg0)->G4double { return a.DistanceToIn(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));
    t.method("DistanceToIn", [](G4Para const* a, const G4ThreeVector & arg0)->G4double { return a->DistanceToIn(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:140:14
    t.method("DistanceToOut", [](G4Para const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"));
    t.method("DistanceToOut", [](G4Para const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"));
    t.method("DistanceToOut", [](G4Para const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"), jlcxx::arg("validNorm"));
    t.method("DistanceToOut", [](G4Para const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3, G4ThreeVector * arg4)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3, arg4); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"), jlcxx::arg("validNorm"), jlcxx::arg("n"));
    t.method("DistanceToOut", [](G4Para const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"));
    t.method("DistanceToOut", [](G4Para const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"));
    t.method("DistanceToOut", [](G4Para const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"), jlcxx::arg("validNorm"));
    t.method("DistanceToOut", [](G4Para const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3, G4ThreeVector * arg4)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3, arg4); }, jlcxx::arg("this"), jlcxx::arg("p"), jlcxx::arg("v"), jlcxx::arg("calcNorm"), jlcxx::arg("validNorm"), jlcxx::arg("n"));

    DEBUG_MSG("Adding wrapper for G4double G4Para::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Para::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:144:14
    t.method("DistanceToOut", [](G4Para const& a, const G4ThreeVector & arg0)->G4double { return a.DistanceToOut(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));
    t.method("DistanceToOut", [](G4Para const* a, const G4ThreeVector & arg0)->G4double { return a->DistanceToOut(arg0); }, jlcxx::arg("this"), jlcxx::arg("p"));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4Para::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4Para::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:146:20
    t.method("GetEntityType", [](G4Para const& a)->G4GeometryType { return a.GetEntityType(); }, jlcxx::arg("this"));
    t.method("GetEntityType", [](G4Para const* a)->G4GeometryType { return a->GetEntityType(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Para::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Para::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:148:19
    t.method("GetPointOnSurface", [](G4Para const& a)->G4ThreeVector { return a.GetPointOnSurface(); }, jlcxx::arg("this"));
    t.method("GetPointOnSurface", [](G4Para const* a)->G4ThreeVector { return a->GetPointOnSurface(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4bool G4Para::IsFaceted() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4Para::IsFaceted()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:150:12
    t.method("IsFaceted", [](G4Para const& a)->G4bool { return a.IsFaceted(); }, jlcxx::arg("this"));
    t.method("IsFaceted", [](G4Para const* a)->G4bool { return a->IsFaceted(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4Para::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4Para::Clone()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:152:15
    t.method("Clone", [](G4Para const& a)->G4VSolid * { return a.Clone(); }, jlcxx::arg("this"));
    t.method("Clone", [](G4Para const* a)->G4VSolid * { return a->Clone(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Para::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Para::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:159:19
    t.method("CreatePolyhedron", [](G4Para const& a)->G4Polyhedron * { return a.CreatePolyhedron(); }, jlcxx::arg("this"));
    t.method("CreatePolyhedron", [](G4Para const* a)->G4Polyhedron * { return a->CreatePolyhedron(); }, jlcxx::arg("this"));


    DEBUG_MSG("Adding wrapper for void G4Para::G4Para(const G4Para &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:166:5
    t.constructor<const G4Para &>(/*finalize=*/jlcxx::finalize_policy::yes, jlcxx::arg("this"), jlcxx::arg("rhs")    );

    DEBUG_MSG("Adding wrapper for G4Para & G4Para::operator=(const G4Para &) (" __HERE__ ")");
    // signature to use in the veto list: G4Para & G4Para::operator=(const G4Para &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4Para.hh:167:13
    t.method("assign", [](G4Para& a, const G4Para & arg0)->G4Para & { return a.operator=(arg0); }, jlcxx::arg("this"), jlcxx::arg("rhs"));
    t.method("assign", [](G4Para* a, const G4Para & arg0)->G4Para & { return a->operator=(arg0); }, jlcxx::arg("this"), jlcxx::arg("rhs"));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Para>> type_;
};
std::shared_ptr<Wrapper> newJlG4Para(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Para(module));
}
