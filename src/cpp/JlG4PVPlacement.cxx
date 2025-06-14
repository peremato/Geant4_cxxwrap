// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4PVPlacement> : std::false_type { };
  template<> struct DefaultConstructible<G4PVPlacement> : std::false_type { };
template<> struct SuperType<G4PVPlacement> { typedef G4VPhysicalVolume type; };
}

// Class generating the wrapper for type G4PVPlacement
// signature to use in the veto file: G4PVPlacement
struct JlG4PVPlacement: public Wrapper {

  JlG4PVPlacement(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4PVPlacement (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:41:7
    jlcxx::TypeWrapper<G4PVPlacement>  t = jlModule.add_type<G4PVPlacement>("G4PVPlacement",
      jlcxx::julia_base_type<G4VPhysicalVolume>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4PVPlacement>>(new jlcxx::TypeWrapper<G4PVPlacement>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4PVPlacement::G4PVPlacement(G4RotationMatrix *, const G4ThreeVector &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int, G4bool) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:45:5
    t.constructor<G4RotationMatrix *, const G4ThreeVector &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int>(/*finalize=*/jlcxx::finalize_policy::no, jlcxx::arg("this"), jlcxx::arg("pRot"), jlcxx::arg("tlate"), jlcxx::arg("pCurrentLogical"), jlcxx::arg("pName"), jlcxx::arg("pMotherLogical"), jlcxx::arg("pMany"), jlcxx::arg("pCopyNo")    );
    t.constructor<G4RotationMatrix *, const G4ThreeVector &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int, G4bool>(/*finalize=*/jlcxx::finalize_policy::no, jlcxx::arg("this"), jlcxx::arg("pRot"), jlcxx::arg("tlate"), jlcxx::arg("pCurrentLogical"), jlcxx::arg("pName"), jlcxx::arg("pMotherLogical"), jlcxx::arg("pMany"), jlcxx::arg("pCopyNo"), jlcxx::arg("pSurfChk")    );


    DEBUG_MSG("Adding wrapper for void G4PVPlacement::G4PVPlacement(const G4Transform3D &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int, G4bool) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:67:5
    t.constructor<const G4Transform3D &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int>(/*finalize=*/jlcxx::finalize_policy::no, jlcxx::arg("this"), jlcxx::arg("Transform3D"), jlcxx::arg("pCurrentLogical"), jlcxx::arg("pName"), jlcxx::arg("pMotherLogical"), jlcxx::arg("pMany"), jlcxx::arg("pCopyNo")    );
    t.constructor<const G4Transform3D &, G4LogicalVolume *, const G4String &, G4LogicalVolume *, G4bool, G4int, G4bool>(/*finalize=*/jlcxx::finalize_policy::no, jlcxx::arg("this"), jlcxx::arg("Transform3D"), jlcxx::arg("pCurrentLogical"), jlcxx::arg("pName"), jlcxx::arg("pMotherLogical"), jlcxx::arg("pMany"), jlcxx::arg("pCopyNo"), jlcxx::arg("pSurfChk")    );


    DEBUG_MSG("Adding wrapper for void G4PVPlacement::G4PVPlacement(G4RotationMatrix *, const G4ThreeVector &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int, G4bool) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:86:5
    t.constructor<G4RotationMatrix *, const G4ThreeVector &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int>(/*finalize=*/jlcxx::finalize_policy::no, jlcxx::arg("this"), jlcxx::arg("pRot"), jlcxx::arg("tlate"), jlcxx::arg("pName"), jlcxx::arg("pLogical"), jlcxx::arg("pMother"), jlcxx::arg("pMany"), jlcxx::arg("pCopyNo")    );
    t.constructor<G4RotationMatrix *, const G4ThreeVector &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int, G4bool>(/*finalize=*/jlcxx::finalize_policy::no, jlcxx::arg("this"), jlcxx::arg("pRot"), jlcxx::arg("tlate"), jlcxx::arg("pName"), jlcxx::arg("pLogical"), jlcxx::arg("pMother"), jlcxx::arg("pMany"), jlcxx::arg("pCopyNo"), jlcxx::arg("pSurfChk")    );


    DEBUG_MSG("Adding wrapper for void G4PVPlacement::G4PVPlacement(const G4Transform3D &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int, G4bool) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:98:5
    t.constructor<const G4Transform3D &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int>(/*finalize=*/jlcxx::finalize_policy::no, jlcxx::arg("this"), jlcxx::arg("Transform3D"), jlcxx::arg("pName"), jlcxx::arg("pLogical"), jlcxx::arg("pMother"), jlcxx::arg("pMany"), jlcxx::arg("pCopyNo")    );
    t.constructor<const G4Transform3D &, const G4String &, G4LogicalVolume *, G4VPhysicalVolume *, G4bool, G4int, G4bool>(/*finalize=*/jlcxx::finalize_policy::no, jlcxx::arg("this"), jlcxx::arg("Transform3D"), jlcxx::arg("pName"), jlcxx::arg("pLogical"), jlcxx::arg("pMother"), jlcxx::arg("pMany"), jlcxx::arg("pCopyNo"), jlcxx::arg("pSurfChk")    );

    DEBUG_MSG("Adding wrapper for G4int G4PVPlacement::GetCopyNo() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4PVPlacement::GetCopyNo()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:111:18
    t.method("GetCopyNo", [](G4PVPlacement const& a)->G4int { return a.GetCopyNo(); }, jlcxx::arg("this"));
    t.method("GetCopyNo", [](G4PVPlacement const* a)->G4int { return a->GetCopyNo(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4PVPlacement::SetCopyNo(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4PVPlacement::SetCopyNo(G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:113:10
    t.method("SetCopyNo", [](G4PVPlacement& a, G4int arg0)->void { a.SetCopyNo(arg0); }, jlcxx::arg("this"), jlcxx::arg("CopyNo"));
    t.method("SetCopyNo", [](G4PVPlacement* a, G4int arg0)->void { a->SetCopyNo(arg0); }, jlcxx::arg("this"), jlcxx::arg("CopyNo"));

    DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::CheckOverlaps(G4int, G4double, G4bool, G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4PVPlacement::CheckOverlaps(G4int, G4double, G4bool, G4int)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:116:12
    t.method("CheckOverlaps", [](G4PVPlacement& a)->G4bool { return a.CheckOverlaps(); }, jlcxx::arg("this"));
    t.method("CheckOverlaps", [](G4PVPlacement& a, G4int arg0)->G4bool { return a.CheckOverlaps(arg0); }, jlcxx::arg("this"), jlcxx::arg("res"));
    t.method("CheckOverlaps", [](G4PVPlacement& a, G4int arg0, G4double arg1)->G4bool { return a.CheckOverlaps(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("res"), jlcxx::arg("tol"));
    t.method("CheckOverlaps", [](G4PVPlacement& a, G4int arg0, G4double arg1, G4bool arg2)->G4bool { return a.CheckOverlaps(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("res"), jlcxx::arg("tol"), jlcxx::arg("verbose"));
    t.method("CheckOverlaps", [](G4PVPlacement& a, G4int arg0, G4double arg1, G4bool arg2, G4int arg3)->G4bool { return a.CheckOverlaps(arg0, arg1, arg2, arg3); }, jlcxx::arg("this"), jlcxx::arg("res"), jlcxx::arg("tol"), jlcxx::arg("verbose"), jlcxx::arg("maxErr"));
    t.method("CheckOverlaps", [](G4PVPlacement* a)->G4bool { return a->CheckOverlaps(); }, jlcxx::arg("this"));
    t.method("CheckOverlaps", [](G4PVPlacement* a, G4int arg0)->G4bool { return a->CheckOverlaps(arg0); }, jlcxx::arg("this"), jlcxx::arg("res"));
    t.method("CheckOverlaps", [](G4PVPlacement* a, G4int arg0, G4double arg1)->G4bool { return a->CheckOverlaps(arg0, arg1); }, jlcxx::arg("this"), jlcxx::arg("res"), jlcxx::arg("tol"));
    t.method("CheckOverlaps", [](G4PVPlacement* a, G4int arg0, G4double arg1, G4bool arg2)->G4bool { return a->CheckOverlaps(arg0, arg1, arg2); }, jlcxx::arg("this"), jlcxx::arg("res"), jlcxx::arg("tol"), jlcxx::arg("verbose"));
    t.method("CheckOverlaps", [](G4PVPlacement* a, G4int arg0, G4double arg1, G4bool arg2, G4int arg3)->G4bool { return a->CheckOverlaps(arg0, arg1, arg2, arg3); }, jlcxx::arg("this"), jlcxx::arg("res"), jlcxx::arg("tol"), jlcxx::arg("verbose"), jlcxx::arg("maxErr"));

    DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::IsMany() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4PVPlacement::IsMany()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:135:12
    t.method("IsMany", [](G4PVPlacement const& a)->G4bool { return a.IsMany(); }, jlcxx::arg("this"));
    t.method("IsMany", [](G4PVPlacement const* a)->G4bool { return a->IsMany(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::IsReplicated() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4PVPlacement::IsReplicated()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:136:12
    t.method("IsReplicated", [](G4PVPlacement const& a)->G4bool { return a.IsReplicated(); }, jlcxx::arg("this"));
    t.method("IsReplicated", [](G4PVPlacement const* a)->G4bool { return a->IsReplicated(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::IsParameterised() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4PVPlacement::IsParameterised()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:137:12
    t.method("IsParameterised", [](G4PVPlacement const& a)->G4bool { return a.IsParameterised(); }, jlcxx::arg("this"));
    t.method("IsParameterised", [](G4PVPlacement const* a)->G4bool { return a->IsParameterised(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4VPVParameterisation * G4PVPlacement::GetParameterisation() (" __HERE__ ")");
    // signature to use in the veto list: G4VPVParameterisation * G4PVPlacement::GetParameterisation()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:138:28
    t.method("GetParameterisation", [](G4PVPlacement const& a)->G4VPVParameterisation * { return a.GetParameterisation(); }, jlcxx::arg("this"));
    t.method("GetParameterisation", [](G4PVPlacement const* a)->G4VPVParameterisation * { return a->GetParameterisation(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for void G4PVPlacement::GetReplicationData(EAxis &, G4int &, G4double &, G4double &, G4bool &) (" __HERE__ ")");
    // signature to use in the veto list: void G4PVPlacement::GetReplicationData(EAxis &, G4int &, G4double &, G4double &, G4bool &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:139:10
    t.method("GetReplicationData", [](G4PVPlacement const& a, EAxis & arg0, G4int & arg1, G4double & arg2, G4double & arg3, G4bool & arg4)->void { a.GetReplicationData(arg0, arg1, arg2, arg3, arg4); }, jlcxx::arg("this"), jlcxx::arg("axis"), jlcxx::arg("nReplicas"), jlcxx::arg("width"), jlcxx::arg("offset"), jlcxx::arg("consuming"));
    t.method("GetReplicationData", [](G4PVPlacement const* a, EAxis & arg0, G4int & arg1, G4double & arg2, G4double & arg3, G4bool & arg4)->void { a->GetReplicationData(arg0, arg1, arg2, arg3, arg4); }, jlcxx::arg("this"), jlcxx::arg("axis"), jlcxx::arg("nReplicas"), jlcxx::arg("width"), jlcxx::arg("offset"), jlcxx::arg("consuming"));

    DEBUG_MSG("Adding wrapper for G4bool G4PVPlacement::IsRegularStructure() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4PVPlacement::IsRegularStructure()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:144:12
    t.method("IsRegularStructure", [](G4PVPlacement const& a)->G4bool { return a.IsRegularStructure(); }, jlcxx::arg("this"));
    t.method("IsRegularStructure", [](G4PVPlacement const* a)->G4bool { return a->IsRegularStructure(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for G4int G4PVPlacement::GetRegularStructureId() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4PVPlacement::GetRegularStructureId()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:145:11
    t.method("GetRegularStructureId", [](G4PVPlacement const& a)->G4int { return a.GetRegularStructureId(); }, jlcxx::arg("this"));
    t.method("GetRegularStructureId", [](G4PVPlacement const* a)->G4int { return a->GetRegularStructureId(); }, jlcxx::arg("this"));

    DEBUG_MSG("Adding wrapper for EVolume G4PVPlacement::VolumeType() (" __HERE__ ")");
    // signature to use in the veto list: EVolume G4PVPlacement::VolumeType()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4PVPlacement.hh:147:13
    t.method("VolumeType", [](G4PVPlacement const& a)->EVolume { return a.VolumeType(); }, jlcxx::arg("this"));
    t.method("VolumeType", [](G4PVPlacement const* a)->EVolume { return a->VolumeType(); }, jlcxx::arg("this"));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4PVPlacement>> type_;
};
std::shared_ptr<Wrapper> newJlG4PVPlacement(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4PVPlacement(module));
}
