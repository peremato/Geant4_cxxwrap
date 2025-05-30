// this file was auto-generated by wrapit v1.3.1-15-g5168a24
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4LogicalBorderSurface> : std::false_type { };
  template<> struct DefaultConstructible<G4LogicalBorderSurface> : std::false_type { };
}

// Class generating the wrapper for type G4LogicalBorderSurface
// signature to use in the veto file: G4LogicalBorderSurface
struct JlG4LogicalBorderSurface: public Wrapper {

  JlG4LogicalBorderSurface(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4LogicalBorderSurface (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:50:7
    jlcxx::TypeWrapper<G4LogicalBorderSurface>  t = jlModule.add_type<G4LogicalBorderSurface>("G4LogicalBorderSurface");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4LogicalBorderSurface>>(new jlcxx::TypeWrapper<G4LogicalBorderSurface>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4LogicalBorderSurface::G4LogicalBorderSurface(const G4String &, G4VPhysicalVolume *, G4VPhysicalVolume *, G4SurfaceProperty *) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:54:5
    t.constructor<const G4String &, G4VPhysicalVolume *, G4VPhysicalVolume *, G4SurfaceProperty *>(/*finalize=*/jlcxx::finalize_policy::no);
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4LogicalBorderSurface::operator==(const G4LogicalBorderSurface &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4LogicalBorderSurface::operator==(const G4LogicalBorderSurface &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:65:12
    t.method("==", static_cast<G4bool (G4LogicalBorderSurface::*)(const G4LogicalBorderSurface &)  const>(&G4LogicalBorderSurface::operator==));

    DEBUG_MSG("Adding wrapper for G4bool G4LogicalBorderSurface::operator!=(const G4LogicalBorderSurface &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4LogicalBorderSurface::operator!=(const G4LogicalBorderSurface &)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:66:12
    t.method("!=", static_cast<G4bool (G4LogicalBorderSurface::*)(const G4LogicalBorderSurface &)  const>(&G4LogicalBorderSurface::operator!=));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for G4LogicalBorderSurface * G4LogicalBorderSurface::GetSurface(const G4VPhysicalVolume *, const G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: G4LogicalBorderSurface * G4LogicalBorderSurface::GetSurface(const G4VPhysicalVolume *, const G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:69:36
    module_.method("G4LogicalBorderSurface!GetSurface", static_cast<G4LogicalBorderSurface * (*)(const G4VPhysicalVolume *, const G4VPhysicalVolume *) >(&G4LogicalBorderSurface::GetSurface));

    DEBUG_MSG("Adding wrapper for void G4LogicalBorderSurface::SetPhysicalVolumes(G4VPhysicalVolume *, G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4LogicalBorderSurface::SetPhysicalVolumes(G4VPhysicalVolume *, G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:71:17
    t.method("SetPhysicalVolumes", static_cast<void (G4LogicalBorderSurface::*)(G4VPhysicalVolume *, G4VPhysicalVolume *) >(&G4LogicalBorderSurface::SetPhysicalVolumes));

    DEBUG_MSG("Adding wrapper for const G4VPhysicalVolume * G4LogicalBorderSurface::GetVolume1() (" __HERE__ ")");
    // signature to use in the veto list: const G4VPhysicalVolume * G4LogicalBorderSurface::GetVolume1()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:73:37
    t.method("GetVolume1", static_cast<const G4VPhysicalVolume * (G4LogicalBorderSurface::*)()  const>(&G4LogicalBorderSurface::GetVolume1));

    DEBUG_MSG("Adding wrapper for const G4VPhysicalVolume * G4LogicalBorderSurface::GetVolume2() (" __HERE__ ")");
    // signature to use in the veto list: const G4VPhysicalVolume * G4LogicalBorderSurface::GetVolume2()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:74:37
    t.method("GetVolume2", static_cast<const G4VPhysicalVolume * (G4LogicalBorderSurface::*)()  const>(&G4LogicalBorderSurface::GetVolume2));

    DEBUG_MSG("Adding wrapper for void G4LogicalBorderSurface::SetVolume1(G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4LogicalBorderSurface::SetVolume1(G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:78:17
    t.method("SetVolume1", static_cast<void (G4LogicalBorderSurface::*)(G4VPhysicalVolume *) >(&G4LogicalBorderSurface::SetVolume1));

    DEBUG_MSG("Adding wrapper for void G4LogicalBorderSurface::SetVolume2(G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4LogicalBorderSurface::SetVolume2(G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:79:17
    t.method("SetVolume2", static_cast<void (G4LogicalBorderSurface::*)(G4VPhysicalVolume *) >(&G4LogicalBorderSurface::SetVolume2));

    DEBUG_MSG("Adding wrapper for void G4LogicalBorderSurface::CleanSurfaceTable() (" __HERE__ ")");
    // signature to use in the veto list: void G4LogicalBorderSurface::CleanSurfaceTable()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:82:17
    module_.method("G4LogicalBorderSurface!CleanSurfaceTable", static_cast<void (*)() >(&G4LogicalBorderSurface::CleanSurfaceTable));

    DEBUG_MSG("Adding wrapper for void G4LogicalBorderSurface::DumpInfo() (" __HERE__ ")");
    // signature to use in the veto list: void G4LogicalBorderSurface::DumpInfo()
    // defined in /Users/mato/.julia/artifacts/c08a070cdc1b892bb33db4924fdac1694e77d3a1/include/Geant4/G4LogicalBorderSurface.hh:85:17
    module_.method("G4LogicalBorderSurface!DumpInfo", static_cast<void (*)() >(&G4LogicalBorderSurface::DumpInfo));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4LogicalBorderSurface>> type_;
};
std::shared_ptr<Wrapper> newJlG4LogicalBorderSurface(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4LogicalBorderSurface(module));
}
