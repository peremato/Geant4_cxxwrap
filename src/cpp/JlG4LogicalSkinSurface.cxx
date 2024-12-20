// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4LogicalSkinSurface> : std::false_type { };
  template<> struct DefaultConstructible<G4LogicalSkinSurface> : std::false_type { };
}

// Class generating the wrapper for type G4LogicalSkinSurface
// signature to use in the veto file: G4LogicalSkinSurface
struct JlG4LogicalSkinSurface: public Wrapper {

  JlG4LogicalSkinSurface(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4LogicalSkinSurface (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:47:7
    jlcxx::TypeWrapper<G4LogicalSkinSurface>  t = jlModule.add_type<G4LogicalSkinSurface>("G4LogicalSkinSurface");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4LogicalSkinSurface>>(new jlcxx::TypeWrapper<G4LogicalSkinSurface>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4LogicalSkinSurface::G4LogicalSkinSurface(const G4String &, G4LogicalVolume *, G4SurfaceProperty *) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:51:5
    t.constructor<const G4String &, G4LogicalVolume *, G4SurfaceProperty *>(/*finalize=*/jlcxx::finalize_policy::no);
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4LogicalSkinSurface::operator==(const G4LogicalSkinSurface &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4LogicalSkinSurface::operator==(const G4LogicalSkinSurface &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:61:12
    t.method("==", static_cast<G4bool (G4LogicalSkinSurface::*)(const G4LogicalSkinSurface &)  const>(&G4LogicalSkinSurface::operator==));

    DEBUG_MSG("Adding wrapper for G4bool G4LogicalSkinSurface::operator!=(const G4LogicalSkinSurface &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4LogicalSkinSurface::operator!=(const G4LogicalSkinSurface &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:62:12
    t.method("!=", static_cast<G4bool (G4LogicalSkinSurface::*)(const G4LogicalSkinSurface &)  const>(&G4LogicalSkinSurface::operator!=));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for G4LogicalSkinSurface * G4LogicalSkinSurface::GetSurface(const G4LogicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: G4LogicalSkinSurface * G4LogicalSkinSurface::GetSurface(const G4LogicalVolume *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:65:34
    module_.method("G4LogicalSkinSurface!GetSurface", static_cast<G4LogicalSkinSurface * (*)(const G4LogicalVolume *) >(&G4LogicalSkinSurface::GetSurface));

    DEBUG_MSG("Adding wrapper for const G4LogicalVolume * G4LogicalSkinSurface::GetLogicalVolume() (" __HERE__ ")");
    // signature to use in the veto list: const G4LogicalVolume * G4LogicalSkinSurface::GetLogicalVolume()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:66:35
    t.method("GetLogicalVolume", static_cast<const G4LogicalVolume * (G4LogicalSkinSurface::*)()  const>(&G4LogicalSkinSurface::GetLogicalVolume));

    DEBUG_MSG("Adding wrapper for void G4LogicalSkinSurface::SetLogicalVolume(G4LogicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4LogicalSkinSurface::SetLogicalVolume(G4LogicalVolume *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:67:18
    t.method("SetLogicalVolume", static_cast<void (G4LogicalSkinSurface::*)(G4LogicalVolume *) >(&G4LogicalSkinSurface::SetLogicalVolume));

    DEBUG_MSG("Adding wrapper for void G4LogicalSkinSurface::CleanSurfaceTable() (" __HERE__ ")");
    // signature to use in the veto list: void G4LogicalSkinSurface::CleanSurfaceTable()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:70:17
    module_.method("G4LogicalSkinSurface!CleanSurfaceTable", static_cast<void (*)() >(&G4LogicalSkinSurface::CleanSurfaceTable));

    DEBUG_MSG("Adding wrapper for const G4LogicalSkinSurfaceTable * G4LogicalSkinSurface::GetSurfaceTable() (" __HERE__ ")");
    // signature to use in the veto list: const G4LogicalSkinSurfaceTable * G4LogicalSkinSurface::GetSurfaceTable()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:71:45
    module_.method("G4LogicalSkinSurface!GetSurfaceTable", static_cast<const G4LogicalSkinSurfaceTable * (*)() >(&G4LogicalSkinSurface::GetSurfaceTable));

    DEBUG_MSG("Adding wrapper for size_t G4LogicalSkinSurface::GetNumberOfSkinSurfaces() (" __HERE__ ")");
    // signature to use in the veto list: size_t G4LogicalSkinSurface::GetNumberOfSkinSurfaces()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:72:24
    module_.method("G4LogicalSkinSurface!GetNumberOfSkinSurfaces", static_cast<size_t (*)() >(&G4LogicalSkinSurface::GetNumberOfSkinSurfaces));

    DEBUG_MSG("Adding wrapper for void G4LogicalSkinSurface::DumpInfo() (" __HERE__ ")");
    // signature to use in the veto list: void G4LogicalSkinSurface::DumpInfo()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4LogicalSkinSurface.hh:73:17
    module_.method("G4LogicalSkinSurface!DumpInfo", static_cast<void (*)() >(&G4LogicalSkinSurface::DumpInfo));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4LogicalSkinSurface>> type_;
};
std::shared_ptr<Wrapper> newJlG4LogicalSkinSurface(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4LogicalSkinSurface(module));
}
