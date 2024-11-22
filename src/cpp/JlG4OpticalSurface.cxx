// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4OpticalSurface> : std::false_type { };
  template<> struct DefaultConstructible<G4OpticalSurface> : std::false_type { };
template<> struct SuperType<G4OpticalSurface> { typedef G4SurfaceProperty type; };
}

// Class generating the wrapper for type G4OpticalSurface
// signature to use in the veto file: G4OpticalSurface
struct JlG4OpticalSurface: public Wrapper {

  JlG4OpticalSurface(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4OpticalSurface (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:119:7
    jlcxx::TypeWrapper<G4OpticalSurface>  t = jlModule.add_type<G4OpticalSurface>("G4OpticalSurface",
      jlcxx::julia_base_type<G4SurfaceProperty>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4OpticalSurface>>(new jlcxx::TypeWrapper<G4OpticalSurface>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::G4OpticalSurface(const G4String &, G4OpticalSurfaceModel, G4OpticalSurfaceFinish, G4SurfaceType, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:123:3
    t.constructor<const G4String &>(/*finalize=*/jlcxx::finalize_policy::no);
    t.constructor<const G4String &, G4OpticalSurfaceModel>(/*finalize=*/jlcxx::finalize_policy::no);
    t.constructor<const G4String &, G4OpticalSurfaceModel, G4OpticalSurfaceFinish>(/*finalize=*/jlcxx::finalize_policy::no);
    t.constructor<const G4String &, G4OpticalSurfaceModel, G4OpticalSurfaceFinish, G4SurfaceType>(/*finalize=*/jlcxx::finalize_policy::no);
    t.constructor<const G4String &, G4OpticalSurfaceModel, G4OpticalSurfaceFinish, G4SurfaceType, G4double>(/*finalize=*/jlcxx::finalize_policy::no);


    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::G4OpticalSurface(const G4OpticalSurface &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:129:3
    t.constructor<const G4OpticalSurface &>(/*finalize=*/jlcxx::finalize_policy::no);

    DEBUG_MSG("Adding wrapper for G4OpticalSurface & G4OpticalSurface::operator=(const G4OpticalSurface &) (" __HERE__ ")");
    // signature to use in the veto list: G4OpticalSurface & G4OpticalSurface::operator=(const G4OpticalSurface &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:130:21
    t.method("assign", static_cast<G4OpticalSurface & (G4OpticalSurface::*)(const G4OpticalSurface &) >(&G4OpticalSurface::operator=));
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4OpticalSurface::operator==(const G4OpticalSurface &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4OpticalSurface::operator==(const G4OpticalSurface &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:132:10
    t.method("==", static_cast<G4bool (G4OpticalSurface::*)(const G4OpticalSurface &)  const>(&G4OpticalSurface::operator==));

    DEBUG_MSG("Adding wrapper for G4bool G4OpticalSurface::operator!=(const G4OpticalSurface &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4OpticalSurface::operator!=(const G4OpticalSurface &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:133:10
    t.method("!=", static_cast<G4bool (G4OpticalSurface::*)(const G4OpticalSurface &)  const>(&G4OpticalSurface::operator!=));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::SetType(const G4SurfaceType &) (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::SetType(const G4SurfaceType &)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:135:8
    t.method("SetType", static_cast<void (G4OpticalSurface::*)(const G4SurfaceType &) >(&G4OpticalSurface::SetType));

    DEBUG_MSG("Adding wrapper for G4OpticalSurfaceFinish G4OpticalSurface::GetFinish() (" __HERE__ ")");
    // signature to use in the veto list: G4OpticalSurfaceFinish G4OpticalSurface::GetFinish()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:138:33
    t.method("GetFinish", static_cast<G4OpticalSurfaceFinish (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetFinish));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::SetFinish(const G4OpticalSurfaceFinish) (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::SetFinish(const G4OpticalSurfaceFinish)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:141:8
    t.method("SetFinish", static_cast<void (G4OpticalSurface::*)(const G4OpticalSurfaceFinish) >(&G4OpticalSurface::SetFinish));

    DEBUG_MSG("Adding wrapper for G4OpticalSurfaceModel G4OpticalSurface::GetModel() (" __HERE__ ")");
    // signature to use in the veto list: G4OpticalSurfaceModel G4OpticalSurface::GetModel()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:144:32
    t.method("GetModel", static_cast<G4OpticalSurfaceModel (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetModel));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::SetModel(const G4OpticalSurfaceModel) (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::SetModel(const G4OpticalSurfaceModel)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:147:15
    t.method("SetModel", static_cast<void (G4OpticalSurface::*)(const G4OpticalSurfaceModel) >(&G4OpticalSurface::SetModel));

    DEBUG_MSG("Adding wrapper for G4double G4OpticalSurface::GetSigmaAlpha() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4OpticalSurface::GetSigmaAlpha()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:150:19
    t.method("GetSigmaAlpha", static_cast<G4double (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetSigmaAlpha));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::SetSigmaAlpha(const G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::SetSigmaAlpha(const G4double)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:153:15
    t.method("SetSigmaAlpha", static_cast<void (G4OpticalSurface::*)(const G4double) >(&G4OpticalSurface::SetSigmaAlpha));

    DEBUG_MSG("Adding wrapper for G4double G4OpticalSurface::GetPolish() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4OpticalSurface::GetPolish()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:156:12
    t.method("GetPolish", static_cast<G4double (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetPolish));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::SetPolish(const G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::SetPolish(const G4double)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:159:15
    t.method("SetPolish", static_cast<void (G4OpticalSurface::*)(const G4double) >(&G4OpticalSurface::SetPolish));

    DEBUG_MSG("Adding wrapper for G4MaterialPropertiesTable * G4OpticalSurface::GetMaterialPropertiesTable() (" __HERE__ ")");
    // signature to use in the veto list: G4MaterialPropertiesTable * G4OpticalSurface::GetMaterialPropertiesTable()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:163:37
    t.method("GetMaterialPropertiesTable", static_cast<G4MaterialPropertiesTable * (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetMaterialPropertiesTable));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::SetMaterialPropertiesTable(G4MaterialPropertiesTable *) (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::SetMaterialPropertiesTable(G4MaterialPropertiesTable *)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:169:15
    t.method("SetMaterialPropertiesTable", static_cast<void (G4OpticalSurface::*)(G4MaterialPropertiesTable *) >(&G4OpticalSurface::SetMaterialPropertiesTable));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::DumpInfo() (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::DumpInfo()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:175:8
    t.method("DumpInfo", static_cast<void (G4OpticalSurface::*)()  const>(&G4OpticalSurface::DumpInfo));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::ReadDataFile() (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::ReadDataFile()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:178:8
    t.method("ReadDataFile", static_cast<void (G4OpticalSurface::*)() >(&G4OpticalSurface::ReadDataFile));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::ReadLUTFile() (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::ReadLUTFile()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:184:8
    t.method("ReadLUTFile", static_cast<void (G4OpticalSurface::*)() >(&G4OpticalSurface::ReadLUTFile));

    DEBUG_MSG("Adding wrapper for G4double G4OpticalSurface::GetAngularDistributionValue(G4int, G4int, G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4OpticalSurface::GetAngularDistributionValue(G4int, G4int, G4int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:187:19
    t.method("GetAngularDistributionValue", static_cast<G4double (G4OpticalSurface::*)(G4int, G4int, G4int) >(&G4OpticalSurface::GetAngularDistributionValue));

    DEBUG_MSG("Adding wrapper for G4double G4OpticalSurface::GetAngularDistributionValueLUT(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4OpticalSurface::GetAngularDistributionValueLUT(G4int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:190:19
    t.method("GetAngularDistributionValueLUT", static_cast<G4double (G4OpticalSurface::*)(G4int) >(&G4OpticalSurface::GetAngularDistributionValueLUT));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::ReadLUTDAVISFile() (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::ReadLUTDAVISFile()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:193:8
    t.method("ReadLUTDAVISFile", static_cast<void (G4OpticalSurface::*)() >(&G4OpticalSurface::ReadLUTDAVISFile));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::ReadReflectivityLUTFile() (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::ReadReflectivityLUTFile()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:196:8
    t.method("ReadReflectivityLUTFile", static_cast<void (G4OpticalSurface::*)() >(&G4OpticalSurface::ReadReflectivityLUTFile));

    DEBUG_MSG("Adding wrapper for G4double G4OpticalSurface::GetReflectivityLUTValue(G4int) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4OpticalSurface::GetReflectivityLUTValue(G4int)
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:199:19
    t.method("GetReflectivityLUTValue", static_cast<G4double (G4OpticalSurface::*)(G4int) >(&G4OpticalSurface::GetReflectivityLUTValue));

    DEBUG_MSG("Adding wrapper for G4int G4OpticalSurface::GetInmax() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4OpticalSurface::GetInmax()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:202:9
    t.method("GetInmax", static_cast<G4int (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetInmax));

    DEBUG_MSG("Adding wrapper for G4int G4OpticalSurface::GetLUTbins() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4OpticalSurface::GetLUTbins()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:205:9
    t.method("GetLUTbins", static_cast<G4int (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetLUTbins));

    DEBUG_MSG("Adding wrapper for G4int G4OpticalSurface::GetRefMax() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4OpticalSurface::GetRefMax()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:208:9
    t.method("GetRefMax", static_cast<G4int (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetRefMax));

    DEBUG_MSG("Adding wrapper for G4int G4OpticalSurface::GetThetaIndexMax() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4OpticalSurface::GetThetaIndexMax()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:210:9
    t.method("GetThetaIndexMax", static_cast<G4int (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetThetaIndexMax));

    DEBUG_MSG("Adding wrapper for G4int G4OpticalSurface::GetPhiIndexMax() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4OpticalSurface::GetPhiIndexMax()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:211:9
    t.method("GetPhiIndexMax", static_cast<G4int (G4OpticalSurface::*)()  const>(&G4OpticalSurface::GetPhiIndexMax));

    DEBUG_MSG("Adding wrapper for void G4OpticalSurface::ReadDichroicFile() (" __HERE__ ")");
    // signature to use in the veto list: void G4OpticalSurface::ReadDichroicFile()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:214:8
    t.method("ReadDichroicFile", static_cast<void (G4OpticalSurface::*)() >(&G4OpticalSurface::ReadDichroicFile));

    DEBUG_MSG("Adding wrapper for G4Physics2DVector * G4OpticalSurface::GetDichroicVector() (" __HERE__ ")");
    // signature to use in the veto list: G4Physics2DVector * G4OpticalSurface::GetDichroicVector()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4OpticalSurface.hh:216:29
    t.method("GetDichroicVector", static_cast<G4Physics2DVector * (G4OpticalSurface::*)() >(&G4OpticalSurface::GetDichroicVector));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4OpticalSurface>> type_;
};
std::shared_ptr<Wrapper> newJlG4OpticalSurface(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4OpticalSurface(module));
}
