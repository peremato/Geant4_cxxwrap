// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Torus> : std::false_type { };
  template<> struct DefaultConstructible<G4Torus> : std::false_type { };
template<> struct SuperType<G4Torus> { typedef G4CSGSolid type; };
}

// Class generating the wrapper for type G4Torus
// signature to use in the veto file: G4Torus
struct JlG4Torus: public Wrapper {

  JlG4Torus(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Torus (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:91:7
    jlcxx::TypeWrapper<G4Torus>  t = jlModule.add_type<G4Torus>("G4Torus",
      jlcxx::julia_base_type<G4CSGSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Torus>>(new jlcxx::TypeWrapper<G4Torus>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4Torus::G4Torus(const G4String &, G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:96:5
    t.constructor<const G4String &, G4double, G4double, G4double, G4double, G4double>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetRmin() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetRmin()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:107:21
    t.method("GetRmin", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetRmin));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetRmax() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetRmax()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:108:21
    t.method("GetRmax", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetRmax));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetRtor() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetRtor()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:109:21
    t.method("GetRtor", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetRtor));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetSPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetSPhi()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:110:21
    t.method("GetSPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetSPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetDPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetDPhi()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:111:21
    t.method("GetDPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetDPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetSinStartPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetSinStartPhi()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:112:21
    t.method("GetSinStartPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetSinStartPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetCosStartPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetCosStartPhi()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:113:21
    t.method("GetCosStartPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetCosStartPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetSinEndPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetSinEndPhi()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:114:21
    t.method("GetSinEndPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetSinEndPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetCosEndPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetCosEndPhi()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:115:21
    t.method("GetCosEndPhi", static_cast<G4double (G4Torus::*)()  const>(&G4Torus::GetCosEndPhi));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:119:21
    t.method("GetCubicVolume", static_cast<G4double (G4Torus::*)() >(&G4Torus::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:120:21
    t.method("GetSurfaceArea", static_cast<G4double (G4Torus::*)() >(&G4Torus::GetSurfaceArea));

    DEBUG_MSG("Adding wrapper for EInside G4Torus::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4Torus::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:122:13
    t.method("Inside", static_cast<EInside (G4Torus::*)(const G4ThreeVector &)  const>(&G4Torus::Inside));

    DEBUG_MSG("Adding wrapper for void G4Torus::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Torus::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:123:10
    t.method("BoundingLimits", static_cast<void (G4Torus::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Torus::BoundingLimits));

    DEBUG_MSG("Adding wrapper for void G4Torus::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Torus::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:128:10
    t.method("ComputeDimensions", static_cast<void (G4Torus::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Torus::ComputeDimensions));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Torus::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Torus::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:131:19
    t.method("SurfaceNormal", static_cast<G4ThreeVector (G4Torus::*)(const G4ThreeVector &)  const>(&G4Torus::SurfaceNormal));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:132:14
    t.method("DistanceToIn", static_cast<G4double (G4Torus::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Torus::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:134:14
    t.method("DistanceToIn", static_cast<G4double (G4Torus::*)(const G4ThreeVector &)  const>(&G4Torus::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Torus::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:135:14
    t.method("DistanceToOut", static_cast<G4double (G4Torus::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Torus::DistanceToOut));
    t.method("DistanceToOut", [](G4Torus const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Torus const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Torus const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); });
    t.method("DistanceToOut", [](G4Torus const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Torus const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Torus const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for G4double G4Torus::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Torus::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:139:14
    t.method("DistanceToOut", static_cast<G4double (G4Torus::*)(const G4ThreeVector &)  const>(&G4Torus::DistanceToOut));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4Torus::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4Torus::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:141:20
    t.method("GetEntityType", static_cast<G4GeometryType (G4Torus::*)()  const>(&G4Torus::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Torus::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Torus::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:143:19
    t.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Torus::*)()  const>(&G4Torus::GetPointOnSurface));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4Torus::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4Torus::Clone()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:145:15
    t.method("Clone", static_cast<G4VSolid * (G4Torus::*)()  const>(&G4Torus::Clone));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Torus::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Torus::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:152:19
    t.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Torus::*)()  const>(&G4Torus::CreatePolyhedron));

    DEBUG_MSG("Adding wrapper for void G4Torus::SetAllParameters(G4double, G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Torus::SetAllParameters(G4double, G4double, G4double, G4double, G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:154:10
    t.method("SetAllParameters", static_cast<void (G4Torus::*)(G4double, G4double, G4double, G4double, G4double) >(&G4Torus::SetAllParameters));


    DEBUG_MSG("Adding wrapper for void G4Torus::G4Torus(const G4Torus &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:162:5
    t.constructor<const G4Torus &>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for G4Torus & G4Torus::operator=(const G4Torus &) (" __HERE__ ")");
    // signature to use in the veto list: G4Torus & G4Torus::operator=(const G4Torus &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Torus.hh:163:14
    t.method("assign", static_cast<G4Torus & (G4Torus::*)(const G4Torus &) >(&G4Torus::operator=));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Torus>> type_;
};
std::shared_ptr<Wrapper> newJlG4Torus(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Torus(module));
}