// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4Orb> : std::false_type { };
  template<> struct DefaultConstructible<G4Orb> : std::false_type { };
template<> struct SuperType<G4Orb> { typedef G4CSGSolid type; };
}

// Class generating the wrapper for type G4Orb
// signature to use in the veto file: G4Orb
struct JlG4Orb: public Wrapper {

  JlG4Orb(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4Orb (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:55:7
    jlcxx::TypeWrapper<G4Orb>  t = jlModule.add_type<G4Orb>("G4Orb",
      jlcxx::julia_base_type<G4CSGSolid>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4Orb>>(new jlcxx::TypeWrapper<G4Orb>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4Orb::G4Orb(const G4String &, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:59:5
    t.constructor<const G4String &, G4double>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for G4double G4Orb::GetRadius() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Orb::GetRadius()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:65:21
    t.method("GetRadius", static_cast<G4double (G4Orb::*)()  const>(&G4Orb::GetRadius));

    DEBUG_MSG("Adding wrapper for G4double G4Orb::GetRadialTolerance() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Orb::GetRadialTolerance()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:66:21
    t.method("GetRadialTolerance", static_cast<G4double (G4Orb::*)()  const>(&G4Orb::GetRadialTolerance));

    DEBUG_MSG("Adding wrapper for void G4Orb::SetRadius(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4Orb::SetRadius(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:68:17
    t.method("SetRadius", static_cast<void (G4Orb::*)(G4double) >(&G4Orb::SetRadius));

    DEBUG_MSG("Adding wrapper for G4double G4Orb::GetCubicVolume() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Orb::GetCubicVolume()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:72:21
    t.method("GetCubicVolume", static_cast<G4double (G4Orb::*)() >(&G4Orb::GetCubicVolume));

    DEBUG_MSG("Adding wrapper for G4double G4Orb::GetSurfaceArea() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Orb::GetSurfaceArea()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:73:21
    t.method("GetSurfaceArea", static_cast<G4double (G4Orb::*)() >(&G4Orb::GetSurfaceArea));

    DEBUG_MSG("Adding wrapper for void G4Orb::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) (" __HERE__ ")");
    // signature to use in the veto list: void G4Orb::ComputeDimensions(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:75:10
    t.method("ComputeDimensions", static_cast<void (G4Orb::*)(G4VPVParameterisation *, const G4int, const G4VPhysicalVolume *) >(&G4Orb::ComputeDimensions));

    DEBUG_MSG("Adding wrapper for void G4Orb::BoundingLimits(G4ThreeVector &, G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4Orb::BoundingLimits(G4ThreeVector &, G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:79:10
    t.method("BoundingLimits", static_cast<void (G4Orb::*)(G4ThreeVector &, G4ThreeVector &)  const>(&G4Orb::BoundingLimits));

    DEBUG_MSG("Adding wrapper for EInside G4Orb::Inside(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: EInside G4Orb::Inside(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:86:13
    t.method("Inside", static_cast<EInside (G4Orb::*)(const G4ThreeVector &)  const>(&G4Orb::Inside));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Orb::SurfaceNormal(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Orb::SurfaceNormal(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:88:19
    t.method("SurfaceNormal", static_cast<G4ThreeVector (G4Orb::*)(const G4ThreeVector &)  const>(&G4Orb::SurfaceNormal));

    DEBUG_MSG("Adding wrapper for G4double G4Orb::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Orb::DistanceToIn(const G4ThreeVector &, const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:90:14
    t.method("DistanceToIn", static_cast<G4double (G4Orb::*)(const G4ThreeVector &, const G4ThreeVector &)  const>(&G4Orb::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Orb::DistanceToIn(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Orb::DistanceToIn(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:93:14
    t.method("DistanceToIn", static_cast<G4double (G4Orb::*)(const G4ThreeVector &)  const>(&G4Orb::DistanceToIn));

    DEBUG_MSG("Adding wrapper for G4double G4Orb::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Orb::DistanceToOut(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:95:14
    t.method("DistanceToOut", static_cast<G4double (G4Orb::*)(const G4ThreeVector &, const G4ThreeVector &, const G4bool, G4bool *, G4ThreeVector *)  const>(&G4Orb::DistanceToOut));
    t.method("DistanceToOut", [](G4Orb const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a.DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Orb const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a.DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Orb const& a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a.DistanceToOut(arg0, arg1, arg2, arg3); });
    t.method("DistanceToOut", [](G4Orb const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1)->G4double { return a->DistanceToOut(arg0, arg1); });
    t.method("DistanceToOut", [](G4Orb const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2)->G4double { return a->DistanceToOut(arg0, arg1, arg2); });
    t.method("DistanceToOut", [](G4Orb const* a, const G4ThreeVector & arg0, const G4ThreeVector & arg1, const G4bool arg2, G4bool * arg3)->G4double { return a->DistanceToOut(arg0, arg1, arg2, arg3); });

    DEBUG_MSG("Adding wrapper for G4double G4Orb::DistanceToOut(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: G4double G4Orb::DistanceToOut(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:101:14
    t.method("DistanceToOut", static_cast<G4double (G4Orb::*)(const G4ThreeVector &)  const>(&G4Orb::DistanceToOut));

    DEBUG_MSG("Adding wrapper for G4GeometryType G4Orb::GetEntityType() (" __HERE__ ")");
    // signature to use in the veto list: G4GeometryType G4Orb::GetEntityType()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:103:20
    t.method("GetEntityType", static_cast<G4GeometryType (G4Orb::*)()  const>(&G4Orb::GetEntityType));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4Orb::GetPointOnSurface() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4Orb::GetPointOnSurface()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:105:19
    t.method("GetPointOnSurface", static_cast<G4ThreeVector (G4Orb::*)()  const>(&G4Orb::GetPointOnSurface));

    DEBUG_MSG("Adding wrapper for G4VSolid * G4Orb::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4VSolid * G4Orb::Clone()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:107:15
    t.method("Clone", static_cast<G4VSolid * (G4Orb::*)()  const>(&G4Orb::Clone));

    DEBUG_MSG("Adding wrapper for G4Polyhedron * G4Orb::CreatePolyhedron() (" __HERE__ ")");
    // signature to use in the veto list: G4Polyhedron * G4Orb::CreatePolyhedron()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:115:19
    t.method("CreatePolyhedron", static_cast<G4Polyhedron * (G4Orb::*)()  const>(&G4Orb::CreatePolyhedron));


    DEBUG_MSG("Adding wrapper for void G4Orb::G4Orb(const G4Orb &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:122:5
    t.constructor<const G4Orb &>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for G4Orb & G4Orb::operator=(const G4Orb &) (" __HERE__ ")");
    // signature to use in the veto list: G4Orb & G4Orb::operator=(const G4Orb &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4Orb.hh:123:12
    t.method("assign", static_cast<G4Orb & (G4Orb::*)(const G4Orb &) >(&G4Orb::operator=));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4Orb>> type_;
};
std::shared_ptr<Wrapper> newJlG4Orb(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4Orb(module));
}
