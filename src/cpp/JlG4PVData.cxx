// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4PVData> : std::false_type { };
  template<> struct DefaultConstructible<G4PVData> : std::false_type { };
}

// Class generating the wrapper for type G4PVData
// signature to use in the veto file: G4PVData
struct JlG4PVData: public Wrapper {

  JlG4PVData(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4PVData (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:55:7
    jlcxx::TypeWrapper<G4PVData>  t = jlModule.add_type<G4PVData>("G4PVData");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4PVData>>(new jlcxx::TypeWrapper<G4PVData>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for void G4PVData::initialize() (" __HERE__ ")");
    // signature to use in the veto list: void G4PVData::initialize()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:65:10
    t.method("initialize", static_cast<void (G4PVData::*)() >(&G4PVData::initialize));

    DEBUG_MSG("Adding frot methods  to provide read access to the field frot (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:71:23
    // signature to use in the veto list: G4PVData::frot
    t.method("frot", [](const G4PVData& a) -> G4RotationMatrix * { return a.frot; });
    t.method("frot", [](G4PVData& a) -> G4RotationMatrix * { return a.frot; });
    t.method("frot", [](const G4PVData* a) -> G4RotationMatrix * { return a->frot; });
    t.method("frot", [](G4PVData* a) -> G4RotationMatrix * { return a->frot; });
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:71:23
    // signature to use in the veto list: G4PVData::frot
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding frot! methods to provide write access to the field frot (" __HERE__ ")");
    t.method("frot!", [](G4PVData& a, G4RotationMatrix * val) -> G4RotationMatrix * { return a.frot = val; });

    DEBUG_MSG("Adding frot! methods to provide write access to the field frot (" __HERE__ ")");
    t.method("frot!", [](G4PVData* a, G4RotationMatrix * val) -> G4RotationMatrix * { return a->frot = val; });

    DEBUG_MSG("Adding tx methods  to provide read access to the field tx (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:72:14
    // signature to use in the veto list: G4PVData::tx
    t.method("tx", [](const G4PVData& a) -> G4double { return a.tx; });
    t.method("tx", [](G4PVData& a) -> G4double { return a.tx; });
    t.method("tx", [](const G4PVData* a) -> G4double { return a->tx; });
    t.method("tx", [](G4PVData* a) -> G4double { return a->tx; });
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:72:14
    // signature to use in the veto list: G4PVData::tx
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding tx! methods to provide write access to the field tx (" __HERE__ ")");
    t.method("tx!", [](G4PVData& a, G4double val) -> G4double { return a.tx = val; });

    DEBUG_MSG("Adding tx! methods to provide write access to the field tx (" __HERE__ ")");
    t.method("tx!", [](G4PVData* a, G4double val) -> G4double { return a->tx = val; });

    DEBUG_MSG("Adding ty methods  to provide read access to the field ty (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:72:23
    // signature to use in the veto list: G4PVData::ty
    t.method("ty", [](const G4PVData& a) -> G4double { return a.ty; });
    t.method("ty", [](G4PVData& a) -> G4double { return a.ty; });
    t.method("ty", [](const G4PVData* a) -> G4double { return a->ty; });
    t.method("ty", [](G4PVData* a) -> G4double { return a->ty; });
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:72:23
    // signature to use in the veto list: G4PVData::ty
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding ty! methods to provide write access to the field ty (" __HERE__ ")");
    t.method("ty!", [](G4PVData& a, G4double val) -> G4double { return a.ty = val; });

    DEBUG_MSG("Adding ty! methods to provide write access to the field ty (" __HERE__ ")");
    t.method("ty!", [](G4PVData* a, G4double val) -> G4double { return a->ty = val; });

    DEBUG_MSG("Adding tz methods  to provide read access to the field tz (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:72:32
    // signature to use in the veto list: G4PVData::tz
    t.method("tz", [](const G4PVData& a) -> G4double { return a.tz; });
    t.method("tz", [](G4PVData& a) -> G4double { return a.tz; });
    t.method("tz", [](const G4PVData* a) -> G4double { return a->tz; });
    t.method("tz", [](G4PVData* a) -> G4double { return a->tz; });
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VPhysicalVolume.hh:72:32
    // signature to use in the veto list: G4PVData::tz
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding tz! methods to provide write access to the field tz (" __HERE__ ")");
    t.method("tz!", [](G4PVData& a, G4double val) -> G4double { return a.tz = val; });

    DEBUG_MSG("Adding tz! methods to provide write access to the field tz (" __HERE__ ")");
    t.method("tz!", [](G4PVData* a, G4double val) -> G4double { return a->tz = val; });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4PVData>> type_;
};
std::shared_ptr<Wrapper> newJlG4PVData(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4PVData(module));
}
