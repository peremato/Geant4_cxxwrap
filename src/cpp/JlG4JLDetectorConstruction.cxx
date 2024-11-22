// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4JLDetectorConstruction> : std::false_type { };
  template<> struct DefaultConstructible<G4JLDetectorConstruction> : std::false_type { };
template<> struct SuperType<G4JLDetectorConstruction> { typedef G4VUserDetectorConstruction type; };
}

// Class generating the wrapper for type G4JLDetectorConstruction
// signature to use in the veto file: G4JLDetectorConstruction
struct JlG4JLDetectorConstruction: public Wrapper {

  JlG4JLDetectorConstruction(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4JLDetectorConstruction (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:79:7
    jlcxx::TypeWrapper<G4JLDetectorConstruction>  t = jlModule.add_type<G4JLDetectorConstruction>("G4JLDetectorConstruction",
      jlcxx::julia_base_type<G4VUserDetectorConstruction>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4JLDetectorConstruction>>(new jlcxx::TypeWrapper<G4JLDetectorConstruction>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4JLDetectorConstruction::G4JLDetectorConstruction(construct_f, void *, sdandf_f, void *) (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:81:5
    t.constructor<construct_f, void *>(/*finalize=*/jlcxx::finalize_policy::no);
    t.constructor<construct_f, void *, sdandf_f>(/*finalize=*/jlcxx::finalize_policy::no);
    t.constructor<construct_f, void *, sdandf_f, void *>(/*finalize=*/jlcxx::finalize_policy::no);

    DEBUG_MSG("Adding wrapper for G4VPhysicalVolume * G4JLDetectorConstruction::Construct() (" __HERE__ ")");
    // signature to use in the veto list: G4VPhysicalVolume * G4JLDetectorConstruction::Construct()
    // defined in ./cpp/Geant4Wrap.h:85:32
    t.method("Construct", static_cast<G4VPhysicalVolume * (G4JLDetectorConstruction::*)() >(&G4JLDetectorConstruction::Construct));

    DEBUG_MSG("Adding wrapper for void G4JLDetectorConstruction::ConstructSDandField() (" __HERE__ ")");
    // signature to use in the veto list: void G4JLDetectorConstruction::ConstructSDandField()
    // defined in ./cpp/Geant4Wrap.h:86:18
    t.method("ConstructSDandField", static_cast<void (G4JLDetectorConstruction::*)() >(&G4JLDetectorConstruction::ConstructSDandField));

    DEBUG_MSG("Adding wrapper for void G4JLDetectorConstruction::SetSensitiveDetector(const G4String &, G4JLSensDet *, G4bool) (" __HERE__ ")");
    // signature to use in the veto list: void G4JLDetectorConstruction::SetSensitiveDetector(const G4String &, G4JLSensDet *, G4bool)
    // defined in ./cpp/Geant4Wrap.h:87:10
    t.method("SetSensitiveDetector", static_cast<void (G4JLDetectorConstruction::*)(const G4String &, G4JLSensDet *, G4bool) >(&G4JLDetectorConstruction::SetSensitiveDetector));
    t.method("SetSensitiveDetector", [](G4JLDetectorConstruction& a, const G4String & arg0, G4JLSensDet * arg1)->void { a.SetSensitiveDetector(arg0, arg1); });
    t.method("SetSensitiveDetector", [](G4JLDetectorConstruction* a, const G4String & arg0, G4JLSensDet * arg1)->void { a->SetSensitiveDetector(arg0, arg1); });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4JLDetectorConstruction>> type_;
};
std::shared_ptr<Wrapper> newJlG4JLDetectorConstruction(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4JLDetectorConstruction(module));
}
