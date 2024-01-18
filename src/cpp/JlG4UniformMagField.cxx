// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4UniformMagField> : std::false_type { };
  template<> struct DefaultConstructible<G4UniformMagField> : std::false_type { };
template<> struct SuperType<G4UniformMagField> { typedef G4MagneticField type; };
}

// Class generating the wrapper for type G4UniformMagField
// signature to use in the veto file: G4UniformMagField
struct JlG4UniformMagField: public Wrapper {

  JlG4UniformMagField(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4UniformMagField (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UniformMagField.hh:41:7
    jlcxx::TypeWrapper<G4UniformMagField>  t = jlModule.add_type<G4UniformMagField>("G4UniformMagField",
      jlcxx::julia_base_type<G4MagneticField>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4UniformMagField>>(new jlcxx::TypeWrapper<G4UniformMagField>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4UniformMagField::G4UniformMagField(const G4ThreeVector &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UniformMagField.hh:45:5
    t.constructor<const G4ThreeVector &>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4UniformMagField::G4UniformMagField(G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UniformMagField.hh:48:5
    t.constructor<G4double, G4double, G4double>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4UniformMagField::G4UniformMagField(const G4UniformMagField &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UniformMagField.hh:54:5
    t.constructor<const G4UniformMagField &>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for G4UniformMagField & G4UniformMagField::operator=(const G4UniformMagField &) (" __HERE__ ")");
    // signature to use in the veto list: G4UniformMagField & G4UniformMagField::operator=(const G4UniformMagField &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UniformMagField.hh:55:24
    t.method("assign", static_cast<G4UniformMagField & (G4UniformMagField::*)(const G4UniformMagField &) >(&G4UniformMagField::operator=));

    DEBUG_MSG("Adding wrapper for void G4UniformMagField::SetFieldValue(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4UniformMagField::SetFieldValue(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UniformMagField.hh:61:10
    t.method("SetFieldValue", static_cast<void (G4UniformMagField::*)(const G4ThreeVector &) >(&G4UniformMagField::SetFieldValue));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4UniformMagField::GetConstantFieldValue() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4UniformMagField::GetConstantFieldValue()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UniformMagField.hh:63:19
    t.method("GetConstantFieldValue", static_cast<G4ThreeVector (G4UniformMagField::*)()  const>(&G4UniformMagField::GetConstantFieldValue));

    DEBUG_MSG("Adding wrapper for G4Field * G4UniformMagField::Clone() (" __HERE__ ")");
    // signature to use in the veto list: G4Field * G4UniformMagField::Clone()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4UniformMagField.hh:66:14
    t.method("Clone", static_cast<G4Field * (G4UniformMagField::*)()  const>(&G4UniformMagField::Clone));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4UniformMagField>> type_;
};
std::shared_ptr<Wrapper> newJlG4UniformMagField(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4UniformMagField(module));
}