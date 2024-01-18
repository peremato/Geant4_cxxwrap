// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4PrimaryParticle> : std::false_type { };
  template<> struct DefaultConstructible<G4PrimaryParticle> : std::false_type { };
}

// Class generating the wrapper for type G4PrimaryParticle
// signature to use in the veto file: G4PrimaryParticle
struct JlG4PrimaryParticle: public Wrapper {

  JlG4PrimaryParticle(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4PrimaryParticle (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:66:7
    jlcxx::TypeWrapper<G4PrimaryParticle>  t = jlModule.add_type<G4PrimaryParticle>("G4PrimaryParticle");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4PrimaryParticle>>(new jlcxx::TypeWrapper<G4PrimaryParticle>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::G4PrimaryParticle(G4int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:71:5
    t.constructor<G4int>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::G4PrimaryParticle(G4int, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:72:5
    t.constructor<G4int, G4double, G4double, G4double>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::G4PrimaryParticle(G4int, G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:73:5
    t.constructor<G4int, G4double, G4double, G4double, G4double>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::G4PrimaryParticle(const G4ParticleDefinition *) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:74:5
    t.constructor<const G4ParticleDefinition *>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::G4PrimaryParticle(const G4ParticleDefinition *, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:75:5
    t.constructor<const G4ParticleDefinition *, G4double, G4double, G4double>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::G4PrimaryParticle(const G4ParticleDefinition *, G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:76:5
    t.constructor<const G4ParticleDefinition *, G4double, G4double, G4double, G4double>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::G4PrimaryParticle(const G4PrimaryParticle &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:85:5
    t.constructor<const G4PrimaryParticle &>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for G4PrimaryParticle & G4PrimaryParticle::operator=(const G4PrimaryParticle &) (" __HERE__ ")");
    // signature to use in the veto list: G4PrimaryParticle & G4PrimaryParticle::operator=(const G4PrimaryParticle &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:86:24
    t.method("assign", static_cast<G4PrimaryParticle & (G4PrimaryParticle::*)(const G4PrimaryParticle &) >(&G4PrimaryParticle::operator=));
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for G4bool G4PrimaryParticle::operator==(const G4PrimaryParticle &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4PrimaryParticle::operator==(const G4PrimaryParticle &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:90:12
    t.method("==", static_cast<G4bool (G4PrimaryParticle::*)(const G4PrimaryParticle &)  const>(&G4PrimaryParticle::operator==));

    DEBUG_MSG("Adding wrapper for G4bool G4PrimaryParticle::operator!=(const G4PrimaryParticle &) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4PrimaryParticle::operator!=(const G4PrimaryParticle &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:91:12
    t.method("!=", static_cast<G4bool (G4PrimaryParticle::*)(const G4PrimaryParticle &)  const>(&G4PrimaryParticle::operator!=));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::Print() (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::Print()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:98:10
    t.method("Print", static_cast<void (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::Print));

    DEBUG_MSG("Adding wrapper for G4int G4PrimaryParticle::GetPDGcode() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4PrimaryParticle::GetPDGcode()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:109:18
    t.method("GetPDGcode", static_cast<G4int (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetPDGcode));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetPDGcode(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetPDGcode(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:110:10
    t.method("SetPDGcode", static_cast<void (G4PrimaryParticle::*)(G4int) >(&G4PrimaryParticle::SetPDGcode));

    DEBUG_MSG("Adding wrapper for G4ParticleDefinition * G4PrimaryParticle::GetG4code() (" __HERE__ ")");
    // signature to use in the veto list: G4ParticleDefinition * G4PrimaryParticle::GetG4code()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:111:34
    t.method("GetG4code", static_cast<G4ParticleDefinition * (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetG4code));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetG4code(const G4ParticleDefinition *) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetG4code(const G4ParticleDefinition *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:112:17
    t.method("SetG4code", static_cast<void (G4PrimaryParticle::*)(const G4ParticleDefinition *) >(&G4PrimaryParticle::SetG4code));

    DEBUG_MSG("Adding wrapper for const G4ParticleDefinition * G4PrimaryParticle::GetParticleDefinition() (" __HERE__ ")");
    // signature to use in the veto list: const G4ParticleDefinition * G4PrimaryParticle::GetParticleDefinition()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:113:40
    t.method("GetParticleDefinition", static_cast<const G4ParticleDefinition * (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetParticleDefinition));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetParticleDefinition(const G4ParticleDefinition *) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetParticleDefinition(const G4ParticleDefinition *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:114:10
    t.method("SetParticleDefinition", static_cast<void (G4PrimaryParticle::*)(const G4ParticleDefinition *) >(&G4PrimaryParticle::SetParticleDefinition));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetMass() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetMass()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:115:21
    t.method("GetMass", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetMass));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetMass(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetMass(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:116:17
    t.method("SetMass", static_cast<void (G4PrimaryParticle::*)(G4double) >(&G4PrimaryParticle::SetMass));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetCharge() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetCharge()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:117:21
    t.method("GetCharge", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetCharge));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetCharge(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetCharge(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:118:17
    t.method("SetCharge", static_cast<void (G4PrimaryParticle::*)(G4double) >(&G4PrimaryParticle::SetCharge));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetKineticEnergy() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetKineticEnergy()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:119:21
    t.method("GetKineticEnergy", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetKineticEnergy));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetKineticEnergy(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetKineticEnergy(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:120:17
    t.method("SetKineticEnergy", static_cast<void (G4PrimaryParticle::*)(G4double) >(&G4PrimaryParticle::SetKineticEnergy));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4PrimaryParticle::GetMomentumDirection() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4PrimaryParticle::GetMomentumDirection()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:121:33
    t.method("GetMomentumDirection", static_cast<const G4ThreeVector & (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetMomentumDirection));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetMomentumDirection(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetMomentumDirection(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:122:17
    t.method("SetMomentumDirection", static_cast<void (G4PrimaryParticle::*)(const G4ThreeVector &) >(&G4PrimaryParticle::SetMomentumDirection));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetTotalMomentum() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetTotalMomentum()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:123:21
    t.method("GetTotalMomentum", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetTotalMomentum));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::Set4Momentum(G4double, G4double, G4double, G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::Set4Momentum(G4double, G4double, G4double, G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:124:10
    t.method("Set4Momentum", static_cast<void (G4PrimaryParticle::*)(G4double, G4double, G4double, G4double) >(&G4PrimaryParticle::Set4Momentum));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetTotalEnergy() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetTotalEnergy()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:125:21
    t.method("GetTotalEnergy", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetTotalEnergy));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetTotalEnergy(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetTotalEnergy(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:126:17
    t.method("SetTotalEnergy", static_cast<void (G4PrimaryParticle::*)(G4double) >(&G4PrimaryParticle::SetTotalEnergy));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4PrimaryParticle::GetMomentum() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4PrimaryParticle::GetMomentum()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:127:26
    t.method("GetMomentum", static_cast<G4ThreeVector (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetMomentum));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetMomentum(G4double, G4double, G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetMomentum(G4double, G4double, G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:128:10
    t.method("SetMomentum", static_cast<void (G4PrimaryParticle::*)(G4double, G4double, G4double) >(&G4PrimaryParticle::SetMomentum));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetPx() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetPx()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:129:21
    t.method("GetPx", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetPx));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetPy() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetPy()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:130:21
    t.method("GetPy", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetPy));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetPz() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetPz()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:131:21
    t.method("GetPz", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetPz));

    DEBUG_MSG("Adding wrapper for G4PrimaryParticle * G4PrimaryParticle::GetNext() (" __HERE__ ")");
    // signature to use in the veto list: G4PrimaryParticle * G4PrimaryParticle::GetNext()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:132:31
    t.method("GetNext", static_cast<G4PrimaryParticle * (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetNext));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetNext(G4PrimaryParticle *) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetNext(G4PrimaryParticle *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:133:17
    t.method("SetNext", static_cast<void (G4PrimaryParticle::*)(G4PrimaryParticle *) >(&G4PrimaryParticle::SetNext));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::ClearNext() (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::ClearNext()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:134:17
    t.method("ClearNext", static_cast<void (G4PrimaryParticle::*)() >(&G4PrimaryParticle::ClearNext));

    DEBUG_MSG("Adding wrapper for G4PrimaryParticle * G4PrimaryParticle::GetDaughter() (" __HERE__ ")");
    // signature to use in the veto list: G4PrimaryParticle * G4PrimaryParticle::GetDaughter()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:135:31
    t.method("GetDaughter", static_cast<G4PrimaryParticle * (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetDaughter));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetDaughter(G4PrimaryParticle *) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetDaughter(G4PrimaryParticle *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:136:17
    t.method("SetDaughter", static_cast<void (G4PrimaryParticle::*)(G4PrimaryParticle *) >(&G4PrimaryParticle::SetDaughter));

    DEBUG_MSG("Adding wrapper for G4int G4PrimaryParticle::GetTrackID() (" __HERE__ ")");
    // signature to use in the veto list: G4int G4PrimaryParticle::GetTrackID()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:137:18
    t.method("GetTrackID", static_cast<G4int (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetTrackID));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetTrackID(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetTrackID(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:138:17
    t.method("SetTrackID", static_cast<void (G4PrimaryParticle::*)(G4int) >(&G4PrimaryParticle::SetTrackID));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4PrimaryParticle::GetPolarization() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4PrimaryParticle::GetPolarization()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:139:26
    t.method("GetPolarization", static_cast<G4ThreeVector (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetPolarization));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetPolarization(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetPolarization(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:140:17
    t.method("SetPolarization", static_cast<void (G4PrimaryParticle::*)(const G4ThreeVector &) >(&G4PrimaryParticle::SetPolarization));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetPolarization(G4double, G4double, G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetPolarization(G4double, G4double, G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:141:17
    t.method("SetPolarization", static_cast<void (G4PrimaryParticle::*)(G4double, G4double, G4double) >(&G4PrimaryParticle::SetPolarization));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetPolX() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetPolX()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:142:21
    t.method("GetPolX", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetPolX));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetPolY() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetPolY()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:143:21
    t.method("GetPolY", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetPolY));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetPolZ() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetPolZ()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:144:21
    t.method("GetPolZ", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetPolZ));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetWeight() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetWeight()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:145:21
    t.method("GetWeight", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetWeight));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetWeight(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetWeight(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:146:17
    t.method("SetWeight", static_cast<void (G4PrimaryParticle::*)(G4double) >(&G4PrimaryParticle::SetWeight));

    DEBUG_MSG("Adding wrapper for G4double G4PrimaryParticle::GetProperTime() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4PrimaryParticle::GetProperTime()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:147:21
    t.method("GetProperTime", static_cast<G4double (G4PrimaryParticle::*)()  const>(&G4PrimaryParticle::GetProperTime));

    DEBUG_MSG("Adding wrapper for void G4PrimaryParticle::SetProperTime(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4PrimaryParticle::SetProperTime(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4PrimaryParticle.hh:148:17
    t.method("SetProperTime", static_cast<void (G4PrimaryParticle::*)(G4double) >(&G4PrimaryParticle::SetProperTime));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4PrimaryParticle>> type_;
};
std::shared_ptr<Wrapper> newJlG4PrimaryParticle(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4PrimaryParticle(module));
}