// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4SPSPosDistribution> : std::false_type { };
  template<> struct DefaultConstructible<G4SPSPosDistribution> : std::false_type { };
}

// Class generating the wrapper for type G4SPSPosDistribution
// signature to use in the veto file: G4SPSPosDistribution
struct JlG4SPSPosDistribution: public Wrapper {

  JlG4SPSPosDistribution(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4SPSPosDistribution (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:59:7
    jlcxx::TypeWrapper<G4SPSPosDistribution>  t = jlModule.add_type<G4SPSPosDistribution>("G4SPSPosDistribution");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4SPSPosDistribution>>(new jlcxx::TypeWrapper<G4SPSPosDistribution>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetPosDisType(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetPosDisType(const G4String &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:71:10
    t.method("SetPosDisType", static_cast<void (G4SPSPosDistribution::*)(const G4String &) >(&G4SPSPosDistribution::SetPosDisType));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetPosDisShape(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetPosDisShape(const G4String &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:75:10
    t.method("SetPosDisShape", static_cast<void (G4SPSPosDistribution::*)(const G4String &) >(&G4SPSPosDistribution::SetPosDisShape));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetCentreCoords(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetCentreCoords(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:80:10
    t.method("SetCentreCoords", static_cast<void (G4SPSPosDistribution::*)(const G4ThreeVector &) >(&G4SPSPosDistribution::SetCentreCoords));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetPosRot1(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetPosRot1(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:83:10
    t.method("SetPosRot1", static_cast<void (G4SPSPosDistribution::*)(const G4ThreeVector &) >(&G4SPSPosDistribution::SetPosRot1));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetPosRot2(const G4ThreeVector &) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetPosRot2(const G4ThreeVector &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:87:10
    t.method("SetPosRot2", static_cast<void (G4SPSPosDistribution::*)(const G4ThreeVector &) >(&G4SPSPosDistribution::SetPosRot2));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetHalfX(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetHalfX(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:92:10
    t.method("SetHalfX", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetHalfX));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetHalfY(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetHalfY(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:95:10
    t.method("SetHalfY", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetHalfY));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetHalfZ(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetHalfZ(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:98:10
    t.method("SetHalfZ", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetHalfZ));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetRadius(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetRadius(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:101:10
    t.method("SetRadius", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetRadius));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetRadius0(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetRadius0(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:104:10
    t.method("SetRadius0", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetRadius0));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetBeamSigmaInR(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetBeamSigmaInR(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:107:10
    t.method("SetBeamSigmaInR", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetBeamSigmaInR));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetBeamSigmaInX(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetBeamSigmaInX(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:110:10
    t.method("SetBeamSigmaInX", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetBeamSigmaInX));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetBeamSigmaInY(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetBeamSigmaInY(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:113:10
    t.method("SetBeamSigmaInY", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetBeamSigmaInY));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetParAlpha(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetParAlpha(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:116:10
    t.method("SetParAlpha", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetParAlpha));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetParTheta(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetParTheta(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:119:10
    t.method("SetParTheta", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetParTheta));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetParPhi(G4double) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetParPhi(G4double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:122:10
    t.method("SetParPhi", static_cast<void (G4SPSPosDistribution::*)(G4double) >(&G4SPSPosDistribution::SetParPhi));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::ConfineSourceToVolume(const G4String &) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::ConfineSourceToVolume(const G4String &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:125:10
    t.method("ConfineSourceToVolume", static_cast<void (G4SPSPosDistribution::*)(const G4String &) >(&G4SPSPosDistribution::ConfineSourceToVolume));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetBiasRndm(G4SPSRandomGenerator *) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetBiasRndm(G4SPSRandomGenerator *)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:128:10
    t.method("SetBiasRndm", static_cast<void (G4SPSPosDistribution::*)(G4SPSRandomGenerator *) >(&G4SPSPosDistribution::SetBiasRndm));

    DEBUG_MSG("Adding wrapper for void G4SPSPosDistribution::SetVerbosity(G4int) (" __HERE__ ")");
    // signature to use in the veto list: void G4SPSPosDistribution::SetVerbosity(G4int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:131:10
    t.method("SetVerbosity", static_cast<void (G4SPSPosDistribution::*)(G4int) >(&G4SPSPosDistribution::SetVerbosity));

    DEBUG_MSG("Adding wrapper for G4ThreeVector G4SPSPosDistribution::GenerateOne() (" __HERE__ ")");
    // signature to use in the veto list: G4ThreeVector G4SPSPosDistribution::GenerateOne()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:134:19
    t.method("GenerateOne", static_cast<G4ThreeVector (G4SPSPosDistribution::*)() >(&G4SPSPosDistribution::GenerateOne));

    DEBUG_MSG("Adding wrapper for const G4String & G4SPSPosDistribution::GetPosDisType() (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4SPSPosDistribution::GetPosDisType()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:137:21
    t.method("GetPosDisType", static_cast<const G4String & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetPosDisType));

    DEBUG_MSG("Adding wrapper for const G4String & G4SPSPosDistribution::GetPosDisShape() (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4SPSPosDistribution::GetPosDisShape()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:138:21
    t.method("GetPosDisShape", static_cast<const G4String & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetPosDisShape));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4SPSPosDistribution::GetCentreCoords() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4SPSPosDistribution::GetCentreCoords()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:139:26
    t.method("GetCentreCoords", static_cast<const G4ThreeVector & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetCentreCoords));

    DEBUG_MSG("Adding wrapper for G4double G4SPSPosDistribution::GetHalfX() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SPSPosDistribution::GetHalfX()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:140:14
    t.method("GetHalfX", static_cast<G4double (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetHalfX));

    DEBUG_MSG("Adding wrapper for G4double G4SPSPosDistribution::GetHalfY() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SPSPosDistribution::GetHalfY()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:141:14
    t.method("GetHalfY", static_cast<G4double (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetHalfY));

    DEBUG_MSG("Adding wrapper for G4double G4SPSPosDistribution::GetHalfZ() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SPSPosDistribution::GetHalfZ()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:142:14
    t.method("GetHalfZ", static_cast<G4double (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetHalfZ));

    DEBUG_MSG("Adding wrapper for G4double G4SPSPosDistribution::GetRadius() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SPSPosDistribution::GetRadius()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:143:14
    t.method("GetRadius", static_cast<G4double (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetRadius));

    DEBUG_MSG("Adding wrapper for G4double G4SPSPosDistribution::GetRadius0() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SPSPosDistribution::GetRadius0()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:144:21
    t.method("GetRadius0", static_cast<G4double (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetRadius0));

    DEBUG_MSG("Adding wrapper for G4double G4SPSPosDistribution::GetParAlpha() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SPSPosDistribution::GetParAlpha()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:145:21
    t.method("GetParAlpha", static_cast<G4double (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetParAlpha));

    DEBUG_MSG("Adding wrapper for G4double G4SPSPosDistribution::GetParTheta() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SPSPosDistribution::GetParTheta()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:146:21
    t.method("GetParTheta", static_cast<G4double (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetParTheta));

    DEBUG_MSG("Adding wrapper for G4double G4SPSPosDistribution::GetParPhi() (" __HERE__ ")");
    // signature to use in the veto list: G4double G4SPSPosDistribution::GetParPhi()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:147:21
    t.method("GetParPhi", static_cast<G4double (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetParPhi));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4SPSPosDistribution::GetRotx() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4SPSPosDistribution::GetRotx()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:148:33
    t.method("GetRotx", static_cast<const G4ThreeVector & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetRotx));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4SPSPosDistribution::GetRoty() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4SPSPosDistribution::GetRoty()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:149:33
    t.method("GetRoty", static_cast<const G4ThreeVector & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetRoty));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4SPSPosDistribution::GetRotz() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4SPSPosDistribution::GetRotz()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:150:33
    t.method("GetRotz", static_cast<const G4ThreeVector & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetRotz));

    DEBUG_MSG("Adding wrapper for G4bool G4SPSPosDistribution::GetConfined() (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4SPSPosDistribution::GetConfined()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:151:19
    t.method("GetConfined", static_cast<G4bool (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetConfined));

    DEBUG_MSG("Adding wrapper for const G4String & G4SPSPosDistribution::GetConfineVolume() (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4SPSPosDistribution::GetConfineVolume()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:152:28
    t.method("GetConfineVolume", static_cast<const G4String & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetConfineVolume));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4SPSPosDistribution::GetSideRefVec1() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4SPSPosDistribution::GetSideRefVec1()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:154:26
    t.method("GetSideRefVec1", static_cast<const G4ThreeVector & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetSideRefVec1));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4SPSPosDistribution::GetSideRefVec2() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4SPSPosDistribution::GetSideRefVec2()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:155:26
    t.method("GetSideRefVec2", static_cast<const G4ThreeVector & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetSideRefVec2));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4SPSPosDistribution::GetSideRefVec3() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4SPSPosDistribution::GetSideRefVec3()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:156:26
    t.method("GetSideRefVec3", static_cast<const G4ThreeVector & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetSideRefVec3));

    DEBUG_MSG("Adding wrapper for const G4String & G4SPSPosDistribution::GetSourcePosType() (" __HERE__ ")");
    // signature to use in the veto list: const G4String & G4SPSPosDistribution::GetSourcePosType()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:157:21
    t.method("GetSourcePosType", static_cast<const G4String & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetSourcePosType));

    DEBUG_MSG("Adding wrapper for const G4ThreeVector & G4SPSPosDistribution::GetParticlePos() (" __HERE__ ")");
    // signature to use in the veto list: const G4ThreeVector & G4SPSPosDistribution::GetParticlePos()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4SPSPosDistribution.hh:158:26
    t.method("GetParticlePos", static_cast<const G4ThreeVector & (G4SPSPosDistribution::*)()  const>(&G4SPSPosDistribution::GetParticlePos));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4SPSPosDistribution>> type_;
};
std::shared_ptr<Wrapper> newJlG4SPSPosDistribution(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4SPSPosDistribution(module));
}