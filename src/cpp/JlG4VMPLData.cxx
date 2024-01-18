// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VMPLData> : std::false_type { };
  template<> struct DefaultConstructible<G4VMPLData> : std::false_type { };
}

// Class generating the wrapper for type G4VMPLData
// signature to use in the veto file: G4VMPLData
struct JlG4VMPLData: public Wrapper {

  JlG4VMPLData(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VMPLData (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VModularPhysicsList.hh:54:7
    jlcxx::TypeWrapper<G4VMPLData>  t = jlModule.add_type<G4VMPLData>("G4VMPLData");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VMPLData>>(new jlcxx::TypeWrapper<G4VMPLData>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);

    DEBUG_MSG("Adding wrapper for void G4VMPLData::initialize() (" __HERE__ ")");
    // signature to use in the veto list: void G4VMPLData::initialize()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VModularPhysicsList.hh:60:10
    t.method("initialize", static_cast<void (G4VMPLData::*)() >(&G4VMPLData::initialize));

    DEBUG_MSG("Adding physicsVector methods  to provide read access to the field physicsVector (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VModularPhysicsList.hh:63:28
    // signature to use in the veto list: G4VMPLData::physicsVector
    t.method("physicsVector", [](const G4VMPLData& a) -> G4VMPLData::G4PhysConstVectorData * { return a.physicsVector; });
    t.method("physicsVector", [](G4VMPLData& a) -> G4VMPLData::G4PhysConstVectorData * { return a.physicsVector; });
    t.method("physicsVector", [](const G4VMPLData* a) -> G4VMPLData::G4PhysConstVectorData * { return a->physicsVector; });
    t.method("physicsVector", [](G4VMPLData* a) -> G4VMPLData::G4PhysConstVectorData * { return a->physicsVector; });
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VModularPhysicsList.hh:63:28
    // signature to use in the veto list: G4VMPLData::physicsVector
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding physicsVector! methods to provide write access to the field physicsVector (" __HERE__ ")");
    t.method("physicsVector!", [](G4VMPLData& a, G4VMPLData::G4PhysConstVectorData * val) -> G4VMPLData::G4PhysConstVectorData * { return a.physicsVector = val; });

    DEBUG_MSG("Adding physicsVector! methods to provide write access to the field physicsVector (" __HERE__ ")");
    t.method("physicsVector!", [](G4VMPLData* a, G4VMPLData::G4PhysConstVectorData * val) -> G4VMPLData::G4PhysConstVectorData * { return a->physicsVector = val; });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VMPLData>> type_;
};
std::shared_ptr<Wrapper> newJlG4VMPLData(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VMPLData(module));
}