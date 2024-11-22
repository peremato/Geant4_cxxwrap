// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4ReplicaData> : std::false_type { };
  template<> struct DefaultConstructible<G4ReplicaData> : std::false_type { };
}

// Class generating the wrapper for type G4ReplicaData
// signature to use in the veto file: G4ReplicaData
struct JlG4ReplicaData: public Wrapper {

  JlG4ReplicaData(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4ReplicaData (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4PVReplica.hh:72:7
    jlcxx::TypeWrapper<G4ReplicaData>  t = jlModule.add_type<G4ReplicaData>("G4ReplicaData");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4ReplicaData>>(new jlcxx::TypeWrapper<G4ReplicaData>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/jlcxx::finalize_policy::yes);

    DEBUG_MSG("Adding wrapper for void G4ReplicaData::initialize() (" __HERE__ ")");
    // signature to use in the veto list: void G4ReplicaData::initialize()
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4PVReplica.hh:81:8
    t.method("initialize", static_cast<void (G4ReplicaData::*)() >(&G4ReplicaData::initialize));

    DEBUG_MSG("Adding fcopyNo methods  to provide read access to the field fcopyNo (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4PVReplica.hh:83:9
    // signature to use in the veto list: G4ReplicaData::fcopyNo
    t.method("fcopyNo", [](const G4ReplicaData& a) -> G4int { return a.fcopyNo; });
    t.method("fcopyNo", [](G4ReplicaData& a) -> G4int { return a.fcopyNo; });
    t.method("fcopyNo", [](const G4ReplicaData* a) -> G4int { return a->fcopyNo; });
    t.method("fcopyNo", [](G4ReplicaData* a) -> G4int { return a->fcopyNo; });
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4PVReplica.hh:83:9
    // signature to use in the veto list: G4ReplicaData::fcopyNo
    // with ! suffix to veto the setter only.
    DEBUG_MSG("Adding fcopyNo! methods to provide write access to the field fcopyNo (" __HERE__ ")");
    t.method("fcopyNo!", [](G4ReplicaData& a, G4int val) -> G4int { return a.fcopyNo = val; });

    DEBUG_MSG("Adding fcopyNo! methods to provide write access to the field fcopyNo (" __HERE__ ")");
    t.method("fcopyNo!", [](G4ReplicaData* a, G4int val) -> G4int { return a->fcopyNo = val; });
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4ReplicaData>> type_;
};
std::shared_ptr<Wrapper> newJlG4ReplicaData(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4ReplicaData(module));
}
