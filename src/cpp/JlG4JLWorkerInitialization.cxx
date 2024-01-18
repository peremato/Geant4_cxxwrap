// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4JLWorkerInitialization> : std::false_type { };
  template<> struct DefaultConstructible<G4JLWorkerInitialization> : std::false_type { };
template<> struct SuperType<G4JLWorkerInitialization> { typedef G4UserWorkerInitialization type; };
}

// Class generating the wrapper for type G4JLWorkerInitialization
// signature to use in the veto file: G4JLWorkerInitialization
struct JlG4JLWorkerInitialization: public Wrapper {

  JlG4JLWorkerInitialization(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4JLWorkerInitialization (" __HERE__ ")");
    // defined in ./cpp/Geant4Wrap.h:127:7
    jlcxx::TypeWrapper<G4JLWorkerInitialization>  t = jlModule.add_type<G4JLWorkerInitialization>("G4JLWorkerInitialization",
      jlcxx::julia_base_type<G4UserWorkerInitialization>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4JLWorkerInitialization>>(new jlcxx::TypeWrapper<G4JLWorkerInitialization>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/false);

    DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerInitialize() (" __HERE__ ")");
    // signature to use in the veto list: void G4JLWorkerInitialization::WorkerInitialize()
    // defined in ./cpp/Geant4Wrap.h:131:18
    t.method("WorkerInitialize", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerInitialize));

    DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerStart() (" __HERE__ ")");
    // signature to use in the veto list: void G4JLWorkerInitialization::WorkerStart()
    // defined in ./cpp/Geant4Wrap.h:132:18
    t.method("WorkerStart", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerStart));

    DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerRunStart() (" __HERE__ ")");
    // signature to use in the veto list: void G4JLWorkerInitialization::WorkerRunStart()
    // defined in ./cpp/Geant4Wrap.h:133:18
    t.method("WorkerRunStart", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerRunStart));

    DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerRunEnd() (" __HERE__ ")");
    // signature to use in the veto list: void G4JLWorkerInitialization::WorkerRunEnd()
    // defined in ./cpp/Geant4Wrap.h:134:18
    t.method("WorkerRunEnd", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerRunEnd));

    DEBUG_MSG("Adding wrapper for void G4JLWorkerInitialization::WorkerStop() (" __HERE__ ")");
    // signature to use in the veto list: void G4JLWorkerInitialization::WorkerStop()
    // defined in ./cpp/Geant4Wrap.h:135:18
    t.method("WorkerStop", static_cast<void (G4JLWorkerInitialization::*)()  const>(&G4JLWorkerInitialization::WorkerStop));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4JLWorkerInitialization>> type_;
};
std::shared_ptr<Wrapper> newJlG4JLWorkerInitialization(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4JLWorkerInitialization(module));
}