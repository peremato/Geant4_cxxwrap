// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VFastSimulationModel> : std::false_type { };
  template<> struct DefaultConstructible<G4VFastSimulationModel> : std::false_type { };
}

// Class generating the wrapper for type G4VFastSimulationModel
// signature to use in the veto file: G4VFastSimulationModel
struct JlG4VFastSimulationModel: public Wrapper {

  JlG4VFastSimulationModel(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VFastSimulationModel (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VFastSimulationModel.hh:59:7
    jlcxx::TypeWrapper<G4VFastSimulationModel>  t = jlModule.add_type<G4VFastSimulationModel>("G4VFastSimulationModel");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VFastSimulationModel>>(new jlcxx::TypeWrapper<G4VFastSimulationModel>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VFastSimulationModel>> type_;
};
std::shared_ptr<Wrapper> newJlG4VFastSimulationModel(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VFastSimulationModel(module));
}
