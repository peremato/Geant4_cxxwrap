// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4ICRU90StoppingData> : std::false_type { };
  template<> struct DefaultConstructible<G4ICRU90StoppingData> : std::false_type { };
}

// Class generating the wrapper for type G4ICRU90StoppingData
// signature to use in the veto file: G4ICRU90StoppingData
struct JlG4ICRU90StoppingData: public Wrapper {

  JlG4ICRU90StoppingData(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4ICRU90StoppingData (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4ICRU90StoppingData.hh:54:7
    jlcxx::TypeWrapper<G4ICRU90StoppingData>  t = jlModule.add_type<G4ICRU90StoppingData>("G4ICRU90StoppingData");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4ICRU90StoppingData>>(new jlcxx::TypeWrapper<G4ICRU90StoppingData>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4ICRU90StoppingData>> type_;
};
std::shared_ptr<Wrapper> newJlG4ICRU90StoppingData(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4ICRU90StoppingData(module));
}
