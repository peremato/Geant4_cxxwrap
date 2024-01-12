// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4VExceptionHandler> : std::false_type { };
  template<> struct DefaultConstructible<G4VExceptionHandler> : std::false_type { };
}

// Class generating the wrapper for type G4VExceptionHandler
// signature to use in the veto file: G4VExceptionHandler
struct JlG4VExceptionHandler: public Wrapper {

  JlG4VExceptionHandler(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4VExceptionHandler (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/G4VExceptionHandler.hh:43:7
    jlcxx::TypeWrapper<G4VExceptionHandler>  t = jlModule.add_type<G4VExceptionHandler>("G4VExceptionHandler");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4VExceptionHandler>>(new jlcxx::TypeWrapper<G4VExceptionHandler>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4VExceptionHandler>> type_;
};
std::shared_ptr<Wrapper> newJlG4VExceptionHandler(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4VExceptionHandler(module));
}
