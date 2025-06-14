// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<G4JLStateDependent> : std::false_type { };
  template<> struct DefaultConstructible<G4JLStateDependent> : std::false_type { };
template<> struct SuperType<G4JLStateDependent> { typedef G4VStateDependent type; };
}

// Class generating the wrapper for type G4JLStateDependent
// signature to use in the veto file: G4JLStateDependent
struct JlG4JLStateDependent: public Wrapper {

  JlG4JLStateDependent(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4JLStateDependent (" __HERE__ ")");
    // defined in cpp/Geant4Wrap.h:267:7
    jlcxx::TypeWrapper<G4JLStateDependent>  t = jlModule.add_type<G4JLStateDependent>("G4JLStateDependent",
      jlcxx::julia_base_type<G4VStateDependent>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<G4JLStateDependent>>(new jlcxx::TypeWrapper<G4JLStateDependent>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;


    DEBUG_MSG("Adding wrapper for void G4JLStateDependent::G4JLStateDependent(notify_f, void *) (" __HERE__ ")");
    // defined in cpp/Geant4Wrap.h:269:5
    t.constructor<notify_f, void *>(/*finalize=*/jlcxx::finalize_policy::no, jlcxx::arg("this"), jlcxx::arg("notify"), jlcxx::arg("data")    );

    DEBUG_MSG("Adding wrapper for G4bool G4JLStateDependent::Notify(G4ApplicationState) (" __HERE__ ")");
    // signature to use in the veto list: G4bool G4JLStateDependent::Notify(G4ApplicationState)
    // defined in cpp/Geant4Wrap.h:270:12
    t.method("Notify", [](G4JLStateDependent& a, G4ApplicationState arg0)->G4bool { return a.Notify(arg0); }, jlcxx::arg("this"), jlcxx::arg("to"));
    t.method("Notify", [](G4JLStateDependent* a, G4ApplicationState arg0)->G4bool { return a->Notify(arg0); }, jlcxx::arg("this"), jlcxx::arg("to"));
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<G4JLStateDependent>> type_;
};
std::shared_ptr<Wrapper> newJlG4JLStateDependent(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4JLStateDependent(module));
}
