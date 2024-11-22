// this file was auto-generated by wrapit 5168a24-dirty
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {

  template<typename K, typename V>
  struct BuildParameterList<G4ParticleTableIterator<K, V>>
  {
    typedef ParameterList<K, V> type;
  };

  template<typename K, typename V> struct IsMirroredType<G4ParticleTableIterator<K, V>> : std::false_type { };
  template<typename K, typename V> struct DefaultConstructible<G4ParticleTableIterator<K, V>> : std::false_type { };
}

// Class generating the wrapper for type G4ParticleTableIterator
// signature to use in the veto file: G4ParticleTableIterator
struct JlG4ParticleTableIterator: public Wrapper {

  JlG4ParticleTableIterator(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type G4ParticleTableIterator (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/04a1f392c53fa9913a6e32dc79e45dcf6f1dd250/include/Geant4/G4ParticleTableIterator.hh:39:7
    jlcxx::TypeWrapper<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>  t =  jlModule.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>("G4ParticleTableIterator");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>>(new jlcxx::TypeWrapper<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>(jlModule, t));
    auto t230_decl_methods = [this]<typename K, typename V> (jlcxx::TypeWrapper<G4ParticleTableIterator<K, V>> wrapped){
      auto module_ = this->module_;
    };
    t.apply<G4ParticleTableIterator<G4String, G4ParticleDefinition *>>(t230_decl_methods);
  }

  void add_methods() const{
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>> type_;
};
std::shared_ptr<Wrapper> newJlG4ParticleTableIterator(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlG4ParticleTableIterator(module));
}
