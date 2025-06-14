// this file was auto-generated by wrapit v1.6.0
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {

  template<typename _T1, typename _T2>
  struct BuildParameterList<std::pair<_T1, _T2>>
  {
    typedef ParameterList<_T1, _T2> type;
  };

  template<typename _T1, typename _T2> struct IsMirroredType<std::pair<_T1, _T2>> : std::false_type { };
  template<typename _T1, typename _T2> struct DefaultConstructible<std::pair<_T1, _T2>> : std::false_type { };
}

// Class generating the wrapper for type std::pair
// signature to use in the veto file: std::pair
struct Jlstd_pair: public Wrapper {

  Jlstd_pair(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type std::pair (" __HERE__ ")");
    // defined in /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1/__utility/pair.h:66:29
    jlcxx::TypeWrapper<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>  t =  jlModule.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>("std!pair");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>>(new jlcxx::TypeWrapper<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>(jlModule, t));
    auto t132_decl_methods = [this]<typename _T1, typename _T2> (jlcxx::TypeWrapper<std::pair<_T1, _T2>> wrapped){
      auto module_ = this->module_;
    };
    t.apply<std::pair<double, bool>>(t132_decl_methods);
  }

  void add_methods() const{
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<jlcxx::Parametric<jlcxx::TypeVar<1>, jlcxx::TypeVar<2>>>> type_;
};
std::shared_ptr<Wrapper> newJlstd_pair(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new Jlstd_pair(module));
}
