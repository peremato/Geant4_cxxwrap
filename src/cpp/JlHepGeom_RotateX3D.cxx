// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<HepGeom::RotateX3D> : std::false_type { };
  template<> struct DefaultConstructible<HepGeom::RotateX3D> : std::false_type { };
template<> struct SuperType<HepGeom::RotateX3D> { typedef HepGeom::Rotate3D type; };
}

// Class generating the wrapper for type HepGeom::RotateX3D
// signature to use in the veto file: HepGeom::RotateX3D
struct JlHepGeom_RotateX3D: public Wrapper {

  JlHepGeom_RotateX3D(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type HepGeom::RotateX3D (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:425:9
    jlcxx::TypeWrapper<HepGeom::RotateX3D>  t = jlModule.add_type<HepGeom::RotateX3D>("HepGeom!RotateX3D",
      jlcxx::julia_base_type<HepGeom::Rotate3D>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<HepGeom::RotateX3D>>(new jlcxx::TypeWrapper<HepGeom::RotateX3D>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void HepGeom::RotateX3D::RotateX3D(double) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:433:5
    t.constructor<double>(/*finalize=*/true);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<HepGeom::RotateX3D>> type_;
};
std::shared_ptr<Wrapper> newJlHepGeom_RotateX3D(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlHepGeom_RotateX3D(module));
}