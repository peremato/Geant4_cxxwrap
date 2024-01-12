// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<HepGeom::Rotate3D> : std::false_type { };
  template<> struct DefaultConstructible<HepGeom::Rotate3D> : std::false_type { };
template<> struct SuperType<HepGeom::Rotate3D> { typedef HepGeom::Transform3D type; };
}

// Class generating the wrapper for type HepGeom::Rotate3D
// signature to use in the veto file: HepGeom::Rotate3D
struct JlHepGeom_Rotate3D: public Wrapper {

  JlHepGeom_Rotate3D(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type HepGeom::Rotate3D (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:370:9
    jlcxx::TypeWrapper<HepGeom::Rotate3D>  t = jlModule.add_type<HepGeom::Rotate3D>("HepGeom!Rotate3D",
      jlcxx::julia_base_type<HepGeom::Transform3D>());
    type_ = std::unique_ptr<jlcxx::TypeWrapper<HepGeom::Rotate3D>>(new jlcxx::TypeWrapper<HepGeom::Rotate3D>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void HepGeom::Rotate3D::Rotate3D(const CLHEP::HepRotation &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:378:12
    t.constructor<const CLHEP::HepRotation &>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void HepGeom::Rotate3D::Rotate3D(double, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:386:5
    t.constructor<double, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void HepGeom::Rotate3D::Rotate3D(double, const HepGeom::Vector3D<double> &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:395:12
    t.constructor<double, const HepGeom::Vector3D<double> &>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void HepGeom::Rotate3D::Rotate3D(const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:405:12
    t.constructor<const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &>(/*finalize=*/true);
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<HepGeom::Rotate3D>> type_;
};
std::shared_ptr<Wrapper> newJlHepGeom_Rotate3D(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlHepGeom_Rotate3D(module));
}
