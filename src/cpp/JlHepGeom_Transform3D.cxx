// this file was auto-generated by wrapit 
#include "Wrapper.h"

#include "jlGeant4.h"
#include "dbg_msg.h"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

namespace jlcxx {
  template<> struct IsMirroredType<HepGeom::Transform3D> : std::false_type { };
  template<> struct DefaultConstructible<HepGeom::Transform3D> : std::false_type { };
}

// Class generating the wrapper for type HepGeom::Transform3D
// signature to use in the veto file: HepGeom::Transform3D
struct JlHepGeom_Transform3D: public Wrapper {

  JlHepGeom_Transform3D(jlcxx::Module& jlModule): Wrapper(jlModule){
    DEBUG_MSG("Adding wrapper for type HepGeom::Transform3D (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:170:9
    jlcxx::TypeWrapper<HepGeom::Transform3D>  t = jlModule.add_type<HepGeom::Transform3D>("HepGeom!Transform3D");
    type_ = std::unique_ptr<jlcxx::TypeWrapper<HepGeom::Transform3D>>(new jlcxx::TypeWrapper<HepGeom::Transform3D>(jlModule, t));
  }

  void add_methods() const{
    auto& t = *type_;
    t.template constructor<>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void HepGeom::Transform3D::Transform3D(const CLHEP::HepRotation &, const CLHEP::Hep3Vector &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:217:12
    t.constructor<const CLHEP::HepRotation &, const CLHEP::Hep3Vector &>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void HepGeom::Transform3D::Transform3D(const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:221:5
    t.constructor<const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &, const HepGeom::Point3D<double> &>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for void HepGeom::Transform3D::Transform3D(const HepGeom::Transform3D &) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:230:5
    t.constructor<const HepGeom::Transform3D &>(/*finalize=*/true);


    DEBUG_MSG("Adding wrapper for HepGeom::Transform3D & HepGeom::Transform3D::operator=(const HepGeom::Transform3D &) (" __HERE__ ")");
    // signature to use in the veto list: HepGeom::Transform3D & HepGeom::Transform3D::operator=(const HepGeom::Transform3D &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:242:19
    t.method("assign", static_cast<HepGeom::Transform3D & (HepGeom::Transform3D::*)(const HepGeom::Transform3D &) >(&HepGeom::Transform3D::operator=));

    module_.set_override_module(jl_base_module);


    DEBUG_MSG("Adding getindex method to wrap const HepGeom::Transform3D::Transform3D_row HepGeom::Transform3D::operator[](int) (" __HERE__ ")");
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:250:34
    t.method("getindex",
      [](HepGeom::Transform3D& a, int i){
      return a[i];
    });

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::operator()(int, int) (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::operator()(int, int)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:253:12
    t.method("paren", static_cast<double (HepGeom::Transform3D::*)(int, int)  const>(&HepGeom::Transform3D::operator()));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::xx() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::xx()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:257:12
    t.method("xx", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::xx));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::xy() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::xy()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:260:12
    t.method("xy", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::xy));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::xz() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::xz()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:263:12
    t.method("xz", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::xz));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::yx() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::yx()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:266:12
    t.method("yx", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::yx));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::yy() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::yy()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:269:12
    t.method("yy", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::yy));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::yz() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::yz()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:272:12
    t.method("yz", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::yz));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::zx() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::zx()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:275:12
    t.method("zx", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::zx));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::zy() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::zy()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:278:12
    t.method("zy", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::zy));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::zz() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::zz()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:281:12
    t.method("zz", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::zz));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::dx() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::dx()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:284:12
    t.method("dx", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::dx));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::dy() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::dy()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:287:12
    t.method("dy", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::dy));

    DEBUG_MSG("Adding wrapper for double HepGeom::Transform3D::dz() (" __HERE__ ")");
    // signature to use in the veto list: double HepGeom::Transform3D::dz()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:290:12
    t.method("dz", static_cast<double (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::dz));

    DEBUG_MSG("Adding wrapper for void HepGeom::Transform3D::setIdentity() (" __HERE__ ")");
    // signature to use in the veto list: void HepGeom::Transform3D::setIdentity()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:294:10
    t.method("setIdentity", static_cast<void (HepGeom::Transform3D::*)() >(&HepGeom::Transform3D::setIdentity));

    DEBUG_MSG("Adding wrapper for HepGeom::Transform3D HepGeom::Transform3D::inverse() (" __HERE__ ")");
    // signature to use in the veto list: HepGeom::Transform3D HepGeom::Transform3D::inverse()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:300:17
    t.method("inverse", static_cast<HepGeom::Transform3D (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::inverse));
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for HepGeom::Transform3D HepGeom::Transform3D::operator*(const HepGeom::Transform3D &) (" __HERE__ ")");
    // signature to use in the veto list: HepGeom::Transform3D HepGeom::Transform3D::operator*(const HepGeom::Transform3D &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:304:17
    t.method("*", static_cast<HepGeom::Transform3D (HepGeom::Transform3D::*)(const HepGeom::Transform3D &)  const>(&HepGeom::Transform3D::operator*));

    module_.unset_override_module();

    DEBUG_MSG("Adding wrapper for void HepGeom::Transform3D::getDecomposition(HepGeom::Scale3D &, HepGeom::Rotate3D &, HepGeom::Translate3D &) (" __HERE__ ")");
    // signature to use in the veto list: void HepGeom::Transform3D::getDecomposition(HepGeom::Scale3D &, HepGeom::Rotate3D &, HepGeom::Translate3D &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:321:10
    t.method("getDecomposition", static_cast<void (HepGeom::Transform3D::*)(HepGeom::Scale3D &, HepGeom::Rotate3D &, HepGeom::Translate3D &)  const>(&HepGeom::Transform3D::getDecomposition));

    DEBUG_MSG("Adding wrapper for bool HepGeom::Transform3D::isNear(const HepGeom::Transform3D &, double) (" __HERE__ ")");
    // signature to use in the veto list: bool HepGeom::Transform3D::isNear(const HepGeom::Transform3D &, double)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:329:10
    t.method("isNear", static_cast<bool (HepGeom::Transform3D::*)(const HepGeom::Transform3D &, double)  const>(&HepGeom::Transform3D::isNear));
    t.method("isNear", [](HepGeom::Transform3D const& a, const HepGeom::Transform3D & arg0)->bool { return a.isNear(arg0); });
    t.method("isNear", [](HepGeom::Transform3D const* a, const HepGeom::Transform3D & arg0)->bool { return a->isNear(arg0); });

    DEBUG_MSG("Adding wrapper for CLHEP::HepRotation HepGeom::Transform3D::getRotation() (" __HERE__ ")");
    // signature to use in the veto list: CLHEP::HepRotation HepGeom::Transform3D::getRotation()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:335:31
    t.method("getRotation", static_cast<CLHEP::HepRotation (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::getRotation));

    DEBUG_MSG("Adding wrapper for CLHEP::Hep3Vector HepGeom::Transform3D::getTranslation() (" __HERE__ ")");
    // signature to use in the veto list: CLHEP::Hep3Vector HepGeom::Transform3D::getTranslation()
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:341:30
    t.method("getTranslation", static_cast<CLHEP::Hep3Vector (HepGeom::Transform3D::*)()  const>(&HepGeom::Transform3D::getTranslation));
    module_.set_override_module(jl_base_module);

    DEBUG_MSG("Adding wrapper for bool HepGeom::Transform3D::operator==(const HepGeom::Transform3D &) (" __HERE__ ")");
    // signature to use in the veto list: bool HepGeom::Transform3D::operator==(const HepGeom::Transform3D &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:345:10
    t.method("==", static_cast<bool (HepGeom::Transform3D::*)(const HepGeom::Transform3D &)  const>(&HepGeom::Transform3D::operator==));

    DEBUG_MSG("Adding wrapper for bool HepGeom::Transform3D::operator!=(const HepGeom::Transform3D &) (" __HERE__ ")");
    // signature to use in the veto list: bool HepGeom::Transform3D::operator!=(const HepGeom::Transform3D &)
    // defined in /Users/mato/.julia/artifacts/4afb5743b029965f72ec5a970d92d5344ce830d2/include/Geant4/CLHEP/Geometry/Transform3D.h:349:10
    t.method("!=", static_cast<bool (HepGeom::Transform3D::*)(const HepGeom::Transform3D &)  const>(&HepGeom::Transform3D::operator!=));

    module_.unset_override_module();
  }

private:
  std::unique_ptr<jlcxx::TypeWrapper<HepGeom::Transform3D>> type_;
};
std::shared_ptr<Wrapper> newJlHepGeom_Transform3D(jlcxx::Module& module){
  return std::shared_ptr<Wrapper>(new JlHepGeom_Transform3D(module));
}
