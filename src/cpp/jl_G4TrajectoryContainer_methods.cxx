// this file was auto-generated by wrapit v0.1.0-54-g4322429
#include <type_traits>
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/functions.hpp"
#include "jlcxx/stl.hpp"

#include "cpp/jlGeant4.h"


#ifdef VERBOSE_IMPORT
#  define DEBUG_MSG(a) std::cerr << a << "\n"
#else
#  define DEBUG_MSG(a)
#endif
#define __HERE__  __FILE__ ":" QUOTE2(__LINE__)
#define QUOTE(arg) #arg
#define QUOTE2(arg) QUOTE(arg)
void add_methods_for_G4TrajectoryContainer(jlcxx::Module& types, jlcxx::TypeWrapper<G4TrajectoryContainer>& t130) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4TrajectoryContainer
   */
  types.set_override_module(jl_base_module);

  DEBUG_MSG("Adding wrapper for G4bool G4TrajectoryContainer::operator==(const G4TrajectoryContainer &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4TrajectoryContainer::operator==(const G4TrajectoryContainer &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TrajectoryContainer.hh:60:12
  t130.method("==", static_cast<G4bool (G4TrajectoryContainer::*)(const G4TrajectoryContainer &)  const>(&G4TrajectoryContainer::operator==));

  DEBUG_MSG("Adding wrapper for G4bool G4TrajectoryContainer::operator!=(const G4TrajectoryContainer &) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4TrajectoryContainer::operator!=(const G4TrajectoryContainer &)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TrajectoryContainer.hh:61:12
  t130.method("!=", static_cast<G4bool (G4TrajectoryContainer::*)(const G4TrajectoryContainer &)  const>(&G4TrajectoryContainer::operator!=));

  types.unset_override_module();

  DEBUG_MSG("Adding wrapper for size_t G4TrajectoryContainer::size() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4TrajectoryContainer::size()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TrajectoryContainer.hh:63:24
  t130.method("size", static_cast<size_t (G4TrajectoryContainer::*)()  const>(&G4TrajectoryContainer::size));

  DEBUG_MSG("Adding wrapper for void G4TrajectoryContainer::push_back(G4VTrajectory *) (" __HERE__ ")");
  // signature to use in the veto list: void G4TrajectoryContainer::push_back(G4VTrajectory *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TrajectoryContainer.hh:64:17
  t130.method("push_back", static_cast<void (G4TrajectoryContainer::*)(G4VTrajectory *) >(&G4TrajectoryContainer::push_back));

  DEBUG_MSG("Adding wrapper for size_t G4TrajectoryContainer::entries() (" __HERE__ ")");
  // signature to use in the veto list: size_t G4TrajectoryContainer::entries()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TrajectoryContainer.hh:65:24
  t130.method("entries", static_cast<size_t (G4TrajectoryContainer::*)()  const>(&G4TrajectoryContainer::entries));

  DEBUG_MSG("Adding wrapper for G4bool G4TrajectoryContainer::insert(G4VTrajectory *) (" __HERE__ ")");
  // signature to use in the veto list: G4bool G4TrajectoryContainer::insert(G4VTrajectory *)
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TrajectoryContainer.hh:66:19
  t130.method("insert", static_cast<G4bool (G4TrajectoryContainer::*)(G4VTrajectory *) >(&G4TrajectoryContainer::insert));

  DEBUG_MSG("Adding wrapper for void G4TrajectoryContainer::clearAndDestroy() (" __HERE__ ")");
  // signature to use in the veto list: void G4TrajectoryContainer::clearAndDestroy()
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TrajectoryContainer.hh:67:17
  t130.method("clearAndDestroy", static_cast<void (G4TrajectoryContainer::*)() >(&G4TrajectoryContainer::clearAndDestroy));
  types.set_override_module(jl_base_module);


  DEBUG_MSG("Adding getindex method to wrap G4VTrajectory * G4TrajectoryContainer::operator[](size_t) (" __HERE__ ")");
  // defined in /Users/mato/.julia/artifacts/9d4b417a98ec8f720b8871baefe87108f864656f/include/Geant4/G4TrajectoryContainer.hh:72:27
  t130.method("getindex",
    [](G4TrajectoryContainer& a, size_t i){
    return a[i];
  });

  /* End of G4TrajectoryContainer class method wrappers
   **********************************************************************/

}
