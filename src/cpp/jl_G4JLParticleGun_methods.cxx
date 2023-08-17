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
void add_methods_for_G4JLParticleGun(jlcxx::Module& types, jlcxx::TypeWrapper<G4JLParticleGun>& t108) {


  /**********************************************************************/
  /* Wrappers for the methods of class G4JLParticleGun
   */

  DEBUG_MSG("Adding wrapper for G4ParticleGun * G4JLParticleGun::GetGun() (" __HERE__ ")");
  // signature to use in the veto list: G4ParticleGun * G4JLParticleGun::GetGun()
  // defined in ./cpp/Geant4Wrap.h:116:18
  t108.method("GetGun", static_cast<G4ParticleGun * (G4JLParticleGun::*)()  const>(&G4JLParticleGun::GetGun));

  DEBUG_MSG("Adding wrapper for void G4JLParticleGun::GeneratePrimaries(G4Event *) (" __HERE__ ")");
  // signature to use in the veto list: void G4JLParticleGun::GeneratePrimaries(G4Event *)
  // defined in ./cpp/Geant4Wrap.h:117:8
  t108.method("GeneratePrimaries", static_cast<void (G4JLParticleGun::*)(G4Event *) >(&G4JLParticleGun::GeneratePrimaries));

  /* End of G4JLParticleGun class method wrappers
   **********************************************************************/

}
