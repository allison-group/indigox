#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/vec_avx.hpp>
#include <indigox/utils/triple.hpp>
#include <indigox/utils/quad.hpp>

#include <vector>

#include <cmath>
#include <iomanip>

#ifndef INDIGOX_ALGORITHM_PERCEPTION_HPP
#define INDIGOX_ALGORITHM_PERCEPTION_HPP

namespace indigox::algorithm {

  struct Perceptatron {
    enum class Settings {
      CreateMissing,
      BoolCount,
      IntCount,
      MinimumDistance,
      RadiusTolerance,
      RealCount
    };
    
    DEFAULT_SETTINGS() {
      ResetSettings();
      SetBool(Settings::CreateMissing, true);
      
      SetReal(Settings::RadiusTolerance, 0.040);
      SetReal(Settings::MinimumDistance, 0.050);
    }
    
    using bond_t = std::pair<Atom, Atom>;
    using angle_t = stdx::triple<Atom>;
    using dihedral_t = stdx::quad<Atom>;
    
    Perceptatron() {
      DefaultSettings();
    }
    
    std::vector<bond_t> PerceiveBonds(Molecule& mol) {
      std::vector<double> radii, distance_check, max_distance;
      AtomicCoordinates& coords = mol.GetAtomicCoordinates();
      int64_t sz = mol.NumAtoms();
      radii.resize(sz + 4, 0.);
      distance_check.resize(sz * sz + 4, 0.);
      max_distance.resize(sz * sz + 4, 0.);
      for (Atom atm : mol.GetAtoms()) radii[atm.GetIndex()] = atm.GetElement().GetCovalentRadius();
      
      // Calculate all distances between atoms and compare with radii sum + tolerance
      Vec4d tol(GetReal(Settings::RadiusTolerance));
      Vec4d zero = Vec4d::ZeroVector();
      for (uint32_t i = 0, k = 0; i < mol.NumAtoms(); ++i) {
        Vec4d xa(coords[i].x);
        Vec4d ya(coords[i].y);
        Vec4d za(coords[i].z);
        Vec4d ra(radii[i]);
        for (uint32_t j = 0; j < mol.NumAtoms(); j += 4) {
          if (j + 3 < i) {
            zero.save(&distance_check[k]);
            zero.save(&max_distance[k]);
          } else {
            Vec4d x(&coords.x_vals()[j]);
            Vec4d y(&coords.y_vals()[j]);
            Vec4d z(&coords.z_vals()[j]);
            Vec4d r(&radii[j]);
            x -= xa;
            y -= ya;
            z -= za;
            r += ra + tol;
            x *= x;
            r *= r;
            x += y * y;
            x += z * z;
            x.save(&distance_check[k]);
            r.save(&max_distance[k]);
          }
          k += (j + 4 < mol.NumAtoms()) ? 4 : mol.NumAtoms() - j;
        }
      }
      
      std::vector<bond_t> bonds;
      double min_dist = GetReal(Settings::MinimumDistance);
      min_dist *= min_dist;
      // A bond is found if distance_check is negative
      for (uint32_t i = 0; i < mol.NumAtoms(); ++i) {
        for (uint32_t j = i + 1; j < mol.NumAtoms(); ++j) {
          uint32_t k = i * mol.NumAtoms() + j;
          if (distance_check[k] < max_distance[k])
            bonds.emplace_back(mol.GetAtoms()[i], mol.GetAtoms()[j]);
          if (distance_check[k] <= min_dist)
            throw std::runtime_error("Atoms too close together");
        }
      }
      
      if (GetBool(Settings::CreateMissing)) {
        for (bond_t bnd : bonds) {
          if (!mol.HasBond(bnd.first, bnd.second)) mol.NewBond(bnd.first, bnd.second);
        }
      }
      return bonds;
    }
    
  };

} // namespace indigox::algorithm

#endif /* INDIGOX_ALGORITHM_PERCEPTION_HPP */
