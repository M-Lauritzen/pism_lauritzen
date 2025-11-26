// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021, 2022, 2023 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _PISMGOWANYIELDSTRESS_H_
#define _PISMGOWANYIELDSTRESS_H_

#include "pism/basalstrength/YieldStress.hh"
#include "pism/basalstrength/MohrCoulombYieldStress.hh"
#include "pism/hydrology/Hydrology.hh"

namespace pism {

//! @brief PISM's default basal yield stress model which applies the
//! Mohr-Coulomb model of deformable, pressurized till.
class GowanYieldStress : public MohrCoulombYieldStress {
public:
  GowanYieldStress(std::shared_ptr<const Grid> g);
  virtual ~GowanYieldStress() = default;

protected:
  void restart_impl(const File &input_file, int record);
  void bootstrap_impl(const File &input_file, const YieldStressInputs &inputs);
  void init_impl(const YieldStressInputs &inputs);

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  DiagnosticList diagnostics_impl() const;

  void update_impl(const YieldStressInputs &inputs, double t, double dt);

  void finish_initialization(const YieldStressInputs &inputs);


  std::shared_ptr<array::Forcing> m_delta;

  hydrology::Hydrology *m_hydrology;  // pointer to hydrology model

  array::Scalar m_effective_pressure;   // effective pressure from hydrology
  array::Scalar m_sliding_mechanism;    // identifies which sliding regime applies
  array::Scalar m_till_cover_local;     // local sediment fraction
  array::Scalar hydro_tauc;             // yield stress computed in hydrology
  array::Scalar tauc_ratio;             // ratio between hydrology and sediment tauc
  array::Scalar m_till_phi;             // till friction angle

  double m_dx, m_dy;

private:
  void till_friction_angle(const array::Scalar &bed_topography,
                           array::Scalar &result);

  void till_friction_angle(const array::Scalar &basal_yield_stress,
                           const array::Scalar &till_water_thickness,
                           const array::Scalar &ice_thickness,
                           const array::CellType &cell_type,
                           array::Scalar &result);
  void set_default_fields();
};

} // end of namespace pism

#endif /* _PISMGOWANYIELDSTRESS_H_ */

