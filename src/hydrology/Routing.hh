// Copyright (C) 2012-2019, 2021, 2022 PISM Authors
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

#ifndef _ROUTING_H_
#define _ROUTING_H_

#include "Hydrology.hh"
#include "pism/util/array/Staggered.hh"

namespace pism {

namespace hydrology {

//! \brief A subglacial hydrology model which assumes water pressure
//! equals overburden pressure.
/*!
  This PISM hydrology model has lateral motion of subglacial water and which
  conserves the water mass.  Further documentation is in [\ref BuelervanPeltDRAFT].

  The water velocity is along the steepest descent route for the hydraulic
  potential.  This potential is (mostly) a function of ice sheet geometry,
  because the water pressure is set to the overburden pressure, a simplified but
  well-established model [\ref Shreve1972].  However, the water layer thickness
  is also a part of the hydraulic potential because it is actually the potential
  of the *top* of the water layer.

  This (essential) model has been used for finding locations of subglacial lakes
  [\ref Siegertetal2007, \ref Livingstoneetal2013].  Subglacial lakes occur
  at local minima of the hydraulic potential.  If water builds up significantly
  (e.g. thickness of 10s of meters or more) then in the model here the resulting
  lakes diffuse instead of becoming infinitely deep.  Thus we avoid delta
  functions of water thickness at the minima of the hydraulic potential in this
  well-posed model.

  This model should generally be tested using static ice geometry first.

  The state space includes both the till water effective thickness \f$W_{till}\f$,
  which is in Hydrology, and the transportable water layer thickness \f$W\f$.

  For more complete modeling where the water pressure is determined by a
  physical model for the opening and closing of cavities, and where the state
  space includes a nontrivial pressure variable, see hydrology::Distributed.

  There is an option `-hydrology_null_strip` `X` which produces a strip of
  `X` km around the edge of the computational domain.  In that strip the water flow
  velocity is set to zero.  The water amount is also reset to zero at the end
  of each time step in this strip (in an accounted way).

  As noted this is the minimal model which has a lateral water flux.  This flux is
  \f[ \mathbf{q} = - K \nabla \psi = \mathbf{V} W - D \nabla W \f]
  where \f$\psi\f$ is the hydraulic potential
  \f[ \psi = P + \rho_w g (b + W). \f]
  The generalized conductivity \f$K\f$ is nontrivial and it generally also
  depends on the water thickness:
  \f[ K = k W^{\alpha-1} |\nabla (P+\rho_w g b)|^{\beta-2}. \f]

  This model contains enough information (enough modeled fields) so that we can compute
  the wall melt generated by dissipating the gravitational potential energy in the moving,
  presumably turbulent, subglacial water. If we suppose that this heat is dissipated
  immediately as melt on the cavity/conduit walls then we get a formula for a wall melt
  contribution. (This is in addition to the `basal_melt_rate_grounded` field coming from
  conserving energy in the flowing ice.) See wall_melt(). At this time the wall melt is
  diagnostic only and does not add to the water amount W; such an addition is generally
  unstable.
*/
class Routing : public Hydrology {
public:

  Routing(std::shared_ptr<const IceGrid> g);
  virtual ~Routing() = default;

  const array::Scalar& subglacial_water_pressure() const;
  const array::Staggered& velocity_staggered() const;

protected:
  virtual void restart_impl(const File &input_file, int record);

  virtual void bootstrap_impl(const File &input_file,
                              const array::Scalar &ice_thickness);

  virtual void init_impl(const array::Scalar &W_till,
                               const array::Scalar &W,
                               const array::Scalar &P);

  virtual void update_impl(double t, double dt, const Inputs& inputs);

  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;
  virtual std::map<std::string, TSDiagnostic::Ptr> ts_diagnostics_impl() const;

  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  double max_timestep_W_diff(double KW_max) const;
  double max_timestep_W_cfl() const;
protected:

  // edge-centered (staggered) advection flux
  array::Staggered1 m_Qstag;

  array::Staggered1 m_Qstag_average;

  // edge-centered (staggered) water velocity
  array::Staggered m_Vstag;

  // edge-centered (staggered) W values (averaged from regular)
  array::Staggered1 m_Wstag;

  // edge-centered (staggered) values of nonlinear conductivity
  array::Staggered1 m_Kstag;

  // work space
  array::Scalar m_Wnew, m_Wtillnew;

  // ghosted temporary storage; modified in compute_conductivity and compute_velocity
  mutable array::Scalar1 m_R;

  double m_dx, m_dy;
  double m_rg;

  array::Scalar1 m_bottom_surface;

  void water_thickness_staggered(const array::Scalar &W,
                                 const array::CellType1 &mask,
                                 array::Staggered &result);

  void compute_conductivity(const array::Staggered &W,
                            const array::Scalar &P,
                            const array::Scalar &bed,
                            array::Staggered &result,
                            double &maxKW) const;

  void compute_velocity(const array::Staggered &W,
                        const array::Scalar &P,
                        const array::Scalar &bed,
                        const array::Staggered &K,
                        const array::Scalar1 *no_model_mask,
                        array::Staggered &result) const;

  void advective_fluxes(const array::Staggered &V,
                        const array::Scalar &W,
                        array::Staggered &result) const;

  void W_change_due_to_flow(double dt,
                            const array::Scalar1    &W,
                            const array::Staggered1 &Wstag,
                            const array::Staggered1 &K,
                            const array::Staggered1 &Q,
                            array::Scalar &result);
  void update_W(double dt,
                const array::Scalar     &surface_input_rate,
                const array::Scalar     &basal_melt_rate,
                const array::Scalar1    &W,
                const array::Staggered1 &Wstag,
                const array::Scalar     &Wtill,
                const array::Scalar     &Wtill_new,
                const array::Staggered1 &K,
                const array::Staggered1 &Q,
                array::Scalar           &W_new);

  void update_Wtill(double dt,
                    const array::Scalar &Wtill,
                    const array::Scalar &surface_input_rate,
                    const array::Scalar &basal_melt_rate,
                    array::Scalar &Wtill_new);

private:
  virtual void initialization_message() const;
};

void wall_melt(const Routing &model,
               const array::Scalar &bed_elevation,
               array::Scalar &result);

} // end of namespace hydrology
} // end of namespace pism

#endif /* _ROUTING_H_ */
