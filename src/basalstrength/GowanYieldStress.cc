// Copyright (C) 2004--2023, 2025 PISM Authors
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

#include <algorithm>
#include <cmath>

#include "pism/basalstrength/GowanYieldStress.hh"
#include "pism/basalstrength/MohrCoulombYieldStress.hh"
#include "pism/basalstrength/MohrCoulombPointwise.hh"

#include "pism/util/Grid.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/coupler/util/options.hh" // ForcingOptions

namespace pism {

//! \file GowanYieldStress.cc  Process model which computes pseudo-plastic yield stress for the subglacial layer.

/*! \file GowanYieldStress.cc
The output variable of this submodel is `tauc`, the pseudo-plastic yield stress
field that is used in the ShallowStressBalance objects.  This quantity is
computed by the Mohr-Coulomb criterion [\ref SchoofTill], but using an empirical
relation between the amount of water in the till and the effective pressure
of the overlying glacier resting on the till [\ref Tulaczyketal2000].

The "dry" strength of the till is a state variable which is private to
the submodel, namely `tillphi`.  Its initialization is nontrivial: either the
`-topg_to_phi`  heuristic is used or inverse modeling can be used.  (In the
latter case `tillphi` can be read-in at the beginning of the run.  Currently
`tillphi` does not evolve during the run.)

The effective pressure is derived from the till (pore) water amount (the effective water
layer thickness). Then the effective pressure is combined with tillphi to compute an
updated `tauc` by the Mohr-Coulomb criterion.

This submodel is inactive in floating areas.
*/


/*!
The pseudo-plastic till basal resistance model is governed by this power law
equation,
    @f[ \tau_b = - \frac{\tau_c}{|\mathbf{U}|^{1-q} U_{\mathtt{th}}^q} \mathbf{U}, @f]
where @f$\tau_b=(\tau_{(b)x},\tau_{(b)y})@f$ is the basal shear stress and
@f$U=(u,v)@f$ is the sliding velocity.

We call the scalar field @f$\tau_c(t,x,y)@f$ the *yield stress* even when
the power @f$q@f$ is not zero; when that power is zero the formula describes
a plastic material with an actual yield stress.  The constant
@f$U_{\mathtt{th}}@f$ is the *threshold speed*, and @f$q@f$ is the *pseudo*
*plasticity exponent*.  The current class computes this yield stress field.
See also IceBasalResistancePlasticLaw::drag().

The strength of the saturated till material, the yield stress, is modeled by a
Mohr-Coulomb relation [\ref Paterson, \ref SchoofStream],
    @f[   \tau_c = c_0 + (\tan \varphi) N_{till}, @f]
where @f$N_{till}@f$ is the effective pressure of the glacier on the mineral
till.

The determination of the till friction angle @f$\varphi(x,y)@f$  is important.
It is assumed in this default model to be a time-independent factor which
describes the strength of the unsaturated "dry" (mineral) till material.  Thus
it is assumed to change more slowly than the till water pressure, and it follows
that it changes more slowly than the yield stress and the basal shear stress.

Option `-topg_to_phi` causes call to topg_to_phi() at the beginning of the run.
This determines the map of @f$\varphi(x,y)@f$.  If this option is note given,
the current method leaves `tillphi` unchanged, and thus either in its
read-in-from-file state or with a default constant value from the config file.
*/
GowanYieldStress::GowanYieldStress(std::shared_ptr<const Grid> grid)
  : MohrCoulombYieldStress(grid),
    m_effective_pressure(grid, "effective_pressure"),
    m_sliding_mechanism(grid, "sliding_mechanism"),
    m_till_cover_local(grid, "till_cover_local"),
    hydro_tauc(grid, "hydro_tauc"),
    tauc_ratio(grid, "tauc_ratio"),
    m_till_phi(grid, "till_phi"),
    m_dx(grid->dx()),
    m_dy(grid->dy())
{

  m_name = "Gowan yield stress model";

  m_till_phi.metadata()
      .long_name("friction angle for till under grounded ice sheet")
      .units("degrees")
      .set_time_independent(true);

  // This model requires till to be present
  {
    std::string hydrology_tillwat_max = "hydrology.tillwat_max";
    bool till_is_present = m_config->get_number(hydrology_tillwat_max) > 0.0;

    if (not till_is_present) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "The Mohr-Coulomb yield stress model cannot be used without till.\n"
                                    "Reset %s or choose a different yield stress model.",
                                    hydrology_tillwat_max.c_str());
    }
  }

  // --- Gowan's additional diagnostic/state fields ---

  m_effective_pressure.metadata()
      .long_name("effective pressure in drainage system")
      .units("Pa");

  m_sliding_mechanism.metadata()
      .long_name("sliding mechanism flag")
      .units("1");

  m_till_cover_local.metadata()
      .long_name("local till cover fraction")
      .units("1")
      .set_time_independent(true);

  hydro_tauc.metadata()
      .long_name("yield stress from hydrology module")
      .units("Pa");

  tauc_ratio.metadata()
      .long_name("ratio of hydrology and sediment yield stress")
      .units("1");
}


void GowanYieldStress::restart_impl(const File &input_file, int record) {
  // read common fields from the base model
  m_basal_yield_stress.read(input_file, record);
  m_till_phi.read(input_file, record);

  // read Gowan-specific diagnostic fields (if present)
  m_effective_pressure.read(input_file, record);
  m_sliding_mechanism.read(input_file, record);
  m_till_cover_local.read(input_file, record);
  hydro_tauc.read(input_file, record);
  tauc_ratio.read(input_file, record);
}

void GowanYieldStress::set_default_fields() {
  // Set default values for fields
  m_effective_pressure.set(0.0);
  m_sliding_mechanism.set(0.0);
  m_till_cover_local.set(m_config->get_number("basal_yield_stress.gowan.till_fraction_coverage")); // sets the fraction of till cover uniformly
  hydro_tauc.set(0.0);
  tauc_ratio.set(0.0);
}

//! Initialize the pseudo-plastic till mechanical model.
void GowanYieldStress::bootstrap_impl(const File &input_file,
                                                const YieldStressInputs &inputs) {

  m_log->message(2, "* Bootstrapping GowanYieldStress ...\n");

  const double high_tauc = m_config->get_number("basal_yield_stress.ice_free_bedrock");

  // Initialize the basal yield stress field
  m_basal_yield_stress.regrid(input_file, io::Default(high_tauc));

  double till_phi_default = m_config->get_number("basal_yield_stress.mohr_coulomb.till_phi_default");
  m_till_phi.regrid(input_file, io::Default(till_phi_default));

  // Initialize extra diagnostic fields with default or neutral values.
  set_default_fields();

  finish_initialization(inputs);
}


void GowanYieldStress::init_impl(const YieldStressInputs &inputs) {
  m_log->message(2, "* Initializing GowanYieldStress ...\n");

  MohrCoulombYieldStress::init_impl(inputs);

  const double high_tauc = m_config->get_number("basal_yield_stress.ice_free_bedrock");
  m_basal_yield_stress.set(high_tauc);

  set_default_fields();

  finish_initialization(inputs);
}

/*!
 * Finish initialization after bootstrapping or initializing using constants.
 */
void GowanYieldStress::finish_initialization(const YieldStressInputs &inputs) {  
  m_log->message(2, "* Finishing initialization for GowanYieldStress ...\n");

  // Default initialization for fields
  set_default_fields();

  // Ensure all dependent quantities are consistent before time stepping
  this->update(inputs, time().current(), 1.0);
}

void GowanYieldStress::define_model_state_impl(const File &output) const {
  // Base model fields
  m_basal_yield_stress.define(output, io::PISM_DOUBLE);
  m_till_phi.define(output, io::PISM_DOUBLE);

  // diagnostic and coupling fields
  m_effective_pressure.define(output, io::PISM_DOUBLE);
  m_sliding_mechanism.define(output, io::PISM_DOUBLE);
  m_till_cover_local.define(output, io::PISM_DOUBLE);
  hydro_tauc.define(output, io::PISM_DOUBLE);
  tauc_ratio.define(output, io::PISM_DOUBLE);
}


void GowanYieldStress::write_model_state_impl(const File &output) const {
  // Base model fields
  m_basal_yield_stress.write(output);
  m_till_phi.write(output);

  // diagnostic and coupling fields
  m_effective_pressure.write(output);
  m_sliding_mechanism.write(output);
  m_till_cover_local.write(output);
  hydro_tauc.write(output);
  tauc_ratio.write(output);
}

//! Update the till yield stress for use in the pseudo-plastic till basal stress
//! model.  See also IceBasalResistancePlasticLaw.
/*!
Updates yield stress  @f$ \tau_c @f$  based on modeled till water layer thickness
from a Hydrology object.  We implement the Mohr-Coulomb criterion allowing
a (typically small) till cohesion  @f$ c_0 @f$
and by expressing the coefficient as the tangent of a till friction angle
 @f$ \varphi @f$ :
    @f[   \tau_c = c_0 + (\tan \varphi) N_{till}. @f]
See [@ref Paterson] table 8.1 regarding values.

The effective pressure on the till is empirically-related
to the amount of water in the till.  We use this formula derived from
[@ref Tulaczyketal2000] and documented in [@ref BuelervanPeltDRAFT]:

@f[ N_{till} = \min\left\{P_o, N_0 \left(\frac{\delta P_o}{N_0}\right)^s 10^{(e_0/C_c) (1 - s)}\right\} @f]

where  @f$ s = W_{till} / W_{till}^{max} @f$,  @f$ W_{till}^{max} @f$ =`hydrology_tillwat_max`,
@f$ \delta @f$ =`basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden`,  @f$ P_o @f$  is the
overburden pressure,  @f$ N_0 @f$ =`basal_yield_stress.mohr_coulomb.till_reference_effective_pressure` is a
reference effective pressure,   @f$ e_0 @f$ =`basal_yield_stress.mohr_coulomb.till_reference_void_ratio` is the void ratio
at the reference effective pressure, and  @f$ C_c @f$ =`basal_yield_stress.mohr_coulomb.till_compressibility_coefficient`
is the coefficient of compressibility of the till.  Constants  @f$ N_0, e_0, C_c @f$  are
found by [@ref Tulaczyketal2000] from laboratory experiments on samples of
till.

If `basal_yield_stress.add_transportable_water` is yes then @f$ s @f$ in the above formula
becomes @f$ s = (W + W_{till}) / W_{till}^{max} @f$,
that is, the water amount is the sum @f$ W+W_{till} @f$.
 */

void GowanYieldStress::update_impl(const YieldStressInputs &inputs,
                                             double t, double dt) {
  (void) t; (void) dt;
  
  if (!inputs.hydrology_flux || !inputs.hydrology_gradient) { // Not sure if this is sufficient
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "The Gowan yield stress model requires hydrology flux and gradient inputs.\n"
                                  "Ensure that a routing hydrology model is active.");
  }

  const bool slippery_grounding_lines = m_config->get_flag("basal_yield_stress.slippery_grounding_lines");

  const double ice_density         = m_config->get_number("constants.ice.density"),
               standard_gravity    = m_config->get_number("constants.standard_gravity"),
               latent_heat         = m_config->get_number("constants.fresh_water.latent_heat_of_fusion"),
               fresh_water_density = m_config->get_number("constants.fresh_water.density"),
               arrhenius_parameter = m_config->get_number("flow_law.isothermal_Glen.ice_softness"),
               Glen_exponent       = m_config->get_number("stress_balance.ssa.Glen_exponent"),
               // tuneable parameters
               high_tauc           = m_config->get_number("basal_yield_stress.ice_free_bedrock")*100, //1e6*100 Pa*100, high_tauc is lower than the usual yield stresses calculated by the MohrCoulomb relation
               delta               = m_config->get_number("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden"),
               rocky_phi_deg       = m_config->get_number("basal_yield_stress.gowan.rocky_phi"),
               seddy_phi_deg       = m_config->get_number("basal_yield_stress.gowan.seddy_phi"), 
               tau_ice_rock        = m_config->get_number("basal_yield_stress.gowan.ice_rock_yield_stress"), //1e5 Pa
               rchanneldistance    = m_config->get_number("basal_yield_stress.gowan.Rothsliberger_distance"), // 12 km
               protrusion_height   = m_config->get_number("basal_yield_stress.gowan.protrusion_height"); // 0.01 m

  // Derived constants:
  const double rocky_angle = rocky_phi_deg * M_PI / 180.0;
  const double seddy_angle = seddy_phi_deg * M_PI / 180.0;
  const double c1 = 1.0 / (ice_density * latent_heat);
  const double c2 = 2.0 * arrhenius_parameter * pow(Glen_exponent,-Glen_exponent);
  const double f = 0.1; // This is hard coded for now
  const double c3 = pow(2.0,(1.0/4.0)) * pow(M_PI+2.0,(1/2)) / (pow(M_PI,(1/4)) * pow(fresh_water_density*f,(1/2)));
  const double alpha = 5.0/4.0;

  MohrCoulombPointwise mc(m_config);

  const array::Vector &Q              = *inputs.hydrology_flux;   // water flux (m²/s)
  const array::Scalar &grad           = *inputs.hydrology_gradient, // |∇ψ| (Pa/m)
                      &W_till         = *inputs.till_water_thickness,
                      &ice_thickness  = inputs.geometry->ice_thickness,
                      &sliding_speed  = *inputs.ice_sliding_speed;

  const auto          &cell_type      = inputs.geometry->cell_type;
  const auto          &bed_topography = inputs.geometry->bed_elevation;
  const auto          &sea_level      = inputs.geometry->sea_level_elevation;

  array::AccessScope scope{&W_till, &Q, &grad, &cell_type, &m_basal_yield_stress, &hydro_tauc, &tauc_ratio,
                             &m_sliding_mechanism, &m_till_phi, &m_effective_pressure, &m_till_cover_local, &bed_topography, &sea_level, &ice_thickness, &sliding_speed};

   {
    
    tauc_ratio.set(0.0);

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      // Ocean & ice-free cases
      if (cell_type.ocean(i, j)) {
        m_basal_yield_stress(i, j) = 0.0;
        m_sliding_mechanism(i, j)  = 0;
        continue;
      }
      if (cell_type.ice_free(i, j)) {
        m_basal_yield_stress(i, j) = high_tauc;
        m_sliding_mechanism(i, j)  = 0;
        continue;
      }

      // Grounded ice
      double slippery_tauc = high_tauc;
      double water = W_till(i, j);

      // Optional "slippery grounding lines" // FIXME: needs to be checked, turned off for now
      // if (slippery_grounding_lines &&
      //     bed_topography(i, j) <= sea_level(i, j) &&
      //     (cell_type.next_to_floating_ice(i, j) || cell_type.next_to_ice_free_ocean(i, j))) {

      //   water = tillwat_max;

      //   double grounding_reduction;
      //   const double z = bed_topography(i, j); // bed elevation (negative below sea level) //What is this? Not continous at z=b=-1000,-2000
      //   if (z >= -1000.0) {
      //     grounding_reduction = z * 1.0e-5 + 0.2;      //  ~0.2 at z=0, down to ~0.19 at z=-10m, etc.
      //   } else if (z >= -2000.0) {
      //     grounding_reduction = z * 9.0e-6 + 0.019;    // from ~0.010 at -1000 to ~0.001 at -2000
      //   } else {
      //     grounding_reduction = 0.001;
      //   }
      //   slippery_tauc = ice_thickness(i, j) * standard_gravity * ice_density * grounding_reduction;
      // }

      // This should be the same as in MohrCoulombPiecewise.cc. Does it have ice/rock cap?
      // Naming should be more consistent across files.
      // Sediment = till
      // Two effective pressures N_till, N_hydro
      // Doesnt use a delta field but a constant delta read from config.

      double P_overburden = ice_density * standard_gravity * ice_thickness(i, j);
      // Use MohrCoulombPiecewise for sediment yield stress:

      double tau_sed = mc.yield_stress(delta, P_overburden, water, m_till_phi(i, j));

      // Cap by ice/rock limit, as in Gowan:
      // if (tau_sed > tau_ice_rock) {
      //   tau_sed = tau_ice_rock;
      // }

      // Mix sediment-covered and rock fractions (hydrostatic rock limit):
      double tau_mix = tau_sed * m_till_cover_local(i, j)
                     + tau_ice_rock * (1.0 - m_till_cover_local(i, j));

      m_basal_yield_stress(i, j) = tau_mix;
      m_sliding_mechanism(i, j)  = 1;  // sediment default

      // Hydrology-based sliding yield stress (only if velocity > 0):
      double yield_stress_hydrology = high_tauc;
      if (sliding_speed(i, j) > 0.0) {

        const double qx = Q(i, j).u;
        const double qy = Q(i, j).v;
        const double gradpsi = grad(i, j);
        const double vel2 = sliding_speed(i, j);

        const double Q_rchannel = std::sqrt(qx * qx + qy * qy) * rchanneldistance / m_dx; // FIXME: only uses dx
        
        // Schoof (2010), Eq. (2):
        double num = c1 * Q_rchannel * gradpsi + vel2 * protrusion_height;
        double denom = c2 * std::pow(c3, -1.0/alpha) * std::pow(Q_rchannel, 1.0/alpha) * std::pow(gradpsi, -1.0 / (2.0 * alpha));
        double N_eff = std::pow(num / denom, 1.0 / Glen_exponent);

        if (!std::isfinite(N_eff)) {
          N_eff = high_tauc;
        }

        m_effective_pressure(i, j) = N_eff;
        // equation 3 in Schoof 2010, used to determine the hydrology system type
        double Qc = vel2 * protrusion_height / (c1 * (alpha-1.0) * gradpsi); // FIXME: should be stored

        // Rocky and sediment areas (Gowan used tan(phi)*N_eff)
        double tau_hydro_sed = N_eff * std::tan(seddy_angle);

        // If sediments are weaker than hydro sliding, sediments take precedence (cap):
        if (tau_sed < tau_hydro_sed) {
          tau_hydro_sed = tau_sed;
        }

        yield_stress_hydrology =
            tau_hydro_sed * m_till_cover_local(i, j)
          + N_eff * std::tan(rocky_angle) * (1.0 - m_till_cover_local(i, j));
      }

      hydro_tauc(i, j) = yield_stress_hydrology;

      // Diagnostic ratio (Gowan wrote m_basal / hydro)
      if (yield_stress_hydrology > 0.0) {
        tauc_ratio(i, j) = m_basal_yield_stress(i, j) / yield_stress_hydrology;
      } else {
        tauc_ratio(i, j) = 0.0;
      }

      // Pick the MINIMUM of (sediment yield stress) and (hydrology yield stress),
      if (yield_stress_hydrology < m_basal_yield_stress(i, j)) {
        m_basal_yield_stress(i, j) = yield_stress_hydrology;
        m_sliding_mechanism(i, j)  = 2;  // hydrology
      }

      // then also apply slippery grounding reduction if smaller.
      // if (slippery_tauc < m_basal_yield_stress(i, j)) {
      //   m_basal_yield_stress(i, j) = slippery_tauc;
      //   m_sliding_mechanism(i, j)  = 3;  // slippery GL
      // }
    }
  }
  
  m_effective_pressure.update_ghosts();
  m_basal_yield_stress.update_ghosts();
}

DiagnosticList GowanYieldStress::diagnostics_impl() const {
  return combine({
      {"effective_pressure", Diagnostic::wrap(m_effective_pressure)},
      {"sliding_mechanism",  Diagnostic::wrap(m_sliding_mechanism)},
      {"till_cover_local",   Diagnostic::wrap(m_till_cover_local)},
      {"hydro_tauc",         Diagnostic::wrap(hydro_tauc)},
      {"tauc_ratio",         Diagnostic::wrap(tauc_ratio)}
    },
    MohrCoulombYieldStress::diagnostics_impl());  // include base diagnostics
}

} // end of namespace pism
