// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev and Ed Bueler
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

#include "PISMStressBalance.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "PISMOcean.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"
#include "Mask.hh"
#include "enthalpyConverter.hh"
#include "PISMConfig.hh"

namespace pism {

StressBalance::StressBalance(IceGrid &g,
                             ShallowStressBalance *sb,
                             SSB_Modifier *ssb_mod,
                             const Config &conf)
  : Component(g, conf), m_stress_balance(sb), m_modifier(ssb_mod) {

  m_basal_melt_rate = NULL;
  m_variables = NULL;

  // allocate the vertical velocity field:
  m_w.create(grid, "wvel_rel", WITHOUT_GHOSTS);
  m_w.set_attrs("diagnostic",
                "vertical velocity of ice, relative to base of ice directly below",
                "m s-1", "");
  m_w.set_time_independent(false);
  m_w.set_glaciological_units("m year-1");
  m_w.write_in_glaciological_units = true;

  m_strain_heating.create(grid, "strain_heating", WITHOUT_GHOSTS);
  m_strain_heating.set_attrs("internal",
                             "rate of strain heating in ice (dissipation heating)",
                             "W m-3", "");
}

StressBalance::~StressBalance() {
  delete m_stress_balance;
  delete m_modifier;
}

//! \brief Initialize the StressBalance object.
void StressBalance::init(Vars &vars) {
  m_variables = &vars;

  m_stress_balance->init(vars);
  m_modifier->init(vars);
}

void StressBalance::set_boundary_conditions(IceModelVec2Int &locations,
                                                      IceModelVec2V &velocities) {
  m_stress_balance->set_boundary_conditions(locations, velocities);
}

//! \brief Set the basal melt rate. (If not NULL, it will be included in the
//! computation of the vertical valocity).
void StressBalance::set_basal_melt_rate(IceModelVec2S *bmr_input) {
  m_basal_melt_rate = bmr_input;
}

//! \brief Performs the shallow stress balance computation.
void StressBalance::update(bool fast, double sea_level,
                                     IceModelVec2S &melange_back_pressure) {
  IceModelVec2V *velocity_2d;
  IceModelVec3  *u, *v;

  // Tell the ShallowStressBalance object about the current sea level:
  m_stress_balance->set_sea_level_elevation(sea_level);

  try {
    grid.profiling.begin("SSB");
    m_stress_balance->update(fast, melange_back_pressure);
    grid.profiling.end("SSB");
    m_stress_balance->get_2D_advective_velocity(velocity_2d);

    grid.profiling.begin("SB modifier");
    m_modifier->update(velocity_2d, fast);
    grid.profiling.end("SB modifier");

    if (fast == false) {

      m_modifier->get_horizontal_3d_velocity(u, v);

      grid.profiling.begin("SB strain heat");
      this->compute_volumetric_strain_heating();
      grid.profiling.end("SB strain heat");

      grid.profiling.begin("SB vert. vel.");
      this->compute_vertical_velocity(u, v, m_basal_melt_rate, m_w);
      grid.profiling.end("SB vert. vel.");
    }
  }
  catch (RuntimeError &e) {
    e.add_context("updating the stress balance");
    throw;
  }
}

void StressBalance::get_2D_advective_velocity(IceModelVec2V* &result) {
  m_stress_balance->get_2D_advective_velocity(result);
}

void StressBalance::get_diffusive_flux(IceModelVec2Stag* &result) {
  m_modifier->get_diffusive_flux(result);
}

void StressBalance::get_max_diffusivity(double &D) {
  m_modifier->get_max_diffusivity(D);
}

void StressBalance::get_3d_velocity(IceModelVec3* &u, IceModelVec3* &v, IceModelVec3* &w_out) {
  m_modifier->get_horizontal_3d_velocity(u, v);
  w_out = &m_w;
}

void StressBalance::get_basal_frictional_heating(IceModelVec2S* &result) {
  m_stress_balance->get_basal_frictional_heating(result);
}

void StressBalance::get_volumetric_strain_heating(IceModelVec3* &result) {
  result = &m_strain_heating;
}

void StressBalance::compute_2D_principal_strain_rates(IceModelVec2V &velocity,
                                                                IceModelVec2Int &mask,
                                                                IceModelVec2 &result) {
  m_stress_balance->compute_2D_principal_strain_rates(velocity, mask, result);
}

void StressBalance::compute_2D_stresses(IceModelVec2V &velocity,
                                                  IceModelVec2Int &mask,
                                                  IceModelVec2 &result) {
  m_stress_balance->compute_2D_stresses(velocity, mask, result);
}

//! Compute vertical velocity using incompressibility of the ice.
/*!
The vertical velocity \f$w(x,y,z,t)\f$ is the velocity *relative to the
location of the base of the ice column*.  That is, the vertical velocity
computed here is identified as \f$\tilde w(x,y,s,t)\f$ in the page
[]@ref vertchange.

Thus \f$w<0\f$ here means that that
that part of the ice is getting closer to the base of the ice, and so on.
The slope of the bed (i.e. relative to the geoid) and/or the motion of the
bed (i.e. from bed deformation) do not affect the vertical velocity.

In fact the following statement is exactly true if the basal melt rate is zero:
the vertical velocity at a point in the ice is positive (negative) if and only
if the average horizontal divergence of the horizontal velocity, in the portion
of the ice column below that point, is negative (positive).
In particular, because \f$z=0\f$ is the location of the base of the ice
always, the only way to have \f$w(x,y,0,t) \ne 0\f$ is to have a basal melt
rate.

Incompressibility itself says
   \f[ \nabla\cdot\mathbf{U} + \frac{\partial w}{\partial z} = 0. \f]
This is immediately equivalent to the integral
   \f[ w(x,y,z,t) = - \int_{b(x,y,t)}^{z} \nabla\cdot\mathbf{U}\,d\zeta
                           + w_b(x,y,t). \f]
Here the value \f$w_b(x,y,t)\f$ is either zero or the negative of the basal melt rate
according to the value of the flag `include_bmr_in_continuity`.

The vertical integral is computed by the trapezoid rule.
 */
void StressBalance::compute_vertical_velocity(IceModelVec3 *u, IceModelVec3 *v,
                                              IceModelVec2S *basal_melt_rate,
                                              IceModelVec3 &result) {
  IceModelVec2Int *mask = m_variables->get_2d_mask("mask");
  MaskQuery m(*mask);

  IceModelVec::AccessList list;
  list.add(*u);
  list.add(*v);
  list.add(result);
  list.add(*mask);

  if (basal_melt_rate) {
    list.add(*basal_melt_rate);
  }

  double *w_ij, *u_ij, *u_w, *u_e, *v_ij, *v_s, *v_n;

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result.getInternalColumn(i,j,&w_ij);

    u->getInternalColumn(i-1,j,&u_w);
    u->getInternalColumn(i,j,  &u_ij);
    u->getInternalColumn(i+1,j,&u_e);

    v->getInternalColumn(i,j-1,&v_s);
    v->getInternalColumn(i,j,  &v_ij);
    v->getInternalColumn(i,j+1,&v_n);

    double west = 1, east = 1, south = 1, north = 1,
      D_x = 0,                // 1/(dx), 1/(2dx), or 0
      D_y = 0;                // 1/(dy), 1/(2dy), or 0

    // Switch between second-order centered differences in the interior and
    // first-order one-sided differences at ice margins.

    // x-derivative
    {
      if ((m.icy(i,j) && m.ice_free(i+1,j)) || (m.ice_free(i,j) && m.icy(i+1,j))) {
        east = 0;
      }
      if ((m.icy(i,j) && m.ice_free(i-1,j)) || (m.ice_free(i,j) && m.icy(i-1,j))) {
        west = 0;
      }

      if (east + west > 0) {
        D_x = 1.0 / (grid.dx() * (east + west));
      } else {
        D_x = 0.0;
      }
    }

    // y-derivative
    {
      if ((m.icy(i,j) && m.ice_free(i,j+1)) || (m.ice_free(i,j) && m.icy(i,j+1))) {
        north = 0;
      }
      if ((m.icy(i,j) && m.ice_free(i,j-1)) || (m.ice_free(i,j) && m.icy(i,j-1))) {
        south = 0;
      }

      if (north + south > 0) {
        D_y = 1.0 / (grid.dy() * (north + south));
      } else {
        D_y = 0.0;
      }
    }

    double u_x = D_x * (west * (u_ij[0] - u_w[0]) + east * (u_e[0] - u_ij[0])),
      v_y = D_y * (south * (v_ij[0] - v_s[0]) + north * (v_n[0] - v_ij[0]));

    // at the base: include the basal melt rate
    if (basal_melt_rate) {
      w_ij[0] = - (*basal_melt_rate)(i,j);
    } else {
      w_ij[0] = 0.0;
    }

    // within the ice and above:
    double old_integrand = u_x + v_y;
    for (unsigned int k = 1; k < grid.Mz(); ++k) {
      u_x = D_x * (west  * (u_ij[k] - u_w[k]) + east  * (u_e[k] - u_ij[k]));
      v_y = D_y * (south * (v_ij[k] - v_s[k]) + north * (v_n[k] - v_ij[k]));
      const double new_integrand = u_x + v_y;

      const double dz = grid.z(k) - grid.z(k-1);

      w_ij[k] = w_ij[k-1] - 0.5 * (new_integrand + old_integrand) * dz;

      old_integrand = new_integrand;
    }
  }
}

/**
 * This function computes \f$D^2\f$ defined by
 *
 * \f[ 2D^2 = D_{ij} D_{ij}\f]
 * or
 * \f[
 * D^2 = \frac{1}{2}\,\left(\frac{1}{2}\,(v_{z})^2 + (v_{y} + u_{x})^2 +
 *       (v_{y})^2 + \frac{1}{2}\,(v_{x} + u_{y})^2 + \frac{1}{2}\,(u_{z})^2 +
 *       (u_{x})^2\right)
 * \f]
 *
 * (note the use of the summation convention). Here \f$D_{ij}\f$ is the
 * strain rate tensor. See
 * StressBalance::compute_volumetric_strain_heating() for details.
 *
 * @param u_x,u_y,u_z partial derivatives of \f$u\f$, the x-component of the ice velocity
 * @param v_x,v_y,v_z partial derivatives of \f$v\f$, the y-component of the ice velocity
 *
 * @return \f$D^2\f$, where \f$D\f$ is defined above.
 */
static inline double D2(double u_x, double u_y, double u_z, double v_x, double v_y, double v_z) {
  return 0.5 * (PetscSqr(u_x + v_y) + u_x*u_x + v_y*v_y + 0.5 * (PetscSqr(u_y + v_x) + u_z*u_z + v_z*v_z));
}

/**
  \brief Computes the volumetric strain heating using horizontal
  velocity.

  Following the notation used in [\ref BBssasliding], let \f$u\f$ be a
  three-dimensional *vector* velocity field. Then the strain rate
  tensor \f$D_{ij}\f$ is defined by

  \f[ D_{ij} = \frac 12 \left(\diff{u_{i}}{x_{j}} + \diff{u_{j}}{x_{i}} \right), \f]

  Where \f$i\f$ and \f$j\f$ range from \f$1\f$ to \f$3\f$.

  The flow law in the viscosity form states

  \f[ \tau_{ij} = 2 \eta D_{ij}, \f]

  and the nonlinear ice viscosity satisfies

  \f[ 2 \eta = B(T) D^{(1/n) - 1}. \f]

  Here \f$D^{2}\f$ is defined by \f$2D^{2} = D_{ij}D_{ij}\f$ (using the
  summation convention) and \f$B(T) = A(T)^{-1/n}\f$ is the ice hardness.

  Now the volumetric strain heating is

  \f[ \Sigma = \sum_{i,j=1}^{3}D_{ij}\tau_{ij} = 2 B(T) D^{(1/n) + 1}. \f]

  We use an *approximation* of \f$D_{ij}\f$ common in shallow ice models:

  - we assume that horizontal derivatives of the vertical velocity are
    much smaller than \f$z\f$ derivatives horizontal velocity
    components \f$u\f$ and \f$v\f$. (We drop \f$w_x\f$ and \f$w_y\f$
    terms in \f$D_{ij}\f$.)

  - we use the incompressibility of ice to approximate \f$w_z\f$:

  \f[ w_z = - (u_x + v_y). \f]

  Requires ghosts of `u` and `v` velocity components and uses the fact
  that `u` and `v` above the ice are filled using constant
  extrapolation.

  Resulting field does not have ghosts.

  Below is the *Maxima* code that produces the expression evaluated by D2().

       derivabbrev : true;
       U : [u, v, w]; X : [x, y, z]; depends(U, X);
       gradef(w, x, 0); gradef(w, y, 0);
       gradef(w, z, -(diff(u, x) + diff(v, y)));
       d[i,j] := 1/2 * (diff(U[i], X[j]) + diff(U[j], X[i]));
       D : genmatrix(d, 3, 3), ratsimp, factor;
       tex('D = D);
       tex('D^2 = 1/2 * mat_trace(D . D));

  @return 0 on success
 */
void StressBalance::compute_volumetric_strain_heating() {
  PetscErrorCode ierr;

  const IceFlowLaw *flow_law = m_stress_balance->get_flow_law();
  EnthalpyConverter &EC = m_stress_balance->get_enthalpy_converter();

  IceModelVec3 *u, *v;
  m_modifier->get_horizontal_3d_velocity(u, v);

  IceModelVec2S *thickness = m_variables->get_2d_scalar("land_ice_thickness");
  IceModelVec3  *enthalpy  = m_variables->get_3d_scalar("enthalpy");

  IceModelVec2Int *mask = m_variables->get_2d_mask("mask");
  MaskQuery m(*mask);

  double
    enhancement_factor = flow_law->enhancement_factor(),
    n = flow_law->exponent(),
    exponent = 0.5 * (1.0 / n + 1.0),
    e_to_a_power = pow(enhancement_factor,-1.0/n);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(*enthalpy);
  list.add(m_strain_heating);
  list.add(*thickness);
  list.add(*u);
  list.add(*v);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double H = (*thickness)(i,j);
    int ks = grid.kBelowHeight(H);
    double
      *u_ij, *u_w, *u_n, *u_e, *u_s,
      *v_ij, *v_w, *v_n, *v_e, *v_s;
    double *Sigma, *E_ij;

    double west = 1, east = 1, south = 1, north = 1,
      D_x = 0,                // 1/(dx), 1/(2dx), or 0
      D_y = 0;                // 1/(dy), 1/(2dy), or 0

    // x-derivative
    {
      if ((m.icy(i,j) && m.ice_free(i+1,j)) || (m.ice_free(i,j) && m.icy(i+1,j))) {
        east = 0;
      }
      if ((m.icy(i,j) && m.ice_free(i-1,j)) || (m.ice_free(i,j) && m.icy(i-1,j))) {
        west = 0;
      }

      if (east + west > 0) {
        D_x = 1.0 / (grid.dx() * (east + west));
      } else {
        D_x = 0.0;
      }
    }

    // y-derivative
    {
      if ((m.icy(i,j) && m.ice_free(i,j+1)) || (m.ice_free(i,j) && m.icy(i,j+1))) {
        north = 0;
      }
      if ((m.icy(i,j) && m.ice_free(i,j-1)) || (m.ice_free(i,j) && m.icy(i,j-1))) {
        south = 0;
      }

      if (north + south > 0) {
        D_y = 1.0 / (grid.dy() * (north + south));
      } else {
        D_y = 0.0;
      }
    }


    u->getInternalColumn(i,     j,     &u_ij);
    u->getInternalColumn(i - 1, j,     &u_w);
    u->getInternalColumn(i + 1, j,     &u_e);
    u->getInternalColumn(i,     j - 1, &u_s);
    u->getInternalColumn(i,     j + 1, &u_n);

    v->getInternalColumn(i,     j,     &v_ij);
    v->getInternalColumn(i - 1, j,     &v_w);
    v->getInternalColumn(i + 1, j,     &v_e);
    v->getInternalColumn(i,     j - 1, &v_s);
    v->getInternalColumn(i,     j + 1, &v_n);

    enthalpy->getInternalColumn(i, j, &E_ij);
    m_strain_heating.getInternalColumn(i, j, &Sigma);

    for (int k = 0; k <= ks; ++k) {
      double dz,
        pressure = EC.getPressureFromDepth(H - grid.z(k)),
        B        = flow_law->hardness_parameter(E_ij[k], pressure);

      double u_z = 0.0, v_z = 0.0,
        u_x = D_x * (west  * (u_ij[k] - u_w[k]) + east  * (u_e[k] - u_ij[k])),
        u_y = D_y * (south * (u_ij[k] - u_s[k]) + north * (u_n[k] - u_ij[k])),
        v_x = D_x * (west  * (v_ij[k] - v_w[k]) + east  * (v_e[k] - v_ij[k])),
        v_y = D_y * (south * (v_ij[k] - v_s[k]) + north * (v_n[k] - v_ij[k]));

      if (k > 0) {
        dz = grid.z(k+1) - grid.z(k-1);
        u_z = (u_ij[k+1] - u_ij[k-1]) / dz;
        v_z = (v_ij[k+1] - v_ij[k-1]) / dz;
      } else {
        // use one-sided differences for u_z and v_z on the bottom level
        dz = grid.z(1) - grid.z(0);
        u_z = (u_ij[1] - u_ij[0]) / dz;
        v_z = (v_ij[1] - v_ij[0]) / dz;
      }

      Sigma[k] = 2.0 * e_to_a_power * B * pow(D2(u_x, u_y, u_z, v_x, v_y, v_z), exponent);
    } // k-loop

    int remaining_levels = grid.Mz() - (ks + 1);
    if (remaining_levels > 0) {
      ierr = PetscMemzero(&Sigma[ks+1],
                          remaining_levels*sizeof(double));
      PISM_PETSC_CHK(ierr, "PetscMemzero");
    }
  }

}

void StressBalance::stdout_report(std::string &result) {
  std::string tmp1, tmp2;

  m_stress_balance->stdout_report(tmp1);

  m_modifier->stdout_report(tmp2);

  result = tmp1 + tmp2;
}

void StressBalance::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                               IO_Type nctype) {

  m_stress_balance->define_variables(vars, nc, nctype);
  m_modifier->define_variables(vars, nc, nctype);
}


void StressBalance::write_variables(const std::set<std::string> &vars, const PIO &nc) {

  m_stress_balance->write_variables(vars, nc);
  m_modifier->write_variables(vars, nc);
}

void StressBalance::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {

  m_stress_balance->add_vars_to_output(keyword, result);
  m_modifier->add_vars_to_output(keyword, result);
}

} // end of namespace pism
