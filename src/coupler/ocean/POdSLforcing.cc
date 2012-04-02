// Copyright (C) 2011, 2012 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include "POdSLforcing.hh"

/// -dSLforcing

POdSLforcing::POdSLforcing(IceGrid &g, const NCConfigVariable &conf, PISMOceanModel* in)
  : PScalarForcing<PISMOceanModel,POModifier>(g, conf, in)
{
  option = "-dSLforcing";
  offset_name = "delta_sea_level";
  offset = new Timeseries(&grid, offset_name, config.get_string("time_dimension_name"));

  offset->set_units("m", "");
  offset->set_dimension_units(grid.time->units(), "");
  offset->set_attr("long_name", "sea level elevation offsets");
}

PetscErrorCode POdSLforcing::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "* Initializing sea level forcing...\n"); CHKERRQ(ierr);

  ierr = init_internal(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode POdSLforcing::sea_level_elevation(PetscReal &result) {
  PetscErrorCode ierr = input_model->sea_level_elevation(result); CHKERRQ(ierr);

  if (offset)
    result += (*offset)(t + 0.5*dt);

  return 0;
}
