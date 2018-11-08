// Copyright (C) 2009-2018 Ricarda Winkelmann, Torsten Albrecht, Constantine Khrulev
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

// Implementation of the atmosphere model using constant-in-time precipitation
// and a cosine yearly cycle for near-surface air temperatures.

// This includes the PIK temperature parameterization.

#include "TemperaturePIK.hh"

#include "pism/geometry/Geometry.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace atmosphere {

TemperaturePIK::TemperaturePIK(IceGrid::ConstPtr g)
  : YearlyCycle(g) {

  auto parameterization = m_config->get_string("atmosphere.pik_temp.parameterization");

  std::map<std::string, Parameterization>
    models = {{"huybrechts_dewolde99", HUYBRECHTS_DEWOLDE99},
              {"era_interim",          ERA_INTERIM},
              {"era_interim_sin",      ERA_INTERIM_SIN},
              {"era_interim_lon",      ERA_INTERIM_LON},
              {"default",              DEFAULT}};

  if (models.find(parameterization) == models.end()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid pik_temp parameterization: %s",
                                  parameterization.c_str());
  }

  m_parameterization = models[parameterization];
}

TemperaturePIK::~TemperaturePIK() {
  // empty
}

void TemperaturePIK::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing PIK atmosphere model with air temperature parameterization based on \n"
                 "  Huybrechts & De Wolde (1999) or multiple regression analysis of ERA INTERIM data...\n"
                 "  Uses a time-independent precipitation field read from a file.");

  m_reference = "Winkelmann et al.";

  auto precip_file = m_config->get_string("atmosphere.pik_temp.file");

  if (not precip_file.empty()) {
    YearlyCycle::init_internal(precip_file,
                               true, /* do regrid */
                               0 /* start (irrelevant) */);
  } else {
    // try to read precipitation from the input (-i) file
    YearlyCycle::init_impl(geometry);
  }

  switch (m_parameterization) {
  case HUYBRECHTS_DEWOLDE99:
    m_log->message(2,
                   "    Parameterization based on: Huybrechts & De Wolde (1999).\n");
    break;
  case ERA_INTERIM:
    m_log->message(2,
                   "    Parameterization based on: multiple regression analysis of ERA INTERIM data.\n");
    break;
  case ERA_INTERIM_SIN:
    m_log->message(2,
                   "    Parameterization based on: multiple regression analysis of ERA INTERIM data with a sin(lat) dependence.\n");
    break;
  case ERA_INTERIM_LON:
    m_log->message(2,
                   "    Parameterization based on: multiple regression analysis of ERA INTERIM data with a cos(lon) dependence.\n");
    break;
  default:
    m_log->message(2,
                   "    Mean annual temperature: as in Martin et al (2011).\n"
                   "    Mean summer temperature: anomaly to the parameterization used by Huybrechts & De Wolde (1999).\n");
  }
}

MaxTimestep TemperaturePIK::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere pik_temp");
}

static double huybrechts_dewolde99_mean_annual(double surface_elevation, double latitude) {
  double gamma_a = surface_elevation < 1500.0 ? -0.005102 : -0.014285;
  return 273.15 + 34.46 + gamma_a * surface_elevation - 0.68775 * latitude * (-1.0);
}

static double huybrechts_dewolde99_mean_summer(double surface_elevation, double latitude) {
  return 273.15 + 16.81 - 0.00692 * surface_elevation - 0.27937 * latitude * (-1.0);
}

/*!
 * Parameterization of mean annual and mean summer near-surface temperature as in
 * Huybrechts & DeWolde (1999)
 */
static void huybrechts_dewolde99(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude;

  IceModelVec::AccessList list{&h, &lat, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    T_ma(i, j) = huybrechts_dewolde99_mean_annual(h(i, j), lat(i, j));
    T_ms(i, j) = huybrechts_dewolde99_mean_summer(h(i, j), lat(i, j));
  }
}

/*!
 * Parametrization based on multiple regression analysis of ERA INTERIM data
 */
static void era_interim(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude;

  IceModelVec::AccessList list{&h, &lat, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    T_ma(i, j) = 273.15 + 29.2 - 0.0082 * h(i, j) - 0.576 * lat(i, j) * (-1.0);
    T_ms(i, j) = 273.15 + 16.5 - 0.0068 * h(i, j) - 0.248 * lat(i, j) * (-1.0);
  }
}

/*!
 * Parametrization based on multiple regression analysis of ERA INTERIM data with sin(lat)
 */
static void era_interim_sin(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude;

  IceModelVec::AccessList list{&h, &lat, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    T_ma(i, j) = 273.15 - 2.0 - 0.0082 * h(i, j) + 18.4 * (sin(3.1415 * lat(i, j) / 180.0) + 0.8910) / (1 - 0.8910);
    T_ms(i, j) = 273.15 + 3.2 - 0.0067 * h(i, j) + 8.3 * (sin(3.1415 * lat(i, j) / 180.0) + 0.8910) / (1 - 0.8910);
  }
}

/*!
 * Parametrization based on multiple regression analysis of ERA INTERIM data with cos(lon)
 */
static void era_interim_lon(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude,
    &lon = geometry.longitude;

  IceModelVec::AccessList list{&h, &lat, &lon, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double hmod = std::max(500.0, h(i, j)); // in the fit for ice free ocean hmod was set to 0
    T_ma(i, j) = 273.15 + 36.81 - 0.00797 * hmod - 0.688 * lat(i, j) * (-1.0) + 2.574 * cos(3.1415 * (lon(i, j) - 110.0) / 180.0);
    T_ms(i, j) = 273.15 + 22.58 - 0.00940 * hmod - 0.234 * lat(i, j) * (-1.0) + 0.828 * cos(3.1415 * (lon(i, j) - 110.0) / 180.0);
  }
}

/*!
 * - annual mean temperature as in Martin et al. (2011)
 * - summer mean temperature computed as an anomaly to Huybrechts & DeWolde (1999)
 */
static void martin_huybrechts_dewolde(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude;

  IceModelVec::AccessList list{&h, &lat, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // mean annual surface temperature as in Martin et al. 2011, Eqn. 2.0.2
    T_ma(i, j) = 273.15 + 30 - 0.0075 * h(i, j) - 0.68775 * lat(i, j) * (-1.0);

    double
      TMA = huybrechts_dewolde99_mean_annual(h(i, j), lat(i, j)),
      TMS = huybrechts_dewolde99_mean_summer(h(i, j), lat(i, j));

    T_ms(i, j) = T_ma(i, j) + (TMS - TMA);
  }
}


/*!
 * Updates mean annual and mean summer (January) near-surface air temperatures. Note that
 * the precipitation rate is time-independent and does not need to be updated.
 */
void TemperaturePIK::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  if (geometry.latitude.metadata().has_attribute("missing_at_bootstrap")) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "latitude variable was missing at bootstrap;\n"
                       "TemperaturePIK atmosphere model depends on latitude and would return nonsense!");
  }

  switch (m_parameterization) {
  case HUYBRECHTS_DEWOLDE99:
    huybrechts_dewolde99(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  case ERA_INTERIM:
    era_interim(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  case ERA_INTERIM_SIN:
    era_interim_sin(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  case ERA_INTERIM_LON:
    era_interim_lon(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  default:
    martin_huybrechts_dewolde(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
  }
}

} // end of namespace atmosphere
} // end of namespace pism
