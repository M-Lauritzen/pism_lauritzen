/* Copyright (C) 2016, 2017, 2018, 2019, 2020, 2022, 2023 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <gsl/gsl_math.h>       // GSL_NAN

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace atmosphere {

array::Scalar::Ptr AtmosphereModel::allocate_temperature(IceGrid::ConstPtr grid) {
  array::Scalar::Ptr result(new array::Scalar(grid, "air_temp"));

  result->set_attrs("climate_forcing", "mean annual near-surface air temperature",
                    "Kelvin", "Kelvin", "", 0);

  return result;
}

array::Scalar::Ptr AtmosphereModel::allocate_precipitation(IceGrid::ConstPtr grid) {
  array::Scalar::Ptr result(new array::Scalar(grid, "precipitation"));
  result->set_attrs("climate_forcing", "precipitation rate",
                    "kg m-2 second-1", "kg m-2 year-1",
                    "precipitation_flux", 0);

  return result;
}

AtmosphereModel::AtmosphereModel(IceGrid::ConstPtr g)
  : Component(g) {
  // empty
}

AtmosphereModel::AtmosphereModel(IceGrid::ConstPtr g,
                                 std::shared_ptr<AtmosphereModel> input)
  :Component(g), m_input_model(input) {
  // empty
}

void AtmosphereModel::init(const Geometry &geometry) {
  this->init_impl(geometry);
}

void AtmosphereModel::update(const Geometry &geometry, double t, double dt) {
  this->update_impl(geometry, t, dt);
}

const array::Scalar& AtmosphereModel::precipitation() const {
  return this->precipitation_impl();
}

const array::Scalar& AtmosphereModel::air_temperature() const {
  return this->air_temperature_impl();
}

void AtmosphereModel::begin_pointwise_access() const {
  this->begin_pointwise_access_impl();
}

void AtmosphereModel::end_pointwise_access() const {
  this->end_pointwise_access_impl();
}

void AtmosphereModel::init_timeseries(const std::vector<double> &ts) const {
  this->init_timeseries_impl(ts);
}

void AtmosphereModel::precip_time_series(int i, int j, std::vector<double> &result) const {
  result.resize(m_ts_times.size());
  this->precip_time_series_impl(i, j, result);
}

void AtmosphereModel::temp_time_series(int i, int j, std::vector<double> &result) const {
  result.resize(m_ts_times.size());
  this->temp_time_series_impl(i, j, result);
}

namespace diagnostics {

/*! @brief Instantaneous near-surface air temperature. */
class AirTemperatureSnapshot : public Diag<AtmosphereModel> {
public:
  AirTemperatureSnapshot(const AtmosphereModel *m)
    : Diag<AtmosphereModel>(m) {

    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "air_temp_snapshot")};

    set_attrs("instantaneous value of the near-surface air temperature",
              "",                 // no standard name
              "Kelvin", "Kelvin", 0);
  }
protected:
  array::Array::Ptr compute_impl() const {

    array::Scalar::Ptr result(new array::Scalar(m_grid, "air_temp_snapshot"));
    result->metadata(0) = m_vars[0];

    std::vector<double> current_time(1, m_grid->ctx()->time()->current());
    std::vector<double> temperature(1, 0.0);

    model->init_timeseries(current_time);

    model->begin_pointwise_access();

    array::AccessScope list(*result);
    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        model->temp_time_series(i, j, temperature);

        (*result)(i, j) = temperature[0];
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();

    model->end_pointwise_access();

    return result;
  }
};

/*! @brief Effective near-surface mean-annual air temperature. */
class AirTemperature : public Diag<AtmosphereModel> {
public:
  AirTemperature(const AtmosphereModel *m)
    : Diag<AtmosphereModel>(m) {

    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "effective_air_temp")};

    set_attrs("effective mean-annual near-surface air temperature", "",
              "Kelvin", "Kelvin", 0);
  }
protected:
  array::Array::Ptr compute_impl() const {

    array::Scalar::Ptr result(new array::Scalar(m_grid, "effective_air_temp"));
    result->metadata(0) = m_vars[0];

    result->copy_from(model->air_temperature());

    return result;
  }
};

/*! @brief Effective precipitation rate (average over time step). */
class Precipitation : public Diag<AtmosphereModel> {
public:
  Precipitation(const AtmosphereModel *m)
    : Diag<AtmosphereModel>(m) {

    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "effective_precipitation")};

    set_attrs("effective precipitation rate",
              "precipitation_flux",
              "kg m-2 second-1", "kg m-2 year-1", 0);
  }
protected:
  array::Array::Ptr compute_impl() const {

    array::Scalar::Ptr result(new array::Scalar(m_grid, "effective_precipitation"));
    result->metadata(0) = m_vars[0];

    result->copy_from(model->precipitation());

    return result;
  }
};

} // end of namespace diagnostics

void AtmosphereModel::update_impl(const Geometry &geometry, double t, double dt) {
  if (m_input_model) {
    m_input_model->update_impl(geometry, t, dt);
  }
}

MaxTimestep AtmosphereModel::max_timestep_impl(double my_t) const {
  if (m_input_model) {
    return m_input_model->max_timestep(my_t);
  }
  return MaxTimestep("atmosphere model");
}

DiagnosticList AtmosphereModel::diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    {"air_temp_snapshot",       Diagnostic::Ptr(new AirTemperatureSnapshot(this))},
    {"effective_air_temp",      Diagnostic::Ptr(new AirTemperature(this))},
    {"effective_precipitation", Diagnostic::Ptr(new Precipitation(this))},
  };

  if (m_input_model) {
    result = combine(result, m_input_model->diagnostics());
  }

  return result;
}

TSDiagnosticList AtmosphereModel::ts_diagnostics_impl() const {
  if (m_input_model) {
    return m_input_model->ts_diagnostics();
  }

  return {};
}

void AtmosphereModel::define_model_state_impl(const File &output) const {
  if (m_input_model) {
    m_input_model->define_model_state(output);
  }
}

void AtmosphereModel::write_model_state_impl(const File &output) const {
  if (m_input_model) {
    m_input_model->write_model_state(output);
  }
}

const array::Scalar& AtmosphereModel::precipitation_impl() const {
  if (m_input_model) {
    return m_input_model->precipitation();
  }
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
}

const array::Scalar& AtmosphereModel::air_temperature_impl() const {
  if (m_input_model) {
    return m_input_model->air_temperature();
  }
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
}

void AtmosphereModel::begin_pointwise_access_impl() const {
  if (m_input_model) {
    m_input_model->begin_pointwise_access();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void AtmosphereModel::end_pointwise_access_impl() const {
  if (m_input_model) {
    m_input_model->end_pointwise_access();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void AtmosphereModel::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  if (m_input_model) {
    m_input_model->temp_time_series(i, j, result);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void AtmosphereModel::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  if (m_input_model) {
    m_input_model->precip_time_series(i, j, result);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void AtmosphereModel::init_timeseries_impl(const std::vector<double> &ts) const {
  if (m_input_model) {
    m_input_model->init_timeseries(ts);
    m_ts_times = ts;
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

} // end of namespace atmosphere
} // end of namespace pism
