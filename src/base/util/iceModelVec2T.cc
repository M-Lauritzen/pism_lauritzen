// Copyright (C) 2009--2013 Constantine Khroulev
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

#include <petsc.h>
#include <algorithm>
#include "iceModelVec2T.hh"
#include "PIO.hh"
#include "pism_const.hh"
#include "PISMTime.hh"
#include "LocalInterpCtx.hh"
#include "IceGrid.hh"

IceModelVec2T::IceModelVec2T() : IceModelVec2S() {
  localp = false;
  da3 = PETSC_NULL;
  v3 = PETSC_NULL;
  array3 = NULL;
  first = -1;
  N = 0;
  n_records = 50;		// just a default
  report_range = false;
  lic = NULL;
}

IceModelVec2T::IceModelVec2T(const IceModelVec2T &other) : IceModelVec2S(other) {
  shallow_copy = true;

  array3      = other.array3;
  da3         = other.da3;
  filename    = other.filename;
  first       = other.first;
  N           = other.N;
  lic         = other.lic;
  localp      = other.localp;
  n_records   = other.n_records;
  time        = other.time;
  time_bounds = other.time_bounds;
  v3          = other.v3;

  m_interp_indices = other.m_interp_indices;
  m_interp_weights = other.m_interp_weights;
}

IceModelVec2T::~IceModelVec2T() {
  if (!shallow_copy) {
    delete lic;
    destroy();
  }
}


//! Sets the number of records to store in memory. Call it before calling create().
void IceModelVec2T::set_n_records(unsigned int my_N) {
  n_records = my_N;
}

unsigned int IceModelVec2T::get_n_records() {
  return n_records;
}

PetscErrorCode IceModelVec2T::create(IceGrid &my_grid, string my_short_name,
                                     bool local, int width) {
  PetscErrorCode ierr;

  if (local) {
    SETERRQ(grid->com, 1, "IceModelVec2T cannot be 'local'");
  }

  ierr = IceModelVec2S::create(my_grid, my_short_name, false, width); CHKERRQ(ierr);

  // initialize the da3 member:
  ierr = grid->get_dm(this->n_records, this->da_stencil_width, da3); CHKERRQ(ierr);

  // allocate the 3D Vec:
  ierr = DMCreateGlobalVector(da3, &v3); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2T::destroy() {
  PetscErrorCode ierr;

  ierr = IceModelVec2S::destroy(); CHKERRQ(ierr);

  if (v3 != PETSC_NULL) {
    ierr = VecDestroy(&v3); CHKERRQ(ierr);
    v3 = PETSC_NULL;
  }

  return 0;
}

PetscErrorCode IceModelVec2T::get_array3(PetscScalar*** &a3) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a3 = (PetscScalar***) array3;
  return 0;
}

PetscErrorCode IceModelVec2T::begin_access() {
  PetscErrorCode ierr;
  if (access_counter == 0) {
    ierr = DMDAVecGetArrayDOF(da3, v3, &array3); CHKERRQ(ierr);
  }

  // this call will increment the access_counter
  ierr = IceModelVec2S::begin_access(); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode IceModelVec2T::end_access() {
  // this call will decrement the access_counter
  PetscErrorCode ierr = IceModelVec2S::end_access(); CHKERRQ(ierr);

  if (access_counter == 0) {
    ierr = DMDAVecRestoreArrayDOF(da3, v3, &array3); CHKERRQ(ierr);
    array3 = PETSC_NULL;
  }

  return 0;
}

PetscErrorCode IceModelVec2T::init(string fname) {
  PetscErrorCode ierr;

  filename = fname;

  // We find the variable in the input file and
  // try to find the corresponding time dimension.

  PIO nc(*grid, "guess_mode");
  string name_found;
  bool exists, found_by_standard_name;
  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
  ierr = nc.inq_var(vars[0].short_name, vars[0].get_string("standard_name"),
                    exists, name_found, found_by_standard_name); CHKERRQ(ierr);
  if (!exists) {
    PetscPrintf(grid->com, "PISM ERROR: can't find %s (%s) in %s.\n",
                vars[0].get_string("long_name").c_str(), vars[0].short_name.c_str(),
                filename.c_str());
    PISMEnd();
  }

  // find the time dimension:
  vector<string> dims;
  ierr = nc.inq_vardims(name_found, dims); CHKERRQ(ierr);
  
  string dimname = "";
  bool time_found = false;
  for (unsigned int i = 0; i < dims.size(); ++i) {
    AxisType dimtype;
    dimname = dims[i];

    ierr = nc.inq_dimtype(dimname, dimtype); CHKERRQ(ierr);

    if (dimtype == T_AXIS) {
      time_found = true;
      break;
    }
  }

  if (time_found) {
    // we're found the time dimension
    NCTimeseries time_dimension;
    time_dimension.init(dimname, dimname, grid->com, grid->rank);

    ierr = time_dimension.set_units(grid->time->units()); CHKERRQ(ierr);
    ierr = time_dimension.read(nc, grid->time->use_reference_date(), time); CHKERRQ(ierr);

    string bounds_name;
    ierr = time_dimension.get_bounds_name(nc, bounds_name);

    if (time.size() > 1) {
      if (!bounds_name.empty()) {
        // read time bounds data from a file
        NCTimeBounds tb;
        tb.init(bounds_name, dimname, grid->com, grid->rank);
        ierr = tb.set_units(time_dimension.get_string("units")); CHKERRQ(ierr);

        ierr = tb.read(nc, grid->time->use_reference_date(), time_bounds); CHKERRQ(ierr);

        // time bounds data overrides the time variable: we make t[j] be the
        // right end-point of the j-th interval
        for (unsigned int k = 0; k < time.size(); ++k)
          time[k] = time_bounds[2*k + 1];
      } else {
        // no time bounds attribute
        PetscPrintf(grid->com,
                    "PISM ERROR: Variable '%s' does not have the time_bounds attribute.\n"
                    "  Cannot use time-dependent forcing data '%s' (%s) without time bounds.\n",
                    dimname.c_str(),  vars[0].get_string("long_name").c_str(), vars[0].short_name.c_str());
        PISMEnd();
      }
    } else {
      // only one time record; set fake time bounds:
      time_bounds.resize(2);
      time_bounds[0] = time[0] - 1;
      time_bounds[1] = time[0] + 1;
    }

  } else {
    // no time dimension; assume that we have only one record and set the time
    // to 0
    time.resize(1);
    time[0] = 0;

    // set fake time bounds:
    time_bounds.resize(2);
    time_bounds[0] = -1;
    time_bounds[1] =  1;
  }
  
  if (!is_increasing(time)) {
    ierr = PetscPrintf(grid->com, "PISM ERROR: times have to be strictly increasing (read from '%s').\n",
		       filename.c_str());
    PISMEnd();
  }

  ierr = get_interp_context(nc, lic); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Read some data to make sure that the interval (my_t, my_t + my_dt) is covered.
PetscErrorCode IceModelVec2T::update(double my_t, double my_dt) {
  PetscErrorCode ierr;
  vector<double>::iterator i, j;
  int m, n, last;

  if (N > 0) {
    last = first + (N - 1);

    // find the interval covered by data held in memory:
    double t0 = time_bounds[first * 2],
      t1 = time_bounds[last * 2 + 1];

    // just return if we have all the data we need:
    if (my_t >= t0 && my_t + my_dt <= t1)
      return 0;
  }

  // i will point to the smallest t so that t >= my_t 
  i = lower_bound(time_bounds.begin(), time_bounds.end(), my_t);

  // j will point to the smallest t so that t >= my_t + my_dt
  j = lower_bound(time_bounds.begin(), time_bounds.end(), my_t + my_dt);

  // some of the the interval (my_t, my_t + my_dt) is outside of the
  // time interval available in the file (on the right)
  // 
  if (i == time_bounds.end()) {
    m = (int)(time.size() - 1);
  } else {
    m = (int)((i - time_bounds.begin() - 1) / 2);
  }

  if (j == time_bounds.end()) {
    n = (int)(time.size() - 1);
  } else {
    n = (int)((j - time_bounds.begin() - 1) / 2);
  }

  // check if all the records necessary to cover this interval fit in the
  // buffer:
  if (n - m + 1 > n_records)
    SETERRQ(grid->com, 1, "IceModelVec2T::update(): timestep is too big");

  ierr = update(m); CHKERRQ(ierr);

  return 0;
}

//! Update by reading at most n_records records from the file.
PetscErrorCode IceModelVec2T::update(int start) {
  PetscErrorCode ierr;
  int time_size = (int)time.size();

  if ((start < 0) || (start >= time_size))
    SETERRQ1(grid->com, 1, "IceModelVec2T::update(int start): start = %d is invalid", start);

  int missing = PetscMin(n_records, time_size - start);

  if (start == first)
    return 0;                   // nothing to do

  int kept = 0;
  if (first >= 0) {
    int last = first + (N - 1);
    if ((N > 0) && (start >= first) && (start <= last)) {
      int discarded = start - first;
      kept = last - start + 1;
      ierr = discard(discarded); CHKERRQ(ierr);
      missing -= kept;
      start += kept;
      first += discarded;
    } else {
      first = start;
    }
  } else {
    first = start;
  }

  if (missing <= 0) return 0;
  
  N = kept + missing;

  if (this->get_n_records() > 1 || getVerbosityLevel() > 4) {
    ierr = verbPrintf(2, grid->com,
                      "  reading \"%s\" into buffer\n"
                      "          (short_name = %s): %d records, time intervals (%s, %s) through (%s, %s)...\n",
                      string_attr("long_name").c_str(), name.c_str(), missing,
                      grid->time->date(time_bounds[start*2]).c_str(),
                      grid->time->date(time_bounds[start*2 + 1]).c_str(),
                      grid->time->date(time_bounds[(start + missing - 1)*2]).c_str(),
                      grid->time->date(time_bounds[(start + missing - 1)*2 + 1]).c_str()); CHKERRQ(ierr);
    report_range = false;
  } else {
    report_range = true;
  }

  PIO nc(*grid, "guess_mode");
  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);

  for (int j = 0; j < missing; ++j) {
    if (lic != NULL) {
      lic->start[0] = start + j;
      lic->report_range = report_range;
    }

    ierr = vars[0].regrid(nc, lic, true, false, 0.0, v); CHKERRQ(ierr);

    ierr = verbPrintf(5, grid->com, " %s: reading entry #%02d, year %f...\n",
		      name.c_str(), start + j, time[start + j]);
    ierr = set_record(kept + j); CHKERRQ(ierr);
  }

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Discard the first N records, shifting the rest of them towards the "beginning".
PetscErrorCode IceModelVec2T::discard(int number) {
  PetscErrorCode ierr;
  PetscScalar **a2, ***a3;

  if (number == 0)
    return 0;

  N -= number;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      for (PetscInt k = 0; k < N; ++k)
	a3[i][j][k] = a3[i][j][k + number];
  ierr = end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);
  
  return 0;
}

//! Sets the record number n to the contents of the (internal) Vec v.
PetscErrorCode IceModelVec2T::set_record(int n) {
  PetscErrorCode ierr;
  PetscScalar **a2, ***a3;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      a3[i][j][n] = a2[i][j];
  ierr = end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! Sets the (internal) Vec v to the contents of the nth record.
PetscErrorCode IceModelVec2T::get_record(int n) {
  PetscErrorCode ierr;
  PetscScalar **a2, ***a3;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i)
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j)
      a2[i][j] = a3[i][j][n];
  ierr = end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Given the time my_t and the current selected time-step my_dt,
//! determines the maximum possible time-step this IceModelVec2T allows.
/*!
  Returns -1 if any time step is OK at my_t.
 */
double IceModelVec2T::max_timestep(double my_t) {
  // only allow going to till the next record
  vector<double>::iterator l = upper_bound(time_bounds.begin(),
                                           time_bounds.end(), my_t);
  if (l != time_bounds.end()) {
    PetscReal tmp = *l - my_t;

    if (tmp > 1)                // never take time-steps shorter than 1 second
      return tmp;
    else if ((l + 1) != time_bounds.end() && (l + 2) != time_bounds.end())
      return *(l + 2) - *l;
    else
      return -1;
  } else
    return -1;

}

//! \brief Get the record that is just before my_t.
PetscErrorCode IceModelVec2T::at_time(double my_t, double period, double reference_time) {
  PetscErrorCode ierr;

  int     last = first + (N - 1);
  double *j,
    *start = &time_bounds[2 * first],
    *end   = &time_bounds[2 * last + 1];

  // "Periodize" time if necessary.
  if (period > 1)
    my_t =  grid->time->mod(my_t - reference_time, period);

  j = lower_bound(start, end, my_t); // binary search

  if (j == end) {
    return get_record(N - 1); // out of range (on the right)
  }

  int i = (int)(j - start);

  if (i == 0) {
    return get_record(0);         // out of range (on the left)
  }

  long int index = ((j - &time_bounds[0]) - 1) / 2;

  ierr = verbPrintf(3, grid->com,
		    "  IceModelVec2T::at_time(%3.5f years): using %s:%s[%d]"
                    " (time bounds: [%3.5f years, %3.5f years]).\n",
                    grid->time->seconds_to_years(my_t),
                    filename.c_str(),
                    vars[0].short_name.c_str(),
                    index,
                    grid->time->seconds_to_years(*(j-1)),
                    grid->time->seconds_to_years(*j)); CHKERRQ(ierr);

  if (i % 2 == 0) {
    PetscPrintf(grid->com,
                "PISM ERROR: time bounds array does not represent continguous time intervals.\n"
                "            (PISM was trying to compute %s at time %s.)\n",
                name.c_str(), grid->time->date(my_t).c_str());
    PISMEnd();
  }

  return get_record((i - 1)/2);
}


/*
 * \brief Use linear interpolation to initialize IceModelVec2T with
 * the value at time `my_t`, assuming that data is periodic with
 * period `period` seconds.
 *
 * \note This method does not check if an update() call is necessary!
 *
 * @param[in] my_t requested time
 * @param[in] period period, in seconds
 * @param[in] reference_time reference time, in seconds
 *
 * @return 0 on success
 */
PetscErrorCode IceModelVec2T::interp(double my_t, double period, double reference_time) {
  PetscErrorCode ierr;

  ierr = init_interpolation(&my_t, 1, period, reference_time); CHKERRQ(ierr);

  PetscScalar **a2, ***a3;

  ierr = get_array(a2); CHKERRQ(ierr);
  ierr = get_array3(a3); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      PetscScalar left = a3[i][j][m_interp_indices[0]],
	right = a3[i][j][m_interp_indices[1]];
      a2[i][j] = left + m_interp_weights[0] * (right - left);
    }
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}


/** 
 * Compute the average value over the time interval `[my_t, my_t + my_dt]`.
 *
 * @param my_t  start of the time interval, in seconds
 * @param my_dt length of the time interval, in seconds
 *
 * @return 0 on success
 */
PetscErrorCode IceModelVec2T::average(double my_t, double my_dt,
				      double period, double reference_time) {
  PetscErrorCode ierr;
  PetscScalar **a2;
  PetscReal dt_years = grid->time->seconds_to_years(my_dt); // *not* time->year(my_dt)

  // Determine the number of small time-steps to use for averaging:
  int M = (int) ceil(52 * (dt_years) + 1); // (52 weeks in a year)
  if (M < 2) M = 2;

  vector<double> ts(M);
  double dt = my_dt / (M - 1);
  for (int k = 0; k < M; k++)
    ts[k] = my_t + k * dt;

  ierr = init_interpolation(&ts[0], M, period, reference_time); CHKERRQ(ierr);

  ierr = get_array(a2);         // calls begin_access()
  for (PetscInt   i = grid->xs; i < grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j < grid->ys+grid->ym; ++j) {
      ierr = average(i, j, a2[i][j]); CHKERRQ(ierr);
    }
  }
  ierr = end_access(); CHKERRQ(ierr);

  return 0;
}

/** 
 * \brief Compute weights for the piecewise-linear interpolation. This is used *both*
 * for time-series and "snapshots".
 *
 * @param ts requested times, in seconds
 * @param ts_length number of requested times (length of the `ts` array)
 * @param period period, in seconds
 * @param indices indices, array of `2*ts_length` integers
 * @param weights interpolation weights, array of `ts_length` doubles
 *
 * @return 0 on success
 */
PetscErrorCode IceModelVec2T::init_interpolation(PetscScalar *ts, unsigned int ts_length,
						 double period, double reference_time) {
  int mcurr = first;
  int last = first + (N - 1);

  m_interp_indices.resize(2*ts_length);
  m_interp_weights.resize(ts_length);

  // Compute "periodized" times if necessary.
  vector<double> times_requested(ts_length);
  if (period > 1) {
    for (unsigned int k = 0; k < ts_length; ++k)
      times_requested[k] = grid->time->mod(ts[k] - reference_time, period);
  } else {
    for (unsigned int k = 0; k < ts_length; ++k)
      times_requested[k] = ts[k];
  }
  
  for (unsigned int k = 0; k < ts_length; ++k) {

    if (k > 0 && times_requested[k] < times_requested[k-1])
      mcurr = first; // reset the mcurr index: times_requested are not increasing!

    // extrapolate on the right:
    if (times_requested[k] >= time[last]) {
      m_interp_indices[2*k + 0] = N - 1;
      m_interp_indices[2*k + 1] = N - 1;
      m_interp_weights[k]       = 0;
      continue;
    }

    int left, right;		// indices of the left and right
				// end-points of the current interval
    PetscReal lambda = 0.0,
      epsilon = 1.0;            // 1 second

    // handle periodicity on the left (using the fact that 'time'
    // contains right end-points of time bounds
    if (times_requested[k] < time[first]) {
      right  = first;
      left   = first;
      lambda = 0.0;

      if (period > 1.0) {	// we do not support periods shorter than 1 second
	while (left < last && time[left] < period - epsilon)
	  left++;

	// periodic data has to have time axis starting from 0
	//
	// Here we make errors in time of up to 1 second to protect from
	// division by zero
	if (time[right] > 1)	
	  lambda = (times_requested[k] - 0) / (time[right] - 0);
      }

    } else {

      while (time[mcurr+1] < times_requested[k]) {
	mcurr++;
      }

      left   = mcurr;
      right  = mcurr + 1;
      lambda = (times_requested[k] - time[left]) / (time[right] - time[left]);
    }

    m_interp_indices[2*k + 0] = left - first;
    m_interp_indices[2*k + 1] = right - first;
    m_interp_weights[k]       = lambda;
  }

  return 0;
}

/** 
 * \brief Compute values of the time-series using precomputed indices
 * and interpolation weights (linear interpolation).
 *
 * @param i,j map-plane grid point
 * @param indices pre-computed indices
 * @param weights pre-computed interpolation weights
 * @param values pointer to an allocated array of `weights.size()` `PetscScalar`
 *
 * @return 0 on success
 */
PetscErrorCode IceModelVec2T::interp(int i, int j,
				     PetscScalar *result) {
  PetscScalar ***a3 = (PetscScalar***) array3;
  unsigned int ts_length = m_interp_weights.size();

  for (unsigned int k = 0; k < ts_length; ++k) {
    PetscScalar left = a3[i][j][m_interp_indices[2*k + 0]],
      right	     = a3[i][j][m_interp_indices[2*k + 1]];

    result[k] = left + m_interp_weights[k] * (right - left);
  }
  
  return 0;
}

//! \brief Finds the average value at i,j over the interval (my_t, my_t +
//! my_dt) using trapezoidal rule.
/*!
  Can (and should) be optimized. Later, though.
 */
PetscErrorCode IceModelVec2T::average(int i, int j,
				      double &result) {
  PetscErrorCode ierr;
  unsigned int M = m_interp_weights.size();
  vector<PetscScalar> values(M);

  ierr = interp(i, j, &values[0]); CHKERRQ(ierr);

  // trapezoidal rule; uses the fact that all 'small' time intervals used here
  // are the same:
  result = 0;
  for (unsigned int k = 0; k < M - 1; ++k)
    result += values[k] + values[k + 1];
  result /= 2*(M - 1);

  return 0;
}


