/* Copyright (C) 2018, 2020 PISM Authors
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

#ifndef PICOGEOMETRY_H
#define PICOGEOMETRY_H

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {

class IceModelVec2S;
class IceModelVec2CellType;

namespace ocean {

void eikonal_equation(IceModelVec2Int &mask);

/*!
 * This class isolates geometric computations performed by the PICO ocean model.
 */
class PicoGeometry : public Component {
public:
  PicoGeometry(IceGrid::ConstPtr grid);
  virtual ~PicoGeometry();

  void init_impl();
  void update(const IceModelVec2S &bed_elevation, const IceModelVec2CellType &cell_type);

  const IceModelVec2Int &continental_shelf_mask() const;
  const IceModelVec2Int &box_mask() const;
  const IceModelVec2Int &ice_shelf_mask() const;
  const IceModelVec2Int &ice_rise_mask() const;
  const IceModelVec2Int &basin_mask() const;

  enum IceRiseMask { OCEAN = 0, RISE = 1, CONTINENTAL = 2, FLOATING = 3 };

private:
  void compute_ice_rises(const IceModelVec2CellType &cell_type, bool exclude_ice_rises, IceModelVec2Int &result);
  void compute_lakes(const IceModelVec2CellType &cell_type, IceModelVec2Int &result);
  void compute_ocean_mask(const IceModelVec2CellType &cell_type, IceModelVec2Int &result);
  void compute_continental_shelf_mask(const IceModelVec2S &bed_elevation, const IceModelVec2Int &ice_rises_mask,
                                      double bed_elevation_threshold, IceModelVec2Int &result);
  void compute_ice_shelf_mask(const IceModelVec2Int &ice_rises_mask, const IceModelVec2Int &lake_mask,
                              IceModelVec2Int &result);

  void split_ice_shelves(const IceModelVec2CellType &cell_type,
                         const IceModelVec2Int &basin_mask,
                         //const std::vector< std::vector <int> > m_n_basin_neighbors,
                         const IceModelVec2Int &shelf_mask);

  void get_basin_neighbors(const IceModelVec2CellType &cell_type,
                           const IceModelVec2Int &basin_mask);
                           //const IceModelVec2Int &basin_mask,
                           //std::vector< std::vector <int> > result);
 
  void compute_distances_cf(const IceModelVec2Int &ocean_mask, const IceModelVec2Int &ice_rises, bool exclude_ice_rises,
                            IceModelVec2Int &dist_cf);
  void compute_distances_gl(const IceModelVec2Int &ocean_mask, const IceModelVec2Int &ice_rises, bool exclude_ice_rises,
                            IceModelVec2Int &dist_gl);

  void compute_box_mask(const IceModelVec2Int &D_gl, const IceModelVec2Int &D_cf, const IceModelVec2Int &shelf_mask,
                        int n_boxes, IceModelVec2Int &result);

  void label_tmp();
  void relabel_by_size(IceModelVec2Int &mask);

  // storage for outputs
  IceModelVec2Int m_continental_shelf;
  IceModelVec2Int m_boxes;
  IceModelVec2Int m_ice_shelves;
  IceModelVec2Int m_basin_mask;

  int m_n_basins, m_n_shelves;

  // storage for intermediate fields
  IceModelVec2Int m_distance_gl;
  IceModelVec2Int m_distance_cf;
  IceModelVec2Int m_ocean_mask;
  IceModelVec2Int m_lake_mask;
  IceModelVec2Int m_ice_rises;

  // temporary storage
  IceModelVec2Int m_tmp;
  std::shared_ptr<petsc::Vec> m_tmp_p0;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* PICOGEOMETRY_H */
