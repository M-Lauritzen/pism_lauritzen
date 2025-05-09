/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2023 PISM Authors
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

#include "pism/frontretreat/util/IcebergRemover.hh"
#include "pism/util/connected_components/label_components.hh"
#include "pism/util/Mask.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {
namespace calving {

IcebergRemover::IcebergRemover(std::shared_ptr<const Grid> g)
    : Component(g), m_iceberg_mask(m_grid, "iceberg_mask") {
}

/**
 * Use PISM's ice cover mask to update ice thickness, removing "icebergs".
 *
 * @param[in,out] pism_mask PISM's ice cover mask
 * @param[in,out] ice_thickness ice thickness
 */
void IcebergRemover::update(const array::Scalar &bc_mask,
                            array::CellType1 &cell_type,
                            array::Scalar &ice_thickness) {
  update_impl(bc_mask, cell_type, ice_thickness);
}

void IcebergRemover::update_impl(const array::Scalar &bc_mask,
                                 array::CellType1 &cell_type,
                                 array::Scalar &ice_thickness) {
  const int
    mask_grounded_ice = 1,
    mask_floating_ice = 2;

  // prepare the mask that will be handed to the connected component
  // labeling code:
  {
    m_iceberg_mask.set(0.0);

    array::AccessScope list{&cell_type, &m_iceberg_mask, &bc_mask};

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.grounded_ice(i,j)) {
        m_iceberg_mask(i,j) = mask_grounded_ice;
      } else if (cell_type.floating_ice(i,j)) {
        m_iceberg_mask(i,j) = mask_floating_ice;
      }
    }

    // Mark icy Dirichlet B.C. cells as "grounded" because we don't want them removed.
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (bc_mask(i, j) > 0.5 and cell_type.icy(i, j)) {
        m_iceberg_mask(i, j) = mask_grounded_ice;
      }
    }
  }

  connected_components::label_isolated(m_iceberg_mask, mask_grounded_ice);

  // correct ice thickness and the cell type mask using the resulting
  // "iceberg" mask:
  {
    array::AccessScope list{&ice_thickness, &cell_type, &m_iceberg_mask, &bc_mask};

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_iceberg_mask(i,j) > 0.5 && bc_mask(i,j) < 0.5) {
        ice_thickness(i,j) = 0.0;
        cell_type(i,j)     = MASK_ICE_FREE_OCEAN;
      }
    }
  }

  // update ghosts of the cell_type and the ice thickness (then surface
  // elevation can be updated redundantly)
  cell_type.update_ghosts();
  ice_thickness.update_ghosts();
}

} // end of namespace calving
} // end of namespace pism
