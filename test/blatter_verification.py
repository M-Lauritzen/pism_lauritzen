#!/usr/bin/env python3
"""This script runs verification tests of the Blatter stress balance solver that use
manufactured solutions from

Tezaur et al (2015) Albany/FELIX: a parallel, scalable and robust, finite element,
first-order Stokes approximation ice sheet solver built for advanced analysis
(doi:10.5194/gmd-8-1197-2015)

"""

from unittest import TestCase

import PISM
import PISM.util

import numpy as np

ctx = PISM.Context()
config = ctx.config

# do not ignore thin ice
config.set_number("geometry.ice_free_thickness_standard", 0.0)
config.set_number("stress_balance.ice_free_thickness_standard", 0.0)

# Set flow law parameters
config.set_string("stress_balance.blatter.flow_law", "isothermal_glen")
config.set_number("stress_balance.blatter.Glen_exponent", 3.0)

# Set constants:
config.set_number("constants.ice.density", 910.0)
config.set_number("constants.standard_gravity", 9.81)

def expt(xs, ys):
    "Compute the convergence rate using a polynomial fit."
    return -np.polyfit(np.log(xs), np.log(ys), 1)[0]

class TestXY(TestCase):
    """2D (x-y) verification test using a manufactured solution.

    u = exp(x) * sin(2 * pi * y)
    v = exp(x) * cos(2 * pi * y)

    on [0, 1] * [0, 1] with Dirichlet BC at x = +- 1, y = +- 1.

    Flat bed, no basal drag, constant ice thickness, isothermal Glen flow law with n == 3.

    See section 4.1 in Tezaur et al.

    The source term needed for the chosen manufactured solution is computed (and the code
    is generated) using SymPy.

    """
    def setUp(self):
        "Set PETSc options"
        opt = PISM.PETSc.Options()

        self.opt = opt

        self.A_old = config.get_number("flow_law.isothermal_Glen.ice_softness")
        config.set_number("flow_law.isothermal_Glen.ice_softness",
                          PISM.util.convert(1e-4, "Pa-3 year-1", "Pa-3 s-1"))

        # the default preconditioner is seems to be ineffective
        opt.setValue("-pc_type", "gamg")
        opt.setValue("-snes_monitor", "")

    def tearDown(self):
        "Clear PETSc options"
        config.set_number("flow_law.isothermal_Glen.ice_softness", self.A_old)

        self.opt.delValue("-pc_type")
        self.opt.delValue("-snes_monitor")

    def inputs(self, N):
        "Allocate stress balance inputs for a given grid size"

        P = PISM.GridParameters(config)

        # Domain: [0,1] * [0, 1] * [0, 1]
        P.Lx = 0.5
        P.Ly = 0.5
        P.x0 = 0.5
        P.y0 = 0.5
        # This vertical grid is used for the enthalpy field only *and* we use an
        # isothermal flow law, so 2 levels is enough.
        P.z = PISM.DoubleVector([0.0, 1.0])
        P.registration = PISM.CELL_CORNER
        P.Mx = int(N)
        P.My = int(N)
        P.ownership_ranges_from_options(ctx.size)

        grid = PISM.IceGrid(ctx.ctx, P)

        geometry = PISM.Geometry(grid)

        enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        # initialize enthalpy (the value used here is irrelevant)
        enthalpy.set(1e5)

        yield_stress = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)
        yield_stress.set(0.0)

        geometry.bed_elevation.set(0.0)
        geometry.ice_thickness.set(1.0)
        geometry.sea_level_elevation.set(0.0)
        geometry.ensure_consistency(0.0)

        return geometry, enthalpy, yield_stress

    def exact_solution(self, grid):
        "Returns an array with the exact solution"
        exact = PISM.IceModelVec2V(grid, "exact", PISM.WITHOUT_GHOSTS)

        with PISM.vec.Access(exact):
            for (i, j) in grid.points():
                x = grid.x(i)
                y = grid.y(j)
                exact[i, j] = PISM.blatter_xy_exact(x, y)

        return exact

    def error_norm(self, N):
        "Return the infinity norm of errors for the u and v components (as a tuple)."
        geometry, enthalpy, yield_stress = self.inputs(N)

        grid = enthalpy.grid()

        u_exact = self.exact_solution(grid)

        # no variation in the Z direction, so it's OK to use 2 vertical levels
        blatter_Mz = 2
        # no coarsening in the Z direction
        n_levels = 0
        coarsening_factor = 1

        model = PISM.BlatterTestXY(grid, blatter_Mz, n_levels, coarsening_factor)

        model.init()

        inputs = PISM.StressBalanceInputs()

        inputs.geometry = geometry
        inputs.basal_yield_stress = yield_stress
        inputs.enthalpy = enthalpy

        # run the solver
        model.update(inputs, True)

        # compute the error
        error = PISM.IceModelVec2V(grid, "error", PISM.WITHOUT_GHOSTS)
        error.copy_from(u_exact)
        error.add(-1.0, model.velocity())

        return error.norm_all(PISM.PETSc.NormType.NORM_INFINITY)

    def test(self):
        "Test that the convergence rate for the XY test is at least quadratic"
        Ns = [11, 21]
        norms = [self.error_norm(N) for N in Ns]

        norms_u = [n[0] for n in norms]
        norms_v = [n[1] for n in norms]

        # Compute the exponent for the convergence rate
        expt_u = expt(Ns, norms_u)
        expt_v = expt(Ns, norms_v)

        print("U component conv. rate: dx^{}".format(expt_u))
        print("V component conv. rate: dx^{}".format(expt_v))

        # The convergence rate should be close to quadratic.
        assert expt_u >= 2.0
        assert expt_v >= 2.0

    def plot(self):
        Ns = [11, 21, 41, 81]
        try:
            self.setUp()
            norms = [self.error_norm(N) for N in Ns]
        finally:
            self.tearDown()

        # the domain is [0, 1]*[0, 1]*[0, 1]
        dxs = 1.0 / (np.array(Ns) - 1)

        norms_u = [x[0] for x in norms]
        norms_v = [x[1] for x in norms]

        p_u = np.polyfit(np.log(dxs), np.log(norms_u), 1)
        fit_u = np.exp(np.polyval(p_u, np.log(dxs)))

        p_v = np.polyfit(np.log(dxs), np.log(norms_v), 1)
        fit_v = np.exp(np.polyval(p_v, np.log(dxs)))

        from bokeh.plotting import figure, show, output_file
        from bokeh.layouts import gridplot

        output_file("test-xy.html")
        f = figure(title = "Blatter-Pattyn stress balance: verification test XY, x-component",
                   x_axis_type="log", y_axis_type="log", x_axis_label="dx", y_axis_label="max error")
        f.scatter(dxs, norms_u)
        f.line(dxs, norms_u, legend_label="max error", line_width=2)
        f.line(dxs, fit_u, line_dash="dashed", line_color="red", line_width=2,
               legend_label="O(dx^{})".format(p_u[0]))

        g = figure(title = "Blatter-Pattyn stress balance: verification test XY, y-component",
                   x_axis_type="log", y_axis_type="log", x_axis_label="dx", y_axis_label="max error")
        g.scatter(dxs, norms_v)
        g.line(dxs, norms_v, legend_label="max error", line_width=2)
        g.line(dxs, fit_v, line_dash="dashed", line_color="red", line_width=2,
               legend_label="O(dx^{})".format(p_v[0]))

        gp = gridplot([[f, g]])

        show(gp)

class TestXZ(TestCase):
    """2D (x-z) verification test using a manufactured solution.

    u = SIA + SSA exact solutions
    v = 0

    With sliding and a constant drag (beta), constant ice thickness, parabolic top surface
    and bed elevations, isothermal Glen flow law with n == 3.

    Domain: x in [-50km, 50km], ice thickness of 1km. The domain extent in the Y direction
    is chosen to set dy = dx with My == 3. This avoids possible element aspect ratio issues.

    See section 4.2 in Tezaur et al.

    - Uses Dirichlet BC at x == +-50km (unlike Tezaur et al who use a stress BC derived
      from the manufactured solution).

    - Uses basal BC which is a sum of the sliding BC and a correction computed using the
      manufactured solution.

    - Uses a stress BC at the top surface derived from the manufactured solution.

    All source terms and BC formulas are from Tezaur et al. (The C code is generated using
    SymPy.)

    This verification test uses the code implementing the basal boundary condition and so
    allows us to check its correctness (at least when beta is a constant).

    """
    def setUp(self):
        "Set PETSc options"

        self.opt = PISM.PETSc.Options()

        self.opts = {"-snes_monitor": "",
                     "-pc_type": "mg"
                     }

        for k, v in self.opts.items():
            self.opt.setValue(k, v)

        # the magnitude of the regularization constant affects the accuracy near the
        # "dome"
        config.set_number("flow_law.Schoof_regularizing_velocity", 1e-5)

        # Set sliding law parameters to make "tauc" equivalent to "beta"
        config.set_flag("basal_resistance.pseudo_plastic.enabled", True)
        config.set_number("basal_resistance.pseudo_plastic.q", 1.0)
        config.set_number("basal_resistance.pseudo_plastic.u_threshold",
                          PISM.util.convert(1.0, "m / s", "m / year"))

        # save old ice softness
        self.A_old = config.get_number("flow_law.isothermal_Glen.ice_softness")

        # set ice softness
        config.set_number("flow_law.isothermal_Glen.ice_softness",
                          PISM.util.convert(1e-16, "Pa-3 year-1", "Pa-3 s-1"))

        self.A = config.get_number("flow_law.isothermal_Glen.ice_softness")

        self.s0    = 2000.0     # m
        self.alpha = 4e-8       # 1/m
        self.H     = 1000.0     # m
        self.beta  = PISM.util.convert(1.0, "kPa year m-1", "Pa s m-1")

    def tearDown(self):
        "Clear PETSc options and configuration parameters"

        config.set_number("flow_law.isothermal_Glen.ice_softness", self.A_old)
        config.set_flag("basal_resistance.pseudo_plastic.enabled", False)
        # restore the default value
        config.set_number("flow_law.Schoof_regularizing_velocity", 1.0)

        for k, v in self.opts.items():
            self.opt.delValue(k)

    def inputs(self, N):
        P = PISM.GridParameters(config)

        # Domain: [-50e3, 50e3] * [-dx, dx] * [0, 1000]
        Lx = 50e3
        Mx = N
        dx = (2 * Lx) / (Mx - 1)

        P.Lx = 50e3
        P.Mx = int(N)
        P.x0 = 0.0

        P.Ly = dx
        P.My = 3
        P.y0 = 0.0

        # this vertical grid is used to store ice enthalpy, not ice velocity, so 2 levels
        # is enough
        P.z = PISM.DoubleVector([0.0, self.H])
        P.registration = PISM.CELL_CORNER
        P.periodicity = PISM.Y_PERIODIC
        P.ownership_ranges_from_options(ctx.size)

        grid = PISM.IceGrid(ctx.ctx, P)

        geometry = PISM.Geometry(grid)

        enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        # initialize enthalpy (the value used here is irrelevant)
        enthalpy.set(1e5)

        yield_stress = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)

        yield_stress.set(self.beta)

        with PISM.vec.Access(geometry.bed_elevation):
            for (i, j) in grid.points():
                x = grid.x(i)
                geometry.bed_elevation[i, j] = self.s0 - self.H - self.alpha * x**2

        geometry.ice_thickness.set(self.H)

        # ensure that all ice is grounded
        geometry.sea_level_elevation.copy_from(geometry.bed_elevation)
        geometry.sea_level_elevation.shift(-1.0)

        geometry.ensure_consistency(0.0)

        return geometry, enthalpy, yield_stress

    def exact_solution(self, grid, bed, Z):
        "Returns an array with the exact solution"
        exact = PISM.IceModelVec3(grid, "exact", PISM.WITHOUT_GHOSTS, Z)
        exact.set_attrs("exact", "x-component of the exact solution", "m / s", "m / year", "", 0)

        rho = config.get_number("constants.ice.density")
        g = config.get_number("constants.standard_gravity")

        u = np.zeros_like(Z)
        with PISM.vec.Access([bed, exact]):
            for (i, j) in grid.points():
                x = grid.x(i)
                for k, z_sigma in enumerate(Z):
                    z = bed[i, j] + self.H * z_sigma
                    u[k] = PISM.blatter_xz_exact(x, z, self.A, rho, g,
                                                 self.s0, self.alpha, self.H, self.beta).u
                exact.set_column(i, j, u)

        return exact

    def error_norm(self, N, n_mg):
        "Return the infinity norm of errors for the u component."
        geometry, enthalpy, yield_stress = self.inputs(N)

        # set the number of multigrid levels
        self.opt.setValue("-pc_mg_levels", n_mg)

        grid = enthalpy.grid()

        blatter_Mz = N
        # do not pad the vertical grid
        n_levels = 0
        coarsening_factor = 4

        model = PISM.BlatterTestXZ(grid, blatter_Mz, n_levels, coarsening_factor)

        model.init()

        inputs = PISM.StressBalanceInputs()

        inputs.geometry = geometry
        inputs.basal_yield_stress = yield_stress
        inputs.enthalpy = enthalpy

        # run the solver
        model.update(inputs, True)

        u_model_z = model.velocity_u_sigma().levels()

        u_model = PISM.IceModelVec3(grid, "u_model", PISM.WITHOUT_GHOSTS, u_model_z)
        u_model.set_attrs("model", "modeled velocity", "m / s", "m / year", "", 0)
        u_model.copy_from(model.velocity_u_sigma())

        u_exact = self.exact_solution(grid, geometry.bed_elevation, u_model_z)

        # compute the error
        u_error = PISM.IceModelVec3(grid, "error", PISM.WITHOUT_GHOSTS, u_model_z)
        u_error.copy_from(u_exact)
        u_error.add(-1.0, model.velocity_u_sigma())

        return u_error.norm(PISM.PETSc.NormType.NORM_INFINITY)

    def test(self):
        "Test that the convergence rate for the XZ test is at least quadratic"

        # refinement path
        Ns = [11, 21]
        # number of MG levels to use for a particular grid size, assuming the coarsening
        # factor of 4
        mg_levels = [1, 2]

        norms = [self.error_norm(N, n_mg) for (N, n_mg) in zip(Ns, mg_levels)]

        # Compute the exponent for the convergence rate
        expt_u = expt(Ns, norms)

        print("U component conv. rate: dx^{}".format(expt_u))

        # The convergence rate should be close to quadratic.
        assert expt_u >= 2.0

    def plot(self):
        Ns = [11, 21, 41, 81]
        mg_levels = [1, 2, 3, 4]

        try:
            self.setUp()
            norms = [self.error_norm(N, n_mg) for (N, n_mg) in zip(Ns, mg_levels)]
        finally:
            self.tearDown()

        # the domain is 100km long and 1km high
        dxs = 100e3 / (np.array(Ns) - 1)

        p = np.polyfit(np.log(dxs), np.log(norms), 1)
        fit = np.exp(np.polyval(p, np.log(dxs)))

        from bokeh.plotting import figure, show, output_file

        output_file("test-xz.html")
        f = figure(title = "Blatter-Pattyn stress balance: verification test XZ",
                   x_axis_type="log", y_axis_type="log", x_axis_label="dx", y_axis_label="max error")
        f.scatter(dxs, norms)
        f.line(dxs, norms, legend_label="max error", line_width=2)
        f.line(dxs, fit, line_dash="dashed", line_color="red", line_width=2,
               legend_label="O(dx^{})".format(p[0]))


        show(f)

class TestCFBC(TestCase):
    """Constant viscosity 2D (x-z) verification test checking the implementation of CFBC.

    u = 1/4 * x * z * g * (rho_w - rho_i)
    v = 0

    on a square 0 <= x <= 1, -1 <= z <= 0 with sea level at z = 0 (fully submerged).

    Dirichet BC at x = 0, periodic in the Y direction, lateral BC at x = 1. No basal drag.

    """
    def setUp(self):
        "Set PETSc options"

        self.H = 1e3

        self.opt = PISM.PETSc.Options()

        self.opts = {"-snes_monitor": "",
                     "-ksp_monitor": "",
                     "-pc_type": "mg"
                     }

        for k, v in self.opts.items():
            self.opt.setValue(k, v)

        # the magnitude of the regularization constant affects the accuracy near the
        # "dome"
        config.set_number("flow_law.Schoof_regularizing_velocity", 1e-5)

        # constant viscocity
        n = 1.0
        config.set_number("stress_balance.blatter.Glen_exponent", n)

        config.set_number("flow_law.isothermal_Glen.ice_softness",
                          PISM.util.convert(1e-3, "Pa-3 year-1", "Pa-3 s-1"))
        A = config.get_number("flow_law.isothermal_Glen.ice_softness")
        self.B = 1.0 / A              # note that n = 1

    def tearDown(self):
        "Clear PETSc options and configuration parameters"

        # restore the default Glen exponent
        config.set_number("stress_balance.blatter.Glen_exponent", 3.0)

        for k, v in self.opts.items():
            self.opt.delValue(k)

    def inputs(self, N):
        P = PISM.GridParameters(config)

        # Domain: [0, 1] * [-dx, dx] * [-1, 0]
        Lx = 0.5
        Mx = N
        dx = (2 * Lx) / (Mx - 1)

        P.Lx = Lx
        P.Mx = int(N)
        P.x0 = Lx

        P.Ly = dx
        P.My = 3
        P.y0 = 0.0

        # this vertical grid is used to store ice enthalpy, not ice velocity, so 2 levels
        # is enough
        P.z = PISM.DoubleVector([0.0, self.H])
        P.registration = PISM.CELL_CORNER
        P.periodicity = PISM.Y_PERIODIC
        P.ownership_ranges_from_options(ctx.size)

        grid = PISM.IceGrid(ctx.ctx, P)

        geometry = PISM.Geometry(grid)

        enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        # initialize enthalpy (the value used here is irrelevant)
        enthalpy.set(1e5)

        yield_stress = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)

        yield_stress.set(0.0)

        geometry.bed_elevation.set(-self.H)
        geometry.ice_thickness.set(self.H)
        geometry.ice_surface_elevation.set(0.0)
        geometry.cell_type.set(PISM.MASK_FLOATING)
        geometry.sea_level_elevation.set(0.0)

        # do *not* call geometry.ensure_consistency(): we want to keep surface elevation
        # at the sea levels

        return geometry, enthalpy, yield_stress

    def exact_solution(self, grid, bed, Z):
        "Returns an array with the exact solution"
        exact = PISM.IceModelVec3(grid, "exact", PISM.WITHOUT_GHOSTS, Z)
        exact.set_attrs("exact", "x-component of the exact solution", "m / s", "m / year", "", 0)

        rho_i = config.get_number("constants.ice.density")
        rho_w = config.get_number("constants.sea_water.density")
        g = config.get_number("constants.standard_gravity")

        u = np.zeros_like(Z)
        with PISM.vec.Access([bed, exact]):
            for (i, j) in grid.points():
                x = grid.x(i)
                for k, z_sigma in enumerate(Z):
                    z = bed[i, j] + self.H * z_sigma
                    u[k] = U = 1 / (4 * self.B) * x**2 * z * g * (rho_w - rho_i)
                exact.set_column(i, j, u)

        return exact

    def error_norm(self, N, n_mg):
        "Return the infinity norm of errors for the u component."
        geometry, enthalpy, yield_stress = self.inputs(N)

        # set the number of multigrid levels
        self.opt.setValue("-pc_mg_levels", n_mg)

        grid = enthalpy.grid()

        blatter_Mz = N
        # do not pad the vertical grid
        n_levels = 0
        coarsening_factor = 4

        model = PISM.BlatterTestCFBC(grid, blatter_Mz, n_levels, coarsening_factor)

        model.init()

        inputs = PISM.StressBalanceInputs()

        inputs.geometry = geometry
        inputs.basal_yield_stress = yield_stress
        inputs.enthalpy = enthalpy

        # run the solver
        model.update(inputs, True)

        u_model_z = model.velocity_u_sigma().levels()

        u_model = PISM.IceModelVec3(grid, "u_model", PISM.WITHOUT_GHOSTS, u_model_z)
        u_model.set_attrs("model", "modeled velocity", "m / s", "m / year", "", 0)
        u_model.copy_from(model.velocity_u_sigma())

        u_model.dump("vel-{}-model.nc".format(N))

        u_exact = self.exact_solution(grid, geometry.bed_elevation, u_model_z)

        u_exact.dump("vel-{}-exact.nc".format(N))

        # compute the error
        u_error = PISM.IceModelVec3(grid, "error", PISM.WITHOUT_GHOSTS, u_model_z)
        u_error.copy_from(u_exact)
        u_error.add(-1.0, model.velocity_u_sigma())

        u_error.dump("vel-{}-error.nc".format(N))

        return u_error.norm(PISM.PETSc.NormType.NORM_INFINITY)

    def test(self):
        "Test that the convergence rate for the XZ test is at least quadratic"

        # refinement path
        Ns = [11, 21]
        # number of MG levels to use for a particular grid size, assuming the coarsening
        # factor of 4
        mg_levels = [1, 2]

        norms = [self.error_norm(N, n_mg) for (N, n_mg) in zip(Ns, mg_levels)]

        # Compute the exponent for the convergence rate
        expt_u = expt(Ns, norms)

        print("U component conv. rate: dx^{}".format(expt_u))

        # The convergence rate should be close to quadratic.
        assert expt_u >= 2.0

    def plot(self):
        Ns = [11, 21, 41, 81]
        mg_levels = [1, 2, 3, 4]

        try:
            self.setUp()
            norms = [self.error_norm(N, n_mg) for (N, n_mg) in zip(Ns, mg_levels)]
        finally:
            self.tearDown()

        # the domain is 100km long and 1km high
        dxs = 100e3 / (np.array(Ns) - 1)

        p = np.polyfit(np.log(dxs), np.log(norms), 1)
        fit = np.exp(np.polyval(p, np.log(dxs)))

        from bokeh.plotting import figure, show, output_file

        output_file("test-xz.html")
        f = figure(title = "Blatter-Pattyn stress balance: verification test XZ",
                   x_axis_type="log", y_axis_type="log", x_axis_label="dx", y_axis_label="max error")
        f.scatter(dxs, norms)
        f.line(dxs, norms, legend_label="max error", line_width=2)
        f.line(dxs, fit, line_dash="dashed", line_color="red", line_width=2,
               legend_label="O(dx^{})".format(p[0]))


        show(f)

if __name__ == "__main__":

    # TestXY().plot()
    # TestXZ().plot()

    t = TestCFBC()
    t.setUp()
    print(t.error_norm(41, 2))
    t.tearDown()
