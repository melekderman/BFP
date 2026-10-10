"""Module for coefficients"""

import mfem.ser as mfem
import numpy as np

__all__ = [
    "gfp2_coefficients",
    "TotalXSCoefficientE",
    "ScatteringXSCoefficientE",
    "StoppingPowerCoefficientE",
    "QCoefficientE",
    "QFuncCoefficient",
    "EDependentCoefficient",
    "XDependentCoefficient",
    "InflowCoefficient",
    "ConstantCoefficient",
    "VectorConstCoefficient",
]


def gfp2_coefficients(sigma_a1, sigma_a2):
    """Convert two elastic angular relaxation cross sections to GFP2 alpha and beta.

    ``sigma_a1`` and ``sigma_a2`` are the first and second relaxation cross
    sections (inverse cm), not absorption. Scalars or broadcastable arrays
    must satisfy ``0 < sigma_a1 < sigma_a2 <= 3*sigma_a1`` and be finite.
    Return ``(alpha, beta)`` as NumPy scalars/arrays: alpha has inverse-length
    units and beta is dimensionless. Invalid values raise ValueError.
    ``sigma_a2 == 3*sigma_a1`` gives the FP limit beta=0.
    """
    sigma_a1 = np.asarray(sigma_a1, dtype=float)
    sigma_a2 = np.asarray(sigma_a2, dtype=float)
    if (
        not np.all(np.isfinite(sigma_a1))
        or not np.all(np.isfinite(sigma_a2))
        or np.any(sigma_a1 <= 0.0)
        or np.any(sigma_a2 <= sigma_a1)
        or np.any(sigma_a2 > 3.0 * sigma_a1)
    ):
        raise ValueError("GFP2 requires 0 < sigma_a1 < sigma_a2 <= 3*sigma_a1.")
    beta = (3.0 * sigma_a1 - sigma_a2) / (6.0 * (sigma_a2 - sigma_a1))
    alpha = 0.5 * sigma_a1 * (1.0 + 2.0 * beta)
    if not np.all(np.isfinite(alpha)) or not np.all(np.isfinite(beta)):
        raise ValueError("GFP2 coefficients must be finite.")
    return alpha, beta


class TotalXSCoefficientE(mfem.PyCoefficient):
    """Energy-dependent macroscopic removal coefficient in inverse cm.

    Args:
        xs_t_data: Scalar constant or values on uniform energy intervals,
            ordered from high to low energy.
        E_start, E_end (float): High and low energy bounds in MeV, with
            E_start > E_end for grouped data; ignored for scalar data.

    Evaluation uses physical energy ``x[1]``. Group ordering is a data-lookup
    convention, not a transformed mesh coordinate. Values close to an internal
    edge use the lower-energy group; unmatched energies return the last value,
    including values outside the supplied range. Supply valid bounds/data;
    this helper does not validate them. In CSDA this is removal, not an
    angular scattering gain. For FP/GFP2 supply genuine absorption separately.
    """

    def __init__(self, xs_t_data, E_start, E_end):
        super(TotalXSCoefficientE, self).__init__()
        self.E_start = E_start
        self.E_end = E_end
        if isinstance(xs_t_data, (int, float)):
            self.constant = True
            self.constant_value = float(xs_t_data)
        else:
            self.constant = False
            self.xs_t_data = xs_t_data
            self.n_groups = len(xs_t_data)
            self.E_bins = np.linspace(E_start, E_end, self.n_groups + 1)

    def EvalValue(self, x):
        """Evaluate the total cross-section at the given energy point.

        Args:
            x (list or array-like): Input point, where x[1] is the energy value.

        Returns:
            float: Corresponding total cross-section value.
        """
        if self.constant:
            return self.constant_value
        E = x[1]
        for i in range(self.n_groups - 1):
            if np.isclose(E, self.E_bins[i + 1]):
                return float(self.xs_t_data[i + 1])
            if E <= self.E_bins[i] and E > self.E_bins[i + 1]:
                return float(self.xs_t_data[i])
        return float(self.xs_t_data[-1])


class ScatteringXSCoefficientE(mfem.PyCoefficient):
    """Energy-group total elastic scattering cross section in inverse cm.

    ``xs_s_data`` is a scalar or values ordered from high to low energy.
    ``E_start``/``E_end`` specify the high/low bounds in MeV. Lookup and scalar
    handling follow :class:`TotalXSCoefficientE`. This helper is not automatically
    coupled to any solver. A total scattering cross section is not the elastic
    transport cross section required by :func:`~bfp.solver.Solve_FP`.
    """

    def __init__(self, xs_s_data, E_start, E_end):
        super(ScatteringXSCoefficientE, self).__init__()
        self.E_start = E_start
        self.E_end = E_end
        if isinstance(xs_s_data, (int, float)):
            self.constant = True
            self.constant_value = float(xs_s_data)
        else:
            self.constant = False
            self.xs_s_data = xs_s_data
            self.n_groups = len(xs_s_data)
            self.E_bins = np.linspace(E_start, E_end, self.n_groups + 1)

    def EvalValue(self, x):
        """Evaluate the scattering cross-section at the given energy point.

        Args:
            x (list or array-like): Input point, where x[1] is the energy value.

        Returns:
            float: Corresponding scattering cross-section value.
        """
        if self.constant:
            return self.constant_value
        E = x[1]
        for i in range(self.n_groups - 1):
            if np.isclose(E, self.E_bins[i + 1]):
                return float(self.xs_s_data[i + 1])
            if E <= self.E_bins[i] and E > self.E_bins[i + 1]:
                return float(self.xs_s_data[i])
        return float(self.xs_s_data[-1])


class StoppingPowerCoefficientE(mfem.PyCoefficient):
    """Positive stopping power S(E) in MeV/cm for the conservative term -d_E(S*psi).

    ``S_data`` is a scalar or group values ordered from high to low energy;
    ``E_start``/``E_end`` are the high/low bounds in MeV. Lookup follows
    :class:`TotalXSCoefficientE` at physical energy ``x[1]``; positivity is
    not enforced by this helper. Grouped values are piecewise constant.
    Transport across discontinuous stopping-power interfaces remains unverified;
    constant and continuous variable-S cases have been verified.
    """

    def __init__(self, S_data, E_start, E_end):
        super(StoppingPowerCoefficientE, self).__init__()
        self.E_start = E_start
        self.E_end = E_end
        if isinstance(S_data, (int, float)):
            self.constant = True
            self.constant_value = float(S_data)
        else:
            self.constant = False
            self.S_data = S_data
            self.n_groups = len(S_data)
            self.E_bins = np.linspace(E_start, E_end, self.n_groups + 1)

    def EvalValue(self, x):
        """Evaluate the stopping power at the given energy point.

        Args:
            x (list or array-like): Input point, where x[1] is the energy value.

        Returns:
            float: Corresponding stopping power value.
        """
        if self.constant:
            return self.constant_value
        E = x[1]
        for i in range(self.n_groups - 1):
            if np.isclose(E, self.E_bins[i + 1]):
                return float(self.S_data[i + 1])
            if E <= self.E_bins[i] and E > self.E_bins[i + 1]:
                return float(self.S_data[i])
        return float(self.S_data[-1])


class QCoefficientE(mfem.PyCoefficient):
    """Isotropic source data Q(E), in angular-flux units per cm.

    ``Q_data`` is a scalar or uniform-group values ordered from high to low
    energy. ``E_start``/``E_end`` are the high/low bounds in MeV. Lookup and
    scalar handling follow :class:`TotalXSCoefficientE`. This is the source
    per angular direction, not an angle-integrated source.
    """

    def __init__(self, Q_data, E_start, E_end):
        super(QCoefficientE, self).__init__()
        self.E_start = E_start
        self.E_end = E_end
        if isinstance(Q_data, (int, float)):
            self.constant = True
            self.constant_value = float(Q_data)
        else:
            self.constant = False
            self.Q_data = Q_data
            self.n_groups = len(Q_data)
            self.E_bins = np.linspace(E_start, E_end, self.n_groups + 1)

    def EvalValue(self, x):
        """Evaluate the energy-dependent source term Q(E) at the given energy point.

        Args:
            x (list or array-like): Input point, where x[1] is the energy value.

        Returns:
            float: Corresponding source term value.
        """
        if self.constant:
            return self.constant_value
        E = x[1]
        for i in range(self.n_groups - 1):
            if np.isclose(E, self.E_bins[i + 1]):
                return float(self.Q_data[i + 1])
            if E <= self.E_bins[i] and E > self.E_bins[i + 1]:
                return float(self.Q_data[i])
        return float(self.Q_data[-1])


class QFuncCoefficient(mfem.PyCoefficient):
    """Manufactured source for constant-S CSDA in physical (x, E).

    ``pn=3,4,5`` selects flux ``a_val+b_val*x``, ``a_val+b_val*E`` or
    ``a_val+b_val*x*E``. The source is ``mu_val*d_x(psi) - S_val*d_E(psi)
    + xs_t_val*psi``, with removal ``xs_t_val`` in inverse cm and stopping
    power ``S_val`` in MeV/cm. Flux coefficients must have corresponding
    units. ``q_const`` is retained for compatibility and is not used.
    For variable S, supply a separate source including the derivative of S.
    """

    def __init__(self, pn, a_val=0, b_val=0, xs_t_val=0, mu_val=0, S_val=0, q_const=0):
        super(QFuncCoefficient, self).__init__()
        self.pn = pn
        self.a_val = a_val
        self.b_val = b_val
        self.xs_t_val = xs_t_val
        self.mu_val = float(mu_val)
        self.S_val = S_val
        self.q_const = q_const

    def EvalValue(self, x):
        if self.pn == 3:
            return self.coeff_pn3(x)
        elif self.pn == 4:
            return self.coeff_pn4(x)
        elif self.pn == 5:
            return self.coeff_pn5(x)

    def coeff_pn3(self, x):
        return float(
            self.mu_val * self.b_val + self.xs_t_val * (self.a_val + self.b_val * x[0])
        )

    def coeff_pn4(self, x):
        return float(
            self.xs_t_val * (self.a_val + self.b_val * x[1]) - self.S_val * self.b_val
        )

    def coeff_pn5(self, x):
        sol = (
            self.mu_val * self.b_val * x[1]
            + self.xs_t_val * (self.a_val + self.b_val * x[0] * x[1])
            - self.S_val * self.b_val * x[0]
        )
        return sol


class EDependentCoefficient(mfem.PyCoefficient):
    """Generic scalar or piecewise-constant function of physical energy x[1].

    ``data`` is a scalar or uniform-group values ordered high to low;
    ``E_start``/``E_end`` are high/low energy bounds in MeV. Lookup, edge
    handling and out-of-range fallback follow :class:`TotalXSCoefficientE`.
    The supplied quantity determines the coefficient units.
    """

    def __init__(self, data, E_start, E_end):
        super(EDependentCoefficient, self).__init__()
        self.E_start = E_start
        self.E_end = E_end
        if isinstance(data, (int, float)):
            self.constant = True
            self.constant_value = float(data)
        else:
            self.constant = False
            self.data = data
            self.n_groups = len(data)
            self.E_bins = np.linspace(E_start, E_end, self.n_groups + 1)

    def EvalValue(self, x):
        """Evaluate the coefficient value at the given energy point.

        Args:
            x (list or array-like): Input point, where x[1] is the energy value.

        Returns:
            float: Corresponding coefficient value.
        """
        if self.constant:
            return self.constant_value
        E = x[1]
        for i in range(self.n_groups - 1):
            if np.isclose(E, self.E_bins[i + 1]):
                return float(self.data[i + 1])
            if E <= self.E_bins[i] and E > self.E_bins[i + 1]:
                return float(self.data[i])
        return float(self.data[-1])


class XDependentCoefficient(mfem.PyCoefficient):
    """Scalar or piecewise-constant function of physical position x[0].

    ``data`` is a scalar constant or values on uniform intervals from
    ``x_start`` to ``x_end`` in increasing position order (cm in examples).
    Bounds are ignored for scalar data. Near an internal edge the right-hand
    interval is selected; unmatched positions return the last value, including
    out-of-range positions. Units are those of the supplied quantity.
    """

    def __init__(self, data, x_start, x_end):
        super(XDependentCoefficient, self).__init__()
        self.x_start = x_start
        self.x_end = x_end
        if isinstance(data, (int, float)):
            self.constant = True
            self.constant_value = float(data)
        else:
            self.constant = False
            self.data = data
            self.n_groups = len(data)
            self.x_bins = np.linspace(x_start, x_end, self.n_groups + 1)

    def EvalValue(self, x):
        """Evaluate the spatial-dependent value at the given spatial point.

        Args:
            x (list or array-like): Input point, where x[0] is the spatial value.

        Returns:
            float: Corresponding spatial-dependent value.
        """
        if self.constant:
            return self.constant_value

        x_val = x[0]
        for i in range(self.n_groups - 1):
            if np.isclose(x_val, self.x_bins[i + 1]):
                return float(self.data[i + 1])
            if self.x_bins[i] <= x_val < self.x_bins[i + 1]:
                return float(self.data[i])
        return float(self.data[-1])


class InflowCoefficient(mfem.PyCoefficient):
    """Incoming angular flux coefficient, including manufactured boundary traces.

    Args:
        inflow (float): Constant incoming angular flux.
        mu (float): Direction cosine.
        xl, xr (bool): If only xl is True, return inflow for mu > 0 and zero
            otherwise. If only xr is True, return inflow for mu < 0 and zero
            otherwise. Both or neither enables the constant for all directions.
        pn (int): For 3, 4 or 5, override the constant/flags with the exact
            manufactured trace a+b*x, a+b*E or a+b*x*E, respectively.
        a_val, b_val (float): Manufactured flux coefficients, with consistent units.

    DG flux integrators select the incoming faces, including the high-energy
    face for S > 0. These flags do not identify the location of a mesh face.
    """

    def __init__(self, inflow, mu, xl=None, xr=None, pn=None, a_val=0, b_val=0):
        super(InflowCoefficient, self).__init__()
        self.inflow = inflow
        self.mu = mu
        self.xl = xl
        self.xr = xr
        self.pn = pn
        self.a_val = a_val
        self.b_val = b_val

    def EvalValue(self, x):
        """Return the manufactured trace or direction-filtered inflow at (x, E)."""

        # MMS traces apply on each inflow face, including E_max.
        if self.pn == 3:
            return self.a_val + self.b_val * x[0]
        if self.pn == 4:
            return self.a_val + self.b_val * x[1]
        if self.pn == 5:
            return self.a_val + self.b_val * x[0] * x[1]

        if self.xl is True and self.xr is not True:
            if self.mu > 0:
                return self.inflow
            else:
                return 0.0

        elif self.xr is True and self.xl is not True:
            if self.mu < 0:
                return self.inflow
            else:
                return 0.0

        else:
            return self.inflow


class ConstantCoefficient(mfem.PyCoefficient):
    """MFEM scalar coefficient returning ``const`` everywhere; units are caller-defined."""

    def __init__(self, const):
        super(ConstantCoefficient, self).__init__()
        self.const = const

    def EvalValue(self, x):
        return self.const


def VectorConstCoefficient(vel1, vel2=None, vel3=None):
    """Return an MFEM VectorConstantCoefficient from one, two or three components.

    ``vel1`` is required; ``vel2`` and ``vel3`` add components in coordinate
    order. For coupled (x, E) transport, streaming uses ``(mu, 0)`` and energy
    convection uses ``(0, -S)``. Components carry the corresponding velocity units.
    """
    if vel2 is None:
        return mfem.VectorConstantCoefficient([vel1])
    elif vel3 is None:
        return mfem.VectorConstantCoefficient([vel1, vel2])
    else:
        return mfem.VectorConstantCoefficient([vel1, vel2, vel3])
