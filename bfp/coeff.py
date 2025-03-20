"""Module for coefficients"""

import mfem.ser as mfem
import numpy as np

__all__ = [
    'TotalXSCoefficientE',
    'ScatteringXSCoefficientE',
    'StoppingPowerCoefficientE',
    'QCoefficientE',
    'QFuncCoefficient',
    'EDependentCoefficient',
    'XDependentCoefficient',
    'InflowCoefficient',
    'ConstantCoefficient',
    'VectorConstCoefficient',
    ]

class TotalXSCoefficientE(mfem.PyCoefficient):
    """Coefficient class for the total macroscopic cross-section Σ_t(E). 
    (This coefficient will be updated for Σ_t(x,E))

    This class provides a piecewise-constant coefficient function based on energy groups.
    The energy range is defined from E_start (higher energy) to E_end (lower energy),
    with energy values sorted in descending order. The corresponding total cross-section
    values in xs_t_data are also arranged from higher to lower energy. This ordering is based
    on the logic of the Continuous Slowing Down Approximation (CSDA) approach, where the energy
    loss of charged particles is treated as a continuous process, naturally leading to a descending
    energy order.

    Attributes:
        E_start (float): Upper bound of the energy range.
        E_end (float): Lower bound of the energy range.
        xs_t_data (array-like or float): Total cross-section values per energy group,
            or a constant value.

    Args:
        xs_t_data (float or list of float): Total cross-section data (constant or group-based).
        E_start (float): Start energy for energy group mapping.
        E_end (float): End energy for energy group mapping.

    Examples:
        Constant cross-section:
            TotalXSCoefficientE(1.0, E_start=1.0, E_end=0.01)

        Energy-dependent cross-section:
            TotalXSCoefficientE(xs_t_data, E_start=1.0, E_end=0.01)
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
            if np.isclose(E, self.E_bins[i+1]):
                return float(self.xs_t_data[i+1])
            if E <= self.E_bins[i] and E > self.E_bins[i+1]:
                return float(self.xs_t_data[i])
        return float(self.xs_t_data[-1])


class ScatteringXSCoefficientE(mfem.PyCoefficient):
    """Coefficient class for the scattering cross-section Σ_s(E).
    (This coefficient will be updated for Σ_s(x,E))

    This class provides a piecewise-constant coefficient function based on energy groups.
    The energy range is defined from E_start (higher energy) to E_end (lower energy),
    with energy values sorted in descending order. The corresponding scattering cross-section
    values in xs_s_data are also arranged from higher to lower energy. This ordering is based
    on the logic of the Continuous Slowing Down Approximation (CSDA) approach, where the energy
    loss of charged particles is treated as a continuous process, naturally leading to a descending
    energy order.

    Attributes:
        E_start (float): Upper bound of the energy range.
        E_end (float): Lower bound of the energy range.
        xs_s_data (array-like or float): Scattering cross-section values per energy group,
            or a constant value.

    Args:
        xs_s_data (float or list of float): Scattering cross-section data (constant or group-based).
        E_start (float): Start energy for energy group mapping.
        E_end (float): End energy for energy group mapping.

    Examples:
        Constant cross-section:
            ScatteringXSCoefficientE(1.0, E_start=1.0, E_end=0.01)

        Energy-dependent cross-section:
            ScatteringXSCoefficientE(xs_s_data, E_start=1.0, E_end=0.01)
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
            if np.isclose(E, self.E_bins[i+1]):
                return float(self.xs_s_data[i+1])
            if E <= self.E_bins[i] and E > self.E_bins[i+1]:
                return float(self.xs_s_data[i])
        return float(self.xs_s_data[-1])


class StoppingPowerCoefficientE(mfem.PyCoefficient):
    """Coefficient class for the stopping power S(E).
    (This coefficient will be updated for S(x,E))

    This class provides a piecewise-constant coefficient function based on energy groups.
    The energy range is defined from E_start (higher energy) to E_end (lower energy),
    with energy values sorted in descending order. The corresponding stopping power
    values in S_data are also arranged from higher to lower energy. This ordering is based
    on the logic of the Continuous Slowing Down Approximation (CSDA) approach, where the energy
    loss of charged particles is treated as a continuous process, naturally leading to a descending
    energy order.

    Attributes:
        E_start (float): Upper bound of the energy range.
        E_end (float): Lower bound of the energy range.
        S_data (array-like or float): Stopping power values per energy group,
            or a constant value.

    Args:
        S_data (float or list of float): Stopping power data (constant or group-based).
        E_start (float): Start energy for energy group mapping.
        E_end (float): End energy for energy group mapping.

    Examples:
        Constant stopping power:
            StoppingPowerCoefficientE(1.0, E_start=1.0, E_end=0.01)

        Energy-dependent stopping power:
            StoppingPowerCoefficientE(S_data, E_start=1.0, E_end=0.01)
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
            if np.isclose(E, self.E_bins[i+1]):
                return float(self.S_data[i+1])
            if E <= self.E_bins[i] and E > self.E_bins[i+1]:
                return float(self.S_data[i])
        return float(self.S_data[-1])


class QCoefficientE(mfem.PyCoefficient):
    """Coefficient class for the energy-dependent source term Q(E). 
    (This coefficient will be updated for Q(x,E))

    This class provides a piecewise-constant coefficient function based on energy groups.
    The energy range is defined from E_start (higher energy) to E_end (lower energy),
    with energy values sorted in descending order. The corresponding source term values
    in Q_data are also arranged from higher to lower energy. This ordering is based on the
    logic of the Continuous Slowing Down Approximation (CSDA) approach, where the energy loss
    of charged particles is treated as a continuous process, naturally leading to a descending
    energy order.

    Attributes:
        E_start (float): Higher energy bound of the energy range.
        E_end (float): Lower energy bound of the energy range.
        Q_data (array-like or float): Energy-dependent source term values per energy group,
            or a constant value.

    Args:
        Q_data (float or list of float): Source term data (constant or group-based).
        E_start (float): Start energy for energy group mapping (higher energy).
        E_end (float): End energy for energy group mapping (lower energy).

    Examples:
        Constant source term:
            QCoefficient(1.0, E_start=1.0, E_end=0.01)

        Energy-dependent source term:
            QCoefficient(Q_data, E_start=1.0, E_end=0.01)
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
            if np.isclose(E, self.E_bins[i+1]):
                return float(self.Q_data[i+1])
            if E <= self.E_bins[i] and E > self.E_bins[i+1]:
                return float(self.Q_data[i])
        return float(self.Q_data[-1])


class QFuncCoefficient(mfem.PyCoefficient):
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
        return float(self.mu_val * self.b_val + self.xs_t_val * (self.a_val + self.b_val * x[0]))

    def coeff_pn4(self, x):
        return float(self.xs_t_val * (self.a_val + self.b_val * x[1]) + self.S_val * self.b_val)

    def coeff_pn5(self, x):
        sol =  self.mu_val * self.b_val * x[1] \
            + self.xs_t_val * (self.a_val + self.b_val * x[0] * x[1]) \
            + self.S_val * self.b_val * x[0]
        return sol


class EDependentCoefficient(mfem.PyCoefficient):
    """Coefficient class for energy-dependent values.

    This class provides a piecewise-constant coefficient function based on energy groups.
    The energy range is defined from E_start (higher energy) to E_end (lower energy),
    with energy values sorted in descending order. The corresponding energy-dependent
    values in data are also arranged from higher to lower energy. This ordering is based
    on the logic of the Continuous Slowing Down Approximation (CSDA) approach, where the energy
    loss of charged particles is treated as a continuous process, naturally leading to a descending
    energy order.

    Attributes:
        E_start (float): Upper bound of the energy range.
        E_end (float): Lower bound of the energy range.
        data (array-like or float): Energy-dependent values per energy group,
            or a constant value.

    Args:
        data (float or list of float): Energy-dependent data (constant or group-based).
        E_start (float): Start energy for energy group mapping.
        E_end (float): End energy for energy group mapping.

    Examples:
        Constant energy-dependent value:
            EDependentCoefficient(1.0, E_start=1.0, E_end=0.01)

        Energy-dependent coefficient:
            EDependentCoefficient(data, E_start=1.0, E_end=0.01)
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
            if np.isclose(E, self.E_bins[i+1]):
                return float(self.data[i+1])
            if E <= self.E_bins[i] and E > self.E_bins[i+1]:
                return float(self.data[i])
        return float(self.data[-1])


class XDependentCoefficient(mfem.PyCoefficient):
    """Coefficient for the spatial-dependent values.

    This class provides a piecewise-constant coefficient function based on spatial cell number.
    The spatial range is defined from x_start (lower bound) to x_end (upper bound),
    with spatial values sorted in ascending order. The corresponding spatial-dependent
    values in data are also arranged from lower to upper according to their grouping.
    
    Attributes:
        x_start (float): Lower bound of the spatial range.
        x_end (float): Upper bound of the spatial range.
        data (array-like or float): Spatial-dependent values per cell, or a constant value.

    Args:
        data (float or list of float): Spatial-dependent data (constant or cell-based).
        x_start (float): Start value for spatial mapping.
        x_end (float): End value for spatial mapping.

    Examples:
        Constant spatial-dependent value:
            XDependentCoefficient(1.0, x_start=0.0, x_end=10.0)

        Spatial-dependent coefficient:
            XDependentCoefficient(data, x_start=0.0, x_end=10.0)
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
            if np.isclose(x_val, self.x_bins[i+1]):
                return float(self.data[i+1])
            if self.x_bins[i] <= x_val < self.x_bins[i+1]:
                return float(self.data[i])
        return float(self.data[-1])


class InflowCoefficient(mfem.PyCoefficient):
    """Coefficient class for the inflow boundary condition.

    This class represents a coefficient for the inflow boundary condition.
    The coefficient returns a specified inflow value based on the angular direction `mu` and
    optional boundary flags. If `xl` is set to True, the coefficient returns the inflow value
    only when `mu < 0` (left boundary). If `xr` is set to True, the coefficient returns the
    inflow value only when `mu > 0` (right boundary). If both flags are True or if neither flag
    is provided (default), the inflow value is returned unconditionally.

    Args:
        inflow (float): The inflow value to be used.
        mu (float): The angular direction (e.g., cosine of the angle) used to determine if the 
            inflow condition is met.
        xl (bool, optional): If True, apply the inflow condition on the left boundary (when mu > 0).
            Defaults to None.
        xr (bool, optional): If True, apply the inflow condition on the right boundary (when mu < 0).
            Defaults to None.
    """

    def __init__(self, inflow, mu, xl=None, xr=None):
        super(InflowCoefficient, self).__init__()
        self.inflow = inflow
        self.mu = mu
        self.xl = xl
        self.xr = xr

    def EvalValue(self, x):
        """Evaluate the inflow coefficient.

        Depending on the specified boundary flags, this method returns the inflow value only if
        the corresponding condition based on `mu` is satisfied. If neither or both boundary flags
        are provided, the inflow value is returned unconditionally.

        Args:
            x (list or array-like): Input point.

        Returns:
            float: The inflow value if the condition is met, otherwise 0.0.
        """

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
    """Coefficient that always returns a constant value.

    This class represents a coefficient that returns the same constant value for any input.
    It is useful in finite element computations where a uniform property or parameter is needed.

    Args:
        const (float): The constant value to be returned.
    """
    def __init__(self, const):
        super(ConstantCoefficient, self).__init__()
        self.const = const

    def EvalValue(self, x):
        return self.const


def VectorConstCoefficient(vel1, vel2=None, vel3=None):
    """
    Create and return a constant velocity vector coefficient using mfem.VectorConstantCoefficient.

    This function instantiates a VectorConstCoefficient that returns a constant velocity vector.
    The dimensionality of the vector is determined by the number of velocity values provided:
      - 1D if only vel1 is provided.
      - 2D if both vel1 and vel2 are provided.
      - 3D if vel1, vel2, and vel3 are provided.

    Args:
        vel1 (float): Constant velocity value for the x-direction.
        vel2 (float, optional): Constant velocity value for the y-direction. Defaults to None.
        vel3 (float, optional): Constant velocity value for the z-direction. Defaults to None.

    Returns:
        VectorConstCoefficient: An instance of VectorConstCoefficient representing the constant velocity vector.
    """
    if vel2 is None:
        return mfem.VectorConstantCoefficient([vel1])
    elif vel3 is None:
        return mfem.VectorConstantCoefficient([vel1, vel2])
    else:
        return mfem.VectorConstantCoefficient([vel1, vel2, vel3])
    

