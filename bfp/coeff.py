"""Module for coefficients"""

import mfem.ser as mfem

__all__ = [
    'TotalXSCoefficient',
    'ScatteringXSCoefficient',
    'StoppingPowerCoefficient',
    'StoppingPowerDerivativeCoefficient',
    'InflowCoefficient',
    'QCoefficient',
    'EnergyDependentCoefficient',
    'VelocityCoefficient'
    ]

class TotalXSCoefficient(mfem.PyCoefficient):
    """Coefficient for the total cross-section Σ_t(E).

    This coefficient maps the normalized energy coordinate (x[1]) in [0, 1] to the
    corresponding total cross-section value. Specifically, a value of y = 0 corresponds to
    E = E_start and y = 1 corresponds to E = E_end. If a constant is provided 
    (e.g., 1), the coefficient returns this constant at all points.

    Examples:
        TotalXSCoefficient(1)
            # Always returns 1.

        TotalXSCoefficient(xs_t_data, E_start=0.0, E_end=1.0)
            # Returns the interpolated value from xs_t_data.
    """

    def __init__(self, xs_t_data, E_start, E_end):
        super(TotalXSCoefficient, self).__init__()
        if isinstance(xs_t_data, (int, float)):
            self.constant = True
            self.constant_value = float(xs_t_data)
        else:
            self.constant = False
            self.xs_t_data = xs_t_data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self, x):
        """Evaluates the total cross-section at a given energy coordinate.

        This method converts the normalized energy coordinate found in x[1] to the actual 
        energy value and determines the corresponding energy group to return the associated 
        total cross-section value. 

        Args:
            x (list or array-like): A coordinate array where x[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The total cross-section value corresponding to the computed energy.
        """
        if self.constant:
            return self.constant_value
        y = x[1]
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.xs_t_data)
        group = min(n_groups - 1,
                    int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.xs_t_data[group])


class ScatteringXSCoefficient(mfem.PyCoefficient):
    """Coefficient for the scattering cross-section Σ_s(E).

    This coefficient maps the normalized energy coordinate (x[1]) in [0, 1] to the
    corresponding cross-section value. Specifically, a value of y = 0 corresponds to
    E = E_start and y = 1 corresponds to E = E_end. If a constant is provided 
    (e.g., 1), the coefficient returns this constant at all points.

    Examples:
        ScatteringXSCoefficient(1)
            # Always returns 1.

        ScatteringXSCoefficient(xs_s_data, E_start=0.0, E_end=1.0)
            # Returns the interpolated value from xs_s_data.
    """

    def __init__(self, xs_s_data, E_start, E_end):
        super(ScatteringXSCoefficient, self).__init__()
        if isinstance(xs_s_data, (int, float)):
            self.constant = True
            self.constant_value = float(xs_s_data)
        else:
            self.constant = False
            self.xs_s_data = xs_s_data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self, x):
        """Evaluates the scattering cross-section at a given normalized energy coordinate.

        This method converts the normalized energy coordinate found in x[1] to the actual 
        energy value and determines the corresponding energy group to return the associated 
        scattering cross-section value. 

        Args:
            x (list or array-like): A coordinate array where x[1] is the normalized energy 
                in the range [0, 1].

        Returns:
            float: The scattering cross-section value for the computed energy.
        """
        if self.constant:
            return self.constant_value
        y = x[1]
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.xs_s_data)
        group = min(n_groups - 1,
                    int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.xs_s_data[group])

class StoppingPowerCoefficient(mfem.PyCoefficient):
    """Coefficient for the stopping power S(E).

    This coefficient maps the normalized energy coordinate (x[1]) in the interval [0, 1]
    to the corresponding stopping power S(E). Specifically, a value of y = 0 corresponds to 
    E = E_start, and a value of y = 1 corresponds to E = E_end. If a constant is provided 
    (e.g., 1), the coefficient returns this constant at all points.

    Examples:
        StoppingPowerCoefficient(1)
            # Always returns 1.

        StoppingPowerCoefficient(S_data, E_start=0.0, E_end=1.0)
            # Returns the interpolated value from S_data.
    """

    def __init__(self, S_data, E_start, E_end):
        super(StoppingPowerCoefficient, self).__init__()
        if isinstance(S_data, (int, float)):
            self.constant = True
            self.constant_value = float(S_data)
        else:
            self.constant = False
            self.S_data = S_data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self, x):
        """Evaluates the stopping power at a given normalized energy coordinate.

        This method converts the normalized energy coordinate found in x[1] to the actual 
        energy value and determines the corresponding energy group to return the associated 
        stopping power value. 

        Args:
            x (list or array-like): A coordinate array where x[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The stopping power value corresponding to the computed energy.
        """
        if self.constant:
            return self.constant_value
        y = x[1]
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.S_data)
        group = min(n_groups - 1,
                    int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.S_data[group])

class StoppingPowerDerivativeCoefficient(mfem.PyCoefficient):
    """Coefficient for the derivative of the stopping power S(E), i.e., S'(E).

    This coefficient maps the normalized energy coordinate (x[1]) in the interval [0, 1]
    to the corresponding derivative of the stopping power S(E). A value of y = 0 corresponds 
    to E = E_start and y = 1 corresponds to E = E_end. If a constant is provided 
    (e.g., 1), the coefficient returns this constant at all points.

    Examples:
        StoppingPowerDerivativeCoefficient(1)
            # Always returns 1.

        StoppingPowerDerivativeCoefficient(dS_data, E_start=0.0, E_end=1.0)
            # Returns the interpolated value from data_array.
    """

    def __init__(self, dS_data, E_start, E_end):
        super(StoppingPowerDerivativeCoefficient, self).__init__()
        if isinstance(dS_data, (int, float)):
            self.constant = True
            self.constant_value = float(dS_data)
        else:
            self.constant = False
            self.dS_data = dS_data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self, x):
        """Evaluates the derivative of the stopping power at a given normalized energy coordinate.

        This method converts the normalized energy coordinate found in x[1] to the actual 
        energy value and determines the corresponding energy group to return the associated 
        derivative of the stopping power value. 

        Args:
            x (list or array-like): A coordinate array where x[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The derivative of the stopping power corresponding to the computed energy.
        """
        if self.constant:
            return self.constant_value
        y = x[1]
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.dS_data)
        group = min(n_groups - 1,
                    int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.dS_data[group])

class InflowCoefficient(mfem.PyCoefficient):
    """Coefficient for the inflow boundary condition.

    Returns a prescribed inflow flux value.

    Attributes:
        inflow_value (float): The constant inflow flux value to be used.
    """

    def __init__(self, inflow_value):
        super(InflowCoefficient, self).__init__()
        self.inflow_value = inflow_value

    def EvalValue(self, x):
        """Evaluates the inflow coefficient at a given coordinate.

        This method returns the given constant inflow flux value, irrespective of the 
        input spatial coordinate.

        Args:
            x (iterable): The spatial coordinate (ignored in this coefficient).

        Returns:
            float: The prescribed inflow flux value.
        """
        return self.inflow_value

class QCoefficient(mfem.PyCoefficient):
    """Coefficient for the energy dependent source term Q(E). ((will be updated for Q(x,E)))

    This coefficient maps the normalized energy coordinate (x[1]) in the interval [0, 1]
    to the corresponding the energy dependent source Q(E). Specifically, a value of y = 0 
    corresponds to E = E_start, and a value of y = 1 corresponds to E = E_end. If a constant 
    is provided 
    (e.g., 1), the coefficient returns this constant at all points.

    Examples:
        StoppingPowerCoefficient(1)
            # Always returns 1.

        StoppingPowerCoefficient(S_data, E_start=0.0, E_end=1.0)
            # Returns the interpolated value from S_data.
    """
    def __init__(self, Q_data, E_start=None, E_end=None):
        super(QCoefficient, self).__init__()
        if isinstance(Q_data, (int, float)):
            self.constant = True
            self.constant_value = float(Q_data)
        else:
            self.constant = False
            self.Q_data = Q_data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self, x):
        """Evaluates the energy dependent source term Q(E) at a given normalized energy coordinate.

        For non-constant Q_data, this method converts the normalized energy coordinate x[1] to the 
        corresponding energy value and returns the Q value for the associated energy group. 
        If Q_data is constant, it always 
        returns this constant value.

        Args:
            x (list or array-like): A coordinate array where x[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The energy dependent source value corresponding to the computed energy.
        """

        if self.constant:
            return self.constant_value

        y = x[1]
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.Q_data)
        group = min(n_groups - 1, int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.Q_data[group])


class EnergyDependentCoefficient(mfem.PyCoefficient):
    """Energy-dependent coefficient using either a constant or a one-dimensional array.

    This coefficient maps a normalized energy coordinate (x[1] in [0, 1])
    to a value by linearly interpolating data over an interval defined by E_start and E_end.
    If a constant is provided (e.g., 1), it is converted to an array (with two points)
    and then processed identically to an array input.

    Examples:
        EnergyDependentCoefficient(1)
            # Always returns 1.

        EnergyDependentCoefficient(data_array, E_start=0.0, E_end=1.0)
            # Returns the interpolated value from data_array.
    """

    def __init__(self, data, E_start=None, E_end=None):
        super(EnergyDependentCoefficient, self).__init__()
        if isinstance(data, (int, float)):
            self.constant = True
            self.constant_value = float(data)
        else:
            self.constant = False
            self.data = data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self, x):
        """Evaluates the coefficient at a given normalized energy coordinate.

        This method converts the normalized energy coordinate (x[1]) to a physical energy
        value, determines the corresponding index (group) in the data array via linear 
        interpolation, and returns the associated coefficient value.

        Args:
            x (list or array-like): A coordinate array where x[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The interpolated coefficient value corresponding to the computed energy.
        """
        if self.constant:
            return self.constant_value

        y = x[1]
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.data)
        group = min(n_groups - 1, int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.data[group])


class VelocityCoefficient(mfem.VectorPyCoefficientBase):
    """Velocity vector coefficient for the transport equation.

    Defines a velocity vector:
        v(x) = [μ, S(E)]

    Attributes:
        mu (float): Scalar value, e.g., discrete ordinate.
        S_coef (Coefficient): Coefficient object used to evaluate S(E) at given points.
    """

    def __init__(self, mu, S_coef):
        mfem.VectorPyCoefficientBase.__init__(self, 2, 0)
        self.mu = mu
        self.S_coef = S_coef

    def _EvalPy(self, V, ip):
        """Evaluates the velocity vector at a given integration point.

        Args:
            V (array_like): Output vector of size 2 to store the velocity.
            ip (IntegrationPoint): Integration point where the evaluation occurs.
        """
        V[0] = self.mu
        V[1] = self.S_coef.EvalValue(ip)