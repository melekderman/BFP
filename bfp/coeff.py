"""Module for coefficients"""

import mfem.ser as mfem

__all__ = [
    'TotalXSCoefficient',
    'ScatteringXSCoefficient',
    'StoppingPowerCoefficient',
    'StoppingPowerDerivativeCoefficient',
    'InflowCoefficientSN',
    'QCoefficient',
    'EnergyDependentCoefficient',
    'XDependentCoefficient',
    'VelocityCoefficientOld',
    'VelocityCoefficient',
    'VelocityCoefficient2',
    'ConstantCoefficient'
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

    def __init__(self, xs_t_data, E_start=None, E_end=None):
        super(TotalXSCoefficient, self).__init__()
        if isinstance(xs_t_data, (int, float)):
            self.constant = True
            self.constant_value = float(xs_t_data)
        else:
            self.constant = False
            self.xs_t_data = xs_t_data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self,ip):
        """Evaluates the total cross-section at a given energy coordinate.

        This method converts the normalized energy coordinate found in x[1] to the actual 
        energy value and determines the corresponding energy group to return the associated 
        total cross-section value. 

        Args:
        ip (list or array-like): A coordinate array where ip[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The total cross-section value corresponding to the computed energy.
        """
        if self.constant:
            return self.constant_value

        y = ip[1]
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

    def __init__(self, xs_s_data, E_start=None, E_end=None):
        super(ScatteringXSCoefficient, self).__init__()
        if isinstance(xs_s_data, (int, float)):
            self.constant = True
            self.constant_value = float(xs_s_data)
        else:
            self.constant = False
            self.xs_s_data = xs_s_data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self, ip):
        """Evaluates the scattering cross-section at a given normalized energy coordinate.

        This method converts the normalized energy coordinate found in x[1] to the actual 
        energy value and determines the corresponding energy group to return the associated 
        scattering cross-section value. 

        Args:
            ip (list or array-like): A coordinate array where ip[1] is the normalized energy 
                in the range [0, 1].

        Returns:
            float: The scattering cross-section value for the computed energy.
        """
        if self.constant:
            return self.constant_value

        y = ip[1] 
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.xs_s_data)
        group = min(n_groups - 1,
                    int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.xs_s_data[group])

class StoppingPowerCoefficient(mfem.PyCoefficient):
    """Coefficient for the stopping power S(E).

    This coefficient maps the normalized energy coordinate (ip[1]) in the interval [0, 1]
    to the corresponding stopping power S(E). Specifically, a value of y = 0 corresponds to 
    E = E_start, and a value of y = 1 corresponds to E = E_end. If a constant is provided 
    (e.g., 1), the coefficient returns this constant at all points.

    Examples:
        StoppingPowerCoefficient(1)
            # Always returns 1.

        StoppingPowerCoefficient(S_data, E_start=0.0, E_end=1.0)
            # Returns the interpolated value from S_data.
    """

    def __init__(self, S_data, E_start=None, E_end=None):
        super(StoppingPowerCoefficient, self).__init__()
        if isinstance(S_data, (int, float)):
            self.constant = True
            self.constant_value = float(S_data)
        else:
            self.constant = False
            self.S_data = S_data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self, ip):
        """Evaluates the stopping power at a given normalized energy coordinate.

        This method converts the normalized energy coordinate found in ip[1] to the actual 
        energy value and determines the corresponding energy group to return the associated 
        stopping power value. 

        Args:
            ip (list or array-like): A coordinate array where ip[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The stopping power value corresponding to the computed energy.
        """
        if self.constant:
            return self.constant_value
        y = ip[1] 
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

    def __init__(self, dS_data, E_start=None, E_end=None):
        super(StoppingPowerDerivativeCoefficient, self).__init__()
        if isinstance(dS_data, (int, float)):
            self.constant = True
            self.constant_value = float(dS_data)
        else:
            self.constant = False
            self.dS_data = dS_data
            self.E_start = E_start
            self.E_end = E_end

    def EvalValue(self, ip):
        """Evaluates the derivative of the stopping power at a given normalized energy coordinate.

        This method converts the normalized energy coordinate found in ip[1] to the actual 
        energy value and determines the corresponding energy group to return the associated 
        derivative of the stopping power value. 

        Args:
            ip (list or array-like): A coordinate array where ip[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The derivative of the stopping power corresponding to the computed energy.
        """
        if self.constant:
            return self.constant_value
        y = ip[1] 
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.dS_data)
        group = min(n_groups - 1,
                    int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.dS_data[group])

class InflowCoefficientSN(mfem.PyCoefficient):
    """Coefficient for an inflow boundary condition in (x,E) space with SN angular dependence.

    This class returns an inflow flux value that depends on the discrete ordinate
    (angular) parameter \(\mu\). On the left boundary:
      - If \(\mu > 0\), the inflow flux is applied (e.g., set to a prescribed value).
      - If \(\mu < 0\), the incoming flux is zero.
    
    This is useful in SN (discrete ordinates) methods where the transport equation is
    solved separately for each angular direction.
    
    Attributes:
        in_flux (float): The constant inflow flux value to be used when \(\mu > 0\).
        mu (float): The discrete ordinate (angular direction) for which this coefficient is used.
    """

    def __init__(self, in_flux, mu):
        """Initializes the SN inflow coefficient.
        
        Args:
            in_flux (float): The constant inflow flux value for \(\mu > 0\).
            mu (float): The discrete angular direction value.
        """
        super(InflowCoefficientSN, self).__init__()
        self.in_flux = in_flux
        self.mu = mu

    def EvalValue(self, ip):
        """Evaluates the inflow coefficient at a given integration point.
        
        In this implementation, the integration point is interpreted in (x,E) space.
        The angular parameter \(\mu\) is stored in the class. The function returns
        the prescribed inflow flux if \(\mu > 0\) and zero otherwise.
        
        Args:
            ip (mfem.IntegrationPoint or list[float]): The integration point coordinates,
                where ip[0] is the spatial coordinate \(x\) and ip[1] is the energy \(E\).
                The angular information is not included in ip.
        
        Returns:
            float: The inflow flux value if \(\mu > 0\); otherwise, 0.0.
        """
        # Since the integration point ip does not include the angular variable,
        # we use the stored self.mu to determine the flux.
        if self.mu > 0:
            return self.in_flux
        else:
            return 0.0

class QCoefficient(mfem.PyCoefficient):
    """Coefficient for the energy dependent source term Q(E). ((will be updated for Q(x,E)))

    This coefficient maps the normalized energy coordinate (ip[1]) in the interval [0, 1]
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

    def EvalValue(self, ip):
        """Evaluates the energy dependent source term Q(E) at a given normalized energy coordinate.

        For non-constant Q_data, this method converts the normalized energy coordinate ip[1] to the 
        corresponding energy value and returns the Q value for the associated energy group. 
        If Q_data is constant, it always 
        returns this constant value.

        Args:
            ip (list or array-like): A coordinate array where ip[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The energy dependent source value corresponding to the computed energy.
        """

        if self.constant:
            return self.constant_value

        y = ip[1] 
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.Q_data)
        group = min(n_groups - 1, int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.Q_data[group])

class EnergyDependentCoefficient(mfem.PyCoefficient):
    """Energy-dependent coefficient using either a constant or a one-dimensional array.

    This coefficient maps a normalized energy coordinate (ip[1] in [0, 1])
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

    def EvalValue(self,ip):
        """Evaluates the coefficient at a given normalized energy coordinate.

        This method converts the normalized energy coordinate (ip[1]) to a physical energy
        value, determines the corresponding index (group) in the data array via linear 
        interpolation, and returns the associated coefficient value.

        Args:
            ip (list or array-like): A coordinate array where ip[1] is the normalized energy 
                (in the range [0, 1]).

        Returns:
            float: The interpolated coefficient value corresponding to the computed energy.
        """
        if self.constant:
            return self.constant_value

        y = ip[1] 
        E = self.E_start + y * (self.E_end - self.E_start)
        n_groups = len(self.data)
        group = min(n_groups - 1, int((E - self.E_start) / (self.E_end - self.E_start) * n_groups))
        return float(self.data[group])
    
class XDependentCoefficient(mfem.PyCoefficient):
    """X-dependent coefficient using either a constant or a one-dimensional array.

    This coefficient maps a normalized x coordinate (ip[0] in [0, 1])
    to a value by linearly interpolating data over an interval defined by x_start and x_end.
    If a constant is provided (e.g., 1), it is converted to an array (with two points)
    and then processed identically to an array input.

    Examples:
        XDependentCoefficient(1)
            # Always returns 1.

        XDependentCoefficient(data_array, x_start=0.0, x_end=1.0)
            # Returns the interpolated value from data_array.
    """

    def __init__(self, data, x_start=None, x_end=None):
        super(XDependentCoefficient, self).__init__()
        if isinstance(data, (int, float)):
            self.constant = True
            self.constant_value = float(data)
        else:
            self.constant = False
            self.data = data
            self.x_start = x_start
            self.x_end = x_end

    def EvalValue(self, ip):
        """Evaluates the coefficient at a given normalized x coordinate.

        This method converts the normalized x coordinate (ip[0]) to a physical x
        value, determines the corresponding index (group) in the data array via linear 
        interpolation, and returns the associated coefficient value.

        Args:
            ip (list or array-like): A coordinate array where ip[0] is the normalized x 
                (in the range [0, 1]).

        Returns:
            float: The interpolated coefficient value corresponding to the computed x.
        """
        if self.constant:
            return self.constant_value

        y = ip[0] 
        x_val = self.x_start + y * (self.x_end - self.x_start)
        n_groups = len(self.data)
        group = min(n_groups - 1, int((x_val - self.x_start) / (self.x_end - self.x_start) * n_groups))
        return float(self.data[group])

class VelocityCoefficientOld(mfem.VectorPyCoefficientBase):
    """Velocity vector coefficient for the transport equation.

    Defines a velocity vector:
        v(x) = [μ, S(E)]

    Attributes:
        mu (float): Scalar value, e.g., discrete ordinate.
        S_coef (Coefficient): Coefficient object used to evaluate S(E) at given points.
    """

    def __init__(self, mu, S_coeff):
        mfem.VectorPyCoefficientBase.__init__(self, 2, 0)
        self.mu = mu
        self.S_coeff = S_coeff

    def _EvalPy(self, ip):
        """Evaluates the velocity vector at a given integration point.

        Args:
            V (array_like): Output vector of size 2 to store the velocity.
            ip (IntegrationPoint): Integration point where the evaluation occurs.
        """
        return [self.mu, self.S_coeff(ip)]

class VelocityCoefficient(mfem.VectorPyCoefficient):
    """
    A simple vector coefficient:
       v(x) = [mu, S(E)]
    that returns a Python list in EvalValue.
    """
    def __init__(self, mu, S_arr, E_start, E_end):
        super(VelocityCoefficient, self).__init__(2)  # 2D vector
        self.mu = mu
        self.S_arr = S_arr
        self.E_start = E_start
        self.E_end   = E_end

    def EvalValue(self, x):
        """
        x is the coordinate array, e.g. [x_coord, E_coord].
        Returns [mu, S(E)].
        """
        comp0 = self.mu
        
        E = self.E_start + x[1]*(self.E_end - self.E_start)
        
        n_groups = len(self.S_arr)
        idx = min(n_groups - 1,
                  int((E - self.E_start)/(self.E_end - self.E_start) * n_groups))
        
        comp1 = self.S_arr[idx]
        
        return [comp0, comp1]
    
class VelocityCoefficient2(mfem.VectorPyCoefficientBase):
    """Velocity vector coefficient for the transport equation.

    Defines a velocity vector:
        v(x) = [μ, S(E)]

    Attributes:
        mu (float): Scalar value, e.g., discrete ordinate.
        S_coef (Coefficient): Coefficient object used to evaluate S(E) at given points.
    """

    def __init__(self, mu, S_coeff):
        mfem.VectorPyCoefficientBase.__init__(self, 2, 0)
        self.mu = mu
        self.S_coeff = S_coeff

    def _EvalPy(self, V, ip):
        """Evaluates the velocity vector at a given integration point.

        Args:
            V (array_like): Output vector of size 2 to store the velocity.
            ip (IntegrationPoint): Integration point where the evaluation occurs.
        """
        return [self.mu, self.S_coeff]

class ConstantCoefficient(mfem.PyCoefficient):
    def __init__(self, const):
        super(ConstantCoefficient, self).__init__()
        self.const = const
    def EvalValue(self, x):
        return self.const
