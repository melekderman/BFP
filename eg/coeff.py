import mfem.ser as mfem
import numpy as np
import h5py


class SimpleCoefficient(mfem.PyCoefficient):
    """
    A simple coefficient class that maps a normalized energy value to a corresponding 
    value in the provided data array.

    The energy groups are determined by splitting the energy range defined by E_start 
    (first/highest energy) and E_end (last/lowest energy) into len(data) groups.
    
    Attributes:
        data (list or numpy.array): A list of coefficient values corresponding to the energy groups.
        E_start (float): The starting (highest) energy value.
        E_end (float): The ending (lowest) energy value.
    """
    
    def __init__(self, data, E_start, E_end):
        super(SimpleCoefficient, self).__init__()
        self.data = data
        self.E_start = E_start
        self.E_end = E_end
        self.nE = len(data)  

    def EvalValue(self, x):
        """
        Evaluates the coefficient value at a given point.

        The energy coordinate is taken from the second component of the input vector x 
        (assumed to be normalized between 0 and 1). This normalized value is used to directly
        select one of the energy groups from the data array.

        Args:
            x (mfem.Vector): The coordinate vector, where x[1] is the normalized energy value.

        Returns:
            float: The coefficient value corresponding to the selected energy group.
        """
        # Compute the energy group index from the normalized energy value.
        group = int(x[1] * (self.nE - 1))
        group = min(self.nE - 1, max(0, group))
        return float(self.data[group])

    


class QFunction(mfem.PyCoefficient):
    """
    Define the source term Q(x,E).
    
    Args:
        q_value (float): Source term value.
    """
    def __init__(self, q_value):
        self.q_value = q_value

    def EvalValue(self, x):
        return self.q_value


class InflowCoeff(mfem.PyCoefficient):
    """
    Define the inflow boundary condition g(x,E) at x = 0.
    For this example, the flux is set to the user-specified value.
    
    Args:
        coeff_value (float): The flux value provided by the user.
    """
    def __init__(self, coeff_value):
        self.coeff_value = coeff_value

    def EvalValue(self, x):
        return self.coeff_value
    
    import mfem.ser as mfem
import math

class AngularFluxCoefficient(mfem.PyCoefficient):
    """
    Angular flux coefficient for discrete ordinates using Gauss-Legendre quadrature.
    
    This coefficient defines the angular dependence of the flux. It maps the angular
    variable (assumed to be x[0] in the integration point, representing the cosine of the angle)
    to a prescribed angular flux value.
    
    For example, a simple angular flux distribution is defined as:
        psi(mu) = 1 + mu,
    where mu is in the interval [-1, 1]. This provides a non-uniform angular distribution.
    
    This coefficient is intended to be used in angular discretization where a Gauss-Legendre
    quadrature rule of order p=6 is applied.
    
    Attributes:
        None.
    """
    
    def EvalValue(self, x):
        """
        Evaluates the angular flux at a given angular coordinate.
        
        Args:
            x (mfem.Vector): The coordinate vector, where x[0] is assumed to be the angular variable (mu).
            
        Returns:
            float: The angular flux value at the given coordinate.
        """
        # Assume x[0] is the angular variable (mu) in the range [-1, 1]
        mu = x[0]
        # Define a simple angular flux distribution, e.g., psi(mu) = 1 + mu
        return 1.0 + mu

'''
# Usage Example:
# Create an instance of AngularFluxCoefficient.
angular_flux_coeff = AngularFluxCoefficient()

# When integrating over the angular variable, use Gauss-Legendre quadrature with p=6.
# For example, to obtain the integration rule for a segment (1D interval):
IR = mfem.IntRules.Get(mfem.Geometry.SEGMENT, 6)
# IR now contains the 6 Gauss-Legendre quadrature points and weights which can be used
# in the assembly of integrals involving the angular variable.

'''