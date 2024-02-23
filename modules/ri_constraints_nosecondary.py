"""
Refractive Index Optimization Constraints

This module defines a set of constraints used in the optimization of refractive index (RI) values 
for atmospheric aerosols. These constraints ensure that the optimized RI values adhere to 
expected physical and empirical boundaries. They are crucial for maintaining the physical 
plausibility and accuracy of aerosol models in atmospheric science.

Each function in this module represents a specific constraint and returns a value indicating 
whether the constraint is satisfied. A positive return value signifies that the constraint is met, 
while a negative value indicates a violation of the constraint.

Functions:
- constraint1: Ensures the GFAS is always greater than all others.

"""
def constraint1(ri_values):
    """
    Ensure that the value of GFAS is always greater than all others.

    Parameters:
    ri_values (list): List of values in the order [GFAS, RESI, SHIP, TRAF, OTHR].
    
    Returns:
    float: The difference between GFAS and the maximum of the others.
    """
    return ri_values[0] - max(ri_values[1:])