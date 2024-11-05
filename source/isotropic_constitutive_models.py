# Routine to store the derivatives of the Helmholtz energy function of 
# different models

import ufl_legacy as ufl

########################################################################
#                             Saint-Venant                             #
########################################################################

# Derivatives of the Helmholtz energy function w.r.t. the invariants of
# the right Cauchy-Green tensor

def dPsi_dI1SaintVenant(I1, I2, parameters):

    mu, lmbda = parameters
    
    return (0.5*mu*(I1-1))+(0.25*lmbda*(I1-3))

def dPsi_dI2SaintVenant(I1, I2, parameters):

    mu, lmbda = parameters
    
    return -0.5*mu

def ddPsi_ddI1SaintVenant(I1, I2, parameters):

    mu, lmbda = parameters
    
    return (0.5*mu)+(0.25*lmbda)

def ddPsi_ddI2SaintVenant(I1, I2, parameters):

    mu, lmbda = parameters
    
    return 0.0

def ddPsi_dI1dI2SaintVenant(I1, I2, parameters):

    mu, lmbda = parameters
    
    return 0.0

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the
# jacobian

def k_volSaintVenant(J, parameters):
    
    return 0.0

def dkvol_dJSaintVenant(J, parameters):
    
    return 0.0

########################################################################
#                             Neo-Hookean                              #
########################################################################

# Derivatives of the Helmholtz energy function w.r.t. the invariants of
# the right Cauchy-Green tensor

def dPsi_dI1CompressibleNeoHookean(I1, I2, parameters):

    mu, lmbda = parameters
    
    return (0.5*mu)

def dPsi_dI2CompressibleNeoHookean(I1, I2, parameters):
    
    return 0.0

def ddPsi_ddI1CompressibleNeoHookean(I1, I2, parameters):
    
    return 0.0

def ddPsi_ddI2CompressibleNeoHookean(I1, I2, parameters):
    
    return 0.0

def ddPsi_dI1dI2CompressibleNeoHookean(I1, I2, parameters):
    
    return 0.0

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the
# jacobian

def k_volCompressibleNeoHookean(J, parameters):

    mu, lmbda = parameters
    
    return 0.5*(-mu+(lmbda*ufl.ln(J)))

def dkvol_dJCompressibleNeoHookean(J, parameters):

    mu, lmbda = parameters
    
    return 0.5*(lmbda/J)

########################################################################
#                            Mooney-Rivlin                             #
########################################################################

# Derivatives of the Helmholtz energy function w.r.t. the isochoric in-
# variants of the right Cauchy-Green tensor

def dPsi_dI1MooneyRivlin(I1_iso, I2_iso, parameters):

    mu_1, mu_2, mu_3 = parameters
    
    return (0.5*mu_1)

def dPsi_dI2MooneyRivlin(I1_iso, I2_iso, parameters):

    mu_1, mu_2, mu_3 = parameters
    
    return (0.5*mu_2)

def ddPsi_ddI1MooneyRivlin(I1_iso, I2_iso, parameters):
    
    return 0.0

def ddPsi_ddI2MooneyRivlin(I1_iso, I2_iso, parameters):
    
    return 0.0

def ddPsi_dI1dI2MooneyRivlin(I1_iso, I2_iso, parameters):
    
    return 0.0

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the
# jacobian

def k_volMooneyRivlin(J, parameters):

    mu_1, mu_2, mu_3 = parameters
    
    return 0.5*(mu_3*J*(J-1))

def dkvol_dJMooneyRivlin(J, parameters):

    mu_1, mu_2, mu_3 = parameters
    
    return 0.5*(mu_3*(J-1))+(0.5*mu_3*J)

########################################################################
#                        Uncoupled Neo-Hookean                         #
########################################################################

# Derivatives of the Helmholtz energy function w.r.t. the invariants of
# the right Cauchy-Green tensor

def dPsi_dI1UncoupledNeoHookean(I1_iso, I2_iso, parameters):

    mu_1, bulk = parameters
    
    return (0.5*mu_1)

def dPsi_dI2UncoupledNeoHookean(I1_iso, I2_iso, parameters):
    
    return 0.0

def ddPsi_ddI1UncoupledNeoHookean(I1_iso, I2_iso, parameters):
    
    return 0.0

def ddPsi_ddI2UncoupledNeoHookean(I1_iso, I2_iso, parameters):
    
    return 0.0

def ddPsi_dI1dI2UncoupledNeoHookean(I1_iso, I2_iso, parameters):
    
    return 0.0

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the
# jacobian

def k_volUncoupledNeoHookean(J, parameters):

    mu_1, bulk = parameters
    
    return 0.5*(bulk*J*(J-1))

def dkvol_dJUncoupledNeoHookean(J, parameters):

    mu_1, bulk = parameters
    
    return (0.5*(bulk*(J-1)))+(0.5*J*bulk)