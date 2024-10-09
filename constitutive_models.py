# Routine to store the derivatives of the Helmholtz energy function of 
# different models

########################################################################
#                             Saint-Venant                             #
########################################################################

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