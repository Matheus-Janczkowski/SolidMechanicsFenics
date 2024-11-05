# Routine to store the derivatives of the Helmholtz energy function of 
# different models

import ufl_legacy as ufl

########################################################################
#                        Holzapfel-Gassner-Ogden                       #
########################################################################

# Defines the first and second derivatives of the anisotropic part of 
# the energy function

def dPsi_dI1PseudoHGO(I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, 
parameters_aniso):
    
    alpha1, alpha2 = parameters_aniso
    
    return (alpha1*(I1_isoPseudo-1)*ufl.exp(alpha2*((I1_isoPseudo-1)**2)
    ))

def dPsi_dI2PseudoHGO(I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, 
parameters_aniso):
    
    alpha1, alpha2 = parameters_aniso
    
    return (alpha1*(I2_isoPseudo-1)*ufl.exp(alpha2*((I2_isoPseudo-1)**2)
    ))

def dPsi_dI3PseudoHGO(I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, 
parameters_aniso):
    
    alpha1, alpha2 = parameters_aniso
    
    return (alpha1*(I3_isoPseudo-1)*ufl.exp(alpha2*((I3_isoPseudo-1)**2)
    ))

def ddPsi_ddI1PseudoHGO(I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, 
parameters_aniso):
    
    alpha1, alpha2 = parameters_aniso
    
    return (alpha1*ufl.exp(alpha2*((I1_isoPseudo-1)**2))*(1+(2*alpha2*((
    I1_isoPseudo-1)**2))))

def ddPsi_ddI2PseudoHGO(I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, 
parameters_aniso):
    
    alpha1, alpha2 = parameters_aniso
    
    return (alpha1*ufl.exp(alpha2*((I2_isoPseudo-1)**2))*(1+(2*alpha2*((
    I2_isoPseudo-1)**2))))

def ddPsi_ddI3PseudoHGO(I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, 
parameters_aniso):
    
    alpha1, alpha2 = parameters_aniso
    
    return (alpha1*ufl.exp(alpha2*((I3_isoPseudo-1)**2))*(1+(2*alpha2*((
    I3_isoPseudo-1)**2))))

def ddPsi_dI1dI2PseudoHGO(I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, 
parameters_aniso):
    
    return 0.0

def ddPsi_dI1dI3PseudoHGO(I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, 
parameters_aniso):
    
    return 0.0

def ddPsi_dI2dI3PseudoHGO(I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, 
parameters_aniso):
    
    return 0.0