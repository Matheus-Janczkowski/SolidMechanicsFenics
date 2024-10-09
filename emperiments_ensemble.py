# Routine to test different models of hyperelasticity using Fenics

import hyperelasticity

import constitutive_models

########################################################################
#                             Saint-Venant                             #
########################################################################

# Saint Venant-Kirchhoff model parameters

E = 200E6

nu = 0.3

mu = E/(2*(1+nu))

lmbda = (E*nu)/((1+nu)*(1-2*nu))

parameters = [mu, lmbda]

# Defines the derivatives of the energy function w.r.t. the invariants 
# of the right Cauchy-Green tensor up to the second order

dPsi_dI1Model = constitutive_models.dPsi_dI1SaintVenant

dPsi_dI2Model = constitutive_models.dPsi_dI2SaintVenant

ddPsi_ddI1Model = constitutive_models.ddPsi_ddI1SaintVenant

ddPsi_ddI2Model = constitutive_models.ddPsi_ddI2SaintVenant

ddPsi_dI1dI2Model = constitutive_models.ddPsi_dI1dI2SaintVenant

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the 
# jacobian

k_volModel = constitutive_models.k_volSaintVenant

dkvol_dJModel = constitutive_models.dkvol_dJSaintVenant

# Runs the simulation

hyperelasticity.simulation_engine(parameters, dPsi_dI1Model, 
dPsi_dI2Model, ddPsi_ddI1Model, ddPsi_ddI2Model, ddPsi_dI1dI2Model, 
k_volModel, dkvol_dJModel)
