# Routine to test different models of hyperelasticity using Fenics

import numpy as np

import hyperelasticity_nonIsochoric

import hyperelasticity_isochoric

import hyperelasticity_isochoricAnisotropic

import isotropic_constitutive_models

import anisotropic_constitutive_models

# Sets the maximum traction in Pa to be applied

T_max = -35E6

########################################################################
#                    Models to be experimented with                    #
########################################################################

# Defines a set of flags to inform which constitutive model to test

flag_SaintVenant = False

flag_compressibleNeoHookean = False

flag_uncoupledMooneyRivlin = False

flag_uncoupledNeoHookean = False

flag_HolzapfelGasserOgden = True

########################################################################
#                             Saint-Venant                             #
########################################################################

print("###############################################################"+
"#########\n#                             Saint-Venant                "+
"             #\n#####################################################"+
"###################\n\n")

# Defines the names of the files to save the solutions of displacement
# and von Mises stress

file_displacementName = "resultados/displacement_saint_venant.pvd" 

file_vonMisesName = "resultados/von_mises_stress_saint_venant.pvd"

# Model parameters

E = 200E6

nu = 0.3

mu = E/(2*(1+nu))

lmbda = (E*nu)/((1+nu)*(1-2*nu))

parameters = [mu, lmbda]

# Defines the derivatives of the energy function w.r.t. the invariants 
# of the right Cauchy-Green tensor up to the second order

dPsi_dI1Model = isotropic_constitutive_models.dPsi_dI1SaintVenant

dPsi_dI2Model = isotropic_constitutive_models.dPsi_dI2SaintVenant

ddPsi_ddI1Model = isotropic_constitutive_models.ddPsi_ddI1SaintVenant

ddPsi_ddI2Model = isotropic_constitutive_models.ddPsi_ddI2SaintVenant

ddPsi_dI1dI2Model = isotropic_constitutive_models.ddPsi_dI1dI2SaintVenant

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the 
# jacobian

k_volModel = isotropic_constitutive_models.k_volSaintVenant

dkvol_dJModel = isotropic_constitutive_models.dkvol_dJSaintVenant

# Runs the simulation

if flag_SaintVenant:

    hyperelasticity_nonIsochoric.simulation_engine(parameters, 
    dPsi_dI1Model, dPsi_dI2Model, ddPsi_ddI1Model, ddPsi_ddI2Model, 
    ddPsi_dI1dI2Model, k_volModel, dkvol_dJModel, file_displacementName, 
    file_vonMisesName, T_max)

else:

    print("Do not evaluate the Saint Venant model\n\n")

########################################################################
#                       Compressible Neo-Hookean                       #
########################################################################

print("###############################################################"+
"#########\n#                       Compressible Neo-Hookean          "+
"             #\n#####################################################"+
"###################\n\n")

# Defines the names of the files to save the solutions of displacement
# and von Mises stress

file_displacementName = "resultados/displacement_compressible_neo_hookean.pvd" 

file_vonMisesName = "resultados/von_mises_stress_compressible_neo_hookean.pvd"

# Model parameters

E = 200E6

nu = 0.3

mu = E/(2*(1+nu))

lmbda = (E*nu)/((1+nu)*(1-2*nu))

parameters = [mu, 0.5*lmbda]

# Defines the derivatives of the energy function w.r.t. the invariants 
# of the right Cauchy-Green tensor up to the second order

dPsi_dI1Model = isotropic_constitutive_models.dPsi_dI1CompressibleNeoHookean

dPsi_dI2Model = isotropic_constitutive_models.dPsi_dI2CompressibleNeoHookean

ddPsi_ddI1Model = isotropic_constitutive_models.ddPsi_ddI1CompressibleNeoHookean

ddPsi_ddI2Model = isotropic_constitutive_models.ddPsi_ddI2CompressibleNeoHookean

ddPsi_dI1dI2Model = isotropic_constitutive_models.ddPsi_dI1dI2CompressibleNeoHookean

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the 
# jacobian

k_volModel = isotropic_constitutive_models.k_volCompressibleNeoHookean

dkvol_dJModel = isotropic_constitutive_models.dkvol_dJCompressibleNeoHookean

# Runs the simulation

if flag_compressibleNeoHookean:

    hyperelasticity_nonIsochoric.simulation_engine(parameters, 
    dPsi_dI1Model, dPsi_dI2Model, ddPsi_ddI1Model, ddPsi_ddI2Model, 
    ddPsi_dI1dI2Model, k_volModel, dkvol_dJModel, file_displacementName, 
    file_vonMisesName, T_max)

else:

    print("Do not evaluate the compressible Neo Hookean model\n\n")

########################################################################
#                       Uncoupled Mooney-Rivlin                        #
########################################################################

print("###############################################################"+
"#########\n#                       Uncoupled Mooney-Rivlin           "+
"             #\n#####################################################"+
"###################\n\n")

# Defines the names of the files to save the solutions of displacement
# and von Mises stress

file_displacementName = "resultados/displacement_mooney_rivlin.pvd" 

file_vonMisesName = "resultados/von_mises_stress_mooney_rivlin.pvd"

# Model parameters

E = 200E6

nu = 0.3

mu_3 = E/(2*(1+nu))

G = (E*nu)/((1+nu)*(1-2*nu))

mu_1 = 0.25*G

mu_2 = 0.25*G

parameters = [mu_1, mu_2, mu_3]

# Defines the derivatives of the energy function w.r.t. the invariants 
# of the right Cauchy-Green tensor up to the second order

dPsi_dI1Model = isotropic_constitutive_models.dPsi_dI1MooneyRivlin

dPsi_dI2Model = isotropic_constitutive_models.dPsi_dI2MooneyRivlin

ddPsi_ddI1Model = isotropic_constitutive_models.ddPsi_ddI1MooneyRivlin

ddPsi_ddI2Model = isotropic_constitutive_models.ddPsi_ddI2MooneyRivlin

ddPsi_dI1dI2Model = isotropic_constitutive_models.ddPsi_dI1dI2MooneyRivlin

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the 
# jacobian

k_volModel = isotropic_constitutive_models.k_volMooneyRivlin

dkvol_dJModel = isotropic_constitutive_models.dkvol_dJMooneyRivlin

# Runs the simulation

if flag_uncoupledMooneyRivlin:

    hyperelasticity_isochoric.simulation_engine(parameters, 
    dPsi_dI1Model, dPsi_dI2Model, ddPsi_ddI1Model, ddPsi_ddI2Model, 
    ddPsi_dI1dI2Model, k_volModel, dkvol_dJModel, file_displacementName, 
    file_vonMisesName, T_max)

else:

    print("Do not evaluate the uncoupled Mooney-Rivlin model\n\n")

########################################################################
#                        Uncoupled Neo-Hookean                         #
########################################################################

print("###############################################################"+
"#########\n#                        Uncoupled Neo-Hookean            "+
"             #\n#####################################################"+
"###################\n\n")

# Defines the names of the files to save the solutions of displacement
# and von Mises stress

file_displacementName = "resultados/displacement_uncoupled_neo_hookean.pvd" 

file_vonMisesName = "resultados/von_mises_stress_uncoupled_neo_hookean.pvd"

# Model parameters

E = 200E6

nu = 0.3

bulk = E/(2*(1+nu))

G = (E*nu)/((1+nu)*(1-2*nu))

mu_1 = 0.5*G

parameters = [mu_1, bulk]

# Defines the derivatives of the energy function w.r.t. the invariants 
# of the right Cauchy-Green tensor up to the second order

dPsi_dI1Model = isotropic_constitutive_models.dPsi_dI1UncoupledNeoHookean

dPsi_dI2Model = isotropic_constitutive_models.dPsi_dI2UncoupledNeoHookean

ddPsi_ddI1Model = isotropic_constitutive_models.ddPsi_ddI1UncoupledNeoHookean

ddPsi_ddI2Model = isotropic_constitutive_models.ddPsi_ddI2UncoupledNeoHookean

ddPsi_dI1dI2Model = isotropic_constitutive_models.ddPsi_dI1dI2UncoupledNeoHookean

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the 
# jacobian

k_volModel = isotropic_constitutive_models.k_volUncoupledNeoHookean

dkvol_dJModel = isotropic_constitutive_models.dkvol_dJUncoupledNeoHookean

# Runs the simulation

if flag_uncoupledNeoHookean:

    hyperelasticity_isochoric.simulation_engine(parameters, 
    dPsi_dI1Model, dPsi_dI2Model, ddPsi_ddI1Model, ddPsi_ddI2Model, 
    ddPsi_dI1dI2Model, k_volModel, dkvol_dJModel, file_displacementName, 
    file_vonMisesName, T_max)

else:

    print("Do not evaluate the uncoupled Neo-Hookean model\n\n")

########################################################################
#                        Holzapfel-Gasser-Ogden                        #
########################################################################

print("###############################################################"+
"#########\n#                        Holzapfel-Gasser-Ogden           "+
"             #\n#####################################################"+
"###################\n\n")

# Defines the names of the files to save the solutions of displacement
# and von Mises stress

file_displacementName = "resultados/displacement_HGO.pvd" 

file_vonMisesName = "resultados/von_mises_stress_HGO.pvd"

# Model parameters

E = 200E6

nu = 0.3

bulk = E/(2*(1+nu))

G = (E*nu)/((1+nu)*(1-2*nu))

mu_1 = 0.5*G

parameters = [mu_1, bulk]

# Anisotropic parameters

alpha1 = 4E7

alpha2 = 1.0

parameters_aniso = [alpha1, alpha2]

# Defines the directions of the fibers. An unitary  vector for each row.
# If there is not 3 directions of fibers, simply put vectors of zeros

fibers_vectors = [[0.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

norm_a = [np.sqrt((fibers_vectors[0][0]**2)+(fibers_vectors[0][1]**2)+(
fibers_vectors[0][2]**2)), np.sqrt((fibers_vectors[1][0]**2)+(
fibers_vectors[1][1]**2)+(fibers_vectors[1][2]**2)), np.sqrt((
fibers_vectors[2][0]**2)+(fibers_vectors[2][1]**2)+(
fibers_vectors[2][2]**2))]

for i in range(3):

    if norm_a[i]>0:

        for j in range(3):

            fibers_vectors[i][j] = fibers_vectors[i][j]/norm_a[i]

# Defines the derivatives of the energy function w.r.t. the invariants 
# of the right Cauchy-Green tensor up to the second order

dPsi_dI1Model = isotropic_constitutive_models.dPsi_dI1UncoupledNeoHookean

dPsi_dI2Model = isotropic_constitutive_models.dPsi_dI2UncoupledNeoHookean

ddPsi_ddI1Model = isotropic_constitutive_models.ddPsi_ddI1UncoupledNeoHookean

ddPsi_ddI2Model = isotropic_constitutive_models.ddPsi_ddI2UncoupledNeoHookean

ddPsi_dI1dI2Model = isotropic_constitutive_models.ddPsi_dI1dI2UncoupledNeoHookean

# Defines the k factor for the volumetric part of the second Piola-
# Kirchhoff stress tensor and the derivative of the former w.r.t. the 
# jacobian

k_volModel = isotropic_constitutive_models.k_volUncoupledNeoHookean

dkvol_dJModel = isotropic_constitutive_models.dkvol_dJUncoupledNeoHookean

# Defines the derivatives of the energy function w.r.t. the pseudo inva-
# riants

dPsi_dI1PseudoModel = anisotropic_constitutive_models.dPsi_dI1PseudoHGO

dPsi_dI2PseudoModel = anisotropic_constitutive_models.dPsi_dI2PseudoHGO

dPsi_dI3PseudoModel = anisotropic_constitutive_models.dPsi_dI3PseudoHGO

ddPsi_ddI1PseudoModel = anisotropic_constitutive_models.ddPsi_ddI1PseudoHGO

ddPsi_ddI2PseudoModel = anisotropic_constitutive_models.ddPsi_ddI2PseudoHGO

ddPsi_ddI3PseudoModel = anisotropic_constitutive_models.ddPsi_ddI3PseudoHGO

ddPsi_dI1dI2PseudoModel = anisotropic_constitutive_models.ddPsi_dI1dI2PseudoHGO

ddPsi_dI1dI3PseudoModel = anisotropic_constitutive_models.ddPsi_dI1dI3PseudoHGO

ddPsi_dI2dI3PseudoModel = anisotropic_constitutive_models.ddPsi_dI2dI3PseudoHGO

# Runs the simulation

if flag_HolzapfelGasserOgden:

    hyperelasticity_isochoricAnisotropic.simulation_engine(parameters, 
    fibers_vectors, dPsi_dI1Model, dPsi_dI2Model, ddPsi_ddI1Model, 
    ddPsi_ddI2Model, ddPsi_dI1dI2Model, dPsi_dI1PseudoModel, 
    dPsi_dI2PseudoModel, dPsi_dI3PseudoModel, ddPsi_ddI1PseudoModel, 
    ddPsi_ddI2PseudoModel, ddPsi_ddI3PseudoModel, 
    ddPsi_dI1dI2PseudoModel, ddPsi_dI1dI3PseudoModel, 
    ddPsi_dI2dI3PseudoModel, parameters_aniso, k_volModel, dkvol_dJModel, 
    file_displacementName, file_vonMisesName, T_max)

else:

    print("Do not evaluate the HGO model\n\n")