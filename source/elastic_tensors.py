from fenics import *

import ufl_legacy as ufl

i,j,k,l = ufl.indices(4)

# Defines a function to perform the operation of double contraction be-
# tween a fourth order tensor and a second order tensor

def double_contraction(C,A):

    return as_tensor(C[i,j,k,l]*A[k,l], (i,j))

# Defines a function to perform the tensor product of a vector by itself

def tensor_vectorByItself(v):

    return as_tensor(v[i]*v[j], (i,j))

# Defines a function to perform the tensor product between to second or-
# der tensors

def tensor_2ndBy2nd(A,B):

    return as_tensor(A[i,j]*B[k,l],(i,j,k,l))

# Defines a function to perform the symmetric part of the tensor product 
# between to second order tensors

def tensor_sym2ndBy2nd(A,B):

    return as_tensor(0.5*((A[i,j]*B[k,l])+(A[k,l]*B[i,j])),(i,j,k,l))

# Defines a function to evaluate the fourth order elastic tensor using
# hyperelastic models that do not use isochoric invariants of the right
# Cauchy-Green tensor

def C_hyperNonIsochoric(C, I1, kvol, J, dkvol_dJ, dPsi_dI2,
ddPsi_ddI1, ddPsi_ddI2, ddPsi_dI1dI2, dimension):
    
    # Initiates the Kronecker's delta
    
    delta = Identity(dimension)
    
    # Calculates the inverse of the right Cauchy-Green tensor, C
    
    C_inv = inv(C)

    # Calculates the fourth order identity tensor
    
    II = as_tensor(delta[i,k]*delta[j,l], (i,j,k,l))

    # Calculates the ball product of the inverse of C (see Holzapfel(2000)
    # page 43 for further information)
    
    C1 = as_tensor((-0.5*((C_inv[i,k]*C_inv[l,j])+(C_inv[i,l]*C_inv[k,j]
    ))),(i,j,k,l))

    # Calculates the tensorial product of the inverse of C by itself
    
    C2 = as_tensor(C_inv[i,j]*C_inv[k,l], (i,j,k,l))

    # Calculates the tensorial product of C by itself
    
    C3 = as_tensor(C[i,j]*C[k,l], (i,j,k,l))

    # Calculates the symmetric part of the tensorial product of the se-
    # cond order identity tensor by C
    
    C4 = as_tensor(0.5*((delta[i,j]*C[k,l])+(delta[k,l]*C[i,j])), 
    (i,j,k,l))

    # Calculates the tensorial product of the second order identity ten-
    # sor by itself
    
    IxI = as_tensor(delta[i,j]*delta[k,l], (i,j,k,l))

    ### Calculas the coefficients that multiply the tensor evaluated a-
    ### bove

    # Calculates the coefficient of the tensorial product of the second
    # order identity tensor by itself
    
    coeff_IxI = (dPsi_dI2+ddPsi_ddI1+(2*I1*ddPsi_dI1dI2)+(I1*I1*
    ddPsi_ddI2))

    # Calculates the coefficient of the symmetric part of the tensorial
    # product of the second order identity tensor by C

    coeff_IxC = -2*(ddPsi_dI1dI2+(I1*ddPsi_ddI2))

    # Sums everything up and returns it

    return ((4*coeff_IxI*IxI)-(4*coeff_IxC*C4)+(4*ddPsi_ddI2*C3)-(4*
    dPsi_dI2*II)+(2*J*dkvol_dJ*C2)+(4*kvol*C1))

# Defines a function to evaluate the fourth order elastic tensor using
# hyperelastic models that do use isochoric invariants of the right 
# Cauchy-Green tensor

def C_hyperIsochoric(C, I1_iso, I2_iso, kvol, J, dkvol_dJ, dPsi_dI1Iso,
dPsi_dI2Iso, ddPsi_ddI1Iso, ddPsi_ddI2Iso, ddPsi_dI1dI2Iso, dimension):
    
    # Initiates the Kronecker's delta
    
    delta = Identity(dimension)
    
    # Calculates the inverse of the right Cauchy-Green tensor, C
    
    C_inv = inv(C)

    # Calculates the fourth order identity tensor
    
    II = as_tensor(delta[i,k]*delta[j,l], (i,j,k,l))

    # Calculates the ball product of the inverse of C (see Holzapfel(2000)
    # page 43 for further information)
    
    Cinv_ballCinv = as_tensor((-0.5*((C_inv[i,k]*C_inv[l,j])+(C_inv[i,l]
    *C_inv[k,j]))),(i,j,k,l))

    # Calculates the tensorial product of the inverse of C by itself
    
    Cinv_xCinv = as_tensor(C_inv[i,j]*C_inv[k,l], (i,j,k,l))

    # Calculates the tensorial product of C by itself
    
    CxC = as_tensor(C[i,j]*C[k,l], (i,j,k,l))

    # Calculates the symmetric part of the tensorial product of the se-
    # cond order identity tensor by C
    
    sym_CxI = as_tensor(0.5*((delta[i,j]*C[k,l])+(delta[k,l]*C[i,j])), 
    (i,j,k,l))

    # Calculates the symmetric part of the tensorial product of the in-
    # verse of C by the second order identity tensor

    sym_CinvxI = as_tensor(0.5*((delta[i,j]*C_inv[k,l])+(delta[k,l]*
    C_inv[i,j])),(i,j,k,l))

    # Calculates the symmetric part of the tensorial product of C by its
    # inverse

    sym_CxCinv = as_tensor(0.5*((C_inv[i,j]*C[k,l])+(C_inv[k,l]*C[i,j])), 
    (i,j,k,l))

    # Calculates the tensorial product of the second order identity ten-
    # sor by itself
    
    IxI = as_tensor(delta[i,j]*delta[k,l], (i,j,k,l))

    ### Calculas the coefficients that multiply the tensor evaluated a-
    ### bove

    J_23 = (J**(-2/3))

    # Calculates the coefficient for the fourth order identity tensor

    coeff_II = -J_23*J_23*dPsi_dI2Iso

    # Calculates the coefficient of the tensorial product of the second
    # order identity tensor by itself
    
    coeff_IxI = J_23*J_23*(dPsi_dI2Iso+ddPsi_ddI1Iso+(2*I1_iso*
    ddPsi_dI1dI2Iso)+(I1_iso*I1_iso*ddPsi_ddI2Iso))

    # Calculates the coefficient of the symmetric part of the tensorial
    # product of the second order identity tensor by the inverse of C

    coeff_symIxCinv = -((2*J_23)/3)*(dPsi_dI1Iso+(2*I1_iso*dPsi_dI2Iso)+
    +(((2*I2_iso)+(I1_iso*I1_iso))*ddPsi_dI1dI2Iso)+(I1_iso*
    ddPsi_ddI1Iso)+(2*I1_iso*I2_iso*ddPsi_ddI2Iso))

    # Calculates the coefficient of the symmetric part of the tensorial
    # product of the second order identity tensor by C

    coeff_symIxC = -(2/(J*J))*(ddPsi_dI1dI2Iso+(I1_iso*ddPsi_ddI2Iso))

    # Calculates the coefficient of the symmetric part of the tensorial
    # product of C by its inverse

    coeff_symCxCinv = ((2*J_23*J_23)/3)*((2*dPsi_dI2Iso)+(I1_iso*
    ddPsi_dI1dI2Iso)+(2*I2_iso*ddPsi_ddI2Iso))

    # Calculates the coefficient of the tensorial product of the inverse
    # of C by itself

    coeff_CinvxCinv = (((1/9)*((I1_iso*dPsi_dI1Iso)+(4*I2_iso*
    dPsi_dI2Iso)+(4*I1_iso*I2_iso*ddPsi_dI1dI2Iso)+(I1_iso*I1_iso*
    ddPsi_ddI1Iso)+(4*I2_iso*I2_iso*ddPsi_ddI2Iso)))+(0.5*J*dkvol_dJ))

    # Calculates the coefficient of the tensorial product of C by itself

    coeff_CxC = -J_23*J_23*J_23*J_23*ddPsi_ddI2Iso

    # Calculates the coefficient of the ball product of the inverse of C
    # by itself

    coeff_CinvBallCinv = (((-1/3)*((I1_iso*dPsi_dI1Iso)+(2*I2_iso*
    dPsi_dI2Iso)))+(kvol))

    # Sums everything up and returns it

    return (4*((coeff_II*II)+(coeff_IxI*IxI)+(coeff_symIxCinv*sym_CinvxI
    )+(coeff_symIxC*sym_CxI)+(coeff_symCxCinv*sym_CxCinv)+(
    coeff_CinvxCinv*Cinv_xCinv)+(coeff_CxC*CxC)+(coeff_CinvBallCinv*
    Cinv_ballCinv)))

# Defines a function to evaluate the fourth order elastic tensor using
# hyperelastic models that do use isochoric invariants of the right 
# Cauchy-Green tensor. Besides, the energy function has a contribution
# for the anisotropic behavior

def C_hyperIsochoricAnisotropic(C, I1_iso, I2_iso, kvol, J, dkvol_dJ, 
dPsi_dI1Iso, dPsi_dI2Iso, ddPsi_ddI1Iso, ddPsi_ddI2Iso, ddPsi_dI1dI2Iso,
dimension, I1_isoPseudo, I2_isoPseudo, I3_isoPseudo, dPsi_dI1IsoPseudo, 
dPsi_dI2IsoPseudo, dPsi_dI3IsoPseudo, ddPsi_ddI1IsoPseudo, 
ddPsi_ddI2IsoPseudo, ddPsi_ddI3IsoPseudo, ddPsi_dI1dI2IsoPseudo, 
ddPsi_dI1dI3IsoPseudo, ddPsi_dI2dI3IsoPseudo, a1xa1, a2xa2, a3xa3):
    
    # Initiates the Kronecker's delta
    
    delta = Identity(dimension)
    
    # Calculates the inverse of the right Cauchy-Green tensor, C
    
    C_inv = inv(C)

    # Calculates the fourth order identity tensor
    
    II = as_tensor(delta[i,k]*delta[j,l], (i,j,k,l))

    # Calculates the ball product of the inverse of C (see Holzapfel(2000)
    # page 43 for further information)
    
    Cinv_ballCinv = as_tensor((-0.5*((C_inv[i,k]*C_inv[l,j])+(C_inv[i,l]
    *C_inv[k,j]))),(i,j,k,l))

    # Calculates the tensorial product of the inverse of C by itself
    
    Cinv_xCinv = as_tensor(C_inv[i,j]*C_inv[k,l], (i,j,k,l))

    # Calculates the tensorial product of C by itself
    
    CxC = as_tensor(C[i,j]*C[k,l], (i,j,k,l))

    # Calculates the symmetric part of the tensorial product of the se-
    # cond order identity tensor by C
    
    sym_CxI = as_tensor(0.5*((delta[i,j]*C[k,l])+(delta[k,l]*C[i,j])), 
    (i,j,k,l))

    # Calculates the symmetric part of the tensorial product of the in-
    # verse of C by the second order identity tensor

    sym_CinvxI = as_tensor(0.5*((delta[i,j]*C_inv[k,l])+(delta[k,l]*
    C_inv[i,j])),(i,j,k,l))

    # Calculates the symmetric part of the tensorial product of C by its
    # inverse

    sym_CxCinv = as_tensor(0.5*((C_inv[i,j]*C[k,l])+(C_inv[k,l]*C[i,j])), 
    (i,j,k,l))

    # Calculates the tensorial product of the second order identity ten-
    # sor by itself
    
    IxI = as_tensor(delta[i,j]*delta[k,l], (i,j,k,l))

    ### Calculas the coefficients that multiply the tensor evaluated a-
    ### bove

    J_23 = (J**(-2/3))

    # Calculates the coefficient for the fourth order identity tensor

    coeff_II = -J_23*J_23*dPsi_dI2Iso

    # Calculates the coefficient of the tensorial product of the second
    # order identity tensor by itself
    
    coeff_IxI = J_23*J_23*(dPsi_dI2Iso+ddPsi_ddI1Iso+(2*I1_iso*
    ddPsi_dI1dI2Iso)+(I1_iso*I1_iso*ddPsi_ddI2Iso))

    # Calculates the coefficient of the symmetric part of the tensorial
    # product of the second order identity tensor by the inverse of C

    coeff_symIxCinv = -((2*J_23)/3)*(dPsi_dI1Iso+(2*I1_iso*dPsi_dI2Iso)+
    +(((2*I2_iso)+(I1_iso*I1_iso))*ddPsi_dI1dI2Iso)+(I1_iso*
    ddPsi_ddI1Iso)+(2*I1_iso*I2_iso*ddPsi_ddI2Iso))

    # Calculates the coefficient of the symmetric part of the tensorial
    # product of the second order identity tensor by C

    coeff_symIxC = -(2/(J*J))*(ddPsi_dI1dI2Iso+(I1_iso*ddPsi_ddI2Iso))

    # Calculates the coefficient of the symmetric part of the tensorial
    # product of C by its inverse

    coeff_symCxCinv = ((2*J_23*J_23)/3)*((2*dPsi_dI2Iso)+(I1_iso*
    ddPsi_dI1dI2Iso)+(2*I2_iso*ddPsi_ddI2Iso))

    # Calculates the coefficient of the tensorial product of the inverse
    # of C by itself

    coeff_CinvxCinv = (((1/9)*((I1_iso*dPsi_dI1Iso)+(4*I2_iso*
    dPsi_dI2Iso)+(4*I1_iso*I2_iso*ddPsi_dI1dI2Iso)+(I1_iso*I1_iso*
    ddPsi_ddI1Iso)+(4*I2_iso*I2_iso*ddPsi_ddI2Iso)))+(0.5*J*dkvol_dJ))

    # Calculates the coefficient of the tensorial product of C by itself

    coeff_CxC = -J_23*J_23*J_23*J_23*ddPsi_ddI2Iso

    # Calculates the coefficient of the ball product of the inverse of C
    # by itself

    coeff_CinvBallCinv = (((-1/3)*((I1_iso*dPsi_dI1Iso)+(2*I2_iso*
    dPsi_dI2Iso)))+(kvol))

    ####################################################################
    #                     Anisotropic contribution                     #
    ####################################################################

    # Calculates the summation IPseudo_i*dPsi_dIPseudo_i

    sum_IidPsidIi = ((I1_isoPseudo*dPsi_dI1IsoPseudo)+(I2_isoPseudo*
    dPsi_dI2IsoPseudo)+(I3_isoPseudo*dPsi_dI3IsoPseudo))

    # Calculates the summation IPseudo_i*IPseudo_j*ddPsi_dIidIj

    sum_IiIjddPsidIidIj = ((I1_isoPseudo*((I1_isoPseudo*
    ddPsi_ddI1IsoPseudo)+(I2_isoPseudo*ddPsi_dI1dI2IsoPseudo)+(
    I3_isoPseudo*ddPsi_dI1dI3IsoPseudo)))+(I2_isoPseudo*((I1_isoPseudo*
    ddPsi_ddI2IsoPseudo)+(I2_isoPseudo*ddPsi_dI1dI2IsoPseudo)+(
    I3_isoPseudo*ddPsi_dI2dI3IsoPseudo)))+(
    I3_isoPseudo*((I1_isoPseudo*ddPsi_ddI3IsoPseudo)+(I2_isoPseudo*
    ddPsi_dI1dI3IsoPseudo)+(I3_isoPseudo*ddPsi_dI2dI3IsoPseudo))))

    # Updates the coefficients of the tensor C_inv x C_inv

    coeff_CinvxCinv += ((1/9)*(sum_IidPsidIi+sum_IiIjddPsidIidIj))

    # Updates the coefficients of the tensor C_inv ball C_inv

    coeff_CinvBallCinv -= ((1/3)*sum_IidPsidIi)

    # Calculates the symmetric part of C_inv x a_i x a_i

    Cinv_xaixai = ((((dPsi_dI1IsoPseudo+(I1_isoPseudo*ddPsi_ddI1IsoPseudo)
    +(I2_isoPseudo*ddPsi_dI1dI2Iso)+(I3_isoPseudo*ddPsi_dI1dI3IsoPseudo)
    )*2*tensor_sym2ndBy2nd(C_inv, a1xa1))+((dPsi_dI2IsoPseudo+(I1_isoPseudo*
    ddPsi_dI1dI2Iso)+(I2_isoPseudo*ddPsi_ddI2IsoPseudo)+(I3_isoPseudo*
    ddPsi_dI2dI3IsoPseudo))*2*tensor_sym2ndBy2nd(C_inv, a2xa2))+((
    dPsi_dI3IsoPseudo+(I1_isoPseudo*ddPsi_dI1dI3IsoPseudo)+(I2_isoPseudo*
    ddPsi_dI2dI3IsoPseudo)+(I3_isoPseudo*ddPsi_ddI3IsoPseudo))*2*
    tensor_sym2ndBy2nd(C_inv, a3xa3)))*(-J_23/3))
    
    # Calculates the tensor a_i x a_i x a_j x a_j

    aixai_xajxaj = (J_23*J_23*((ddPsi_ddI1IsoPseudo*tensor_2ndBy2nd(
    a1xa1,a1xa1))+(ddPsi_ddI2IsoPseudo*tensor_2ndBy2nd(a2xa2,a2xa2))+(
    ddPsi_ddI3IsoPseudo*tensor_2ndBy2nd(a3xa3,a3xa3))+(2*ddPsi_dI1dI2Iso
    *tensor_sym2ndBy2nd(a1xa1,a2xa2))+(2*ddPsi_dI1dI3IsoPseudo*
    tensor_sym2ndBy2nd(a1xa1,a3xa3))+(2*ddPsi_dI2dI3IsoPseudo*
    tensor_sym2ndBy2nd(a2xa2,a3xa3))))

    # Calculates the tensor 

    # Sums everything up and returns it

    return (4*((coeff_II*II)+(coeff_IxI*IxI)+(coeff_symIxCinv*sym_CinvxI
    )+(coeff_symIxC*sym_CxI)+(coeff_symCxCinv*sym_CxCinv)+(
    coeff_CinvxCinv*Cinv_xCinv)+(coeff_CxC*CxC)+(coeff_CinvBallCinv*
    Cinv_ballCinv)+aixai_xajxaj+Cinv_xaixai))