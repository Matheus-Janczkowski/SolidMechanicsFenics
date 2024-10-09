from dolfin import *

import ufl_legacy as ufl

i,j,k,l = ufl.indices(4)

import elastic_tensors

import newton_raphson_methods

# Defines a function to be the engine of the simulation

def simulation_engine(parameters, dPsi_dI1Model, dPsi_dI2Model,
ddPsi_ddI1Model, ddPsi_ddI2Model, ddPsi_dI1dI2Model, k_volModel, 
dkvol_dJModel):

    ####################################################################
    #                               Mesh                               #
    ####################################################################

    # Defines the degree of the interpolation polynomial 

    polynomial_degree = 1

    element_type = 'P'

    # Defines the mesh

    z_max = 2.0

    mesh = BoxMesh(Point(0,0,0), Point(1,1,z_max), 16,16,48)

    # Defines the solution space using polynomial interpolation functions

    V = VectorFunctionSpace(mesh, element_type, polynomial_degree)

    ####################################################################
    #                      Newton-Raphson settings                     #
    ####################################################################

    # Sets the final time

    final_time = 1.0

    # Sets the number of pseudotime points for loading

    n_pseudoTimePoints = 5

    # Sets the step in pseudo time

    delta_t = final_time/n_pseudoTimePoints

    # Defines the maximum number of iterations

    n_maxIter = 40

    # Defines a pseudotime variable to control the incremental process

    t = 0.0

    # Defines the norm of the residual to be used as stopping criterion
    # to the Newton-Raphson method

    residual_tol = 1E-4

    ####################################################################
    #                  Dirichlet boundary conditions                   #
    ####################################################################

    # Defines domain subregions

    fixed_edge = CompiledSubDomain("near(x[0], side, tol) && near(x[1]"+
    ", side, tol) && near(x[2], side, tol) && on_boundary", side=0.0, 
    tol=1e-2)

    below =  CompiledSubDomain("near(x[2], side) && on_boundary", side=
    0.0)

    above = CompiledSubDomain("near(x[2], side) && on_boundary", side=
    z_max)

    # Defines the expressions for the Dirichlet boundary conditions, 
    # i.e., their values for the displacement field

    expression_below = Expression("0.0", degree=polynomial_degree)

    expression_fixedEdge = Expression(("0.0", "0.0", "0.0"), degree=
    polynomial_degree)

    # Defines two objects of Dirichlet boundary conditions for the dis-
    # placement field. Constraints the z direction only. For the fixed 
    # node, applies pointwise method, for it is a point

    bc_DirichletBelowDisplacement = DirichletBC(V.sub(2), 
    expression_below, below)

    bc_DirichletEdgeDisplacement = DirichletBC(V, expression_fixedEdge,
    fixed_edge, method="pointwise")

    def BCs_dirichletDisplacement(t):
        
        return [bc_DirichletBelowDisplacement, 
        bc_DirichletEdgeDisplacement]

    # The update direction must be null in the regions where Dirichlet 
    # boundary conditions are applied to the displacement field

    bc_DirichletBelowDirection = DirichletBC(V.sub(2), Constant("0.0"), 
    below)

    bc_DirichletEdgeDirection = DirichletBC(V, Constant(("0.0", "0.0", 
    "0.0")), fixed_edge, method="pointwise")

    BCs_dirichletDirection = [bc_DirichletBelowDirection, 
    bc_DirichletEdgeDirection]

    ####################################################################
    #                   Neumann boundary conditions                    #
    ####################################################################

    # Creates a mesh function for the boundaries to apply Neumann bound-
    # ary conditions. The dimension of the integration over the boundary 
    # is one less than the domain integration. Besides, the last argu-
    # ment, 0, is to initialize the mesh with this value; so it substi-
    # tutes the command set_all(0)

    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 
    0)

    # Tags the labels for the subregions in the boundary

    above.mark(boundaries, 1)

    # Creates the surface area differential for the above facet

    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    # Defines the maximum traction

    T_max = 1000000.0

    # Defines the boundary traction vector in the reference configura-
    # tion

    T = Expression(("0.0", "0.0", "(t/t_final)*T_max"), t=0.0, T_max=
    T_max, t_final=final_time, degree=0)

    ####################################################################
    #               Definition of solutions and variation              #
    ####################################################################

    # Defines the function for the update direction in the Newton-Raphson 
    # scheme

    Delta_u = TrialFunction(V)

    # Defines the variation

    delta_u = TestFunction(V)

    # Defines the displacement field

    u = Function(V)

    # Defines the body forces vector in the reference configuration

    B = Expression(("0.0", "0.0", "0.0"), t=0.0, degree=0)

    ####################################################################
    #                        Tensors definition                        #
    ####################################################################

    # Recovers the problem's dimension

    dimension = Delta_u.geometric_dimension()

    # Initializes the second order identity tensor

    I = Identity(dimension)

    # Calculates the deformation gradient tensor in terms of the solution
    # u

    F = grad(u)+I

    # Calculates the right Cauchy-Green tensor

    C = (F.T)*F

    # Calculates the first two invariants of the right Cauchy-Green ten-
    # sor

    I1 = tr(C)

    I2 = 0.5*((tr(C)**2)-tr(C*C))

    # Calculates the jacobian of the kinematic mapping from the refe-
    # rence to the deformed configurations

    J = det(F)

    ####################################################################
    #                      Constitutive modelling                      #
    ####################################################################

    # Defines the derivatives of the energy function w.r.t. the invari- 
    # ants of the right Cauchy-Green tensor up to the second order

    dPsi_dI1 = dPsi_dI1Model(I1, I2, parameters)

    dPsi_dI2 = dPsi_dI2Model(I1, I2, parameters)

    ddPsi_ddI1 = ddPsi_ddI1Model(I1, I2, parameters)

    ddPsi_ddI2 = ddPsi_ddI2Model(I1, I2, parameters)

    ddPsi_dI1dI2 = ddPsi_dI1dI2Model(I1, I2, parameters)

    # Defines the k factor for the volumetric part of the second Piola-
    # Kirchhoff stress tensor and the derivative of the former w.r.t. 
    # the jacobian

    k_vol = k_volModel(J, parameters)

    dkvol_dJ = dkvol_dJModel(J, parameters)

    # Calculates the second Piola-Kirchhoff stress tensor

    S = (2*(((dPsi_dI1+(I1*dPsi_dI2))*I)-(dPsi_dI2*C)+(k_vol*inv(C))))

    # Defines the Cauchy stress tensor by the push-forward of the second 
    # Piola-Kirchhoff stress tensor

    sigma = (1/J)*F*S*F.T

    # Evaluates the deviatoric parcel of the Cauchy stress tensor

    dev_sigma = sigma-((1/3)*(tr(sigma))*I)

    # Defines a function to evaluate the von Mises stress

    von_mises = sqrt((3/2)*inner(dev_sigma,dev_sigma))

    # Defines the fourth order elastic tensor

    C_elastic = elastic_tensors.C_hyperNonIsochoric(C, I1, k_vol, J,
    dkvol_dJ, dPsi_dI2, ddPsi_ddI1, ddPsi_ddI2, ddPsi_dI1dI2, dimension)

    ####################################################################
    #             Variational form and solution procedure              #
    ####################################################################

    # Defines the linearized variational form in terms of the update di-
    # rection, a==L. Note that a is constituted by the directional deri-
    # vative of the variation of the power/work; while L is formed by 
    # the residual itself. This scheme states the Newton-Raphson method 
    # into a variational form

    L = (-inner((F.T)*grad(delta_u),S)*dx+(dot(B,delta_u)*dx)+(dot(T,
    delta_u)*ds(1)))

    # Analytical derivative of the fourth order elastic tensor

    a = ((inner((F.T)*grad(delta_u),elastic_tensors.double_contraction(
    C_elastic,sym((F.T)*grad(Delta_u))))*dx))+(inner(grad(delta_u)*S,
    grad(Delta_u))*dx)

    # Numerical derivative of the residue entire. Use it to validate the
    # analytical derivative only

    #a = -derivative(L, u, Delta_u)

    Delta_u = Function(V)

    # Initializes the solution plotting file

    file_displacement = File("displacement.pvd")

    # Initializes the solution file for the von Mises stress

    V_vonMises = FunctionSpace(mesh, element_type, polynomial_degree)

    file_vonMises = File("von_misesStress.pvd")

    # Performs Newton-Raphson method and evaluate displacement and von 
    # Mises stress

    newton_raphson_methods.Newton_RaphsonVonMises(u, a, L, T, B, 
    BCs_dirichletDisplacement, BCs_dirichletDirection, t, delta_t, 
    n_maxIter, n_pseudoTimePoints, residual_tol, V_vonMises, von_mises, 
    file_displacement, file_vonMises)