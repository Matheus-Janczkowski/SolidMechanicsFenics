# Routine to store Newton-Raphson methods

from dolfin import *

from petsc4py import PETSc

# Defines a function to perform the Newton-Raphson method for displa-
# cement of a hyperelastic media. Evalautes the von Mises stress too

def Newton_RaphsonVonMises(u, a, L, T, B, BCs_dirichletDisplacement,
BCs_dirichletDirection, t, delta_t, n_maxIter, n_pseudoTimePoints,
residual_tol, V_vonMises, von_mises, file_displacement, file_vonMises):

    # Iterates through the pseudotime points

    for j in range(n_pseudoTimePoints):

        # Updates time

        t += delta_t

        print("\n#####################################################"+
        "####\n#           Initiates the", j+1,"pseudotime point      "+
        "      #\n####################################################"+
        "#####\n")

        print("t =", t)

        # Updates the traction force

        T.t = t

        # Updates the body vector force

        B.t = t

        # Applies boundary conditions to the displacement field

        BCs_dirichletU = BCs_dirichletDisplacement(t)

        for CC in BCs_dirichletU:

            CC.apply(u.vector())

        # Iterates through the loop of the Newton-Raphson method

        for i in range(n_maxIter):

            print("\nInitiates iteration", i+1, "of the Newton-Raphson"+
            " method\n")

            # Solves the linearized variational form
            
            Kt, R = assemble_system(a, L, BCs_dirichletDirection)

            ### PETSc solver

            # Gets the tangent matrix

            Kt_matrix = as_backend_type(Kt).mat()

            vectorDelta_u, R_vector = Kt_matrix.getVecs()

            # Gets the residual vector

            R_vector.setValues(range(R.size()), R)

            # Sets the Krylov linear solver. Uses conjugate gradient 
            # ('cg')

            krylov_linearSolver = PETSc.KSP().create()

            krylov_linearSolver.setType('cg')

            krylov_linearSolver.setOperators(Kt_matrix)

            krylov_linearSolver.solve(R_vector, vectorDelta_u)

            # Evaluates the norm of the residual

            norm_residual = norm(R, "l2")

            print("Residual norm is", norm_residual)

            # Verifies stop criterion

            if norm_residual<residual_tol:

                print("Residual criterion met at iteration", i)

                break

            # Updates the displacement field
            
            u.vector()[:] = u.vector()[:]+vectorDelta_u

        print("\n\n")

        # Saves the solution into the VTK file

        file_displacement << u

        # Calculates the von Mises stress

        von_mises = project(von_mises, V_vonMises)

        file_vonMises << von_mises