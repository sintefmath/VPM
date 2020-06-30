#pragma once

#include <vector>
#include <memory>

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines
     petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets
     petscksp.h - Krylov subspace methods
     petscviewer.h - viewers
     petscpc.h  - preconditioners
*/
#include <petscksp.h>

#include <petscdm.h>
#include <petscdmda.h>

extern PetscErrorCode ComputeMatrix(KSP,Mat,Mat,void*);
extern PetscErrorCode ComputeRHS(KSP,Vec,void*);

/**
  * Types of boundary conditions supported
  */
enum Poisson2D_BoundaryConditionType {
    POISSON2D_BC_DIRICHLET=0,
    POISSON2D_BC_NEUMANN=1
};

struct PoissonSolver2D_BC
{
    Poisson2D_BoundaryConditionType north;
    Poisson2D_BoundaryConditionType east;
    Poisson2D_BoundaryConditionType south;
    Poisson2D_BoundaryConditionType west;

    PoissonSolver2D_BC( )
    {
        north = POISSON2D_BC_NEUMANN;
        east = POISSON2D_BC_NEUMANN;
        south = POISSON2D_BC_NEUMANN;
        west = POISSON2D_BC_NEUMANN;
    }

    PoissonSolver2D_BC(
            Poisson2D_BoundaryConditionType north_,
            Poisson2D_BoundaryConditionType east_,
            Poisson2D_BoundaryConditionType south_,
            Poisson2D_BoundaryConditionType west_
            )
        :
            north( north_ ),
            east( east_ ),
            south( south_ ),
            west( west_ )
    {}

};


typedef struct {
    PoissonSolver2D_BC bc_type;
    std::vector<double> rhs;
    PetscInt nx;
    PetscInt ny;
    PetscScalar dx;
    PetscScalar dy;
} UserContext;

namespace ModulesHotel
{
    class PoissonSolver2D
    {
    public:


        PoissonSolver2D(
                //int argc, char** argv,
                const PoissonSolver2D_BC boundaryconditions=PoissonSolver2D_BC()
                );
        ~PoissonSolver2D();

        /*
         * @params f, rhs of A*x = f, where A is Laplace operator
         * returns x
         */
        std::vector<double> solve(
                const std::vector<double> & f,
                const int nx,
                const int ny,
                const double dx,
                const double dy
                );
    protected:

        /*
         * for unit test only
         * returns the 2norm of the residual
         */
        double solveResidual(
                const std::vector<double> & f,
                const int nx,
                const int ny,
                const double dx,
                const double dy
                );


    private:

        bool m_initialized;
        Vec m_x;
        KSP m_ksp; /* linear solver context */
        PetscInt m_N;
        std::vector<int> m_indices;

        DM m_da;
        UserContext m_user;

        void init(const int nx, const int ny);

    };

}
