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
enum Poisson3D_BoundaryConditionType {
    POISSON3D_BC_DIRICHLET=0,
    POISSON3D_BC_NEUMANN=1
};

struct PoissonSolver3D_BC
{
    Poisson3D_BoundaryConditionType north;
    Poisson3D_BoundaryConditionType east;
    Poisson3D_BoundaryConditionType south;
    Poisson3D_BoundaryConditionType west;
    Poisson3D_BoundaryConditionType up;
    Poisson3D_BoundaryConditionType down;

    PoissonSolver3D_BC( )
    {
        north = POISSON3D_BC_NEUMANN;
        east = POISSON3D_BC_NEUMANN;
        south = POISSON3D_BC_NEUMANN;
        west = POISSON3D_BC_NEUMANN;
        up = POISSON3D_BC_NEUMANN;
        down = POISSON3D_BC_NEUMANN;
    }

    PoissonSolver3D_BC(
            Poisson3D_BoundaryConditionType north_,
            Poisson3D_BoundaryConditionType east_,
            Poisson3D_BoundaryConditionType south_,
            Poisson3D_BoundaryConditionType west_,
            Poisson3D_BoundaryConditionType up_,
            Poisson3D_BoundaryConditionType down_
            )
        :
            north( north_ ),
            east( east_ ),
            south( south_ ),
            west( west_ ),
            up( up_ ),
            down( down_ )
    {}

};


typedef struct {
    PoissonSolver3D_BC bc_type;
    std::vector<double> rhs;
    PetscInt nx;
    PetscInt ny;
    PetscInt nz;
    PetscScalar dx;
    PetscScalar dy;
    PetscScalar dz;
} UserContext;

namespace ModulesHotel
{
    class PoissonSolver3D
    {
    public:


        PoissonSolver3D(
                //int argc, char** argv,
                const PoissonSolver3D_BC boundaryconditions=PoissonSolver3D_BC()
                );
        ~PoissonSolver3D();

        /*
         * @params f, rhs of A*x = f, where A is Laplace operator
         * returns x
         */
        std::vector<double> solve(
                const std::vector<double> & f,
                const int nx,
                const int ny,
                const int nz,
                const double dx,
                const double dy,
                const double dz
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
                const int nz,
                const double dx,
                const double dy,
                const double dz
                );


    private:

        bool m_initialized;
        Vec m_x;
        KSP m_ksp; /* linear solver context */
        PetscInt m_N;
        std::vector<int> m_indices;

        DM m_da;
        UserContext m_user;

        void init(const int nx, const int ny, const int nz);

    };

}
