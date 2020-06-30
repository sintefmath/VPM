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

extern PetscErrorCode ComputeMatrix_DIFFEQ2D(KSP,Mat,Mat,void*);
extern PetscErrorCode ComputeRHS_DIFFEQ2D(KSP,Vec,void*);

/**
  * Types of boundary conditions supported
  */
enum DiffusionEquation2D_BoundaryConditionType {
    DIFFEQ2D_BC_DIRICHLET=0,
    DIFFEQ2D_BC_NEUMANN=1
};

struct DiffusionEquationSolver2D_BC
{
    DiffusionEquation2D_BoundaryConditionType north;
    DiffusionEquation2D_BoundaryConditionType east;
    DiffusionEquation2D_BoundaryConditionType south;
    DiffusionEquation2D_BoundaryConditionType west;

    DiffusionEquationSolver2D_BC( )
    {
        north = DIFFEQ2D_BC_NEUMANN;
        east = DIFFEQ2D_BC_NEUMANN;
        south = DIFFEQ2D_BC_NEUMANN;
        west = DIFFEQ2D_BC_NEUMANN;
    }

    DiffusionEquationSolver2D_BC(
            DiffusionEquation2D_BoundaryConditionType north_,
            DiffusionEquation2D_BoundaryConditionType east_,
            DiffusionEquation2D_BoundaryConditionType south_,
            DiffusionEquation2D_BoundaryConditionType west_
            )
        :
            north( north_ ),
            east( east_ ),
            south( south_ ),
            west( west_ )
    {}

};


typedef struct {
    DiffusionEquationSolver2D_BC bc_type;
    std::vector<double> rhs;
    PetscInt nx;
    PetscInt ny;
    PetscScalar dx;
    PetscScalar dy;
    PetscScalar lambda;
} UserContext;

namespace ModulesHotel
{
    class DiffusionEquationSolver2D
    {
    public:


        DiffusionEquationSolver2D(
                //int argc, char** argv,
                const DiffusionEquationSolver2D_BC boundaryconditions=DiffusionEquationSolver2D_BC()
                );
        ~DiffusionEquationSolver2D();

        /*
         */
        void solve(
                std::vector<double> & f,
                const int nx,
                const int ny,
                const double dx,
                const double dy,
                const double nu,
                const double dt
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
                const double dy,
                const double nu,
                const double dt
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
