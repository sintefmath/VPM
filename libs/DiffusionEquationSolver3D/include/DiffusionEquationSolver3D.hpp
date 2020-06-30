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

extern PetscErrorCode ComputeMatrix_DIFFEQ3D(KSP,Mat,Mat,void*);
extern PetscErrorCode ComputeRHS_DIFFEQ3D(KSP,Vec,void*);

/**
  * Types of boundary conditions supported
  */
enum DiffusionEquation3D_BoundaryConditionType {
    DIFFEQ3D_BC_DIRICHLET=0,
    DIFFEQ3D_BC_NEUMANN=1
};

struct DiffusionEquationSolver3D_BC
{
    DiffusionEquation3D_BoundaryConditionType north;
    DiffusionEquation3D_BoundaryConditionType east;
    DiffusionEquation3D_BoundaryConditionType south;
    DiffusionEquation3D_BoundaryConditionType west;
    DiffusionEquation3D_BoundaryConditionType up;
    DiffusionEquation3D_BoundaryConditionType down;

    DiffusionEquationSolver3D_BC( )
    {
        north = DIFFEQ3D_BC_NEUMANN;
        east = DIFFEQ3D_BC_NEUMANN;
        south = DIFFEQ3D_BC_NEUMANN;
        west = DIFFEQ3D_BC_NEUMANN;
        up = DIFFEQ3D_BC_NEUMANN;
        down = DIFFEQ3D_BC_NEUMANN;
    }

    DiffusionEquationSolver3D_BC(
            DiffusionEquation3D_BoundaryConditionType north_,
            DiffusionEquation3D_BoundaryConditionType east_,
            DiffusionEquation3D_BoundaryConditionType south_,
            DiffusionEquation3D_BoundaryConditionType west_,
            DiffusionEquation3D_BoundaryConditionType up_,
            DiffusionEquation3D_BoundaryConditionType down_
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
    DiffusionEquationSolver3D_BC bc_type;
    std::vector<double> rhs;
    PetscInt nx;
    PetscInt ny;
    PetscInt nz;
    PetscScalar dx;
    PetscScalar dy;
    PetscScalar dz;
    PetscScalar lambda;
} UserContext;

namespace ModulesHotel
{
    class DiffusionEquationSolver3D
    {
    public:


        DiffusionEquationSolver3D(
                //int argc, char** argv,
                const DiffusionEquationSolver3D_BC boundaryconditions=DiffusionEquationSolver3D_BC()
                );
        ~DiffusionEquationSolver3D();

        /*
         */
        void solve(
                std::vector<double> & f,
                const int nx,
                const int ny,
                const int nz,
                const double dx,
                const double dy,
                const double dz,
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
                const int nz,
                const double dx,
                const double dy,
                const double dz,
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

        void init(const int nx, const int ny, const int nz);

    };

}
