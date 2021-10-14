#include "DiffusionEquationSolver2D.hpp"

#include <assert.h>
#include <cstring>
#include <iostream>
#include <string>
#include <algorithm>
#define MY_CHKERRQ(ierr) do { \
    PetscErrorCode ierr__ = (ierr); \
    if (PetscUnlikely(ierr__)) { \
        auto error = PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr__,PETSC_ERROR_REPEAT," "); \
        std::stringstream ss; \
        ss << "Petsc error at " << __FILE__ << ":" << __LINE__ << std::endl; \
        throw std::runtime_error(ss.str()); \
    } \
} while (0)

namespace ModulesHotel
{
    /*
     * just wrappers for error checking PETSc
     */
    namespace {
        int chkptscwrp(PetscErrorCode e)
        {
            CHKERRQ(e);
            return 0;
        }
        void chkptsc(PetscErrorCode e)
        {
            chkptscwrp(e);
        }
        char *convert(const std::string & s)
        {
               char *pc = new char[s.size()+1];
                  std::strcpy(pc, s.c_str());
                     return pc;
        }

    }

    DiffusionEquationSolver2D::DiffusionEquationSolver2D(
            //int argc, char** argv,
            const DiffusionEquationSolver2D_BC boundaryconditions
            )
        :
            m_initialized(false),
            m_N(0)
    {
        m_user.bc_type = DiffusionEquationSolver2D_BC(
                boundaryconditions.north,
                boundaryconditions.east,
                boundaryconditions.south,
                boundaryconditions.west
                );
    }

    void DiffusionEquationSolver2D::init(const int nx, const int ny)
    {

        if (m_initialized)
            return;

        std::vector<std::string> options;
            options.push_back("DiffusionEquation2d");
            options.push_back("-da_grid_x"); options.push_back(std::to_string(nx));
            options.push_back("-da_grid_y"); options.push_back(std::to_string(ny));
            options.push_back("-pc_type"); options.push_back("mg");
            options.push_back("-pc_mg_levels"); options.push_back("1");
            options.push_back("-mg_levels_0_pc_type"); options.push_back("ilu");
            options.push_back("-mg_levels_0_pc_factor_levels"); options.push_back("1");
        std::vector<char*>  vc;
        std::transform(options.begin(), options.end(), std::back_inserter(vc), convert);

        int argc = options.size();
        std::cout << "DiffusionEquationSolver2D arguments:" << std::endl;
        for (int i = 1; i < argc; i++) {
            std::cout << options[i] << " ";
        }
        std::cout << std::endl;

        char help[] = "Solves the DiffusionEquation equation in 2D.\n\n";

        char ** argv = &vc[0];

        PetscInitialize(&argc,&argv,(char*)0,help);
        MY_CHKERRQ(KSPCreate(PETSC_COMM_WORLD,&m_ksp));
        MY_CHKERRQ(VecCreate(PETSC_COMM_WORLD,&m_x));

        DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,nx,ny,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&m_da);
        MY_CHKERRQ(DMSetFromOptions(m_da));
        MY_CHKERRQ(DMSetUp(m_da));
        MY_CHKERRQ(KSPSetDM(m_ksp,(DM)m_da));

        //DMSetApplicationContext(m_da,&m_user);
        //KSPSetComputeRHS(m_ksp,ComputeRHS_DIFFEQ2D,&m_user);
        //KSPSetComputeOperators(m_ksp,ComputeMatrix_DIFFEQ2D,&m_user);
        //KSPSetFromOptions(m_ksp);

        m_N = nx*ny;
        m_indices.clear();
        for (int i=0; i<m_N; i++)
        {
            m_indices.push_back(i);
        }

        m_initialized = true;

    }

    DiffusionEquationSolver2D::~DiffusionEquationSolver2D()
    {
        if (!m_initialized)
            return;
        MY_CHKERRQ(DMDestroy(&m_da));
        MY_CHKERRQ(KSPDestroy(&m_ksp));
        PetscFinalize();
    };

    void DiffusionEquationSolver2D::solve(
            std::vector<double> & f,
            const int nx,
            const int ny,
            const double dx,
            const double dy,
            const double nu,
            const double dt
            )
    {

        if (f.size() != nx*ny)
        {
            std::cerr<<"---> DiffusionEquationSolver2D: inconsistent sizes";
        }

        init(nx, ny);

        m_user.rhs.clear();
        m_user.rhs.resize(m_N, 0);
        double dx_sq = dx*dy;
        for (int i=0; i<m_N; i++)
        {
            m_user.rhs[i] = f[i];
        }
        m_user.nx = nx;
        m_user.ny = ny;
        m_user.dx = dx;
        m_user.dy = dy;
        m_user.lambda = nu*dt/(dx*dx);


        DMSetApplicationContext(m_da,&m_user);

        KSPSetComputeRHS(m_ksp,ComputeRHS_DIFFEQ2D,&m_user);
        KSPSetComputeOperators(m_ksp,ComputeMatrix_DIFFEQ2D,&m_user);
        KSPSetFromOptions(m_ksp);
        MY_CHKERRQ(KSPSolve(m_ksp,NULL,NULL));
        MY_CHKERRQ(KSPGetSolution(m_ksp,&m_x));
        //MY_CHKERRQ(VecView(m_x,PETSC_VIEWER_STDOUT_WORLD));

        MY_CHKERRQ(VecGetValues(m_x,m_N,&m_indices[0],&f[0]));

    }

    double DiffusionEquationSolver2D::solveResidual(
            const std::vector<double> & f,
            const int nx,
            const int ny,
            const double dx,
            const double dy,
            const double nu,
            const double dt
            )
    {
    //    solve(f,nx,ny,dx,dy, nu, dt);
    //    //
    //    PetscReal norm;
    //    Vec b,r;
    //    Mat J;
    //    KSPGetRhs(m_ksp,&b);
    //    KSPGetOperators(m_ksp,NULL,&J);
    //    VecDuplicate(b,&r);

    //    MatMult(J,m_x,r);
    //    VecAXPY(r,-1.0,b);
    //    VecNorm(r,NORM_2,&norm);

    //    PetscPrintf(PETSC_COMM_WORLD,"Residual norm %g\n",(double)norm);

    //    return (double)norm;
        return 0;
    }
}

PetscErrorCode ComputeRHS_DIFFEQ2D(KSP ksp,Vec b,void *ctx)
{
    UserContext    *user = (UserContext*)ctx;
    PetscInt       i,j,M,N,xm,ym,xs,ys;
    PetscScalar    **array;
    DM             da;

    KSPGetDM(ksp,&da);
    DMDAGetInfo(da, 0, &M, &N, 0,0,0,0,0,0,0,0,0,0);

    DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0); /* Fine grid */
    //printf(" M N: %d %d; xm ym: %d %d; xs ys: %d %d\n",M,N,xm,ym,xs,ys);
    DMDAVecGetArray(da, b, &array);
    for (j=ys; j<ys+ym; j++)
    {
        for (i=xs; i<xs+xm; i++)
        {
            array[j][i] = user->rhs[i+user->nx*j];
        }
    }
    DMDAVecRestoreArray(da, b, &array);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    /* force right hand side to be consistent for singular matrix */
    /* note this is really a hack, normally the model would provide you with a consistent right handside */
    if (user->bc_type.north == DIFFEQ2D_BC_NEUMANN) {
        MatNullSpace nullspace;

        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
        MatNullSpaceRemove(nullspace,b);
        MatNullSpaceDestroy(&nullspace);
    }
    //VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    return(0);
}

PetscErrorCode ComputeMatrix_DIFFEQ2D(KSP ksp,Mat J, Mat jac,void *ctx)
{
    UserContext    *user = (UserContext*)ctx;
    PetscInt       i, j, M, N, xm, ym, xs, ys, num, numi, numj;
    PetscScalar    v[5], DydDx, DxdDy, lambda;
    MatStencil     row, col[5];
    DM             da;

    KSPGetDM(ksp,&da);
    DMDAGetInfo(da,0,&M,&N,0,0,0,0,0,0,0,0,0,0);
    DxdDy = user->dx/user->dy;
    DydDx = user->dy/user->dx;
    lambda = user->lambda;
    DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);
    for (j=ys; j<ys+ym; j++)
    {
        for (i=xs; i<xs+xm; i++)
        {
            row.i = i; row.j = j;

            if (i==0 || j==0 || i==M-1 || j==N-1)
            {
                if (user->bc_type.north != DIFFEQ2D_BC_NEUMANN) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"Only Neumann boundary conditions are supported !\n");
                else if (user->bc_type.north == DIFFEQ2D_BC_NEUMANN)
                {
                    num=0; numi=0; numj=0;
                    if (j!=0)
                    {
                        v[num] = -lambda;
                        col[num].i = i;
                        col[num].j = j-1;
                        num++; numj++;
                    }
                    if (i!=0)
                    {
                        v[num] = -lambda;
                        col[num].i = i-1;
                        col[num].j = j;
                        num++; numi++;
                    }
                    if (i!=M-1)
                    {
                        v[num] = -lambda;
                        col[num].i = i+1;
                        col[num].j = j;
                        num++; numi++;
                    }
                    if (j!=N-1)
                    {
                        v[num] = -lambda;
                        col[num].i = i;
                        col[num].j = j+1;
                        num++; numj++;
                    }
                    v[num] = 1+(PetscReal)(numi+numj)*lambda;
                    col[num].i = i;
                    col[num].j = j;
                    num++;
                    MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);
                }
            }
            else
            {
                v[0] = -lambda;              col[0].i = i;   col[0].j = j-1;
                v[1] = -lambda;              col[1].i = i-1; col[1].j = j;
                v[2] = 1+4.0*lambda;         col[2].i = i;   col[2].j = j;
                v[3] = -lambda;              col[3].i = i+1; col[3].j = j;
                v[4] = -lambda;              col[4].i = i;   col[4].j = j+1;
                MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);
            }
        }
    }
    MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
    if (user->bc_type.north == DIFFEQ2D_BC_NEUMANN)
    {
        MatNullSpace nullspace;

        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
        MatSetNullSpace(J,nullspace);
        MatNullSpaceDestroy(&nullspace);
    }
    //MatView(jac,PETSC_VIEWER_STDOUT_WORLD);
    return(0);
}
