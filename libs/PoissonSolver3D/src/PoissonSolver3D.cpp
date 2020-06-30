#include "PoissonSolver3D.hpp"

#include <assert.h>
#include <cstring>
#include <iostream>
#include <string>
#include <algorithm>

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

    PoissonSolver3D::PoissonSolver3D(
            //int argc, char** argv,
            const PoissonSolver3D_BC boundaryconditions
            )
        :
            m_initialized(false),
            m_N(0)
    {
        m_user.bc_type = PoissonSolver3D_BC(
                boundaryconditions.north,
                boundaryconditions.east,
                boundaryconditions.south,
                boundaryconditions.west,
                boundaryconditions.up,
                boundaryconditions.down
                );
    }

    void PoissonSolver3D::init(const int nx, const int ny, const int nz)
    {

        if (m_initialized)
            return;

        std::vector<std::string> options;
            options.push_back("Split_Advection3d");
            options.push_back("-da_grid_x"); options.push_back(std::to_string(nx));
            options.push_back("-da_grid_y"); options.push_back(std::to_string(ny));
            options.push_back("-da_grid_z"); options.push_back(std::to_string(nz));
            options.push_back("-pc_type"); options.push_back("mg");
            options.push_back("-pc_mg_levels"); options.push_back("1");
            options.push_back("-mg_levels_0_pc_type"); options.push_back("ilu");
            options.push_back("-mg_levels_0_pc_factor_levels"); options.push_back("1");
            //options.push_back("-ksp_monitor");// -ksp_view

        std::vector<char*>  vc;
        std::transform(options.begin(), options.end(), std::back_inserter(vc), convert);

        int argc = options.size();
        std::cout << "PoissonSolver3D arguments:" << std::endl;
        for (int i = 1; i < argc; i++) {
            std::cout << options[i] << " ";
        }
        std::cout << std::endl;

        char help[] = "Solves the Poisson equation in 3D.\n\n";

        char ** argv = &vc[0];

        PetscInitialize(&argc,&argv,(char*)0,help);
        chkptsc(KSPCreate(PETSC_COMM_WORLD,&m_ksp));
        chkptsc(VecCreate(PETSC_COMM_WORLD,&m_x));

        DMDACreate3d(PETSC_COMM_WORLD,
                DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                DMDA_STENCIL_STAR,
                -12,-12,-12,
                PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                1,1,
                0,0,0,
                &m_da);
        KSPSetDM(m_ksp,(DM)m_da);

        //DMSetApplicationContext(m_da,&m_user);
        //KSPSetComputeRHS(m_ksp,ComputeRHS,&m_user);
        //KSPSetComputeOperators(m_ksp,ComputeMatrix,&m_user);
        //KSPSetFromOptions(m_ksp);
        //KSPSetInitialGuessNonzero(m_ksp, PETSC_TRUE);
        //KSPSetInitialGuessKnoll(m_ksp, PETSC_TRUE);

        m_N = nx*ny*nz;
        m_indices.clear();
        for (int i=0; i<m_N; i++)
        {
            m_indices.push_back(i);
        }

        m_initialized = true;

    }

    PoissonSolver3D::~PoissonSolver3D()
    {
        chkptsc(DMDestroy(&m_da));
        chkptsc(KSPDestroy(&m_ksp));
        PetscFinalize();
    };

    std::vector<double> PoissonSolver3D::solve(
            const std::vector<double> & f,
            const int nx,
            const int ny,
            const int nz,
            const double dx,
            const double dy,
            const double dz
            )
    {

        if (f.size() != nx*ny*nz)
        {
            std::cerr<<"---> PoissonSolver3D: inconsistent sizes";
            return f;
        }

        init(nx, ny, nz);

        m_user.rhs.clear();
        m_user.rhs.resize(m_N, 0);
        double dx_sq = dx*dy*dz;
        for (int i=0; i<m_N; i++)
        {
            m_user.rhs[i] = f[i]*dx_sq;
        }
        m_user.nx = nx;
        m_user.ny = ny;
        m_user.nz = nz;
        m_user.dx = dx;
        m_user.dy = dy;
        m_user.dz = dz;

        KSPSetComputeRHS(m_ksp,ComputeRHS,&m_user);
        DMSetApplicationContext(m_da,&m_user);
        KSPSetComputeOperators(m_ksp,ComputeMatrix,&m_user);
        KSPSetFromOptions(m_ksp);

        chkptsc(KSPSolve(m_ksp,NULL,NULL));

//        KSPConvergedReason reason;
//        PetscInt its;
//        KSPGetConvergedReason(m_ksp,&reason);
//          if (reason==KSP_DIVERGED_INDEFINITE_PC) {
//  printf("\nDivergence because of indefinite preconditioner;\n");
//  printf("Run the executable again but with -pc_factor_shift_positive_definite option.\n");
//} else if (reason<0) {
//  printf("\nOther kind of divergence: this should not happen.\n");
//} else {
//  KSPGetIterationNumber(m_ksp,&its);
//  printf("\nConvergence in %d iterations.\n",(int)its);
//}
//printf("\n");

        chkptsc(KSPGetSolution(m_ksp,&m_x));
        //chkptsc(VecView(m_x,PETSC_VIEWER_STDOUT_WORLD));

        std::vector<double> ret(m_N);
        chkptsc(VecGetValues(m_x,m_N,&m_indices[0],&ret[0]));

        return ret;

    }

    double PoissonSolver3D::solveResidual(
            const std::vector<double> & f,
            const int nx,
            const int ny,
            const int nz,
            const double dx,
            const double dy,
            const double dz
            )
    {
        solve(f,nx,ny,ny,dx,dy,dz);
        //
        PetscReal norm;
        Vec b,r;
        Mat J;
        KSPGetRhs(m_ksp,&b);
        KSPGetOperators(m_ksp,NULL,&J);
        VecDuplicate(b,&r);

        MatMult(J,m_x,r);

        VecAXPY(r,-1.0,b);
        VecNorm(r,NORM_2,&norm);

        PetscPrintf(PETSC_COMM_WORLD,"Residual norm %g\n",(double)norm);

        return (double)norm;
    }
}

PetscErrorCode ComputeRHS(KSP ksp,Vec b,void *ctx)
{
    UserContext    *user = (UserContext*)ctx;
    PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
    PetscScalar    ***array;
    DM             da;

    KSPGetDM(ksp,&da);
    DMDAGetInfo(da, 0, &mx, &my, &mz, 0,0,0,0,0,0,0,0,0);

    DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);
    DMDAVecGetArray(da, b, &array);
    for (k=zs; k<zs+zm; k++)
    {
        for (j=ys; j<ys+ym; j++)
        {
            for (i=xs; i<xs+xm; i++)
            {
                array[k][j][i] = user->rhs[i+user->nx*(j+user->ny*k)];
            }
        }
    }
    DMDAVecRestoreArray(da, b, &array);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    /* force right hand side to be consistent for singular matrix */
    /* note this is really a hack, normally the model would provide you with a consistent right handside */
    if (user->bc_type.north == POISSON3D_BC_NEUMANN) {
        MatNullSpace nullspace;

        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
        MatNullSpaceRemove(nullspace,b);
        MatNullSpaceDestroy(&nullspace);
    }
    //VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    return(0);
}

PetscErrorCode ComputeMatrix(KSP ksp,Mat J, Mat jac,void *ctx)
{
    UserContext    *user = (UserContext*)ctx;
    PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,num, numi, numj, numk;
    PetscScalar    v[7],DyDzdDx,DxDzdDy,DxDydDz;
    MatStencil     row, col[7];
    DM             da;
    PetscErrorCode ierr;

    KSPGetDM(ksp,&da);
    DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);
    DyDzdDx = user->dy*user->dz/user->dx;
    DxDzdDy = user->dx*user->dz/user->dy;
    DxDydDz = user->dx*user->dy/user->dz;
    DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);

    for (k=zs; k<zs+zm; k++)
    {
        for (j=ys; j<ys+ym; j++)
        {
            for (i=xs; i<xs+xm; i++)
            {
                row.i = i; row.j = j; row.k = k;
                if (i==0 || j==0 || k==0 || i==mx-1 || j==my-1 || k==mz-1)
                {
                    num = 0; numi=0; numj=0; numk=0;
                    if (k!=0)
                    {
                        v[num]     = -DxDydDz;
                        col[num].i = i;
                        col[num].j = j;
                        col[num].k = k-1;
                        num++; numk++;
                    }
                    if (j!=0)
                    {
                        v[num]     = -DxDzdDy;
                        col[num].i = i;
                        col[num].j = j-1;
                        col[num].k = k;
                        num++; numj++;
                    }
                    if (i!=0)
                    {
                        v[num]     = -DyDzdDx;
                        col[num].i = i-1;
                        col[num].j = j;
                        col[num].k = k;
                        num++; numi++;
                    }
                    if (i!=mx-1)
                    {
                        v[num]     = -DyDzdDx;
                        col[num].i = i+1;
                        col[num].j = j;
                        col[num].k = k;
                        num++; numi++;
                    }
                    if (j!=my-1)
                    {
                        v[num]     = -DxDzdDy;
                        col[num].i = i;
                        col[num].j = j+1;
                        col[num].k = k;
                        num++; numj++;
                    }
                    if (k!=mz-1)
                    {
                        v[num]     = -DxDydDz;
                        col[num].i = i;
                        col[num].j = j;
                        col[num].k = k+1;
                        num++; numk++;
                    }
                    v[num]     = (PetscReal)(numk)*DxDydDz + (PetscReal)(numj)*DxDzdDy + (PetscReal)(numi)*DyDzdDx;
                    col[num].i = i;   col[num].j = j;   col[num].k = k;
                    num++;
                    ierr = MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);CHKERRQ(ierr);
                }
                else
                {
                    v[0] = -DxDydDz;                          col[0].i = i;   col[0].j = j;   col[0].k = k-1;
                    v[1] = -DxDzdDy;                          col[1].i = i;   col[1].j = j-1; col[1].k = k;
                    v[2] = -DyDzdDx;                          col[2].i = i-1; col[2].j = j;   col[2].k = k;
                    v[3] = 2.0*(DyDzdDx + DxDzdDy + DxDydDz); col[3].i = i;   col[3].j = j;   col[3].k = k;
                    v[4] = -DyDzdDx;                          col[4].i = i+1; col[4].j = j;   col[4].k = k;
                    v[5] = -DxDzdDy;                          col[5].i = i;   col[5].j = j+1; col[5].k = k;
                    v[6] = -DxDydDz;                          col[6].i = i;   col[6].j = j;   col[6].k = k+1;
                    ierr = MatSetValuesStencil(jac,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);
                }
            }
        }
    }

    MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
    if (user->bc_type.north == POISSON3D_BC_NEUMANN)
    {
        MatNullSpace nullspace;

        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);
        MatSetNullSpace(J,nullspace);
        MatNullSpaceDestroy(&nullspace);
    }
    //MatView(jac,PETSC_VIEWER_STDOUT_WORLD);
    return(0);
}
