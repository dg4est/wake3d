//
//  visit_interface.c
//  cdriver
//
//  Created by Zhi Yang on 1/16/17.
//


#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "driver.h"
#include "alloc.h"

extern driver_t *d;

typedef struct basis
{
    int ndim;
    int nrefine;
    int nfrac;
    double *frac;
    double *dfracdr, *dfracds, *dfracdt;
    double *r;
    double *phi;
    double *dphi;
} basis_t;

void visit_load_dynamic_library() {
  visit_t *v = d->visit;
  char *filename = v->visit_so_file;
  char *error;

  if(v->visit_flag == 0) return;

  /* get group handle */
  v->group_handle = dlopen (filename, RTLD_LAZY);
  if (!v->group_handle) {
      if(d->rank == 0) fputs(dlerror(),stderr);
      v->visit_flag = 0;
      return;
  }

  /* load each function */
  *(void **)(&v->visit_extract_ini) = dlsym(v->group_handle, "visit_extract_ini");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&v->visit_extract_this) = dlsym(v->group_handle, "visit_extract_this");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}
}

void visit_close_dynamic_library(){
  if(d->visit->group_handle) dlclose(d->visit->group_handle);
}

void refine_mesh(flow_solver_t *f, visit_t *v)
{
  int pmax = 10; /*max p*/
  basis_t *bp;
  int dimfrac = (pmax+1)*(pmax+1)*(pmax+1); /*param [in] dimfrac large integer (set inside tioga to be 11^3=1331 for a p=10 hex) used for allocation frac if necessary  */
  double frac[dimfrac];  /*param [out] frac basis evaluated at point \a rst_in length \a nfrac */
  int nfrac_p[dimfrac+1];
  // adjust below array for the refinement for each p=1...10


  //int refinement[] = {1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; /*size of pmax+1*/
  int refinement[] = {1, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4}; /*size of pmax+1*/
  double xyz;    /*param [in]  xyz physical location (not used anymore since \a rst_in is always passed in)*/
  int nfrac;     /*param [out] nfrac number of modes in the basis used for interpolation */
  int index;     /*param [out] index index to access solution array index is base 1 */
  double rst[]={0.,0.,0.}; /*param [in]  rst_in natural coordinates input to evaulate the basis at  (-1<r,s,t<1, length 3) */
//  double r[refinement[pmax+1]+1];
//  double phi[refinement[pmax+1]+1];
//  double dphi[refinement[pmax+1]+1];
  double drdx, dsdy, dtdz;
  int i,j,k,im,i1,p,ndim,ih,ih1,n,m;
  int ii,jj,kk,ind1,ind2;
  //double strain[3][3],omega[3][3];

  bp = my_alloc(basis_t, pmax);

  for(p=0; p<pmax; p++)
  {

    int p1 = p+1;
    int nmode = p1*p1*p1;
    int ndim = refinement[p]+1;

    bp[p].ndim = ndim;
    bp[p].nfrac = nmode;
    bp[p].nrefine = ndim*ndim*ndim;

    bp[p].frac    = my_alloc(double, nmode*bp[p].nrefine);
    bp[p].dfracdr = my_alloc(double, nmode*bp[p].nrefine);
    bp[p].dfracds = my_alloc(double, nmode*bp[p].nrefine);
    bp[p].dfracdt = my_alloc(double, nmode*bp[p].nrefine);

    bp[p].r    = my_alloc(double, ndim);
    bp[p].phi  = my_alloc(double, p1*ndim);
    bp[p].dphi = my_alloc(double, p1*ndim);

    for (i=0; i<ndim; i++)
    {
        bp[p].r[i] = -1. + 2.*((double) i)/((double) (ndim-1));
    }

    (f->basis_and_derivative)(&p, &ndim, bp[p].r, bp[p].phi, bp[p].dphi);

    ind1=-1;
    for(ii=0;ii<ndim;ii++)
      for(jj=0;jj<ndim;jj++)
        for(kk=0;kk<ndim;kk++)
        {
          ind1++;
          ind2=-1;
          for(k=0;k<p1;k++)
            for(j=0;j<p1;j++)
              for(i=0;i<p1;i++)
              {
                ind2++;
                bp[p].frac[ind1*nmode+ind2]    = bp[p].phi[ii*p1+i] *bp[p].phi[jj*p1+j] *bp[p].phi[kk*p1+k];
                bp[p].dfracdr[ind1*nmode+ind2] = bp[p].dphi[ii*p1+i]*bp[p].phi[jj*p1+j] *bp[p].phi[kk*p1+k];
                bp[p].dfracds[ind1*nmode+ind2] = bp[p].phi[ii*p1+i] *bp[p].dphi[jj*p1+j]*bp[p].phi[kk*p1+k];
                bp[p].dfracdt[ind1*nmode+ind2] = bp[p].phi[ii*p1+i] *bp[p].phi[jj*p1+j] *bp[p].dphi[kk*p1+k];
              }
        }
  }


// map the number of modes to p
  for (i=0; i<dimfrac+1; i++)
    nfrac_p[i] = 0;

  for (i=1; i<pmax+1; i++)
  {
    nfrac_p[(i+1)*(i+1)*(i+1)] = i;
  }

  v->nnode = 0;
  v->nhex  = 0;
  for (i=0; i<f->nhex; i++)
  {
//     if(f->iblank_cell[i] == 1)
   {
     i1 = i+1;
     f->create_donor_frac(&i1, &xyz, &nfrac, &index, frac, rst, &dimfrac);
     p = nfrac_p[nfrac];
     if (p >= 0)
     {
       ndim = refinement[p];
       v->nnode = v->nnode + (ndim+1)*(ndim+1)*(ndim+1);
       v->nhex  = v->nhex  + ndim*ndim*ndim;
     }
     else
     {
       printf("ERROR nfract: %d", nfrac);
       return;
     }
    }
  }

  if(v->hex) my_free(v->hex);
  if(v->coord) my_free(v->coord);
  if(v->iblank) my_free(v->iblank);
  v->hex    = my_alloc(int,    v->nhex*8);
  v->coord  = my_alloc(double, v->nnode*3);
  v->iblank = my_alloc(int,    v->nnode);

  if(v->var)
  {
    for(i=0; i<v->nvar; i++)
    {
      if(v->var[i]) free(v->var[i]);
    }
    free(v->var);
  }
  v->var = (double**)malloc(v->nvar * sizeof(double*));
  for(i=0; i<v->nvar; i++)
  {
     v->var[i] = (double*)malloc(v->nnode *sizeof(double));
  }

  int inode = 0;
  int ihex  = 0;
  for (ih=0; ih<f->nhex; ih++)
  {
//    if(f->iblank_cell[ih] == 1)
    {
      ih1=ih+1;
      rst[0]=0.; rst[1]=0.; rst[2]=0.;
      (f->create_donor_frac)(&ih1, &xyz, &nfrac, &index, frac, rst, &dimfrac);
      index--;
      int p = nfrac_p[nfrac];
      int ndim = refinement[p];
      int iplus1 = (ndim+1)*(ndim+1);
      int jplus1 =  ndim+1;
      int kplus1 =       1;

      int n0 = f->ndc8[ih*8+0]-1;
      int n6 = f->ndc8[ih*8+6]-1;
      double x0 = f->xgeom[n0*3+0];
      double y0 = f->xgeom[n0*3+1];
      double z0 = f->xgeom[n0*3+2];
      double x6 = f->xgeom[n6*3+0];
      double y6 = f->xgeom[n6*3+1];
      double z6 = f->xgeom[n6*3+2];

      drdx = 2./(x6-x0);
      dsdy = 2./(y6-y0);
      dtdz = 2./(z6-z0);

      int in=inode;
      ind1 = -1;
      for (i=0; i<ndim+1; i++)
      {
        for (j=0; j<ndim+1; j++)
        {
          for (k=0; k<ndim+1; k++)
          {

            ind1++;
            v->iblank[in] = f->iblank_cell[ih];

            v->coord[in*3+0] = x0 + 0.5*(bp[p].r[i]+1.)*(x6-x0);
            v->coord[in*3+1] = y0 + 0.5*(bp[p].r[j]+1.)*(y6-y0);
            v->coord[in*3+2] = z0 + 0.5*(bp[p].r[k]+1.)*(z6-z0);

            v->var[0][in] = 0.;
            v->var[1][in] = 0.;
            v->var[2][in] = 0.;
            v->var[3][in] = 0.;
            v->var[4][in] = 0.;

            v->var[5][in] = v->iblank[i];

            double dqdx[5] = {0.,0.,0.,0.,0.};
            double dqdy[5] = {0.,0.,0.,0.,0.};
            double dqdz[5] = {0.,0.,0.,0.,0.};

            for (im=0; im<nfrac; im++)
            {
              for(n=0;n<5;++n){
                v->var[n][in] +=bp[p].frac[ind1*nfrac+im]*f->soln[index+im*5+n];
                dqdx[n] += drdx*bp[p].dfracdr[ind1*nfrac+im]*f->soln[index+im*5+n];
                dqdy[n] += dsdy*bp[p].dfracds[ind1*nfrac+im]*f->soln[index+im*5+n];
                dqdz[n] += dtdz*bp[p].dfracdt[ind1*nfrac+im]*f->soln[index+im*5+n];
              }
            }
            // rho inverse for below
            double rho1 = 1.0/v->var[0][in];
            // velocity magnitude
            v->var[5][in] = sqrt(v->var[1][in] * v->var[1][in] +
                                 v->var[2][in] * v->var[2][in] +
                                 v->var[3][in] * v->var[3][in]) * rho1;
            // velocity gradients
            double dudx = (dqdx[1] - v->var[1][in]*rho1*dqdx[0])*rho1;
            double dvdx = (dqdx[2] - v->var[2][in]*rho1*dqdx[0])*rho1;
            double dwdx = (dqdx[3] - v->var[3][in]*rho1*dqdx[0])*rho1;

            double dudy = (dqdy[1] - v->var[1][in]*rho1*dqdy[0])*rho1;
            double dvdy = (dqdy[2] - v->var[2][in]*rho1*dqdy[0])*rho1;
            double dwdy = (dqdy[3] - v->var[3][in]*rho1*dqdy[0])*rho1;

            double dudz = (dqdz[1] - v->var[1][in]*rho1*dqdz[0])*rho1;
            double dvdz = (dqdz[2] - v->var[2][in]*rho1*dqdz[0])*rho1;
            double dwdz = (dqdz[3] - v->var[3][in]*rho1*dqdz[0])*rho1;

            /*
            strain[0][0] = dudx;
            strain[0][1] = 0.5*(dudy + dvdx);
            strain[0][2] = 0.5*(dudz + dwdx);
            strain[1][0] = 0.5*(dudy + dvdx);
            strain[1][1] = dvdy;
            strain[1][2] = 0.5*(dvdz + dwdy);
            strain[2][0] = 0.5*(dudz + dwdx);
            strain[2][1] = 0.5*(dvdz + dwdy);
            strain[2][2] = dwdz;

            omega[0][0] = 0.0;
            omega[0][1] = 0.5*(dudy - dvdx);
            omega[0][2] = 0.5*(dudz - dwdx);
            omega[1][0] = 0.5*(dvdx - dudy);
            omega[1][1] = 0.0;
            omega[1][2] = 0.5*(dvdz - dwdy);
            omega[2][0] = 0.5*(dwdx - dudz);
            omega[2][1] = 0.5*(dwdy - dvdz);
            omega[2][2] = 0.0;

            // calculate frobenius norm
            double snorm = 0.0;
            double onorm = 0.0;
            for(n=0;n<3;n++){
              for(m=0;m<3;m++){
                snorm += strain[n][m]*strain[n][m];
                onorm += omega[n][m]*omega[n][m];
              }
            }

            //snorm = sqrt(snorm);
            //onorm = sqrt(onorm);
            */

            // vorticity magnitude
            v->var[6][in] = sqrt((dwdy - dvdz)*(dwdy - dvdz) + (dudz - dwdx)*(dudz - dwdx) + (dvdx - dudy)*(dvdx - dudy));

            double grad_vel[3][3];
            grad_vel[0][0] = dudx; grad_vel[0][1] = dvdx; grad_vel[0][2] = dwdx;
            grad_vel[1][0] = dudy; grad_vel[1][1] = dvdy; grad_vel[1][2] = dwdy;
            grad_vel[2][0] = dudz; grad_vel[2][1] = dvdz; grad_vel[2][2] = dwdz;

            double S,S_tot;
            double R,R_tot;
            S_tot = R_tot = 0.0;
            for (n = 0; n < 3; n++) {
              for (m = 0; m < 3; m++) {
                S = 0.5*(grad_vel[n][m] + grad_vel[m][n]); /* shear rate */
                R = 0.5*(grad_vel[n][m] - grad_vel[m][n]); /* rotation rate */
                S_tot += S*S;
                R_tot += R*R;
              }
            }

            // Q-criterion
            v->var[7][in] = 0.5*(R_tot/(S_tot + 1.0E-14) - 1.0); //0.5*Q_loc; //0.5*(onorm/snorm - 1.0);

            in++;
          }
        }
      }

      int ic = ihex;
      for (i=0; i<ndim; i++)
      {
        for(j=0; j<ndim; j++)
        {
          for(k=0; k<ndim; k++)
          {
            int n0 = inode + i*(ndim+1)*(ndim+1)+j*(ndim+1)+k +1;
            v->hex[ic*8+0] = n0;
            v->hex[ic*8+1] = n0 + iplus1;
            v->hex[ic*8+2] = n0 + iplus1 + jplus1;
            v->hex[ic*8+3] = n0          + jplus1;
            v->hex[ic*8+4] = n0                   + kplus1;
            v->hex[ic*8+5] = n0 + iplus1          + kplus1;
            v->hex[ic*8+6] = n0 + iplus1 + jplus1 + kplus1;
            v->hex[ic*8+7] = n0          + jplus1 + kplus1;

            ic++;
          }
        }
      }

      inode = in;
      ihex  = ic;
    }
  }

  if(inode != v->nnode)
  {
    printf("ERROR nnode: %d %d\n", inode, v->nnode);
  }
  if(ihex != v->nhex)
  {
    printf("ERROR nhex: %d %d\n", ihex, v->nhex);
  }

  if(bp)
  {
      for(p=0; p<pmax; p++)
      {
          if(bp[p].frac) my_free(bp[p].frac);
          if(bp[p].dfracdr) my_free(bp[p].dfracdr);
          if(bp[p].dfracds) my_free(bp[p].dfracds);
          if(bp[p].dfracdt) my_free(bp[p].dfracdt);
          if(bp[p].r) my_free(bp[p].r);
          if(bp[p].phi) my_free(bp[p].phi);
          if(bp[p].dphi) my_free(bp[p].dphi);
      }
      my_free(bp);
  }

}

void compute_vm(int nnode, double *soln, double *vm)
{
    int i;
    double rho,rhou,rhov,rhow;
    for(i=0; i<nnode; i++)
    {
        rho  = soln[i*5];
        rhou = soln[i*5+1];
        rhov = soln[i*5+2];
        rhow = soln[i*5+3];
        vm[i] = sqrt(rhou*rhou+rhov*rhov+rhow*rhow)/rho;
    }
}

void visit_init()
{

  if(d->visit->visit_flag == 0) return;

  int i;

  visit_t *v = d->visit;
  flow_solver_t *f = d->flow;

  v->nvar = 8;
  v->varName = (char**)malloc(v->nvar * sizeof(char*));
  for(i=0; i<v->nvar; i++)
    v->varName[i] = (char*)malloc(30 * sizeof(char));


  strcpy(v->varName[0], "rho");
  strcpy(v->varName[1], "rhou");
  strcpy(v->varName[2], "rhov");
  strcpy(v->varName[3], "rhow");
  strcpy(v->varName[4], "rhoe");
  strcpy(v->varName[5], "vm");
  strcpy(v->varName[6], "omega");
  strcpy(v->varName[7], "qcriterion");

//printf("point 1\n");
  if(d->high_order_group_flag == 0)
  {
// for nsu3d
    int i;

    v->nnode  = f->nnode;
    v->ntetra = f->ntet;
    v->npyr   = f->npyramid;
    v->nprism = f->nprism;
    v->nhex   = f->nhex;
    v->tetra  = f->ndc4;
    v->pyr    = f->ndc5;
    v->prism  = f->ndc6;
    v->hex    = f->ndc8;
    v->coord  = f->xgeom_translated;
    v->iblank = f->iblank;

    v->vm = my_alloc(double,f->nnode);
    compute_vm(f->nnode,f->soln,v->vm);

    v->omega = my_alloc(double,f->nnode);
    v->qcriterion = my_alloc(double,f->nnode);
    for(i=0; i<v->nnode; i++){
      v->omega[i] = 0.0;
      v->qcriterion[i]=0.0;
    }

    v->narray = v->nvar-5+1;
    v->mvar = (int*)malloc(v->narray * sizeof(int));
    v->mvar[0] = 5;
    v->mvar[1] = 1;
    v->mvar[2] = 1;
    v->mvar[3] = 1;

    v->var = (double**)malloc(v->narray * sizeof(double*));
    v->var[0] = f->soln;
    v->var[1] = v->vm;
    v->var[2] = v->omega;
    v->var[3] = v->qcriterion;

    v->node_flag = (double*)malloc(v->nnode *sizeof(double));
    for(i=0; i<v->nnode; i++)
        v->node_flag[i] = 1.e10;

    for(i=0; i<f->nwbc; i++)
        v->node_flag[f->iwbcnode[i]-1] = 1.;

  }
  else
  {
    v->narray = v->nvar;
    v->mvar = (int*)malloc(v->narray * sizeof(int));
    for(i=0; i<v->narray; i++)
      v->mvar[i] = 1;
    refine_mesh(f, v);
    v->ntetra = 0;
    v->npyr   = 0;
    v->nprism = 0;
    v->tetra = my_alloc(int, v->ntetra*4);
    v->pyr   = my_alloc(int, v->npyr*5);
    v->prism = my_alloc(int, v->nprism*6);

    v->node_flag = (double*)malloc(v->nnode *sizeof(double));
    for(i=0; i<v->nnode; i++)
        v->node_flag[i] = 1.e10;

  }

  v->visit_extract_ini(v->nnode, v->ntetra, v->npyr, v->nprism, v->nhex,
                                 v->tetra,  v->pyr,  v->prism,  v->hex,
                       v->coord, v->iblank, v->nvar, v->varName, v->narray, v->mvar, v->var, v->node_flag,
                       d->group_comm);
}

void visit_extract(int icyc)
{
  if(d->visit->visit_flag == 0) return;

  int i;
  visit_t *v = d->visit;
  flow_solver_t *f = d->flow;

  if(d->high_order_group_flag == 0)
  {
// for nsu3d
    if(v->vm) my_free(v->vm);
    v->vm = my_alloc(double,f->nnode);
    if(v->omega) my_free(v->omega);
    v->omega = my_alloc(double,f->nnode);
    if(v->qcriterion) my_free(v->qcriterion);
    v->qcriterion = my_alloc(double,f->nnode);
    for(i=0; i<v->nnode; i++){
      v->omega[i] = 0.0;
      v->qcriterion[i] = 0.0;
    }

    v->nnode  = f->nnode;
    v->ntetra = f->ntet;
    v->npyr   = f->npyramid;
    v->nprism = f->nprism;
    v->nhex   = f->nhex;
    v->tetra  = f->ndc4;
    v->pyr    = f->ndc5;
    v->prism  = f->ndc6;
    v->hex    = f->ndc8;
    v->coord  = f->xgeom_translated;
    v->iblank = f->iblank;

    compute_vm(f->nnode,f->soln,v->vm);

    v->var[0] = f->soln;
    v->var[1] = v->vm;
    v->var[2] = v->omega;
    v->var[3] = v->qcriterion;

  }
  else
  {
    refine_mesh(f, v);
    if(v->node_flag) free(v->node_flag);
    v->node_flag = (double*)malloc(v->nnode *sizeof(double));
    int i;
    for(i=0; i<v->nnode; i++)
        v->node_flag[i] = 1.e10;
  }
  double t = icyc;
  v->visit_extract_this(icyc,t,
                        v->nnode, v->ntetra, v->npyr, v->nprism, v->nhex,
                                  v->tetra,  v->pyr,  v->prism,  v->hex,
                        v->coord, v->iblank, v->nvar, v->varName, v->narray, v->mvar, v->var, v->node_flag);
}


