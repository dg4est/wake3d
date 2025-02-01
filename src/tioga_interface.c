//
//  tioga_interface.c
//  cdriver
//
//  Created by Michael Brazell on 1/6/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//


#include <dlfcn.h>
#include <stdio.h>
#include "driver.h"
#include "alloc.h"

/* this is to match the node connectivities when they are base 1 */
#define BASE 1

extern driver_t *d;

void tioga_load_dynamic_library() {

  char *error;
  tioga_t *t = d->tioga;
  char *filename = t->tioga_so_file;

  /* get group handle */
  t->group_handle = dlopen (filename, RTLD_LAZY);
  if (!t->group_handle) {fputs (dlerror(), stderr);}

  /* load each function */
  *(void **)(&t->dataupdate) = dlsym(t->group_handle, "tioga_dataupdate_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->tioga_delete) = dlsym(t->group_handle, "tioga_delete_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->getdonorcount) = dlsym(t->group_handle, "tioga_getdonorcount_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->getdonorinfo) = dlsym(t->group_handle, "tioga_getdonorinfo_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->init) = dlsym(t->group_handle, "tioga_init_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->performconnectivity) = dlsym(t->group_handle, "tioga_performconnectivity_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->performconnectivity_highorder) = dlsym(t->group_handle, "tioga_performconnectivity_highorder_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->preprocess_grids) = dlsym(t->group_handle, "tioga_preprocess_grids_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->registergrid_data) = dlsym(t->group_handle, "tioga_registergrid_data_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->register_composite_body) = dlsym(t->group_handle, "tioga_register_composite_body_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->set_highorder_callback) = dlsym(t->group_handle,"tioga_set_highorder_callback_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->setcelliblank) = dlsym(t->group_handle, "tioga_setcelliblank_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->setresolutions) = dlsym(t->group_handle, "tioga_setresolutions_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->setsymmetry) = dlsym(t->group_handle, "tioga_setsymmetry_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->writeoutputfiles) = dlsym(t->group_handle, "tioga_writeoutputfiles_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->set_p4est) = dlsym(t->group_handle, "tioga_set_p4est_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->set_p4est_search_callback) = dlsym(t->group_handle, "tioga_set_p4est_search_callback_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&t->setcomposite) = dlsym(t->group_handle, "tioga_setcomposite_");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}
}

void tioga_set_p4est(){
  // only for off-body setting (sets resolution 1000x)
  if(d->off_body_group_flag) d->tioga->set_p4est();
}

void tioga_close_dynamic_library(){
  if(d->tioga->group_handle) dlclose(d->tioga->group_handle);
}

void tioga_init(MPI_Comm mpi_comm_in){
  d->tioga->init(mpi_comm_in);
}

void amr_bounding_box(double *amrxlo, double *amrxhi, double *amrdx){
  double xlo[3];
  double xhi[3];
  double dx;
  int i,j,k;

  flow_solver_t *f = d->flow;

  xlo[0]=xlo[1]=xlo[2]=1e15;
  xhi[0]=xhi[1]=xhi[2]=-1e15;
  dx = 0.0;

  if(d->off_body_group_flag){
    for(i=0;i<f->nhex;i++){
      for(j=0;j<8;j++){
        int n = f->ndc8[8*i+j]-BASE;
        double x[3];
        for(k=0;k<3;k++){
          x[k] = f->xgeom_translated[3*n+k];
          if(x[k] < xlo[k]) xlo[k] = x[k];
          if(x[k] > xhi[k]) xhi[k] = x[k];
        }
      }
      int n1 = f->ndc8[8*i+0]-BASE;
      int n2 = f->ndc8[8*i+1]-BASE;
      double x1 = f->xgeom[3*n1+0];
      double x2 = f->xgeom[3*n2+0];

      if(x2-x1 > dx) dx = x2-x1;

    }
  }

  /* all reduce xlo/hi so that all procs know the amr bounding box */
  MPI_Allreduce(xlo,amrxlo,3,MPI_DOUBLE,MPI_MIN,d->mpicomm);
  MPI_Allreduce(xhi,amrxhi,3,MPI_DOUBLE,MPI_MAX,d->mpicomm);
  MPI_Allreduce(&dx,amrdx,1,MPI_DOUBLE,MPI_MAX,d->mpicomm);
}

void atm_amr_intersection_compress(){

  if(d->atm_compressed) {
    if(d->atm_group_flag){
      d->flow->nhex = d->atm_nhex;
    }
    return;
  }

  int i,j,icell;
  flow_solver_t *f = d->flow;

  double xlo[3];
  double xhi[3];
  double dx,overlap;
  amr_bounding_box(xlo, xhi, &dx);

  if(d->rank==0) printf("amr dx %f \n",dx);

  overlap = dx*1.1;

  if(d->atm_group_flag){

    if(f->ntet>0) printf("error assumed atm code does not have tets\n");
    if(f->nprism>0) printf("error assumed atm code does not have prism\n");
    if(f->npyramid>0) printf("error assumed atm code does not have pyramids\n");

    d->atm_nhex = 0;

    for(i=0;i<f->nhex;i++){

      icell = 0;

      for(j=0;j<8;j++){
        int n = f->ndc8[8*i+j]-BASE;
        double x[3];
        x[0] = f->xgeom_translated[3*n+0];
        x[1] = f->xgeom_translated[3*n+1];
        x[2] = f->xgeom_translated[3*n+2];

        double eps = 1.0e-6;

        if(x[0] > xlo[0] - overlap && x[0] < xlo[0] + overlap) icell = 1;
        if(x[0] > xhi[0] - overlap && x[0] < xhi[0] + overlap) icell = 1;

        // if(x[1] > xlo[1] - eps && x[1] < xlo[1] + overlap) icell = 1;
        if(x[1] > xhi[1] - overlap && x[1] < xhi[1] + overlap) icell = 1;

        if(x[2] > xlo[2] - overlap && x[2] < xlo[2] + overlap) icell = 1;
        if(x[2] > xhi[2] - overlap && x[2] < xhi[2] + overlap) icell = 1;

      }

      /* compress hex list */
      if(icell){
        for(j=0;j<8;j++){
          f->ndc8[8*d->atm_nhex+j] = f->ndc8[8*i+j];
        }
        ++d->atm_nhex;
      }
    }

    if(d->atm_nhex < d->flow->nhex) printf("proc %d nhex before %d and after %d\n", d->rank,d->flow->nhex,d->atm_nhex);

    d->flow->nhex = d->atm_nhex;

  }

  d->atm_compressed = 1;
}


void set_resolutions(int user_spec_res){

  int i;
  flow_solver_t *f = d->flow;

  if(f->node_res) my_free(f->node_res);
  if(f->cell_res) my_free(f->cell_res);
  f->node_res = NULL;
  f->cell_res = NULL;

  if(!d->atm_group_flag &&  user_spec_res ||
      d->atm_group_flag && !user_spec_res ||
      d->flow->receptor_only_flag){

    double BIGVALUE;
    if(user_spec_res || d->flow->receptor_only_flag){
      BIGVALUE = 1.0e15; // must match Tioga's internal BIGVALUE = 1.0e15
    } else {
      BIGVALUE = 1.0e9;
    }

    f->node_res = my_alloc(double,f->nnode);
    int ncell = f->ntet+f->nprism+f->npyramid+f->nhex;
    f->cell_res = my_alloc(double,ncell);

    for(i=0;i<f->nnode;i++){
      //f->node_res[i] = (double) user_spec_res;
      f->node_res[i] = BIGVALUE;
    }

    for(i=0;i<ncell;i++){
      //f->cell_res[i] = (double) user_spec_res;
      f->cell_res[i] = BIGVALUE;
    }
  }

}


void tioga_register_and_connect(int user_spec_res){
  flow_solver_t *f = d->flow;
  int temp = f->nwbc;
  int i;

  if(user_spec_res) f->nwbc = 0;

  if(!user_spec_res && d->atm_mode_flag){
    atm_amr_intersection_compress();
  }

  /* register unstructured or structured data with tioga */
  tioga_registergrid_data();

  f->nwbc = temp;

  /* if a high order group then give tioga iblank_cell pointer and call back functions */
  tioga_setcelliblank();
  tioga_set_highorder_callback();

  if(!user_spec_res){
    tioga_set_p4est();
    setup_p4est_maps(f->xgeom_translated,&f->nnode);
    d->tioga->set_p4est_search_callback(search_p4est,check_intersect_p4est);
  }

  /* change resolutions for initialization or for atm code */
  set_resolutions(user_spec_res);
  d->tioga->setresolutions(f->node_res,f->cell_res);

  /* preprocess and connect */
  d->tioga->preprocess_grids();

  d->tioga->performconnectivity();

  flow_iblank_flipper();

  /* if any groups high order then everyone calls */
  tioga_performconnectivity_highorder();

  /* turn off all nodes in atm to speed up dataupdate */
  if(d->atm_group_flag){
    for(i=0;i<f->nnode;i++){
      //if(f->iblank[i] == -1) f->iblank[i] = 0;
    }
  }
}


void tioga_registergrid_data(){
  flow_solver_t *f = d->flow;

  /* every proc has to register an unstructured grid even if there is nothing in it */
  d->tioga->registergrid_data(&d->global_body_tag,
                              &f->nnode,
                               f->xgeom_translated,
                               f->iblank,
                              &f->nwbc,
                              &f->nobc,
                               f->iwbcnode,
                               f->iobcnode,
                              &f->ncell_types,
                              &f->kstride4,
                              &f->ntet,
                               f->ndc4,
                              &f->kstride5,
                              &f->npyramid,
                               f->ndc5,
                              &f->kstride6,
                              &f->nprism,
                               f->ndc6,
                              &f->kstride8,
                              &f->nhex,
                               f->ndc8);
}


void tioga_dataupdate(){
  d->tioga->dataupdate(d->flow->soln,&d->nfield,"row");
}

void tioga_set_highorder_callback(){
  if(d->high_order_group_flag){
    d->tioga->set_highorder_callback(d->flow->count_receptor_nodes,
                                     d->flow->create_receptor_nodes,
                                     d->flow->donor_inclusion_test,
                                     d->flow->create_donor_frac,
                                     d->flow->convert_to_receptor_coefficients);
  }
}

void tioga_setcelliblank(){
  if(d->high_order_group_flag) {
    d->tioga->setcelliblank(d->flow->iblank_cell);
  }
}

void tioga_performconnectivity_highorder(){
  if(d->high_order_mode_flag){
    d->tioga->performconnectivity_highorder();
  }
}

