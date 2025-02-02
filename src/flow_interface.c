//
//  flow_interface.c
//  cdriver
//
//  Created by Michael Brazell on 1/6/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#include "driver.h"
#include <dlfcn.h>
#include <stdio.h>
#include "alloc.h"

extern driver_t *d;

void flow_load_dynamic_library() {

  char *error;
  flow_solver_t *f   = d->flow;
  void *group_handle = f->group_handle;
  char *filename     = f->solver_so_file;

  /* get group handle */
  group_handle = dlopen (filename, RTLD_LAZY);
  if (!group_handle) {fputs (dlerror(), stderr);}

  /* load each function */

  /* these are required */
  *(void **)(&f->initialize_group_mpi) = dlsym(group_handle,"driver_interface_initialize_group_mpi");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}


  *(void **)(&f->initialize) = dlsym(group_handle, "driver_interface_initialize");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}


  *(void **)(&f->get_data) = dlsym(group_handle, "driver_interface_get_data");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}


  *(void **)(&f->set_pointers) = dlsym(group_handle, "driver_interface_set_pointers");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}


  *(void **)(&f->steady_update) = dlsym(group_handle, "driver_interface_steady_update");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}


  *(void **)(&f->unsteady_update) = dlsym(group_handle, "driver_interface_unsteady_update");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}


  *(void **)(&f->output_solution) = dlsym(group_handle, "driver_interface_output_solution");
  if ((error = dlerror()) != NULL) {fputs(error, stderr);}


  /* these routines are not required and are initialized with null
   they are pointed to functions in .so if they exist
   to use these check for null first */


  *(void **)(&f->regrid) = NULL;
  *(void **)(&f->regrid) = dlsym(group_handle, "driver_interface_regrid");

  *(void **)(&f->igbp_regrid) = NULL;
  *(void **)(&f->igbp_regrid) = dlsym(group_handle, "driver_interface_igbp_regrid");

  *(void **)(&f->real2ghost) = NULL;
  *(void **)(&f->real2ghost) = dlsym(group_handle, "driver_interface_real2ghost");

  *(void **)(&f->get_igbp_data) = NULL;
  *(void **)(&f->get_igbp_data) = dlsym(group_handle, "driver_interface_get_igbp_data");

  *(void **)(&f->iblank_check) = NULL;
  *(void **)(&f->iblank_check) = dlsym(group_handle, "driver_interface_iblank_check");

  *(void **)(&f->point_inclusion) = NULL;
  *(void **)(&f->point_inclusion) = dlsym(group_handle, "driver_interface_point_inclusion");

  *(void **)(&f->bounding_box_intersection) = NULL;
  *(void **)(&f->bounding_box_intersection) = dlsym(group_handle, "driver_interface_bounding_box_intersection");

  *(void **)(&f->iblank_flipper) = NULL;
  *(void **)(&f->iblank_flipper) = dlsym(group_handle, "driver_interface_iblank_flipper");

  *(void **)(&f->basis_and_derivative) = NULL;
  *(void **)(&f->basis_and_derivative) = dlsym(group_handle, "driver_interface_basis_and_derivative");

  *(void **)(&f->extract_points) = NULL;
  *(void **)(&f->extract_points) = dlsym(group_handle, "driver_interface_extract_points");
}

void flow_close_dynamic_library(){
  if(d->flow->group_handle) dlclose(d->flow->group_handle);
}

void flow_iblank_flipper(){
  if((d->flow->iblank_flipper)!=NULL) {
    d->flow->iblank_flipper();
  }
}

void flow_set_p4est_flag(){
  d->p4est_group_flag = 0;
  if((d->flow->bounding_box_intersection)!=NULL){
    d->p4est_group_flag = 1;
  }
}

void flow_initialize_group_mpi(MPI_Comm group_comm){
  MPI_Fint new_group_comm = MPI_Comm_c2f(group_comm);
  (*d->flow->initialize_group_mpi)(&new_group_comm);
}

void flow_set_pointers(){

  flow_solver_t *f = d->flow;

  f->set_pointers(&f->soln,
                  &f->xgeom,
                  &f->iblank,
                  &f->iwbcnode,
                  &f->iobcnode,
                  &f->ndc4,
                  &f->ndc5,
                  &f->ndc6,
                  &f->ndc8,
                  &f->iblank_cell,
                  &f->count_receptor_nodes,
                  &f->create_receptor_nodes,
                  &f->donor_inclusion_test,
                  &f->create_donor_frac,
                  &f->convert_to_receptor_coefficients);

  f->get_data(&f->body_tag,
              &f->nnode,
              &f->nwbc,
              &f->nobc,
              &f->ntet,
              &f->npyramid,
              &f->nprism,
              &f->nhex);

  /* translate the geometry using translation from input file */
  if(f->xgeom_translated) my_free(f->xgeom_translated);
  f->xgeom_translated = my_alloc(double,3*f->nnode);

  int i,j;
  for(i=0;i<f->nnode;++i){
    for(j=0;j<3;j++){
      f->xgeom_translated[3*i+j] = f->xgeom[3*i+j] + f->translation[j];
    }
  }

}

void flow_output_solution(int t){

  double t1 = MPI_Wtime();
  d->flow->output_solution(&t);
  double t2 = MPI_Wtime();

  if(d->group_master_flag) printf("[wake3d] Group: %d output solution time (sec): %e\n",d->group,t2-t1);

};


/* these routines are not required and are initialized with null
 they are pointed to functions in .so if they exist
 to use these check for null first */

void flow_real2ghost() {
  if((d->flow->real2ghost)!=NULL) {
    d->flow->real2ghost();
  }
}

void flow_regrid() {
  if((d->flow->regrid)!=NULL) {
    d->flow->regrid();
  }
}

void flow_get_igbp_data(){
  flow_solver_t *f = d->flow;
  if((f->get_igbp_data)!=NULL) {
    f->get_igbp_data(&f->igbp->n_lcl, &f->igbp->index, &f->igbp->dx_lcl);
  }
}

void flow_obc_regrid() {
  flow_solver_t *f = d->flow;
  if((f->igbp_regrid)!=NULL){
    (*f->igbp_regrid)(f->obc->nobc,f->obc->x,f->obc->y,f->obc->z,f->obc->dx);
  }

}

void flow_igbp_regrid() {
  flow_solver_t *f = d->flow;
  if((f->igbp_regrid)!=NULL){
    (*f->igbp_regrid)(f->igbp->n,f->igbp->x,f->igbp->y,f->igbp->z,f->igbp->dx);
  }
}

void flow_iblank_check() {
  flow_solver_t *f = d->flow;
  if((f->iblank_check)!=NULL){
    (*f->iblank_check)();
  }
}


void point_inclusion(int *npoint, double *x, int *cell_id) {
  if((d->flow->point_inclusion)!=NULL) {
    d->flow->point_inclusion(npoint,x,cell_id);
  } else {
    int i;
    for(i=0;i<*npoint;++i){
      cell_id[i] = -1;
    }
  }
}

void bounding_box_intersection(int *nbb, double *bb_xlo, double *bb_xhi, int *bb_flag) {
  if((d->flow->bounding_box_intersection)!=NULL){
    d->flow->bounding_box_intersection(nbb,bb_xlo,bb_xhi,bb_flag);
  } else {
    int i;
    for(i=0;i<*nbb;i++){
      bb_flag[i] = 0;
    }
  }
}



