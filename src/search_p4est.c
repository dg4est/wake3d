//
//  p4est_search.c
//  cdriver
//
//  Created by Michael Brazell on 1/11/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#include "driver.h"
#include "alloc.h"

extern driver_t *d;

void setup_p4est_maps(double *xsearch,int *numpts){
  double *bb_xlo,*bb_xhi;
  int *counts,*displacements;
  int *bb_flags;
  int i,j;

  /* use for counting number of bounding boxes and pts later */
  counts = my_alloc(int,d->num_procs);
  displacements = my_alloc(int,d->num_procs);

  int npts = *numpts;
  int near_body = (d->p4est_group_flag) ? 0:1;

  /* check if you are near body and contain points */
  int bb_flag = (near_body && npts > 0) ? 1:0;

  /* gather all procs that have bounding boxes */
  MPI_Allgather(&bb_flag,1,MPI_INT,counts,1,MPI_INT,d->mpicomm);

  /* count number of bounding boxes */
  int num_bb = 0;
  for(i = 0; i < d->num_procs; ++i) num_bb += counts[i];

  /* if no bounding boxes then return */
  if(num_bb == 0) return;

  /* stretch counts by 3 for sending bounding box coordinates */
  for(i = 0;i < d->num_procs; ++i) counts[i] *= 3;

  /* find the displacements when stacked for allgatherv later */
  displacements[0] = 0;
  for(i = 1; i < d->num_procs; i++){
    displacements[i] = counts[i-1] + displacements[i-1];
  }

  /* find bounding box coordinates */
  bb_xlo = my_alloc(double,3*num_bb);
  bb_xhi = my_alloc(double,3*num_bb);
  double bb_xlo_loc[3];
  double bb_xhi_loc[3];

  if(near_body && npts > 0){
    double xlo[3];
    double xhi[3];

    for(j = 0; j < 3; ++j) xlo[j] = xsearch[j];
    for(j = 0; j < 3; ++j) xhi[j] = xsearch[j];

    for(i = 0; i < npts; ++i){
      for(j = 0; j < 3; ++j){
        if(xsearch[3*i+j] < xlo[j]) xlo[j] = xsearch[3*i+j];
        if(xsearch[3*i+j] > xhi[j]) xhi[j] = xsearch[3*i+j];
      }
    }

    for(j = 0; j < 3; ++j) bb_xlo_loc[j] = xlo[j];
    for(j = 0; j < 3; ++j) bb_xhi_loc[j] = xhi[j];
  }

  /* all processors will have all of the bounding boxes */
  /* would it be more efficient for all near-body procs to gather and then scatter to off-body procs? */
  MPI_Allgatherv(bb_xlo_loc, 3*bb_flag, MPI_DOUBLE, bb_xlo, counts, displacements, MPI_DOUBLE, d->mpicomm);
  MPI_Allgatherv(bb_xhi_loc, 3*bb_flag, MPI_DOUBLE, bb_xhi, counts, displacements, MPI_DOUBLE, d->mpicomm);

  /* identify which bounding boxes intersect with p4est, bb_flags = 0 if no intersection and 1 if there is an intersection */
  bb_flags = my_alloc(int,num_bb);
  bounding_box_intersection(&num_bb, bb_xlo, bb_xhi, bb_flags);

  /* convert bb flags to a recvmap and put back in counts for now */
  int ind = 0;
  for(i = 0; i < d->num_procs; ++i){
    if(counts[i]) counts[i] = bb_flags[ind++];
  }

  /* find the processor maps */
  if(d->recvmap == NULL) d->recvmap = my_alloc(int,d->num_procs);
  if(d->sendmap == NULL) d->sendmap = my_alloc(int,d->num_procs);

  for(i = 0; i < d->num_procs; ++i) d->recvmap[i] = counts[i];
  for(i = 0; i < d->num_procs; ++i) d->sendmap[i] = 0;

  /* for future scalability this could probably be swapped out for a mpi-3 window get/put one sided comm */
  MPI_Alltoall(d->recvmap,1,MPI_INT,d->sendmap,1,MPI_INT,d->mpicomm);

  my_free(bb_xlo);
  my_free(bb_xhi);
  my_free(bb_flags);
  my_free(counts);
  my_free(displacements);
}

void check_intersect_p4est(int *proc, int* overlap){
  if(d->recvmap){
    *overlap = d->recvmap[*proc]+d->sendmap[*proc];
  } else {
    *overlap = 1;
  }
}

void search_p4est(double *xsearch,int *process_id,int *cell_id,int *numpts){
    MPI_Request *request;
    MPI_Status *status;
    int *num_nb_pts,*displacements;
    int *sendmap,*recvmap;
    int *gather_cell_id;
    int *nb_cell_id;
    double *nb_pts;

    int nsend,nrecv;
    int tag,irnum;
    int i,j;

    int near_body = (d->p4est_group_flag) ? 0:1;

    /* initialize the cell id's and process id's */
    int npts = *numpts;
    for(j = 0; j < npts; j++) cell_id[j] = process_id[j] = -1;

    /* setup the communication maps by sending bounding boxes around */
    setup_p4est_maps(xsearch,numpts);

    /* slide the mapping, insert rank, and count number of send/recv */
    nsend = nrecv = 0;
    for(i = 0; i < d->num_procs; i++){
        if(d->sendmap[i]) nsend++;
        if(d->recvmap[i]) nrecv++;
    }
    sendmap = my_alloc(int,nsend);
    recvmap = my_alloc(int,nrecv);

    nsend = nrecv = 0;
    for(i = 0; i < d->num_procs; i++){
        if(d->sendmap[i]) sendmap[nsend++] = i;
        if(d->recvmap[i]) recvmap[nrecv++] = i;
    }

    /* send the number of near-body points to all off-body cores and store them in counts */
    request = my_alloc(MPI_Request,(nsend+nrecv));
    status = my_alloc(MPI_Status,(nsend+nrecv));
    displacements = my_alloc(int,nrecv);
    num_nb_pts = my_alloc(int,nrecv);

    /* off-body receives number of near-body points */
    /* near-body sends pts to off-body */
    tag = 1; irnum = 0;
    for(i=0;i<nrecv;i++) MPI_Irecv(&num_nb_pts[i],1,MPI_INT,recvmap[i],tag,d->mpicomm,&request[irnum++]);
    for(i=0;i<nsend;i++) MPI_Isend(&npts,         1,MPI_INT,sendmap[i],tag,d->mpicomm,&request[irnum++]);
    MPI_Waitall(irnum,request,status);

#if 0
    if(near_body) {
      for(i = 0; i < nsend; i++) printf("near body rank: %d npts: %d send: %d\n",d->rank,npts,sendmap[i]);
    }
    if(!near_body) {
      for(i = 0; i < nrecv; i++) printf("off body rank: %d npts %d recv: %d\n",d->rank,num_nb_npts[i],recvmap[i]);
    }
#endif

    /* count total number of near body points collected on off body */
    int total_nb_pts = 0;
    for(i = 0; i < nrecv; i++) total_nb_pts += num_nb_pts[i];

    /* save displacements too for recv/send */
    displacements[0] = 0;
    for(i = 1; i < nrecv; i++) displacements[i] = num_nb_pts[i-1] + displacements[i-1];

    /* store all near-body points on the off-body */
    int  num_total_nb_pts = (total_nb_pts > 0) ? (3*total_nb_pts):1;
    nb_pts = my_alloc(double,num_total_nb_pts);

    tag = 1; irnum = 0;
    for(i=0;i<nrecv;i++) MPI_Irecv(&nb_pts[3*displacements[i]],3*num_nb_pts[i],MPI_DOUBLE,recvmap[i],tag,d->mpicomm,&request[irnum++]);
    for(i=0;i<nsend;i++) MPI_Isend(xsearch,                    3*npts,         MPI_DOUBLE,sendmap[i],tag,d->mpicomm,&request[irnum++]);
    MPI_Waitall(irnum,request,status);

    /* search if any of the pts are in the octree on this processor */
    int num_nb_cell = (total_nb_pts) ? (total_nb_pts):1;
    nb_cell_id = my_alloc(int,num_nb_cell);

    /* call a p4est octree search on a distributed octree that searches for points and if found returns the cell id */
    point_inclusion(&total_nb_pts,nb_pts,nb_cell_id);

    /* send all found pts cell ids back to near-bodies */
    /* fixme no need to send pts not found but it is easier than finding all the sizes again */
    int num_gather_cell = (nsend) ? (nsend*npts):1;
    gather_cell_id = my_alloc(int,num_gather_cell);

    /* reversing the sends and receives here */
    /* receive multiple sets of possible cell ids */
    /* send the cell ids back to the near-bodies and stack them */
    tag = 1; irnum = 0;
    for(i=0;i<nsend;i++) MPI_Irecv(&gather_cell_id[i*npts],      npts,         MPI_INT,sendmap[i],tag,d->mpicomm,&request[irnum++]);
    for(i=0;i<nrecv;i++) MPI_Isend(&nb_cell_id[displacements[i]],num_nb_pts[i],MPI_INT,recvmap[i],tag,d->mpicomm,&request[irnum++]);
    MPI_Waitall(irnum,request,status);

    /* search over multiple sets of cell id's for a match and store process id */
    /* switch the i and j loop? remember to remove break if you do this */
    if(near_body){
        for(j = 0; j < npts; j++){
            for(i = 0; i < nsend; i++){
                if(gather_cell_id[npts*i+j] > -1){
                    cell_id[j] = gather_cell_id[npts*i+j];
                    process_id[j] = sendmap[i];
                    break;
                }
            }
        }
    }

    if(near_body){
        for(j = 0; j < npts; j++){
            if(cell_id[j] == -1) printf("[driver] oh no pt not found in p4est search rank: %d pt: %d\n",d->rank,j);
        }
    }

    /* deallocate memory */
    my_free(displacements);
    my_free(num_nb_pts);
    my_free(request);
    my_free(status);
    my_free(sendmap);
    my_free(recvmap);
    my_free(nb_pts);
    my_free(nb_cell_id);
    my_free(gather_cell_id);
}