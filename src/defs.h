//
//  defs.h
//  cdriver
//
//  Created by Michael Brazell on 1/6/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#ifndef defs_h
#define defs_h

#define buff_size 1024
#include <mpi.h>
#include <sys/stat.h>

/* outer boundary condition points for tagging on the off-body */
typedef struct obc {
  int nobc;   /**< number of obc points */
  double *x;  /**< x location of points */
  double *y;  /**< y location of points */
  double *z;  /**< z location of points */
  double *dx; /**< mesh size at points */
} obc_t;

/** intergrid boundary points for tagging on the off-body */
typedef struct igbp {

  /** variables that come from get igbp */
  int n_lcl;       /**< number of igbp on proc */
  int *index;      /**< index of igbp node on proc */
  double *dx_lcl;  /**< mesh size pointer on proc */

  /** variables that get accumulated and copied to all procs */
  int n;           /**< number of igbp on all procs */
  double *x;       /**< x location of igbp copied on all procs */
  double *y;       /**< y location of igbp copied on all procs */
  double *z;       /**< z location of igbp copied on all procs */
  double *dx;      /**< mesh size of igbp copied on all procs */

} igbp_t;

/* a struct to store flow solver function pointers and mesh */
typedef struct flow_solver {

  char solver_so_file[buff_size];
  char solver_so_name[buff_size];
  void *group_handle;
  int receptor_only_flag;

  /* mandatory callback functions */
  void (*initialize_group_mpi)(int *group_comm);
  void (*initialize)(void);
  void (*unsteady_update)(int *ncyc);
  void (*steady_update)(int *ncyc);
  void (*set_pointers)(double**,double**,int**,int**,int**,int**,int**,int**,int**,int**,void**,void**,void**,
                        void(**create_donor_frac)(int*, double*, int*, int*, double*, double*, int*),void**);
  void (*get_data)(int*,int*,int*,int*,int*,int*,int*,int*);
  void (*output_solution)(int *t);

  /* these routines are not required and are initialized with null
   they are pointed to functions in .so if they exist
   to use these check for null first */
  void (*regrid)(void);
  void (*get_igbp_data)(int *, int **, double **);
  void (*igbp_regrid)(int,double*,double*,double*,double*);
  void (*real2ghost)(void);
  void (*iblank_check)(void);
  void (*point_inclusion)(int *npoint, double *x, int *cell_id);
  void (*bounding_box_intersection)(int *nbb, double *bb_xlo, double *bb_xhi, int *bb_flag);
  void (*iblank_flipper)(void);
  void (*basis_and_derivative)(int *pdegree,int *numpnts,double *xi,double *phi, double *dphidxi);
  void (*extract_points)(int *npts,double *pts,int *nvar,char **var_names,double **extract_data);

  /* pointers to data inside flow solvers */
  double* soln;
  double* xgeom;
  int* iblank;
  int* iwbcnode;
  int* iobcnode;
  int* ndc4;
  int* ndc5;
  int* ndc6;
  int* ndc8;
  int* iblank_cell;
  double translation[3];
  double *xgeom_translated;

  /* node and cell resolution to force tioga to connect */
  double* node_res;
  double* cell_res;

  /* function pointers */
  void (*count_receptor_nodes);
  void (*create_receptor_nodes);
  void (*donor_inclusion_test);
  void (*create_donor_frac)(int*, double*, int*, int*, double*, double*, int*);
  void (*convert_to_receptor_coefficients);

  /* static data that may need to get updated */
  int nnode;
  int nwbc;
  int nobc;
  int ntet;
  int nprism;
  int npyramid;
  int nhex;
  int body_tag;

  int ncell_types;
  int kstride4;
  int kstride5;
  int kstride6;
  int kstride8;

  igbp_t *igbp;
  obc_t *obc;

} flow_solver_t;

/* tioga struct to store function pointers and files */
typedef struct tioga {

  /** names of shared object files */
  char tioga_so_file[buff_size]; /**< filename for the shared object file for tioga */
  void *group_handle; /**< handle name for the shared object save for closing dynamic library later */

  /** tioga mandatory function pointers */
  void (*init)(MPI_Comm); /**< initializes tioga and sets mpi communicator */
  void (*registergrid_data)(int*,int*,double*,int*,int*,int*,int*,int*,int*,...); /**< sends the body tag, nodes, and connectivities to tioga */
  void (*register_composite_body)(int *compid,int *bodyids,int *nbodies,double *searchTol);
  void (*preprocess_grids)(void); /**< performs precossing steps in tioga */
  void (*performconnectivity)(void); /**< performs connecitivity with all grids */
  void (*dataupdate)(double *q,int *nvar,char *itype); /**< passes the solution between grids */

  /** tioga mandatory function pointers for high-order */
  void (*setcelliblank)(int*); /**< gives the location of a the cell iblank array */
  void (*set_highorder_callback)(void*,void*,void*,void*,void*); /**< gives tioga access to high-order callback functions */
  void (*performconnectivity_highorder)(void); /**< performs a high order connectivity */

  /** function pointers for p4est search callback function for tioga */
  void (*set_p4est_search_callback)(void*,void*); /**< gives two new callback functions for tioga to use in p4est */
  void (*set_p4est)(void); /**< sets p4est groups for tioga */

  /** other tioga function pointers that could be used */
  void (*tioga_delete)(); /**< deletes a tioga object */
  void (*getdonorcount)();
  void (*getdonorinfo)();
  void (*setresolutions)(double *nres, double *cres); /**< allows a user specified resolution to be used */
  void (*setsymmetry)();
  void (*setcomposite)(int *ncomp); /**< set the number of composite bodies */
  void (*writeoutputfiles)();

} tioga_t;

typedef struct visit {

  int visit_flag;

  /** names of shared object files */
  char visit_so_file[buff_size];
  void *group_handle;

  int nnode;
  int ntetra,npyr,nprism,nhex;
  int *tetra,*pyr,*prism,*hex;
  double *coord;
//  double *soln;
  double *vm;
  double *omega;
  double *qcriterion;
//  double *soln_tur;
  int *iblank;

  int nvar, narray;
  int *mvar;
  char **varName;
  double **var;

  double *node_flag;
  void (*visit_extract_ini)(int, int, int, int, int, int*, int*, int*, int*,
                            double*, int*, int, char**, int, int*, double**, double*, MPI_Comm);
  void (*visit_extract_this)(int , double,
                             int, int, int, int, int, int*, int*, int*, int*,
                             double*, int*, int, char**, int, int*, double**, double*);

} visit_t;


/* the main global storage struct that holds everything including flow solvers and tioga */
typedef struct driver {

  char input_file[buff_size]; /** input file name */
  time_t input_file_mod_time; /** input file modification time */

  /** input file variables */
  int number_groups;          /**< total number of groups */
  int number_solver_so_files; /**< total number of solvers */
  int regrid_interval;        /**< interval between time steps to perform a regrid */
  int plot_interval;          /**< interval between time steps to output plots and restart files */
  int restart_counter;        /**< restarts the driver enter 0 for no restart and a positive integer for restart */
  int atm_group;              /**< identifies the atmospheric solver group since it looks like a node based finite volume solver */
  int number_time_steps;      /**< number of time steps to perform in the simulation */
  int off_body_group;         /**< identifies the off body group */
  int ncyc;                   /**< near-body sub-iteration count */

  /** group identification for all procs */
  int *group_solver_id;       /**< solver id for each group */
  int *group_num_procs;       /**< number of procs for each group */

  /** local group identification */
  int group_range[1][3];      /**< procs range for this group */
  int group;                  /**< my local group number indexing from 0 */
  int group_base1;            /**< my local group number indexing from 1 */
  int global_body_tag;        /**< my global mesh body tag number indexing from 1 used for replicate in nsu3d */

  /** mpi stuff for all procs */
  MPI_Comm mpicomm;
  int rank;                   /**< mpi rank across all cores */
  int rank_master_flag;       /**< if master rank then flag=1,else=0 */
  int num_procs;              /**< total number of processors */

  /** mpi stuff for this group */
  int group_rank;             /**< group rank */
  int group_master_flag;      /**< if group master rank then flag=1,else=0 */
  MPI_Group mpi_group;        /**< a mpi group that goes along with the group number */
  MPI_Comm group_comm;        /**< mpi communicator for the group */

  /** high order code identification */
  int high_order_group_flag;  /**< high order group flag: not high order=0, high order=1 */
  int high_order_mode_flag;   /**< if no groups are high order then don't call some tioga routines, if any groups high order then mode=1 else =0 */

  /** p4est code identification */
  int p4est_group_flag;       /**< if a group is a p4est codem then 1, else 0 */

  /** off body code identification */
  int off_body_group_flag;    /**< if a group is an off body group then = 1, else 0 */
  int off_body_mode_flag;     /**< if any group is an off body group then mode = 1, else = 0 */

//  MPI_Comm near_body_comm;
//  MPI_Comm off_body_comm;
//  int num_near_body_procs;
//  int num_off_body_procs;

  /** atmospheric code identification */
  int atm_group_flag;         /**< atm group flag: not atm=0, atm=1 */
  int atm_mode_flag;          /**< if any groups are atm then mode flag = 1, else 0 */
  int atm_compressed;
  int atm_nhex;

  /** number of field variables this is equal to 5 and is hard coded in many places */
  int nfield;                 /** the number of field variables */

  int *recvmap;
  int *sendmap;

  /** pointers to more data structs */
  flow_solver_t *flow;        /**< a pointer to an array of flow solver data structs  */
  tioga_t *tioga;             /**< a pointer to an array of tioga data structs */
  visit_t *visit;             /**< a pointer to an array of visit data structs */

  int spin_test;
  int increase_ncyc_flag;

} driver_t;

#endif /* defs_h */
