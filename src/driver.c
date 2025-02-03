//
//  driver.c
//  cdriver
//
//  Created by Michael Brazell on 1/4/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

/* header files */
#include "driver.h"
#include "input.h"
#include "alloc.h"

/* system header files */
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>

extern driver_t *d;

void driver_go(){
    double a,b;
    int i,ind;
    double t1,t2;

    ind = d->restart_counter;
    if(d->rank_master_flag){LINEBREAK};

    /* output solution if initial condition */
    if(!d->restart_counter) {
        flow_output_solution(ind);
        visit_extract(ind);
    }
    ind++;

    MPI_Barrier(d->mpicomm);
    if(d->rank_master_flag) {LINEBREAKSET; fflush(stdout);}

    /* perform time steps */
    for (i = 1; i <= d->number_time_steps; ++ind,++i) {
        t1 = MPI_Wtime();
        if (d->rank_master_flag) {
            printf(" [wake3d]   Time Step: %7d  --sec--\n",i);
            fflush(stdout);
        }

        /* display memory */
        //if(d->group_master_flag) memory_usage(d->rank,i,1,1,1);

        /* check input file changes */
        check_input_file();

        /* =================== */
        /* perform a time step */
        /* =================== */
        driver_time_step();

        /* ================= */
        /* flow solvers plot */
        /* ================= */
        if (d->visit->visit_flag) {
            a = MPI_Wtime();
                visit_extract(ind);
            b = MPI_Wtime();
            if(d->group_master_flag) printf("[Libsim] GID: %d took %f secs\n",d->group,b-a);
        }
        if(i%d->plot_interval == 0) flow_output_solution(ind);

        /* ==================================== */
        /* gather, regrid, and reconnect meshes */
        /* ==================================== */
        if (i%d->regrid_interval == 0) {
            driver_igbp_regrid(1); //efficient mode
            // flow_regrid();
            // driver_obc_regrid();
            // driver_igbp_regrid(0);
        }
        MPI_Barrier(d->mpicomm);
        t2 = MPI_Wtime();
        if(d->rank_master_flag){LINEBREAK};
        if(d->rank_master_flag){printf(" Step: %d, cpu time (sec): %f\n",i,t2-t1);}
        if(d->rank_master_flag){LINEBREAKSET};
    }
    flow_output_solution(ind);

    /* int j;
     * int numfringe = 0;
     * for(j=0;j<d->flow->nnode;j++){
     *   if(d->flow->iblank[j] == -1) ++numfringe;
     * }
     * if(numfringe) printf("proc %d has this many fringe points %d\n",d->rank,numfringe);
     */
}

void spin_test(){
    int ind = 0;
    int i;

    flow_output_solution(ind++);

    /* perform time steps */
    for (i = 1; i <= d->number_time_steps; ++ind,++i) {
        if(d->rank_master_flag) printf(" [wake3d] Time step: %d\n",i);

        driver_time_step_nearbody_move();
        if(i%d->plot_interval == 0) flow_output_solution(ind);

        /* gather, regrid, and reconnect meshes */
        if(i%d->regrid_interval == 0){
           driver_igbp_regrid(1);
        }
    }
    flow_output_solution(ind);
}

void check_input_file(){
    /* check input file changes */
    if (file_is_modified(d->input_file,d->input_file_mod_time)) {
        if(d->rank_master_flag) printf("\n[wake3d]\033[1;33m INPUT FILE CHANGED! Rereading...\033[0m\n");
        read_input_file_options();
    }
}

long mem_usage(){
    struct rusage usage;
    int i;

    getrusage(RUSAGE_SELF, &usage);

    /* convert to MB */
    long mem = (long) ((double) usage.ru_maxrss)/1024.0;

    /* gather all memory usage to proc 0 */
    long memall[d->num_procs];
    MPI_Gather(&mem,1,MPI_LONG,&memall,1,MPI_LONG,0,d->mpicomm);

    if (d->rank == 0) {
        FILE *fp;
        char filename[] = "memusage.dat";

        fp = fopen(filename,"a");
        for(i = 0; i < d->num_procs; i++) fprintf(fp,"%ld ",memall[i]);
        fprintf(fp,"\n");
        fclose(fp);
    }
    return mem;
}

void check_point(int tag){
    int i;

    MPI_Barrier(d->mpicomm);
    for (i = 0; i < d->num_procs; ++i) {
        if(d->rank == i) printf("[wake3d] rank: %d made it here %d\n",d->rank,tag);

        MPI_Barrier(d->mpicomm);
    }
}

void driver_time_step_nearbody_move(){
    int ncyc = 0;

    flow_real2ghost();
    tioga_dataupdate();
    MPI_Barrier(d->mpicomm);

    if (!d->off_body_group_flag) {
        d->flow->unsteady_update(&ncyc);
    }
}

void driver_time_step(){
    static int increase_ncyc_h = 0;
    double t1,t2,wtime,max_wtime;
    int ncyc;

    flow_real2ghost();
    tioga_dataupdate();
    MPI_Barrier(d->mpicomm);

    t1 = MPI_Wtime(); /* start timer */
        ncyc = d->ncyc + increase_ncyc_h;
        d->flow->unsteady_update(&ncyc);
      //d->flow->steady_update(&ncyc);
    t2 = MPI_Wtime(); /* stop timer */

    wtime = t2-t1;
    MPI_Reduce(&wtime,&max_wtime,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,d->group_comm);
    if(d->group_rank==0) {printf("   Grp[%4d][%8s] SOLVER: %f\n",d->group,d->flow->solver_so_name,max_wtime); fflush(stdout);}

    /* check if ncyc is modified on fly based on wall-time */
    if(!d->increase_ncyc_flag) return;

    /* ========================================================== */
    /* Warning: Assumes all near-body groups do the same work.    */
    /* Fix: Check all groups and have each group modify its ncyc. */
    /* ========================================================== */

    /* find near body and off body times */
    double nbt,obt;
    if(d->off_body_group_flag){
        obt = t2-t1;
        nbt = 0.0;
    } else {
        obt = 0.0;
        nbt = t2-t1;
    }

    /* find the max near body and off body time */
    double maxobt,maxnbt;
    MPI_Allreduce(&obt,&maxobt,1,MPI_DOUBLE,MPI_MAX,d->mpicomm);
    MPI_Allreduce(&nbt,&maxnbt,1,MPI_DOUBLE,MPI_MAX,d->mpicomm);

    /* ensure that everything is positive */
    if(maxnbt < 0.0) return;
    if(maxobt < 0.0) return;

    /* if near-body is taking less time than off body set flag */
    int tflag = (maxnbt < 0.9*maxobt) ? 1:0;

    /* all reduce the flag and increment ncyc */
    int increase_ncyc = 0;
    MPI_Allreduce(&tflag,&increase_ncyc,1,MPI_INT,MPI_MAX,d->mpicomm);
    increase_ncyc_h += increase_ncyc;
}

void driver_obc_regrid(){
    MPI_Barrier(d->mpicomm);
    double t1 = MPI_Wtime();    /* gather obc */
        driver_gather_obc();
        MPI_Barrier(d->mpicomm);
    double t2 = MPI_Wtime();    /* regrid */
        flow_obc_regrid();
        MPI_Barrier(d->mpicomm);
    double t3 = MPI_Wtime();    /* register + connect grids */
        flow_set_pointers();
        tioga_register_and_connect(0);
        MPI_Barrier(d->mpicomm);
    double t4 = MPI_Wtime();

    if (d->rank_master_flag) {
        LINEBREAK;
        printf(" [wake3d] REGRID/OVERSET TIMINGS (sec):\n");
        printf("          IGBP:    %f\n",t2-t1);
        printf("          REGRID:  %f\n",t3-t2);
        printf("          OVERSET: %f\n",t4-t3);
    }
}

void driver_igbp_regrid(int efficient){
    /* ======================================================================= */
    /* Non-efficient gather is needed to initialize a tree structure in dg4est */
    /* Use non-efficient on first igbp call and efficient on later ones        */
    /* ======================================================================= */
    MPI_Barrier(d->mpicomm);
    double t1 = MPI_Wtime();    /* gather obc */
        if (efficient) {
            /* could use the last one but doing this for safety */
            setup_p4est_maps(d->flow->xgeom_translated,&d->flow->nnode);
            driver_gather_igbp_efficient();
        } else {
            driver_gather_igbp();
        }
        MPI_Barrier(d->mpicomm);
    double t2 = MPI_Wtime();    /* regrid */
        flow_igbp_regrid();
        MPI_Barrier(d->mpicomm);
    double t3 = MPI_Wtime();    /* register + connect grids */
        flow_set_pointers();
        tioga_register_and_connect(0);
        MPI_Barrier(d->mpicomm);
    double t4 = MPI_Wtime();

    if (d->rank_master_flag) {
        LINEBREAK;
        printf(" [wake3d] REGRID/OVERSET TIMINGS (sec):\n");
        printf("          IGBP:    %e\n",t2-t1);
        printf("          REGRID:  %e\n",t3-t2);
        printf("          OVERSET: %e\n",t4-t3);
    }
}

void display_optional_inputs(){
    if(d->rank == 0){
        printf("+==========================================+\n");
        printf(" WAKE3D Optional Inputs:                    \n");
        printf(" ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾                     \n");
        printf(" regrid_interval: %d\n",d->regrid_interval);
        printf(" plot_interval: %d\n",d->plot_interval);
        printf(" number_time_steps: %d\n",d->number_time_steps);
        printf(" ncyc: %d\n",d->ncyc);
        printf(" variable_ncyc_flag: %d\n",d->increase_ncyc_flag);
        printf("+==========================================+\n");
        printf("\n");
    }
}

void read_input_file(int argc, char **argv){
    char default_filename[] = "input.wake3d";
    char *filename = (argc < 2) ? default_filename:argv[1];

    char keyword[buff_size];
    char cwd[buff_size];
    struct stat file_stat;
    int group_temp;
    int sum;
    int err;
    int i;

    /* check if file exists */
    if (access(filename, F_OK ) == -1) {
        printf("[wake3d] Input file does not exist.\n");
        printf("         Usage: ./wake3d.mpi <input.file>\n");
        exit(1);
    }

    if (getcwd(cwd,sizeof(cwd)) != NULL) {
        char *pwd = trimwhitespace(cwd);

        /* save input file name with full path*/
        strcpy(d->input_file,pwd);
        strcat(d->input_file,"/");
        strcat(d->input_file,basename(filename));
    }

    /* file state for checking runtime modifications */
    stat(filename,&file_stat);
    d->input_file_mod_time = file_stat.st_mtime;

    /* mandatory inputs */
    err = find_keyword_integer(filename, "number_groups:", &d->number_groups, 1);
    err = find_keyword_integer(filename, "number_solver_so_files:", &d->number_solver_so_files, 1);

    if (d->number_groups < 1 && d->rank == 0) {
        printf("\033[1;31m"
               "[wake3d] ERROR input error need at least one group %d\n"
               "\033[0m\n",d->number_groups);
    }

    if (d->number_solver_so_files < 1 && d->rank == 0) {
        printf("\033[1;31m"
               "[wake3d] ERROR input error need at least one solver %d\n"
               "\033[0m\n",d->number_solver_so_files);
    }

    /* these are defaults for non-mandatory inputs */
    d->regrid_interval = 10000000;
    d->plot_interval = 10000000;
    d->restart_counter = 0;
    d->atm_group = -1;
    d->number_time_steps = 10000000;
    d->off_body_group = -1;
    d->flow->translation[0] = 0.0;
    d->flow->translation[1] = 0.0;
    d->flow->translation[2] = 0.0;
    d->flow->receptor_only_flag = 0;
    d->ncyc = 25;
    d->increase_ncyc_flag = 0;

    d->spin_test = 0;

    /* non-mandatory inputs */
    find_keyword_integer(filename, "regrid_interval:",  &d->regrid_interval,  0);
    find_keyword_integer(filename, "plot_interval:",    &d->plot_interval,    0);
    find_keyword_integer(filename, "restart_counter:",  &d->restart_counter,  0);
    find_keyword_integer(filename, "atm_group:",        &d->atm_group,        0);
    find_keyword_integer(filename, "number_time_steps:",&d->number_time_steps,0);
    find_keyword_integer(filename, "off_body_group:",   &d->off_body_group,   0);
    find_keyword_integer(filename, "ncyc:",             &d->ncyc,             0);
    find_keyword_integer(filename, "spin_test:",        &d->spin_test,        0);
    find_keyword_integer(filename, "variable_ncyc_flag:",&d->increase_ncyc_flag,0);

    /* check for nonsense inputs */
    if(d->number_time_steps < 0) d->number_time_steps = 10000000;
    if(d->regrid_interval   < 0) d->regrid_interval = d->number_time_steps+1;
    if(d->plot_interval     < 0) d->plot_interval = d->number_time_steps+1;
    if(d->restart_counter   < 0) d->restart_counter = 0;

    d->group_solver_id = my_alloc(int,d->number_groups);
    d->group_num_procs = my_alloc(int,d->number_groups);

    for (i = 0; i < d->number_groups; i++) {
        sprintf(keyword, "group%d:",i);
        find_keyword_two_integers(filename, keyword,
                                  &d->group_solver_id[i],
                                  &d->group_num_procs[i],
                                  1);
    }

    sum = 0;
    for(i = 0; i < d->number_groups; ++i) sum += d->group_num_procs[i];

    if (sum != d->num_procs && d->rank == 0) {
        printf("\033[1;31m"
               "[wake3d] Warning: total number of MPI Ranks %d does not match sum of groups ranks %d\n"
               "\033[0m\n",d->num_procs,sum);
    }

    sum = 0;
    for (i = 0; i < d->number_groups; i++) {
        if (d->rank >= sum && d->rank < sum+d->group_num_procs[i]) {
          group_temp=i;
          break;
        }
        sum += d->group_num_procs[i];
    }

    sprintf(keyword, "translation%d:", group_temp);
    find_keyword_three_doubles(filename,keyword,
                               &d->flow->translation[0],
                               &d->flow->translation[1],
                               &d->flow->translation[2],
                               0);

    sprintf(keyword,"receptor_only%d:",group_temp);
    find_keyword_integer(filename,keyword,&d->flow->receptor_only_flag,1);

    sprintf(keyword, "so%d:", d->group_solver_id[group_temp]);
    err = find_keyword_string(filename, keyword, d->flow->solver_so_file, 1);
    err = find_keyword_string(filename, "tioga_so_file:",d->tioga->tioga_so_file, 1);
    err = find_keyword_string(filename, "visit_so_file:",d->visit->visit_so_file, 0);

    extractLibraryName(d->flow->solver_so_file,d->flow->solver_so_name);

    d->visit->visit_flag = (err) ? 0:1;
    if(err && d->rank == 0) printf("[wake3d] Warning: visit shared library not found\n");

    /* write out optional inputs to screen */
    display_optional_inputs();
}

void read_input_file_options(){
    char *filename = d->input_file;
    struct stat file_stat;

    /* non-mandatory inputs */
    find_keyword_integer(filename, "regrid_interval:",   &d->regrid_interval,  0);
    find_keyword_integer(filename, "plot_interval:",     &d->plot_interval,    0);
    find_keyword_integer(filename, "number_time_steps:", &d->number_time_steps,0);
    find_keyword_integer(filename, "ncyc:",              &d->ncyc,             0);
    find_keyword_integer(filename, "variable_ncyc_flag:",&d->increase_ncyc_flag,0);

    /* write out optional inputs to screen */
    display_optional_inputs();

    /* check safe inputs */
    if(d->plot_interval     <= 0) d->plot_interval = INT_MAX;
    if(d->regrid_interval   <= 0) d->regrid_interval = INT_MAX;
    if(d->number_time_steps <= 0) d->number_time_steps = INT_MAX;

    /* file state for checking runtime modifications */
    stat(filename,&file_stat);
    d->input_file_mod_time = file_stat.st_mtime;
}

/* probably should not touch anything below this point */
void driver_initialize(int argc, char **argv){
    /* allocate one driver */
    d = my_alloc(driver_t,1);

    /* allocate one flow solver struct but could allocate more here */
    d->flow = my_alloc(flow_solver_t,1);
    d->flow->igbp = my_alloc(igbp_t,1);
    d->flow->obc = my_alloc(obc_t,1);

    /* allocate one tioga struct but could allocate more here */
    d->tioga = my_alloc(tioga_t,1);

    /* allocate visit struct */
    d->visit = my_alloc(visit_t,1);

    MPI_Init(&argc,&argv);

    d->mpicomm = MPI_COMM_WORLD;
    MPI_Comm_rank(d->mpicomm, &d->rank);
    MPI_Comm_size(d->mpicomm, &d->num_procs);

    read_input_file(argc, argv);

    /* use input file to create all the groups and comms */
    create_groups();

    /* identify groups and set flags for easy id later */
    set_flags();

    /* initialize all variables in flow solver and tioga */
    initialize_variables();

    /* load the function pointers from the dynamic libraries */
    flow_load_dynamic_library();
    tioga_load_dynamic_library();
    visit_load_dynamic_library();

    /* move into the group folders, this assumes the folders exist already */
    move_into_group_directory(d->group);

    /* initialize flow solver mpi */
    flow_initialize_group_mpi(d->group_comm);

    /* initialize tioga mpi */
    d->tioga->init(d->mpicomm);

    /* TODO: Implement into input file */
    /* initialize composite bodies */
    //hack
//    {
//        int ncomp = 1;
//        int nbodytags = 2;
//        int compbodytag = 0;
//        d->tioga->setcomposite(&ncomp);

//        int bodytags[2];
//        bodytags[0] = 1;
//        bodytags[1] = 2;
//        double searchTol = 1.0E-6;
//        d->tioga->register_composite_body(&compbodytag,bodytags,&nbodytags,&searchTol);
//    }

    /* initialize flow solver */
    d->flow->initialize();

    /* get pointers from flow solver */
    flow_set_pointers();

    /* now that body tags are known find unique global body tag */
    create_unique_body_tag();

    /* this checks the pointers to detect if high order and p4est */
    set_high_order_flags();

    /* check if flow solver contains bounding box function to set p4est flag */
    flow_set_p4est_flag();

    /* initialize visit */
    visit_init();

    /* this gets the pointers for igbp points */
    flow_get_igbp_data();

    /* if not a restart and atm mode then use set resolution to initialize flow */
    if(!d->restart_counter && d->atm_mode_flag){
        tioga_register_and_connect(1);
        flow_real2ghost();
        tioga_dataupdate();
    }

    /* choose one of these three make sure driver go matches this */
    /* Do NOT use efficient: required for restart combined with turning off off-body flag */
    int efficient = (d->restart_counter > 0) ? 1:0;
    driver_igbp_regrid(efficient);
    //driver_obc_regrid();
    //tioga_register_and_connect();
}

void create_unique_body_tag(){
   flow_solver_t *f = d->flow;

    if (f->body_tag < 1) {
        printf("ahh error this should not be called unless flow set pointers (get data) is called or your flow solver body tag is not base 1\n");
    }

    int max_body_tag = 0;
    MPI_Allreduce(&f->body_tag,&max_body_tag,1,MPI_INT,MPI_MAX,d->mpicomm);

    if (max_body_tag > 20 && d->rank_master_flag) {
        printf("warning you are using a lot of body tags and this will overallocate the hole map in tioga, either use less replicates or modify tioga \n");
    }

    /* this is unique but not contiguous if every group has different amounts of body tags */
    /* this is not a huge issue unless max_body_tag is big */
    d->global_body_tag = d->group*max_body_tag + f->body_tag;

    /* this is unique and does not require the allreduce however tioga will definitely overallocate in holemap
      so do not use this one unless tioga is modified */
    //d->global_body_tag = d->number_groups*(f->body_tag-1) + d->group + 1;

    /* could also do a gather and stack the global body tag in order but will require more communication than an allreduce */
}

void create_groups(){
    MPI_Group orig_group;
    int i,sum;

    /* find group proc range */
    sum = 0;
    for(i=0;i<d->number_groups;i++){
        if (d->rank >= sum && d->rank < sum+d->group_num_procs[i]) {
            d->group_range[0][0] = sum;
            d->group_range[0][1] = sum+d->group_num_procs[i]-1;
            d->group_range[0][2] = 1;
            d->group = i;
        }
        sum+=d->group_num_procs[i];
    }

    /* Extract the original group handle */
    MPI_Comm_group(d->mpicomm, &orig_group);

    /* Divide tasks into distinct groups based upon group range */
    MPI_Group_range_incl(orig_group, 1, d->group_range, &d->mpi_group);

    /* create new communicator for group */
    MPI_Comm_create(d->mpicomm, d->mpi_group, &d->group_comm);

    MPI_Comm_rank(d->group_comm, &d->group_rank);
    return;

//  fixme i think this is buggy
//  /* identify near and off body groups */
//  int group_range[d->number_groups][3];
//  int num_near_body;
//  num_near_body = 0;
//  sum = 0;
//
//  for(i=0;i<d->number_groups;i++){
//    if(i != d->off_body_group) {
//      group_range[num_near_body][0]=sum;
//      group_range[num_near_body][1]=sum+d->group_num_procs[i]-1;
//      group_range[num_near_body][2]=1;
//      num_near_body = num_near_body+1;
//    }
//    sum+=d->group_num_procs[i];
//  }
//
//  MPI_Group near_body_group;
//  MPI_Group off_body_group;
//
//  /* Divide tasks near body group based upon near body group range */
//  MPI_Group_range_incl(orig_group, num_near_body, group_range, &near_body_group);
//  MPI_Comm_create(d->mpicomm, near_body_group, &d->near_body_comm);
//  MPI_Comm_size(d->near_body_comm, &d->num_near_body_procs);
//
//  /* Divide tasks into off body group based upon excluding near body group range */
//  MPI_Group_range_excl(orig_group, num_near_body, group_range, &off_body_group);
//  MPI_Comm_create(d->mpicomm, off_body_group, &d->off_body_comm);
//  MPI_Comm_size(d->off_body_comm, &d->num_off_body_procs);
}

void initialize_variables(){
    flow_solver_t *f = d->flow;

    d->nfield = 5;

    f->group_handle = NULL;
    f->soln = NULL;
    f->xgeom = NULL;
    f->iblank = NULL;
    f->iwbcnode = NULL;
    f->iobcnode = NULL;
    f->ndc4 = NULL;
    f->ndc5 = NULL;
    f->ndc6 = NULL;
    f->ndc8 = NULL;
    f->iblank_cell = NULL;
    f->xgeom_translated = NULL;

    f->nnode = 0;
    f->nwbc = 0;
    f->nobc = 0;
    f->ntet = 0;
    f->nprism = 0;
    f->npyramid = 0;
    f->nhex = 0;

    f->node_res = NULL;
    f->cell_res = NULL;

    f->ncell_types = 4;
    f->kstride4 = 4;
    f->kstride5 = 5;
    f->kstride6 = 6;
    f->kstride8 = 8;

    f->body_tag = 0;
    d->global_body_tag = 0;

    f->obc->nobc = 0;
    f->obc->x = NULL;
    f->obc->y = NULL;
    f->obc->z = NULL;
    f->obc->dx = NULL;

    f->igbp->n = 0;
    f->igbp->n_lcl = 0;
    f->igbp->index = NULL;
    f->igbp->dx_lcl = NULL;
    f->igbp->x = NULL;
    f->igbp->y = NULL;
    f->igbp->z = NULL;
    f->igbp->dx = NULL;

    tioga_t *t = d->tioga;
    t->group_handle = NULL;

    d->recvmap = NULL;
    d->sendmap = NULL;

    d->visit->group_handle = NULL;
}

void set_flags(){
    d->rank_master_flag = (d->rank == 0) ? 1:0;
    d->group_master_flag = (d->group_rank == 0) ? 1:0;

    d->atm_compressed = 0;
    d->atm_group_flag = (d->atm_group == d->group) ? 1:0;
    MPI_Allreduce(&d->atm_group_flag,&d->atm_mode_flag,1,MPI_INT,MPI_MAX,d->mpicomm);

    d->off_body_group_flag = (d->off_body_group == d->group) ? 1:0;
    MPI_Allreduce(&d->off_body_group_flag,&d->off_body_mode_flag,1,MPI_INT,MPI_MAX,d->mpicomm);

    /* set to zero and check flow_set_p4est_flag() */
    d->p4est_group_flag = 0;

    /* set to zero and check iblank_cell array for null later */
    d->high_order_group_flag = 0;
    d->high_order_mode_flag = 0;
}

void move_into_group_directory(int g){
    char string[1000];
    int ret;

    snprintf(string, sizeof(string), "group%d", g);
    ret = chdir(string);
    (void) ret;
}

void set_high_order_flags(){
    /* check if iblank_cell is null or not */
    d->high_order_group_flag = (&d->flow->iblank_cell[0]) ? 1:0;

    /* if any group is high order then we are in high order mode */
    MPI_Allreduce(&d->high_order_group_flag,&d->high_order_mode_flag,1,MPI_INT,MPI_MAX,d->mpicomm);
}

void driver_gather_igbp(){
    int i;
    int* counts;
    int* displacements;
    double *x,*y,*z,*dx;
    igbp_t *igbp;
    igbp = d->flow->igbp;

    x  = my_alloc(double,igbp->n_lcl);
    y  = my_alloc(double,igbp->n_lcl);
    z  = my_alloc(double,igbp->n_lcl);
    dx = my_alloc(double,igbp->n_lcl);

    counts = my_alloc(int,d->num_procs);
    displacements = my_alloc(int,d->num_procs);

    MPI_Allgather(&igbp->n_lcl,1,MPI_INT,counts,1,MPI_INT,d->mpicomm);

    displacements[0] = 0;
    for (i = 1; i < d->num_procs; i++) {
        displacements[i] = counts[i-1] + displacements[i-1];
    }

    igbp->n = 0;
    for(i = 0; i < d->num_procs; i++) igbp->n += counts[i];

    if(igbp->x)  my_free(igbp->x);
    if(igbp->y)  my_free(igbp->y);
    if(igbp->z)  my_free(igbp->z);
    if(igbp->dx) my_free(igbp->dx);

    igbp->x = my_alloc(double,igbp->n);
    igbp->y = my_alloc(double,igbp->n);
    igbp->z = my_alloc(double,igbp->n);
    igbp->dx = my_alloc(double,igbp->n);

    for (i = 0; i < igbp->n_lcl; i++) {
        x[i] = d->flow->xgeom_translated[3*(igbp->index[i]-1)];
        y[i] = d->flow->xgeom_translated[3*(igbp->index[i]-1)+1];
        z[i] = d->flow->xgeom_translated[3*(igbp->index[i]-1)+2];
        dx[i] = igbp->dx_lcl[i];
    }

    MPI_Allgatherv(x, igbp->n_lcl,MPI_DOUBLE,igbp->x, counts,displacements,MPI_DOUBLE,d->mpicomm);
    MPI_Allgatherv(y, igbp->n_lcl,MPI_DOUBLE,igbp->y, counts,displacements,MPI_DOUBLE,d->mpicomm);
    MPI_Allgatherv(z, igbp->n_lcl,MPI_DOUBLE,igbp->z, counts,displacements,MPI_DOUBLE,d->mpicomm);
    MPI_Allgatherv(dx,igbp->n_lcl,MPI_DOUBLE,igbp->dx,counts,displacements,MPI_DOUBLE,d->mpicomm);

#if 0
  if(d->rank == 0){
    char filename[] = "igbp_points.dat";
    FILE *fp = fopen(filename,"w");
    for(i = 0;i < igbp->n; i++) fprintf(fp,"%d %e %e %e %e\n",i,igbp->x[i],igbp->y[i],igbp->z[i],igbp->dx[i]);
    fclose(fp);
  }
#endif
    my_free(x);
    my_free(y);
    my_free(z);
    my_free(dx);
    my_free(counts);
    my_free(displacements);
}

void driver_gather_igbp_efficient(){
    MPI_Request *request;
    MPI_Status *status;
    int *num_igbp_pts,*displacements;
    int *sendmap,*recvmap;
    int nsend,nrecv;
    int i;

    if(!d->recvmap) {
        driver_gather_igbp();
        return;
    }
    igbp_t *igbp = d->flow->igbp;

    /* slide the mapping, insert rank, and count number of send/recv */
    nsend = nrecv = 0;
    for (i = 0; i < d->num_procs; i++){
        if(d->sendmap[i]) nsend++;
        if(d->recvmap[i]) nrecv++;
    }
    sendmap = my_alloc(int,nsend);
    recvmap = my_alloc(int,nrecv);

    nsend = nrecv = 0;
    for (i = 0; i < d->num_procs; i++) {
        if(d->sendmap[i]) sendmap[nsend++] = i;
        if(d->recvmap[i]) recvmap[nrecv++] = i;
    }
    num_igbp_pts = my_alloc(int,nrecv);
    displacements = my_alloc(int,nrecv);

    /* send the number of near-body points to all off-body cores and store them in counts */
    request = my_alloc(MPI_Request, 4*(nsend+nrecv));
    status = my_alloc(MPI_Status, 4*(nsend+nrecv));

    int tag = 1;
    int irnum = 0;

    /* off-body receives number of near-body points */
    /* near-body sends npts to off-body */
    for(i = 0; i < nrecv; i++) MPI_Irecv(&num_igbp_pts[i],1,MPI_INT,recvmap[i],tag,d->mpicomm,&request[irnum++]);
    for(i = 0; i < nsend; i++) MPI_Isend(&igbp->n_lcl,    1,MPI_INT,sendmap[i],tag,d->mpicomm,&request[irnum++]);
    MPI_Waitall(irnum,request,status);

    /* count total number of near body points collected on off body */
    igbp->n = 0;
    for(i = 0; i < nrecv; i++) igbp->n += num_igbp_pts[i];
    if(nsend) igbp->n = igbp->n_lcl;

    /* save displacements too for recv/send */
    if(nrecv) displacements[0] = 0;
    for(i = 1; i < nrecv; i++) displacements[i] = num_igbp_pts[i-1] + displacements[i-1];

    /* store all near-body points on the off-body */
    if(igbp->x)  my_free(igbp->x);
    if(igbp->y)  my_free(igbp->y);
    if(igbp->z)  my_free(igbp->z);
    if(igbp->dx) my_free(igbp->dx);

    igbp->x  = my_alloc(double,igbp->n);
    igbp->y  = my_alloc(double,igbp->n);
    igbp->z  = my_alloc(double,igbp->n);
    igbp->dx = my_alloc(double,igbp->n);

    if(nsend){
        for(i = 0; i < igbp->n_lcl; i++){
            igbp->x[i] = d->flow->xgeom_translated[3*(igbp->index[i]-1)+0];
            igbp->y[i] = d->flow->xgeom_translated[3*(igbp->index[i]-1)+1];
            igbp->z[i] = d->flow->xgeom_translated[3*(igbp->index[i]-1)+2];
            igbp->dx[i] = igbp->dx_lcl[i];
        }
    }

    int tagx = 10;
    int tagy = 11;
    int tagz = 12;
    int tagh = 13;
    irnum = 0;
    for(i = 0; i < nrecv; i++){
        MPI_Irecv(&igbp->x[ displacements[i]],num_igbp_pts[i],MPI_DOUBLE,recvmap[i],tagx,d->mpicomm,&request[irnum++]);
        MPI_Irecv(&igbp->y[ displacements[i]],num_igbp_pts[i],MPI_DOUBLE,recvmap[i],tagy,d->mpicomm,&request[irnum++]);
        MPI_Irecv(&igbp->z[ displacements[i]],num_igbp_pts[i],MPI_DOUBLE,recvmap[i],tagz,d->mpicomm,&request[irnum++]);
        MPI_Irecv(&igbp->dx[displacements[i]],num_igbp_pts[i],MPI_DOUBLE,recvmap[i],tagh,d->mpicomm,&request[irnum++]);
    }
    for(i = 0; i < nsend; i++){
        MPI_Isend(igbp->x, igbp->n_lcl,MPI_DOUBLE,sendmap[i],tagx,d->mpicomm,&request[irnum++]);
        MPI_Isend(igbp->y, igbp->n_lcl,MPI_DOUBLE,sendmap[i],tagy,d->mpicomm,&request[irnum++]);
        MPI_Isend(igbp->z, igbp->n_lcl,MPI_DOUBLE,sendmap[i],tagz,d->mpicomm,&request[irnum++]);
        MPI_Isend(igbp->dx,igbp->n_lcl,MPI_DOUBLE,sendmap[i],tagh,d->mpicomm,&request[irnum++]);
    }
    MPI_Waitall(irnum,request,status);

    /* deallocate memory */
    my_free(displacements);
    my_free(num_igbp_pts);
    my_free(sendmap);
    my_free(recvmap);
    my_free(request);
    my_free(status);
}

void driver_gather_obc(){
    int i;
    int* counts;
    int* displacements;
    double *x,*y,*z,*dx;
    flow_solver_t *f = d->flow;
    obc_t *obc = f->obc;

    counts = my_alloc(int,d->num_procs);
    displacements = my_alloc(int,d->num_procs);

    int nobc = 0;
    if(!d->off_body_group_flag) nobc = d->flow->nobc;
    MPI_Allgather(&nobc,1,MPI_INT,counts,1,MPI_INT,d->mpicomm);

    x  = my_alloc(double,nobc);
    y  = my_alloc(double,nobc);
    z  = my_alloc(double,nobc);
    dx = my_alloc(double,nobc);

    displacements[0] = 0;
    for (i = 1; i < d->num_procs; i++) {
        displacements[i] = counts[i-1] + displacements[i-1];
    }

    // check if nobc changes size in future using tioga fixme
    obc->nobc = 0;
    for (i = 0; i < d->num_procs; i++) {
        obc->nobc += counts[i];
    }

    if(obc->x)  my_free(obc->x);
    if(obc->y)  my_free(obc->y);
    if(obc->z)  my_free(obc->z);
    if(obc->dx) my_free(obc->dx);

    obc->x  = my_alloc(double,obc->nobc);
    obc->y  = my_alloc(double,obc->nobc);
    obc->z  = my_alloc(double,obc->nobc);
    obc->dx = my_alloc(double,obc->nobc);

    for (i = 0; i < nobc; i++) {
        x[i] = f->xgeom_translated[3*(f->iobcnode[i]-1)];
        y[i] = f->xgeom_translated[3*(f->iobcnode[i]-1)+1];
        z[i] = f->xgeom_translated[3*(f->iobcnode[i]-1)+2];
        dx[i] = 1.0;
    }

    MPI_Allgatherv(x,nobc,MPI_DOUBLE,obc->x,counts,displacements,MPI_DOUBLE,d->mpicomm);
    MPI_Allgatherv(y,nobc,MPI_DOUBLE,obc->y,counts,displacements,MPI_DOUBLE,d->mpicomm);
    MPI_Allgatherv(z,nobc,MPI_DOUBLE,obc->z,counts,displacements,MPI_DOUBLE,d->mpicomm);
    MPI_Allgatherv(dx,nobc,MPI_DOUBLE,obc->dx,counts,displacements,MPI_DOUBLE,d->mpicomm);

    my_free(x);
    my_free(y);
    my_free(z);
    my_free(dx);
    my_free(counts);
    my_free(displacements);
}


void driver_finalize(){
    flow_close_dynamic_library();
    tioga_close_dynamic_library();
    visit_close_dynamic_library();

    if(d->flow->cell_res) my_free(d->flow->cell_res);
    if(d->flow->node_res) my_free(d->flow->node_res);
    if(d->flow->obc->x) my_free(d->flow->obc->x);
    if(d->flow->obc->y) my_free(d->flow->obc->y);
    if(d->flow->obc->z) my_free(d->flow->obc->z);
    if(d->flow->obc->dx) my_free(d->flow->obc->dx);
    if(d->flow->obc) my_free(d->flow->obc);

    if(d->flow->igbp->x) my_free(d->flow->igbp->x);
    if(d->flow->igbp->y) my_free(d->flow->igbp->y);
    if(d->flow->igbp->z) my_free(d->flow->igbp->z);
    if(d->flow->igbp->dx) my_free(d->flow->igbp->dx);
    if(d->flow->igbp) my_free(d->flow->igbp);
    if(d->flow->xgeom_translated) my_free(d->flow->xgeom_translated);

    if(d->flow) my_free(d->flow);
    if(d->tioga) my_free(d->tioga);
    if(d->visit) my_free(d->visit);

    if(d->group_solver_id) my_free(d->group_solver_id);
    if(d->group_num_procs) my_free(d->group_num_procs);
    if(d->recvmap) my_free(d->recvmap);
    if(d->sendmap) my_free(d->sendmap);

    /* save these before deleting driver struct d */
    int rank = d->rank;
    MPI_Comm mpicomm = d->mpicomm;

    if(d) my_free(d);

    int max_alloc_count;
    int alloc_count = get_alloc_count();

    //if(alloc_count > 0) printf("[wake3d] Memory leak %d on proc %d\n",alloc_count,rank);
    //if(alloc_count < 0) printf("[wake3d] Alloc count is negative %d on proc %d\n",alloc_count,rank);

    alloc_count = abs(alloc_count);
    MPI_Reduce(&alloc_count,&max_alloc_count,1,MPI_INT,MPI_MAX,0,mpicomm);
    if(max_alloc_count == 0 && rank == 0) printf("[wake3d] YAY no memory leaks!\n");

    MPI_Finalize();
}
