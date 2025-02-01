//
//  driver.h
//  cdriver
//
//  Created by Michael Brazell on 1/4/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#ifndef driver_h
#define driver_h

//#include <mpi.h>
#include "input.h"
#include "defs.h"
#include <malloc.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>

#define LINEBREAK {int e; printf("+"); for(e = 0; e < 60; ++e) printf("="); printf("+\n");}

static inline
double memory_usage(int mpi_rank,int timestep,int display){

    /* get malloc info structure */
    struct mallinfo my_mallinfo = mallinfo();

    /*total memory reserved by the system for malloc currently */
    double reserved_mem = my_mallinfo.arena;

    /* get all the memory currently allocated to user by malloc, etc. */
    double used_mem = my_mallinfo.hblkhd
                    + my_mallinfo.usmblks
                    + my_mallinfo.uordblks;

    /* get memory not currently allocated to user but malloc controls */
    double free_mem = my_mallinfo.fsmblks
                    + my_mallinfo.fordblks;

    /* get number of items currently allocated */
    /* double number_allocated = my_mallinfo.ordblks + my_mallinfo.smblks; */

    /* Print out concise malloc info line */
    if(display && mpi_rank == 0){
        printf("Step[%d]: %f MB(%.0f) malloc: %f MB reserved (%.0f unused)\n",
            timestep,
            used_mem / (1024.0 * 1024.0),
            used_mem,
            reserved_mem / (1024.0 * 1024.0),
            free_mem);

 	if(mpi_rank == 0){
            FILE *fp;
	    char filename[] = "new_memusage.dat";
    	    fp=fopen(filename,"a+");
            fprintf(fp,"Step[%d]: %f MB(%.0f) malloc: %f MB reserved (%.0f unused)\n",
	            timestep,used_mem / (1024.0 * 1024.0),used_mem,reserved_mem / (1024.0 * 1024.0), free_mem);
            fclose(fp);
  	}

    }
    return used_mem;
}

static inline
int file_is_modified(const char *path,time_t oldMTime){
    struct stat file_stat;
    int err = stat(path, &file_stat);
    if (err != 0) {
        perror(" [file_is_modified] stat");
        exit(errno);
    }
    return file_stat.st_mtime > oldMTime;
}

static inline
char *trimwhitespace(char *str){
    char *end;

    // Trim leading space
    while(isspace((unsigned char)*str)) str++;

    if(*str == 0) return str;  // All spaces?

    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;

    // Write new null terminator character
    end[1] = '\0';
    return str;
}

/* main driver routines */
void driver_initialize(int argc, char **argv);
void driver_go();
void driver_finalize();

/* time stepping routines */
void driver_time_step();
void driver_gather_obc();
void driver_gather_igbp();
void driver_gather_igbp_efficient();
void driver_igbp_regrid(int efficient);
void driver_obc_regrid();

/* driver utilities */
void display_optional_inputs();
void read_input_file(int argc, char **argv);
void read_input_file_options();
void move_into_group_directory(int);
void create_groups();
void set_flags();
void initialize_variables();
void create_unique_body_tag();

void set_high_order_flags();
void check_point(int tag);
void search_p4est(double *xsearch,int *process_id, int *cell_id, int *npts);
void check_intersect_p4est(int *proc, int* overlap);
void setup_p4est_maps(double *xsearch, int *numpts);
long mem_usage();

/* flow solver routines */
void flow_load_dynamic_library();
void flow_close_dynamic_library();

void flow_initialize_group_mpi(MPI_Comm group_comm);
void flow_set_pointers();
void flow_output_solution(int t);

/* visit routines */
void visit_load_dynamic_library();
void visit_close_dynamic_library();

void visit_init();
void visit_extract(int icyc);

/* these routines are not required and are initialized with null
 they are pointed to functions in .so if they exist
 to use these check for null first */
void flow_regrid();
void flow_obc_regrid();
void flow_get_igbp_data();
void flow_igbp_regrid();
void flow_real2ghost();
void flow_iblank_check();
void flow_iblank_flipper();
void flow_set_p4est_flag();

void point_inclusion(int *npoint, double *x, int *cell_id);
void bounding_box_intersection(int *nbb, double *bb_xlo, double *bb_xhi, int *bb_flag);


/* tioga wrapper routines */
void tioga_load_dynamic_library();
void tioga_close_dynamic_library();
void tioga_init(MPI_Comm mpi_comm_in);
void tioga_register_and_connect(int user_spec_res);
void tioga_dataupdate();

void tioga_registergrid_data();
void tioga_set_highorder_callback();
void tioga_setcelliblank();
void tioga_performconnectivity_highorder();
void tioga_set_p4est();

void driver_time_step_nearbody_move();
void spin_test();
void check_input_file();

#endif /* driver_h */
