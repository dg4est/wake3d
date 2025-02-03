//
//  driver.h
//  cdriver
//
//  Created by Michael Brazell on 1/4/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#ifndef driver_h
#define driver_h

/* header files */
#include "input.h"
#include "defs.h"
#include "colors.h"

/* system header files */
#include <sys/stat.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#ifndef __APPLE__
#  include <malloc.h>
#endif

/* Check if glibc is at least 2.33 to use mallinfo2() */
#if __GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33)
#  define MALLINFO_T mallinfo2
#else
#  define MALLINFO_T mallinfo
#endif

#define LINEBREAK {int e; printf(" "); for(e = 0; e < 60; ++e) printf("-"); printf(" \n");}
#define LINEBREAKSET {int e; printf("+"); for(e = 0; e < 60; ++e) printf("="); printf("+\n");}

#if 1
static inline
size_t memory_usage(int mpi_rank,int timestep,int display,int details,int write2file){
#ifndef __APPLE__
    const double B2GB = 1.0/(1024.0 * 1024.0 * 1024.0);
    const double B2MB = 1.0/(1024.0 * 1024.0);
    struct MALLINFO_T mi = MALLINFO_T();

    /* get all the memory currently allocated to user by malloc, etc. */
    size_t used_mem = (size_t) mi.hblkhd
                    + (size_t) mi.usmblks
                    + (size_t) mi.uordblks;

    /* print out concise malloc info line */
    if (display) {
        printf(GREEN);
        printf("+======================= MEMORY ALLOCATION ===========================+\n");
        printf(BIWHITE " USAGE Rank[%d] Step[%d]:" GREEN " %f GiB : %f MiB : %zu Bytes  \n",mpi_rank,timestep,used_mem*B2GB,used_mem*B2MB,used_mem);
        if (details) {
            printf(COLOR_OFF " --------------------------------------------------------------------- \n");
            printf(BIWHITE " FIELDS                                    BYTES        MiB\n" GREEN);
            printf(" Total non-mmapped bytes       (arena): %12zu\t%f\n", (size_t)mi.arena   ,mi.arena   *B2MB);
            printf(" # of free chunks            (ordblks): %12zu    \n", (size_t)mi.ordblks                  );
            printf(" # of free fastbin blocks     (smblks): %12zu    \n", (size_t)mi.smblks                   );
            printf(" # of mapped regions           (hblks): %12zu    \n", (size_t)mi.hblks                    );
            printf(" Bytes in mapped regions      (hblkhd): %12zu\t%f\n", (size_t)mi.hblkhd,  mi.hblkhd  *B2MB);
            printf(" Max. total allocated space  (usmblks): %12zu\t%f\n", (size_t)mi.usmblks, mi.usmblks *B2MB);
            printf(" Free bytes held in fastbins (fsmblks): %12zu\t%f\n", (size_t)mi.fsmblks, mi.fsmblks *B2MB);
            printf(" Total allocated space      (uordblks): %12zu\t%f\n", (size_t)mi.uordblks,mi.uordblks*B2MB);
            printf(" Total free space           (fordblks): %12zu\t%f\n", (size_t)mi.fordblks,mi.fordblks*B2MB);
            printf(" Topmost releasable block   (keepcost): %12zu\t%f\n", (size_t)mi.keepcost,mi.keepcost*B2MB);
        }
        printf("+=====================================================================+\n");
        printf(COLOR_OFF);
        fflush(stdout);
    }

    /* write out concise malloc info line */
    if (write2file) {
        char filename[128];
        snprintf(filename, sizeof(filename), "wake3d.mem.%d.txt", mpi_rank);

        FILE *fp = fopen(filename,"a+");
        fprintf(fp,"%6d\t%f GiB : %f MiB : %zu Bytes\n",timestep,used_mem*B2GB,used_mem*B2MB,used_mem);
        fclose(fp);
    }
    return used_mem;
#else
    return 0.0;
#endif
}
#endif

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

// Custom function to find the last occurrence of a substring (like strrchr but for strings)
static
const char* strrstr(const char *haystack, const char *needle) {
    const char *result = NULL;
    const char *p = strstr(haystack, needle);
    while (p) {
        result = p;
        p = strstr(p + 1, needle); // Keep searching for later occurrences
    }
    return result;
}

static
void extractLibraryName(const char *path, char *output) {
    const char *libPos = strrstr(path, "lib"); // Find last occurrence of "lib"
    if (!libPos) {
        output[0] = '\0'; // Return empty string if "lib" is not found
        return;
    }

    libPos += 3; // Move past "lib"

    const char *soPos = strrstr(libPos, ".so"); // Find last occurrence of ".so"
    if (!soPos || soPos <= libPos) {
        output[0] = '\0'; // Return empty string if ".so" not found or misplaced
        return;
    }

    size_t len = soPos - libPos; // Length of extracted name
    if(len >= buff_size) len = buff_size - 1; // Ensure no overflow

    strncpy(output, libPos, len);
    output[len] = '\0'; // Null-terminate the string
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
