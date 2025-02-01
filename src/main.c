//
//  main.c
//  cdriver
//
//  Created by Michael Brazell on 1/4/17.
//  Copyright © 2017 Michael J. Brazell. All rights reserved.
//

#include <stdio.h>
#include "driver.h"
#include "defs.h"

/* global storage pointer */
driver_t *d;

int main(int argc, char **argv) {
    driver_initialize(argc, argv);
        (d->spin_test) ? spin_test() : driver_go();
    driver_finalize();
    return 0;
}
