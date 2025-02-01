#!/bin/bash
# This script moves all files and finds the next restart file.

setup_dg4est () {
  cp originals/input.dg4est.v2.initial .
  cp originals/input.dg4est.v2.restart .

  unlink input.dg4est.v2
  ln -s input.dg4est.v2.initial input.dg4est.v2
}

setup_nsu3d () {
   cp originals/input.nsu3d.initial .
   cp originals/input.nsu3d.restart .
   cp originals/mesh_file.replicate.initial .
   cp originals/mesh_file.replicate.restart .

    unlink input.nsu3d
    ln -s input.nsu3d.initial input.nsu3d

    unlink mesh_file.replicate
    ln -s mesh_file.replicate.initial mesh_file.replicate
}

# --------------------------------------------------------------------------- #
# ============ #
# parse inputs
# ============ #
echo " "
echo "Setting up initial simulation directories..."

# ======= #
# Group 0 #
# ======= #
cd group0
setup_dg4est
cd ..

# ======= #
# Group 1 #
# ======= #
cd group1
setup_nsu3d
cd ..

# ======= #
# Group 2 #
# ======= #
cd group2
setup_nsu3d
cd ..

# ======= #
# cdriver #
# ======= #
sed -i "s/restart_counter:.*/restart_counter: 0000000/g" input.driver

echo "0" > restart_count.txt
