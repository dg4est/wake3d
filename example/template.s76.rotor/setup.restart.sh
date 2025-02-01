#!/bin/bash
# This script moves all files and finds the next restart file.

update_dg4est () {
  unlink input.dg4est.v2
  ln -s input.dg4est.v2.restart input.dg4est.v2

  sed -i "/^restart_file/ s/WRK$OLDINDEX/WRK$NEWINDEX/" input.dg4est.v2.restart
  sed -i "/^restart_file/ s/$PREVID/$LASTID/" input.dg4est.v2.restart
}

update_nsu3d () {
    unlink input.nsu3d
    ln -s input.nsu3d.restart input.nsu3d

    unlink mesh_file.replicate
    ln -s mesh_file.replicate.restart mesh_file.replicate

    sed -i "/^WRK/ s/$OLDINDEX/$NEWINDEX/" input.nsu3d.restart
    sed -i "/^WRK/ s/$NSU_PREVID/$NSU_LASTID/" input.nsu3d.restart
}

# --------------------------------------------------------------------------- #

#====================== #
# ENTRY POINT TO SCRIPT #
# ===================== #

# parse restart_count.txt
OLDINDEX=`sed -n 1p restart_count.txt`
NEWINDEX=$[$OLDINDEX+1]
echo "Old WRK Index: " $OLDINDEX
echo "New WRK Index: " $NEWINDEX

# --------------------------------------------------------------------------- #
echo "Moving work directories..."
mv group0/WRK group0/WRK$NEWINDEX
mv group1/WRK group1/WRK$NEWINDEX
mv group2/WRK group2/WRK$NEWINDEX

# --------------------------------------------------------------------------- # 
PREVID=`sed -n 34p group0/input.dg4est.v2 | sort -V | grep -o '[0-9]*' | tail -1`
G0WRK="group0/WRK$NEWINDEX/checkpoint/."
LASTID=`ls $G0WRK | sort -V | grep -o '[0-9:]*' | tail -1`
NSU_PREVID=${PREVID:1}
NSU_LASTID=${LASTID:1}

echo "The previous id: " $PREVID
echo "The last restart file available:" $LASTID

# --------------------------------------------------------------------------- # 
echo " "
echo "Changing input files..."
# ======= #
# Group 0 #
# ======= #
cd group0
update_dg4est
cd ..

# ======= #
# Group 1 #
# ======= #
cd group1
update_nsu3d
cd ..

# ======= #
# Group 2 #
# ======= #
cd group2
update_nsu3d
cd ..

# ======= #
# cdriver #
# ======= #
sed -i "/^restart_counter/ s/$PREVID/$LASTID/" input.driver

# ====================== #
# Update restart counter #
# ====================== #
echo $NEWINDEX > restart_count.txt

echo "Complete"
