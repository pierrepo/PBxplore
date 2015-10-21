#! /bin/bash

# 2013 - P. Poulain, A. G. de Brevern

# exit script at first error
set -e

function pause(){
    read -r -s -n1 -p "Press any key to continue."
    echo
}


echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Demo script for PBxplore: multiple conformation clustering            |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

pause

# create and move into the demo directory
DATA_PATH=$(PBdata)
DEMO_PATH=demo2_clusters_tmp
INPUT_FILES=(psi_md_traj_all.PB.fasta)
mkdir -p $DEMO_PATH
for input_file in ${INPUT_FILES[@]}
do
    cp ${DATA_PATH}/${input_file} $DEMO_PATH
done

cd $DEMO_PATH

echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Cluster structures                                                    |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "Produce with 3 clusters (--clusters option)"
echo "PBclust -f psi_md_traj_all.PB.fasta -o psi_md_traj_all_3 --clusters 3"
pause
PBclust -f psi_md_traj_all.PB.fasta -o psi_md_traj_all_3 --clusters 3

echo  -e "\n"
echo "Compare all sequences against the first one (--compare option)"
echo "PBclust -f psi_md_traj_all.PB.fasta -o psi_md_traj_all --compare"
pause
PBclust -f psi_md_traj_all.PB.fasta -o psi_md_traj_all --compare


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Demo completed!                                                       |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"
echo
echo "Look at *.PB.* files in the demo2_clusters_tmp directory."
pwd
ls -lh
echo "Do not forget to delete demo2_clusters_tmp directory when you will be done with this demo."

