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
echo "|  Demo script for PBxplore: multiple conformation analysis              |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

pause

# create and move into the demo directory
if [ ! -d "demo2_statistics_tmp" ];
then
    mkdir demo2_statistics_tmp
fi
cd demo2_statistics_tmp
cp ../demo2/* ./

echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignment: many conformations of the same protein     |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "Test with the PSI domain of the human beta3 integrin"
echo "90 conformations issued from 1 molecular dynamics simulation"
echo "PBassign -p psi_md_traj_1.pdb -o psi_md_traj_1"
pause
PBassign -p psi_md_traj_1.pdb -o psi_md_traj_1

echo  -e "\n"
echo "Test with the PSI domain of the human beta3 integrin"
echo "270 conformations issued from 3 molecular dynamics simulation"
echo "PBassign -p psi_md_traj_1.pdb -p psi_md_traj_2.pdb -p psi_md_traj_3.pdb -o psi_md_traj_all"
pause
PBassign -p psi_md_traj_1.pdb -p psi_md_traj_2.pdb -p psi_md_traj_3.pdb -o psi_md_traj_all

echo  -e "\n"
echo "Test with the Barstar protein"
echo "51 conformations issued from 1 molecular dynamics simulation"
echo "PBassign -x barstar_md_traj.xtc -g barstar_md_traj.gro -o barstar_md_traj"
pause
PBassign -x barstar_md_traj.xtc -g barstar_md_traj.gro -o barstar_md_traj || echo "PBassign failed; is MDAnalysis installed?"


echo  -e "\n"
echo "#---------------------------------------------------------------------------#"
echo "|                                                                           |"
echo "|  Compute the frequency of PBs at each position along the protein sequence |"
echo "|                                                                           |"
echo "#---------------------------------------------------------------------------#"

echo  -e "\n"
echo "with one input file"
echo "PBcount -f psi_md_traj_1.PB.fasta -o psi_md_traj_1"
pause
PBcount -f psi_md_traj_1.PB.fasta -o psi_md_traj_1

echo  -e "\n"
echo "with several input files"
echo "PBcount -f psi_md_traj_1.PB.fasta -f psi_md_traj_2.PB.fasta -f psi_md_traj_3.PB.fasta -o psi_md_traj_all"
pause
PBcount -f psi_md_traj_1.PB.fasta -f psi_md_traj_2.PB.fasta -f psi_md_traj_3.PB.fasta -o psi_md_traj_all

echo  -e "\n"
echo "with one input file and the --first-residue option"
echo "PBcount -f psi_md_traj_1.PB.fasta -o psi_md_traj_1_shifted --first-residue 20"
pause
PBcount -f psi_md_traj_1.PB.fasta -o psi_md_traj_1_shifted --first-residue 20


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Generate distribution of PBs, Neq and logo representation of PBs      |"
echo "|  along protein sequence                                                |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "PBstat -f psi_md_traj.PB.count -o psi_md_traj --map --neq --logo"
pause
PBstat -f psi_md_traj_all.PB.count -o psi_md_traj --map --neq --logo || echo 'The command failed, is weblogo installed?'

echo  -e "\n"
echo "Change the file format of the images"
echo "PBstat -f psi_md_traj.PB.count -o psi_md_traj --map --neq --logo --image-format pdf"
pause
PBstat -f psi_md_traj_all.PB.count -o psi_md_traj --map --neq --logo --image-format pdf || echo 'The command failed, is weblogo installed?'

echo  -e "\n"
echo "Define a residue frame (--residue-min and --residue-max options)"
echo "PBstat -f psi_md_traj.PB.count -o psi_md_traj --map --neq --logo --residue-min 10 --residue-max 30"
pause
PBstat -f psi_md_traj_all.PB.count -o psi_md_traj --map --neq --logo --residue-min 10 --residue-max 30 || echo 'The command failed, is weblogo installed?'


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Demo completed!                                                       |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"
echo
echo "Look at *.PB.* files in the demo2_statistics_tmp directory."
pwd
ls -lh
echo "Do not forget to delete demo2_statistics_tmp directory when you will be done with this demo."

