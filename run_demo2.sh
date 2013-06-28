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
mkdir demo2_tmp
cd demo2_tmp
cp ../demo2/* ./

echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignment: many conformations of the same protein     |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "Test with the 225 conformations of the beta3 integrin"
echo "../PBassign.py -p md_traj.pdb -o md_traj"
pause
../PBassign.py -p md_traj.pdb -o md_traj


echo "#---------------------------------------------------------------------------#"
echo "|                                                                           |"
echo "|  Compute the frequency of PBs at each position along the protein sequence |"
echo "|                                                                           |"
echo "#---------------------------------------------------------------------------#"

echo  -e "\n"
echo "with one input file"
echo "../PBcount.py -f md_traj.PB.fasta -o md_traj"
pause
../PBcount.py -f md_traj.PB.fasta -o md_traj


echo  -e "\n"
echo "with several input files"
echo "../PBcount.py -f md_traj_1.PB.fasta -f md_traj_2.PB.fasta -f md_traj_3.PB.fasta -o test_output"
pause
../PBcount.py -f md_traj_1.PB.fasta -f md_traj_2.PB.fasta -f md_traj_3.PB.fasta -o test_output



echo  -e "\n"
echo "with one input file and the --residue-shift option"
echo "../PBcount.py -f md_traj.PB.fasta -o md_traj2 --residue-shift 20"
pause
../PBcount.py -f md_traj.PB.fasta -o md_traj2 --residue-shift 20


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Generate distribution of PBs, Neq and logo representation of PBs      |"
echo "|  along protein sequence                                                |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "../PBstat.py -f md_traj.PB.count -o md_traj --map --neq --logo"
pause
../PBstat.py -f md_traj.PB.count -o md_traj --map --neq --logo


echo  -e "\n"
echo "Test with --residue-min and --residue-max options"
echo "../PBstat.py -f md_traj.PB.count -o md_traj --map --neq --logo --residue-min 10 --residue-max 30"
pause
../PBstat.py -f md_traj.PB.count -o md_traj --map --neq --logo --residue-max 10 --residue-min 30


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Demo completed!                                                       |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"
echo 
echo "Look at *.PB.* files in the demo directory."
pwd
ls -lh 

