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


echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Make statistics (Neq)                                                 |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "../PBstat.py -f md_traj.PB.count -o md_traj"
pause
../PBstat.py -f md_traj.PB.count -o md_traj


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Make statistics with several PB assignment fasta files options        |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"


echo  -e "\n"
echo "Test with several 'PB.count' files"
echo "../PBstat.py -f md_traj_1.PB.count -f md_traj_2.PB.count -f md_traj_3.PB.count -o test_output"
pause
../PBstat.py -f md_traj_1.PB.count -f md_traj_2.PB.count -f md_traj_3.PB.count -o test_output


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Making statitics: playing with Neq graphic                            |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"


echo  -e "\n"
echo "Test with --neq-residue-lower, --neq-residue-upper, and --neq-residue-shift options"
echo "../PBstat.py -f md_traj.PB.count -o md_traj --neq-residue-lower=10 --neq-residue-upper=30 --neq-residue-shift=5"
pause
../PBstat.py -f md_traj.PB.count -o md_traj --neq-residue-lower=10 --neq-residue-upper=30 --neq-residue-shift=5


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

