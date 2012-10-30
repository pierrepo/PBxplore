#! /bin/bash

# 2012 - P. Poulain, A. G. de Brevern

echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  demo script for PBxplore                                              |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

read -n1 -r -p "Press any key to continue."

# create and move into the demo directory
mkdir demo
cd demo
cp ../test/static/*pdb ./

echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignation: static structures                         |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "Test with 3ICH.pdb: RX structure, one chain"
echo "../PBassign.py -f 3ICH.pdb -o 3ICH"
read -n1 -r -p "Press any key to continue."
../PBassign.py -p 3ICH.pdb -o 3ICH


echo  -e "\n"
echo "Test with 1AY7.pdb: RX structure, complex with two chains"
echo "../PBassign.py -f 1AY7.pdb -o 1AY7"
read -n1 -r -p "Press any key to continue."
../PBassign.py -p 1AY7.pdb -o 1AY7


echo  -e "\n"
echo "Test with 2LFU.pdb: RMN structure, 10 models"
echo "../PBassign.py -f 2LFU.pdb -o 2LFU"
read -n1 -r -p "Press any key to continue."
../PBassign.py -p 2LFU.pdb -o 2LFU


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignation: static structures --phipsi option         |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"


echo  -e "\n"
echo "Test with 3ICH.pdb: RX structure, one chain"
echo "../PBassign.py -f 3ICH.pdb -o 3ICH --phipsi"
read -n1 -r -p "Press any key to continue."
../PBassign.py -p 3ICH.pdb -o 3ICH --phipsi


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignation: static structures --flat option           |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"


echo  -e "\n"
echo "Test with 2LFU.pdb: RMN structure, 10 models"
echo "../PBassign.py -f 2LFU.pdb -o 2LFU --flat"
read -n1 -r -p "Press any key to continue."
../PBassign.py -p 2LFU.pdb -o 2LFU --flat


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignation: several static structures                 |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "Test with several PDB files"
echo "../PBassign.py -d ./ -o all"
read -n1 -r -p "Press any key to continue."
../PBassign.py -p 3ICH.pdb -p 2LFU.pdb -o several

echo  -e "\n"
echo "Test with all PDB files from a directory"
echo "../PBassign.py -d ./ -o all"
read -n1 -r -p "Press any key to continue."
../PBassign.py -p ./ -o all


