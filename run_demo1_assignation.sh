#! /bin/bash

# 2012 - P. Poulain, A. G. de Brevern

# exit script at first error
set -e

function pause(){
    read -r -s -n1 -p "Press any key to continue."
    echo
}


echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  demo script for PBxplore                                              |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

pause

DATA_PATH=$(PBdata)
DEMO_PATH=demo1_assignation_tmp
INPUT_FILES=(3ICH.pdb
             1AY7.pdb
             2LFU.pdb)
mkdir -p $DEMO_PATH
for input_file in ${INPUT_FILES[@]}
do
    cp ${DATA_PATH}/${input_file} $DEMO_PATH
done

cd $DEMO_PATH

echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignment: static structures                          |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "Test with 3ICH.pdb: RX structure, one chain"
echo "PBassign -p 3ICH.pdb -o 3ICH"
pause
PBassign -p 3ICH.pdb -o 3ICH


echo  -e "\n"
echo "Test with 1AY7.pdb: RX structure, complex with two chains"
echo "PBassign -p 1AY7.pdb -o 1AY7"
pause
PBassign -p 1AY7.pdb -o 1AY7


echo  -e "\n"
echo "Test with 2LFU.pdb: RMN structure, 10 models"
echo "PBassign -p 2LFU.pdb -o 2LFU"
pause
PBassign -p 2LFU.pdb -o 2LFU


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignment: static structures --phipsi option          |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"


echo  -e "\n"
echo "Test with 3ICH.pdb: RX structure, one chain"
echo "PBassign -p 3ICH.pdb -o 3ICH --phipsi"
pause
PBassign -p 3ICH.pdb -o 3ICH --phipsi


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignment: static structures --flat option            |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"


echo  -e "\n"
echo "Test with 2LFU.pdb: RMN structure, 10 models"
echo "PBassign -p 2LFU.pdb -o 2LFU --flat"
pause
PBassign -p 2LFU.pdb -o 2LFU --flat


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Protein Blocks assignment: several static structures                  |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"

echo  -e "\n"
echo "Test with several PDB files"
echo "PBassign -p 3ICH.pdb -p 2LFU.pdb -o several"
pause
PBassign -p 3ICH.pdb -p 2LFU.pdb -o several

echo  -e "\n"
echo "Test with all PDB files from a directory"
echo "PBassign -p ./ -o all"
pause
PBassign -p ./ -o all


echo  -e "\n"
echo "#------------------------------------------------------------------------#"
echo "|                                                                        |"
echo "|  Demo completed!                                                       |"
echo "|                                                                        |"
echo "#------------------------------------------------------------------------#"
echo
echo "Look at *.PB.* files in the demo1_assignation_tmp directory."
pwd
ls -lh
echo "Do not forget to delete demo1_assignation_tmp directory when you will be done with this demo."

