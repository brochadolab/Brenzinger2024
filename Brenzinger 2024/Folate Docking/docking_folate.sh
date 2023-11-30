#!bin/bash

#required:
## 2D sdf files from PubChem stored in a "ligands/" folder
## PDB file of the protein ; generated from Alphafold or retrieved from PDB
## config file for the starting orientation and location of the ligand in the active site

touch __logger

for file in ligands/*.sdf; do

    echo "picking and starting $file"

    #check availability of the ligand file
    [ -e "$file" ] || continue

    #if available, convert to 3D pdb file
    obabel ${file} -O ${file}.pdb --gen3D
    
    #convert the pdb to pdbqt file
    obabel ${file}.pdb -O ${file}.pdbqt
    
    #perform molecular docking with autodock vina
    vina --ligand ${file}.pdbqt --config config.txt --log ${file}.log
    
    grep 'Output\|0.000' ${file}.log >> __logger


    echo "done and dusted with $file"

    
done

cat __logger
