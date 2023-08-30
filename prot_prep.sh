#!/bin/bash

# Script to clean PDB structures and generate CG MARTINI representation with
# dssp

# Set paths
INPUT_DIR="/path/to/your/input/dir"
OUTPUT_DIR="/path/to/your/output/dir"

# Define the base name for the protein; change this for other proteins
BASE_NAME="1bvl"

# Construct file paths based on the base name
PROTEIN_PDB="${BASE_NAME}.pdb"
CLEAN_PDB="${BASE_NAME}_clean.pdb"
CG_PDB="${BASE_NAME}_CG.pdb"
CG_TOP="${BASE_NAME}.top"

# Use PyMOL to remove hetatms
pymol -c -d "load ${INPUT_DIR}/${PROTEIN_PDB}; remove hetatm; save ${INPUT_DIR}/${CLEAN_PDB}"

# Use martinize to generate the coarse-grained representation
martinize2 -f ${INPUT_DIR}/${CLEAN_PDB} -o ${OUTPUT_DIR}/${CG_TOP} -x ${OUTPUT_DIR}/${CG_PDB} -dssp /usr/local/bin/mkdssp -ff martini22

echo "Processing done!"

