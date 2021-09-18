#!/bin/sh
#SBATCH --ntasks-per-node=40 --nodes=1
#SBATCH -t 8:00:00
#SBATCH -A gsienkf
##SBATCH -q debug
#SBATCH -J generate_fv3mesh
#SBATCH -o generate_fv3mesh.out
#SBATCH -e generate_fv3mesh.err
source ~/bin/condapy
python generate_fv3mesh.py 48 &
python generate_fv3mesh.py 96 &
python generate_fv3mesh.py 192 &
python generate_fv3mesh.py 384 &
python generate_fv3mesh.py 768 &
python generate_fv3mesh.py 1152 &
python generate_fv3mesh.py 3072 &
wait
