#!/bin/sh
#SBATCH --job-name=bs-solctra-multinode
#SBATCH --partition=nu-wide
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=1:00:00
#SBATCH --chdir /home/djimenez/bs-solctra-vec-omp-mpi/resultados

export ID=${SLURM_JOB_ID}



cd $SLURM_SUBMIT_DIR
rm -Rf results_${ID}


##############Load needed modules##################################


#These modules are needed just for testing and HPC measuring
#module load valgrind/3.12.0-mpich-3.2
#source /home/djimenez/intel/advisor_2019/advixe-vars.sh
#source /opt/intel/vtune_amplifier_2019.1.0.579888/amplxe-vars.sh


#These modules are essential for BS-SOLCTRA execution
module load gcc/7.2.0
module load mpich/3.2.1-gcc-7.2.0 


####################################################################




#########Set Environment Variables for parallel execution###########

export OMP_NUM_THREADS=64
export OMP_SCHEDULE=dynamic
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

###################################################################



###############Execution commands####################################

#valgrind --tool=memcheck --leak-check=yes ./bs-solctra-debug -length 128 -particles /home/djimenez/bs-solctra/input.txt -id $ID -resource /home/djimenez/bs-solctra/resources/ -mode 1 -magnetic_prof 1 100 0 -print_type 0

#advixe-cl --collect=survey -project-dir=/home/djimenez/bs-solctra/solctra-dev/advisor_results_omp -- ./bs-solctra-vec-omp -length 128 -particles /home/djimenez/bs-solctra/input.txt -id $ID -resource /home/djimenez/bs-solctra/resources/ -mode 1 -magnetic_prof 1 100 0 2 -print_type 0

#amplxe-cl -collect uarch-exploration -result-dir /home/djimenez/bs-solctra/solctra-dev/vtune_hw1024_256 -knob pmu-collection-mode=detailed ./bs-solctra-vec-omp -length 1024 -particles /home/djimenez/bs-solctra/solctra-dev/input_1000.txt -id $ID -resource /home/djimenez/bs-solctra/resources/ -mode 1 -magnetic_prof 0 100 0 2 -print_type 1
 
#mpiexec -n 1 ./bs-solctra-multinode -length 9998 -particles ../input_36k.txt -id $ID -resource ../resources/ -mode 1 -magnetic_prof 0 100 0 2 -print_type 1

#mpiexec -n 8 -bind-to core:64 ./bs-solctra-multinode -length 1024 -particles ../input_1000.txt -id $ID -resource ../resources/ -mode 1 -magnetic_prof 0 100 0 2 -print_type 1

mpiexec -n 8 ./resultados/bs-solctra-multinode -length 1024 -particles input_1000.txt -id $ID -resource resources/ -mode 1 -magnetic_prof 0 100 0 2 -print_type 1
######################################################################

