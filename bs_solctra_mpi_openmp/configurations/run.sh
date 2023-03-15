ID=100

mpiexec -np 4\
    /home/casch/Devel/bs-solctra/bs-solctra-implementations/bs_solctra_mpi_openmp/build/main/bs_solctra\
    -length 1024\
    -particles ../data/input_1000.txt\
    -id $ID\
    -resource ../data/resources/\
    -steps 10\
    -mode 1\
    -magnetic_prof 0 100 0 2\
    -print_type 1