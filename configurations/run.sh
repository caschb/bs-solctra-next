ID=$1

mpiexec -np 4\
    ../build/main/bs_solctra\
    -length 1024\
    -particles ../data/input_1000.txt\
    -id $ID\
    -resource ../data/resources/\
    -steps 100\
    -mode 1\
    -magnetic_prof 0 100 0 2\
    -print_type 1