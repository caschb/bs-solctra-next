ID=$1

mpiexec -np 4\
    ../build/main/bs_solctra\
    --length 127\
    --particles ../data/input_1000.txt\
    --job-id $ID\
    --resource ../data/resources/\
    --steps 256\
    --mode 1\
    --magnetic-profile 0\
    --num-points 100\
    --phi-angle 0\
    --dimension 2\
    --debug 1