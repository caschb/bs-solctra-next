//
//
// Created by lchavarr on 4/19/16.
//

#include <cassert>
#include "utils.h"
#include <cstdlib>
#include <sys/time.h>
#include <sys/stat.h>
#include <sstream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#ifndef __INTEL_COMPILER
#include <cstdlib>

void* _mmm_malloc(size_t size, size_t alignment)
{
    return aligned_alloc(alignment, size);
}
void _mmm_free(void* pointer)
{
    free(pointer);
}
#endif

void loadCartesianFile(double* x, double* y, double* z, const int length, const std::string& path)
{
    FILE* file_buff;
    //Open file
    file_buff = fopen(path.c_str(), "r");
    if (file_buff == nullptr)
    {
        printf("Error al abrir archivo \n");
   	exit(0);
    }
    else
    {
        double localX, localY, localZ;
        printf("Loading %s with length=%d\n", path.c_str(), length);
        for (int point = 0; point < length; point++)
        {
            fscanf(file_buff, "%le %le %le", &localX, &localY, &localZ);
            x[point] = localX;
            y[point] = localY;
            z[point] = localZ;
        }
        fclose(file_buff);
    }
}

void loadParticleFile(Particle* particles, const int length, const std::string& path)
{
    FILE* file_buff;
    //Open file
    file_buff = fopen(path.c_str(), "r");
    if (file_buff == nullptr)
    {
        printf("Error al abrir archivo \n");
	exit(0);
    }
    else
    {
        double localX, localY, localZ;
        printf("Loading %s with length=%d\n", path.c_str(), length);
        for (int point = 0; point < length; point++)
        {
            fscanf(file_buff, "%le %le %le", &localX, &localY, &localZ);
            particles[point].x = localX;
            particles[point].y = localY;
            particles[point].z = localZ;
        }
        fclose(file_buff);
    }
}


double randomGenerator(const double min, const double max,const int seedValue){
    //std::random_device rd; //Used to obtain a seed for the random number
    std::mt19937 gen(seedValue); 
    std::uniform_real_distribution<double> distribution(min, max);
    double result = distribution(gen);
    return result;    
}


/*double randomGeneratorExponential(const double maxNumber){
    std::random_device rd; //Used to obtain a seed for the random number
    std::mt19937 gen(rd()); 
    std::exponential_distribution<double> distribution(1.0);
    double result = distribution(gen);
    if(result>1.0){
        double factor = gen.max()/maxNumber;
        result = result/factor;   
    }else{ result = result*maxNumber;}
    return result;    
}*/


void initializeParticles(Particle * particles, const int length, const int seedValue){
    double radius;
    double minorRadius;
    double toroidalAngle;
    double poloidalAngle;
    double toroidalRadians;
    double poloidalRadians;
    printf("Initializing with seed %d\n",seedValue);
    for (int point = 0; point < length; point++)
    {
        radius = randomGenerator(0.0,0.4,seedValue);
        toroidalAngle = randomGenerator(0.0,360.0,seedValue);
        poloidalAngle = randomGenerator(0.0,360.0,seedValue);
        
        toroidalRadians = toroidalAngle*PI/180.0;
        poloidalRadians = poloidalAngle*PI/180.0;
        minorRadius = MINOR_RADIUS*radius;
        particles[point].x = (MAJOR_RADIUS+(minorRadius)*cos(poloidalRadians))*cos(toroidalRadians);
        particles[point].y = (MAJOR_RADIUS+(minorRadius)*cos(poloidalRadians))*sin(toroidalRadians);
        particles[point].z = (minorRadius)*sin(poloidalRadians);
    } 
    
    
}




double getCurrentTime()
{
    struct timeval tod;
    gettimeofday(&tod, nullptr);
    return static_cast<double>(tod.tv_sec) + static_cast<double>(tod.tv_usec) * 1.0e-6;
}

void createDirectoryIfNotExists(const std::string& path)
{
    if(!directoryExists(path))
    {
        const int error = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (0 > error)
        {
            printf("Error creating directory!n");
            exit(1);
        }
    }
}

bool directoryExists(const std::string& path)
{
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
    {
        return false;
    }
    else if (info.st_mode & S_IFDIR)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::string getZeroPadded(const int num)
{
    std::ostringstream convert;
    convert << num;
    std::string value = convert.str();
    const size_t numSize = value.size();
    int numberZeros = 5; //Number of 0s wanted in file name: particle00000...X.txt
    for(size_t i = 0 ; i < numberZeros - numSize ; ++i)
    {
        value = "0" + value;
    }
    return value;
}

void setupParticleType(){
    const int nitems = 3;
    int blocklenghts[3]={1,1,1};
    MPI_Datatype memberTypes[3] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    MPI_Aint offsets[3];

    offsets[0] = offsetof(struct Particle, x);
    offsets[1] = offsetof(struct Particle, y);
    offsets[2] = offsetof(struct Particle, z);

    MPI_Type_create_struct(nitems, blocklenghts,offsets,memberTypes,&mpi_particle);
    MPI_Type_commit(&mpi_particle);
}
