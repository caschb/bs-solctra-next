//
// utils.h
// Created by Diego Jimenez
// Basic simulation structures and function prototypes declarations
//

#ifndef SOLCTRA_UTILS_H
#define SOLCTRA_UTILS_H


#define PI      3.141592654
#define miu     1.2566e-06
#define I       -4350
#define MINOR_RADIUS 0.0944165
#define MAJOR_RADIUS 0.2381
#define ALIGNMENT_SIZE 64
#ifdef KNL
#define GRADES_PER_PAGE ALIGNMENT_SIZE  * KNL / sizeof(double)
#else
#define GRADES_PER_PAGE ALIGNMENT_SIZE / sizeof(double)
#endif
#define TOTAL_OF_GRADES 360
#define TOTAL_OF_GRADES_PADDED 384
#define TOTAL_OF_COILS 12

#include <sstream>
#include <cstdio>
#include <string>
#include "mpi.h"

extern MPI_Datatype mpi_particle;

struct cartesian
{
    double x, y, z;
    void print()
    {
        printf("X=[%.17g]. Y=[%.17g]. Z=[%.17g].\n", x, y, z);
    }
};


struct Particle
{
    double x, y, z;
};


struct Coil
{
    double* x;
    double* y;
    double* z;

};

struct GlobalData
{
    Coil coils;
    Coil e_roof;
    double* leng_segment;
};

#ifndef __INTEL_COMPILER

void* _mmm_malloc(size_t size, size_t alignment);
void _mmm_free(void* pointer);

#define nullptr NULL

#endif

void setupParticleType();

void loadCartesianFile(double* x, double* y, double* z, const int length, const std::string& path);
void loadParticleFile(Particle* particles, const int length, const std::string& path);
double getCurrentTime();
void createDirectoryIfNotExists(const std::string& path);
bool directoryExists(const std::string& path);
std::string getZeroPadded(const int num);
double randomGenerator(const double min, const double max, const int seedValue);
//double randomGeneratorExponential(const double lambda);
void initializeParticles(Particle * particles, const int length, const int seedValue);


#endif //SOLCTRA_UTILS_H
