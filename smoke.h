#ifndef SMOKE_H
#define SMOKE_H
#include "data.h"
void simulation_step(
    float **pressure,
    float **pressure_buffer,

    float **density,
    float **density_buffer,

    float **b,
    float ***forces,

    float ***velocity,
    float ***velocity_buffer,

    Data data);


#endif