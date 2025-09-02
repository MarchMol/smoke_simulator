#ifndef SMOKE_H
#define SMOKE_H
#include "data.h"

#ifdef _OPENMP
#include <omp.h>
static inline int choose_threads(const Data *data) {
    long long cells = (long long)data->x * (long long)data->y;
    int max_t = omp_get_max_threads();
    if (cells < 20000)   return 1;
    if (cells < 40000)   return max_t > 2 ? 2 : 1;
    if (cells < 200000)  return max_t > 4 ? 4 : max_t;
    return max_t;
}
#else
// En versión secuencial, siempre solo se llama a un hilo
static inline int choose_threads(const Data *data) {
    return 1;
}
#endif

void simulation_step(
    float **pressure,
    float **pressure_buffer,
    float **density,
    float **density_buffer,
    float **b,
    float ***forces,
    float ***velocity,
    float ***velocity_buffer,
    Data *data);

#endif
