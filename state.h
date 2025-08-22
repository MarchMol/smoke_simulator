#ifndef STATE_H
#define STATE_H

extern float **pressure;
extern float **pressure_buffer;

extern float **density;
extern float **density_buffer;

extern float **b;
extern float ***forces;

extern float ***velocity;
extern float ***velocity_buffer;

void allocate_arrays(int X, int Y);
void free_arrays(int X, int Y);

#endif