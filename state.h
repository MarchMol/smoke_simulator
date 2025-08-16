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
extern int X, Y;

void allocate_arrays(void);
void free_arrays(void);

#endif