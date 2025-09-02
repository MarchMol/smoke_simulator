#include <stdlib.h>
#include "smoke.h"
#include "data.h"
#include "state.h"

// Variables globales
float **pressure;
float **pressure_buffer;

float **density;
float **density_buffer;

float **b;
float ***forces;

float ***velocity;
float ***velocity_buffer;

// Bloques contiguos de memoria
static float *pressure_mem, *pressure_buffer_mem;
static float *density_mem, *density_buffer_mem;
static float *b_mem;
static float *forces_mem;
static float *velocity_mem, *velocity_buffer_mem;

void allocate_arrays(int X, int Y){
    // ====== VELOCITY (X x Y x 2) ======
    velocity = malloc(X * sizeof(float**));
    velocity_buffer = malloc(X * sizeof(float**));

    velocity_mem = malloc(X * Y * 2 * sizeof(float));
    velocity_buffer_mem = malloc(X * Y * 2 * sizeof(float));

    for(int i = 0; i < X; i++){
        velocity[i] = malloc(Y * sizeof(float*));
        velocity_buffer[i] = malloc(Y * sizeof(float*));
        for(int j = 0; j < Y; j++){
            velocity[i][j] = &velocity_mem[(i*Y + j) * 2];
            velocity_buffer[i][j] = &velocity_buffer_mem[(i*Y + j) * 2];
        }
    }

    // ====== FORCES (X x Y x 2) ======
    forces = malloc(X * sizeof(float**));
    forces_mem = malloc(X * Y * 2 * sizeof(float));

    for(int i = 0; i < X; i++){
        forces[i] = malloc(Y * sizeof(float*));
        for(int j = 0; j < Y; j++){
            forces[i][j] = &forces_mem[(i*Y + j) * 2];
        }
    }

    // ====== PRESSURE (X x Y) ======
    pressure = malloc(X * sizeof(float*));
    pressure_buffer = malloc(X * sizeof(float*));
    pressure_mem = malloc(X * Y * sizeof(float));
    pressure_buffer_mem = malloc(X * Y * sizeof(float));

    for(int i = 0; i < X; i++){
        pressure[i] = &pressure_mem[i*Y];
        pressure_buffer[i] = &pressure_buffer_mem[i*Y];
    }

    // ====== DENSITY (X x Y) ======
    density = malloc(X * sizeof(float*));
    density_buffer = malloc(X * sizeof(float*));
    density_mem = malloc(X * Y * sizeof(float));
    density_buffer_mem = malloc(X * Y * sizeof(float));

    for(int i = 0; i < X; i++){
        density[i] = &density_mem[i*Y];
        density_buffer[i] = &density_buffer_mem[i*Y];
    }

    // ====== B (X x Y) ======
    b = malloc(X * sizeof(float*));
    b_mem = malloc(X * Y * sizeof(float));

    for(int i = 0; i < X; i++){
        b[i] = &b_mem[i*Y];
    }

    // ====== Inicializar en 0 ======
    for (int i = 0; i < X*Y; i++){
        pressure_mem[i] = 0.0f;
        pressure_buffer_mem[i] = 0.0f;
        density_mem[i] = 0.0f;
        density_buffer_mem[i] = 0.0f;
        b_mem[i] = 0.0f;
    }
    for (int i = 0; i < X*Y*2; i++){
        velocity_mem[i] = 0.0f;
        velocity_buffer_mem[i] = 0.0f;
        forces_mem[i] = 0.0f;
    }
}

void free_arrays(int X, int Y){
    // Liberar punteros secundarios
    for(int i = 0; i < X; i++){
        free(velocity[i]);
        free(velocity_buffer[i]);
        free(forces[i]);
    }
    free(velocity);
    free(velocity_buffer);
    free(forces);

    free(pressure);
    free(pressure_buffer);
    free(density);
    free(density_buffer);
    free(b);

    // Liberar bloques de memoria contigua
    free(pressure_mem);
    free(pressure_buffer_mem);
    free(density_mem);
    free(density_buffer_mem);
    free(b_mem);
    free(velocity_mem);
    free(velocity_buffer_mem);
    free(forces_mem);
}
