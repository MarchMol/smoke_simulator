#include <stdlib.h>
#include "smoke.h"
#include "data.h"
#include "state.h"

int X = 300;
int Y = 300;

float **pressure;
float **pressure_buffer;

float **density;
float **density_buffer;

float **b;
float ***forces;

float ***velocity;
float ***velocity_buffer;

void allocate_arrays(){
    // Velocity allocation
    velocity = malloc(X * sizeof(float**));
    velocity_buffer = malloc(X * sizeof(float**));
    forces = malloc(X * sizeof(float**));
    for(int i = 0; i<X; i++){
        velocity[i] = malloc(Y * sizeof(float*));
        velocity_buffer[i] = malloc(Y * sizeof(float*));
        forces[i] = malloc(Y * sizeof(float*));
        for(int j= 0; j<Y; j++){
            velocity[i][j] = malloc(2 * sizeof(float));
            velocity_buffer[i][j] = malloc(2 * sizeof(float));
            forces[i][j] = malloc(2 * sizeof(float));
        }
    }
    // Init
    for (int i = 0; i < X; i++) {
        for (int j = 0; j < Y; j++) {
            forces[i][j][0] = 0.0f;
            forces[i][j][1] = 0.0f;
            velocity[i][j][0] = 0.0f;
            velocity[i][j][1] = 0.0f;
            velocity_buffer[i][j][0] = 0.0f;
            velocity_buffer[i][j][1] = 0.0f;
        }
    }

    // pressure allocation
    pressure = malloc(X * sizeof(float*));
    pressure_buffer = malloc(X * sizeof(float*));
    density = malloc(X * sizeof(float*));
    density_buffer = malloc(X * sizeof(float*));
    b = malloc(X * sizeof(float*));
    for(int i = 0; i<X; i++){
        pressure[i] = malloc(Y * sizeof(float));
        pressure_buffer[i] = malloc(Y * sizeof(float));
        density[i] = malloc(Y * sizeof(float));
        density_buffer[i] = malloc(Y * sizeof(float));
        b[i] = malloc(Y * sizeof(float));
    }

    // Init
        for (int i = 0; i < X; i++) {
        for (int j = 0; j < Y; j++) {
            pressure[i][j] = 0.0f;
            pressure_buffer[i][j] = 0.0f;
            density[i][j] = 0.0f;
            density_buffer[i][j] = 0.0f;
            b[i][j] = 0.0f;
        }
    }
}

void free_arrays(){
    // Free velocity
    for(int i = 0; i < X; i++) {
        for(int j = 0; j < Y; j++) {
            free(forces[i][j]);
            free(velocity[i][j]);
            free(velocity_buffer[i][j]);
        }
        free(forces[i]);
        free(velocity[i]);
        free(velocity_buffer[i]);
    }
    free(forces);
    free(velocity);
    free(velocity_buffer);

    // Free other
    for(int i = 0; i < X; i++){
        free(pressure[i]);
        free(pressure_buffer[i]);
        free(density[i]);
        free(density_buffer[i]);
        free(b[i]);
    }
    free(pressure);
    free(pressure_buffer);
    free(density);
    free(density_buffer);
    free(b);
}