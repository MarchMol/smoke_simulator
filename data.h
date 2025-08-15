#ifndef DATA_H
#define DATA_H

typedef struct{
    // Grid parameters
    float x;
    float y;
    float h;

    // Sim parameters
    float dt;
    int jacobi_iter;
    
    // Smoke Constants
    float viscosity;
    float scalar_diffusion;
    float fluid_density;

    // Other constants
    float buoyancy_coeff;

} Data;

#endif