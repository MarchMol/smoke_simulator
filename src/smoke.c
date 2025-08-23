#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <math.h>
#include "smoke.h"
#include "data.h"
#include "state.h"

// ---------------- Memory ----------------------- //

void copy_vec(float ***target, float ***origin, Data *data){
    int x = data->x;
    int y = data->y;
    for(int i = 0; i<x; i++){
        for(int j = 0; j<y; j++){
            target[i][j][0] = origin[i][j][0];
            target[i][j][1] = origin[i][j][1];
        }
    }
}

void copy_scalar(float **target, float **origin, Data *data){
    int x = data->x;
    int y = data->y;
    for(int i = 0; i<x; i++){
        for(int j = 0; j<y; j++){
            target[i][j] = origin[i][j];
        }
    }
}

void print_all(
    int stage,
    float ***force,
    float ***velocity,
    float **density,
    float **pressure,
    Data *data
){
    int x = data->x;
    int y = data->y;

    printf("At stage = %d\n", stage);
    printf("\tForce\n");
    for(int i = 0; i<x; i++){
        for(int j = 0; j<y; j++){
            printf("(%.2f, %.2f) ",force[i][j][0],force[i][j][1]);
        }
        printf("\n");
    }
    printf("\tVelocity\n");
    for(int i = 0; i<x; i++){
        for(int j = 0; j<y; j++){
            printf("(%.2f, %.2f) ",velocity[i][j][0],velocity[i][j][1]);
        }
        printf("\n");
    }
    printf("\tDensity\n");
    for(int i = 0; i<x; i++){
        for(int j = 0; j<y; j++){
            printf("%.2f ",density[i][j]);
        }
        printf("\n");
    }
    printf("\tPressure\n");
    for(int i = 0; i<x; i++){
        for(int j = 0; j<y; j++){
            printf("%.2f ",pressure[i][j]);
        }
        printf("\n");
    }
}
// --------------- Simulation ----------------- //

void update_forces(
    float ***velocity,
    float ***forces,
    float ***buffer_vec,
    float **buffer_scalar,
    Data *data
){
    int x = data->x;
    int y = data->y;
    float h = data->h;
    float w = 0.0f;
    float dt = data->dt;
    // 1. Calculte Vorticity index
    for(int i = 1; i<x-1; i++){
        for(int j = 1; j<y-1; j++){
            //      w = (∂U / ∂x)  ​− (∂U / ∂y)​
            //      w = ∇U / 2*h
            w = 
            (velocity[i+1][j][0] - velocity[i-1][j][0])/(2*h)
             -
            (velocity[i][j+1][1] - velocity[i][j-1][1])/(2*h);

            // |w| = abs(w)
            buffer_scalar[i][j] = fabsf(w);
            // (active region check, to avoid doing later)
        }
    }
    // ∇|w|
    float wx = 0.0f;
    float wy = 0.0f;
    float mag_w = 0.0f;

    float fx = 0.0f;
    float fy = 0.0f;

    // 2. Calculate ∇|w| and resulting force
    int alpha = data->buoyancy_coeff;
    float eps_h = data->conf_strenght * h;
    for(int i = 1; i<x-1; i++){
        for(int j = 1; j<y-1; j++){
            // Calc ∇|w|
            wx = (
                buffer_scalar[i+1][j] - buffer_scalar[i-1][j]
            )/(2*h);
            wy = (
                buffer_scalar[i][j+1] - buffer_scalar[i][j-1]
            )/(2*h);
            // Normalization
            mag_w = sqrtf(wx*wx + wy*wy);
            if(mag_w <= 0) mag_w = 1e-12;

            // force: F = eps * h (N x w)
            fx = eps_h* wx/mag_w * buffer_scalar[i][j];
            fy = -(eps_h * wy/mag_w * buffer_scalar[i][j]) +  alpha*density[i][j];

            forces[i][j][0] = fx;
            // fy = - eps * h (Nx X w)
            forces[i][j][1] = fy;

            // 3. Update Velocities
            velocity[i][j][0]+=dt*fx;
            velocity[i][j][1]+=dt*fy;
        }
    }
}

void density_steps(
    float ***velocity,
    float **density,
    float **density_buffer,
    Data *data
){
    // Saved referenced data
    int x = data->x;
    int y = data->y;
    float dt = data->dt;
    float h = data->h;
    float scalar_diffusion = data->scalar_diffusion;
    float jacobi_iter = data->jacobi_iter;

    // 1. Advection
    
    copy_scalar(density_buffer, density, data);
    for(int i = 1; i<=x-1; i++){
        for(int j = 1; j<=y-1; j++){
            float x_prev = i - (dt * velocity[i][j][0])/h;
            float y_prev = j - (dt * velocity[i][j][1])/h;

            // Clamp
            if (x_prev < 0.5f) x_prev = 0.5f;
            if (x_prev > x - 1.5f) x_prev = x - 1.5f;
            if (y_prev < 0.5f) y_prev = 0.5f;
            if (y_prev > y - 1.5f) y_prev = y - 1.5f;

            // Saving for later
            int i0 = (int)x_prev;
            int j0 = (int)y_prev;
            float sx = x_prev - i0;
            float sy = y_prev - j0;

            // Index check
            if(i0 < 0) i0 = 0;
            if(j0 < 0) j0 = 0;
            if(i0 > x-2) i0 = x-2;
            if(j0 > y-2) j0 = y-2;
            // Bilinear interpolation
            float density_interp = 
                (1-sx) * (1-sy) * density_buffer[i0  ][j0  ] + 
                sx     * (1-sy) * density_buffer[i0+1][j0  ] + 
                (1-sx) * sy     * density_buffer[i0  ][j0+1] + 
                sx     * sy     * density_buffer[i0+1][j0+1]
            ;
            // Saving into complete array
            density[i][j] = density_interp;
        }  
    }

    // 3. Diffuse
    copy_scalar(density_buffer, density, data);
    float alpha = (scalar_diffusion*dt)/(h*h);
    for(int jac =0; jac<jacobi_iter;jac++){
        for(int i = 1; i<x-1; i++){
            for(int j = 1; j<y-1; j++){
                density[i][j] = 
                    (density_buffer[i][j] + alpha *
                    (density_buffer[i-1][j]+density_buffer[i+1][j]+density_buffer[i][j-1]+density_buffer[i][j+1]))/
                    (1 + 4*alpha)
                ;
            }
        }
        for(int i = 0; i < x; i++){
            density[i][0] = density[i][1];
            density[i][y-1] = density[i][y-2];
        }
        for(int j = 0; j < y; j++){
            density[0][j] = density[1][j];
            density[x-1][j] = density[x-2][j];
        }
        copy_scalar(density_buffer, density, data);
    }
    
}


void velocity_steps(
    float ***velocity,
    float ***velocity_buffer,
    Data *data
){
    // printf("Entered!\n");
    // Dimention save
    int x = data->x;
    int y = data->y;

    // Constant save
    float dt = data->dt;
    float h = data->h;
    float viscosity = data->viscosity;
    float jacobi_iter = data->jacobi_iter;

    // 1. Forces
    copy_vec(velocity_buffer, velocity, data);
    // printf("completed U^f\n");

    // 2. Advection
    //      velocity is now U^a
    for(int i = 0; i<x; i++){
        for(int j = 0; j<y; j++){
            float x_prev = i - dt * velocity_buffer[i][j][0]/ h;
            float y_prev = j - dt * velocity_buffer[i][j][1]/h;
            // Clamp
            if (x_prev < 0.5f) x_prev = 0.5f;
            if (x_prev > x - 1.5f) x_prev = x - 1.5f;
            if (y_prev < 0.5f) y_prev = 0.5f;
            if (y_prev > y - 1.5f) y_prev = y - 1.5f;

            // Saving for later
            int i0 = (int)x_prev; // 2
            int j0 = (int)y_prev; // 2
            float sx = x_prev - i0;  // 0
            float sy = y_prev - j0; // 0.39

            if(i0 < 0) i0 = 0;
            if(j0 < 0) j0 = 0;
            if(i0 > x-2) i0 = x-2;
            if(j0 > y-2) j0 = y-2;

            // Bilinear interpolation
            float hor_interp = 
                (1-sx)*(1-sy) * velocity_buffer[i0][j0][0] +
                sx*(1-sy) * velocity_buffer[i0+1][j0][0] +
                (1-sx)*sy * velocity_buffer[i0][j0+1][0] +
                sx*sy * velocity_buffer[i0+1][j0+1][0]
            ;
            float vers_interp = 
                (1-sx)*(1-sy) * velocity_buffer[i0][j0][1] +
                sx*(1-sy) * velocity_buffer[i0+1][j0][1] +
                (1-sx)*sy * velocity_buffer[i0][j0+1][1] +
                sx*sy * velocity_buffer[i0+1][j0+1][1]
            ;
            // Saving into complete array
            velocity[i][j][0] = hor_interp;
            velocity[i][j][1] = vers_interp;
            
        }  
    }
    
    
    // printf("completed U^a\n");
    // 3. Diffuse | Jacobi iterate
    //      buffer is now prev, and velocity is U^*
    copy_vec(velocity_buffer, velocity, data);
    float alpha = (viscosity*dt)/(h*h);
    for(int jac = 0; jac<jacobi_iter; jac++){
        for(int i = 1; i<x-1; i++){
            for(int j = 1; j<y-1; j++){
                velocity[i][j][0] = 
                    (velocity_buffer[i][j][0] + alpha *
                    (velocity_buffer[i-1][j][0] + velocity_buffer[i+1][j][0] + velocity_buffer[i][j-1][0] + velocity_buffer[i][j+1][0]))/
                    (1 + 4*alpha)
                ;
                velocity[i][j][1] = 
                    (velocity_buffer[i][j][1] + alpha *
                    (velocity_buffer[i-1][j][1]+velocity_buffer[i+1][j][1] + velocity_buffer[i][j-1][1]+velocity_buffer[i][j+1][1]))/
                    (1 + 4*alpha)
                ;
            }
        }
        // Setting boundary conditions to 0
        for (int i = 1; i < x-1; i++) {
            velocity[i][0][0] = velocity[i][1][0];
            velocity[i][0][1] = velocity[i][1][1];
            velocity[i][y-1][0] = velocity[i][y-2][0];
            velocity[i][y-1][1] = velocity[i][y-2][1];
        }
        for (int j = 1; j < y-1; j++) {
            velocity[0][j][0] = velocity[1][j][0];
            velocity[0][j][1] = velocity[1][j][1];
            velocity[x-1][j][0] = velocity[x-2][j][0];
            velocity[x-1][j][1] = velocity[x-2][j][1];
        }
        copy_vec(velocity_buffer, velocity, data);
    }
}

void pressure_projection(
    float **pressure,
    float **pressure_buffer,

    float **b,
    float **density,    // rou_star
    float ***velocity,  // u_star
    Data *data
){
    // Save referenced data
    int x = data->x;
    int y = data->y;
    float dt = data->dt;
    float h = data->h;
    float jacobi_iter = data->jacobi_iter;
    // Calculate b (solution)
    for(int i = 1; i<x-1;i++){
        for(int j = 1; j<y-1; j++){
            // b = (ρ/Δt) * ∇u
            b[i][j] = 
            (density[i][j]/dt) * // (ρ/Δt)
            ((velocity[i+1][j][0] - velocity[i-1][j][0])/(2*h) + 
            (velocity[i][j+1][1] - velocity[i][j-1][1])/(2*h))
            ;
        }  
    }

    // Poisson pressure equations
    copy_scalar(pressure_buffer, pressure, data);
    for(int jac =0;jac<jacobi_iter;jac++){
        for(int i = 1; i<x-1; i++){
            for(int j = 1; j<y-1; j++){
                // P^n+1 = (1/4)*
                pressure[i][j] = 
                    (
                    pressure_buffer[i-1][j]+pressure_buffer[i+1][j]+
                    pressure_buffer[i][j-1] + pressure_buffer[i][j+1] -
                    (h*h*b[i][j]))/(4.0f)
                ;
            }
        }
        for(int i = 0; i < x; i++) {
            pressure[i][0] = 0;
            pressure[i][y-1] = 0;
        }
        for(int j = 0; j < y; j++) {
            pressure[0][j]= 0;
            pressure[x-1][j]= 0;
        }
        copy_scalar(pressure_buffer, pressure, data);
    }
}


void correct_velocity(
    float ***velocity,
    float **density,
    float **pressure,
    Data *data
){
    // Save referenced data
    int x = data->x;
    int y = data->y;
    float dt = data->dt;
    float h = data->h;
    for(int i = 1; i<x-1; i++){
        for(int j = 1; j<y-1; j++){
            velocity[i][j][0] -= 
                (dt/(2.0f*h)) *
                (pressure[i+1][j] - pressure[i-1][j])
            ;

            velocity[i][j][1] -=
                (dt/(2*h)) *
                (pressure[i][j+1] - pressure[i][j-1])
            ;

        }
    }
    
}

// -------------- Main function -----------------//
void simulation_step(
    float **pressure,
    float **pressure_buffer,

    float **density,
    float **density_buffer,

    float **b,
    float ***forces,

    float ***velocity,
    float ***velocity_buffer,

    Data *data

){
    update_forces(
        velocity,
        forces,
        velocity_buffer,
        density_buffer,
        data
    );
    density_steps(
        velocity,
        density,
        density_buffer,
        data
    );
    velocity_steps(
        velocity,
        velocity_buffer,
        data
    );
    pressure_projection(
        pressure,
        pressure_buffer,
        b,
        density,
        velocity,
        data
    );
    correct_velocity(
        velocity,
        density,
        pressure,
        data
    );
}