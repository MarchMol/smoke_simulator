#include <GLFW/glfw3.h>
#define _USE_MATH_DEFINES 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <windows.h>
#include <time.h>
#include "state.h"
#include "data.h"
#include "smoke.h"

// Display
#define DISPLAY_WIDTH 900
#define DISPLAY_HEIGHT 600

// Simulation
#define GRID_WIDTH 200
#define GRID_HEIGHT 200

// Initial conditions

FILE *stream;

typedef struct {
    float r, g, b;
} Color;

float lerp(float a, float b, float t) {
    return a + t * (b - a);
}

Color colormap(float t) {
    Color c;
    if (t < 0.33f) {
        // Blue → Green
        float local = t / 0.33f;
        c.r = lerp(0.0f, 0.0f, local);
        c.g = lerp(0.0f, 1.0f, local);
        c.b = lerp(1.0f, 0.0f, local);
    }
    else if (t < 0.66f) {
        // Green → Yellow
        float local = (t - 0.33f) / 0.33f;
        c.r = lerp(0.0f, 1.0f, local);
        c.g = lerp(1.0f, 1.0f, local);
        c.b = lerp(0.0f, 0.0f, local);
    }
    else {
        // Yellow → Red
        float local = (t - 0.66f) / 0.34f;
        c.r = lerp(1.0f, 1.0f, local);
        c.g = lerp(1.0f, 0.0f, local);
        c.b = lerp(0.0f, 0.0f, local);
    }
    return c;
}

float rand_float(float min, float max) {
    return min + ((float)rand() / RAND_MAX) * (max - min);
}



void write_to_window(
    float **density,
    unsigned char *framebuffer
){
    double x_ratio = GRID_WIDTH / DISPLAY_WIDTH;
    double y_ratio = GRID_HEIGHT / DISPLAY_HEIGHT;
    for(int i = 0; i < DISPLAY_WIDTH; i++){
        for(int j = 0; j < DISPLAY_HEIGHT; j++){
            // Round index
            int x_orig = (int)round( i * GRID_WIDTH / DISPLAY_WIDTH);
            int y_orig = (int)round( j * GRID_HEIGHT / DISPLAY_HEIGHT);
            // Get density
            float density_val = density[x_orig][y_orig];
            if(density_val > 1.0f) density_val = 1.0f;
            if(density_val < 0.0f) density_val = 0.0f;

            // Shader (lerp)
            // Write to buffer
            int fb_index = (j * DISPLAY_WIDTH + i) * 3;
            framebuffer[fb_index + 0] = (unsigned char)(density_val* 255.0f);
            framebuffer[fb_index + 1] = (unsigned char)(density_val* 255.0f);
            framebuffer[fb_index + 2] = (unsigned char)(density_val* 255.0f);
        }
    }
}

int render(
    GLFWwindow* window,
    unsigned char *framebuffer,
    int emission_area,
    float emission_rate,
    Data data
){

    int counter = 0;
    int min_x = emission_area;
    int max_x = GRID_WIDTH - emission_area;
    int min_y = emission_area;
    int max_y = GRID_HEIGHT - emission_area;

    glfwMakeContextCurrent(window);
    while (!glfwWindowShouldClose(window)) {
        clock_t start = clock();

        // Clear buffer
        glClear(GL_COLOR_BUFFER_BIT);

       // Add source of fluid
        // if(counter % 20 == 1){
        if(counter == 0){
            for( int q = 0; q < 10; q++){
                int x_s = min_x + rand() % ( max_x - min_x);
                int y_s = min_y + rand() % ( max_y - min_y);


                double theta = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
                float dir_x = (float)cos(theta);
                float dir_y = (float)sin(theta);

                for(int a_i = -emission_area; a_i< emission_area; a_i++){
                    for(int a_j = -emission_area; a_j< emission_area; a_j++){
                          if (a_i*a_i + a_j*a_j <= emission_area*emission_area) {
                            density[x_s + a_i][y_s + a_j]+= emission_rate;
                            velocity[x_s + a_i][y_s + a_j][0] += 150.0f * dir_x;
                            velocity[x_s + a_i][y_s + a_j][1] += 150.0f* dir_y;
                          }
                    }
                }
            }
        }

        // Simulate next step
        simulation_step(pressure, pressure_buffer, density, density_buffer, b, forces, velocity, velocity_buffer, data);
        // Write to framebuffer
        write_to_window(density, framebuffer);
        
        // Framebuffer to window
        glDrawPixels(DISPLAY_WIDTH, DISPLAY_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, framebuffer);
        glfwSwapBuffers(window);
        glfwPollEvents();

        // fps
        counter++;
        clock_t end = clock();
        double elapsed_sec = (double)(end - start) / CLOCKS_PER_SEC;
        if(counter %5 == 0){
            clock_t end = clock();
            double elapsed_sec = (double)(end - start) / CLOCKS_PER_SEC;
            printf("FPS %.3f\n",1.0f/elapsed_sec);
        }
    }

    glfwTerminate();
    return 0;
}

int main() {

    // Display Init
    if (!glfwInit())
        return -1;

    GLFWwindow* window = glfwCreateWindow(
        DISPLAY_WIDTH, DISPLAY_HEIGHT, "Fluid SImulation", NULL, NULL
    );

    if (!window) {
        glfwTerminate();
        return -1;
    }
    unsigned char *framebuffer = malloc(DISPLAY_WIDTH * DISPLAY_HEIGHT * 3);


    // Simulation Init
    allocate_arrays(GRID_WIDTH, GRID_HEIGHT);
    Data data;
    data.x = GRID_WIDTH;
    data.y = GRID_HEIGHT;
    data.h = 10.0f;
    data.dt = 0.05f;
    data.jacobi_iter = 100;
    data.viscosity = 0.01f;
    data.scalar_diffusion = 0.01f;
    data.buoyancy_coeff = 1.0f;
    data.conf_strenght = 1.0f;
    int emission_area = GRID_WIDTH/10;
    
    float emission_rate =0.8f;
    render(window, framebuffer, emission_area, emission_rate, data);
    free_arrays(GRID_WIDTH, GRID_HEIGHT);

    return 0;

}
