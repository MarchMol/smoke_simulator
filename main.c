#include <GLFW/glfw3.h>
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
#define WIDTH 200
#define HEIGHT 200

void set_pixel(unsigned char *fb, int width, int x, int y,
               float **density) {
    float b = density[x][y];
    b = powf(b, 0.5f);
    if (b < 0.0f) b = 0.0f;
    if (b > 1.0f) b = 1.0f;
    unsigned char value = (unsigned char)(b * 255.0f);
       int index = (y * width + x) * 3;
    fb[index + 0] = value; // R
    fb[index + 1] = value; // G
    fb[index + 2] = value; // B
}
float rand_float(float min, float max) {
    return min + ((float)rand() / RAND_MAX) * (max - min);
}


int main() {
    if (!glfwInit())
        return -1;

    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Smoke sim", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Set up
    unsigned char *framebuffer = malloc(WIDTH * HEIGHT * 3);
    allocate_arrays();
    Data data;
    data.x = WIDTH;
    data.y = HEIGHT;
    data.h = 1.0f;
    data.dt = 0.001;
    data.jacobi_iter = 40;
    data.viscosity = 5.0f;
    data.scalar_diffusion =100;
    data.buoyancy_coeff = 1.0f;
    int source_x = HEIGHT/2;
    int source_y = WIDTH/2;
    int emission_area = 70;
    float emission_rate = 0.001f;


    srand((unsigned int)time(NULL)); // seed once at start

    glfwMakeContextCurrent(window);
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        for(int a_i = -emission_area; a_i< emission_area; a_i++){
            for(int a_j = -emission_area; a_j< emission_area; a_j++){
                density[source_x + a_i][source_y + a_j]+=emission_rate;
                if(density[source_x + a_i][source_y + a_j] > 1.0f) density[source_x][source_y] = 1.0f;
                // velocity[source_x + a_i][source_y + a_j][0] += 0.1f;
                // velocity[source_x + a_i][source_y + a_j][1] += 0.1f;
            }
        }
        for(int i = 0; i<WIDTH;i++){
            for(int j = 0; j<HEIGHT;j++){
                float factor = rand_float(-0.1f, 0.1f);
                velocity[i][j][0] += factor;
                velocity[i][j][1] += factor;
            }
        }
        
        simulation_step(pressure, pressure_buffer, density, density_buffer, b, buoyancy, velocity, velocity_buffer, data);

        for(int i = 0; i<HEIGHT; i++){
            for(int j = 0; j<WIDTH; j++){
                set_pixel(framebuffer, WIDTH, i, j, density);
            }
        }        
        
        glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, framebuffer);
        glfwSwapBuffers(window);
        glfwPollEvents();
        Sleep(500);
        printf("Done!\n");
    }

    glfwTerminate();
    free_arrays();

    return 0;

}
