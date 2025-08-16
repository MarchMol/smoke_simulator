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
#define WIDTH 300
#define HEIGHT 300
FILE *stream;

typedef struct {
    float r, g, b;
} Color;

float lerp(float a, float b, float t) {
    return a + t * (b - a);
}

Color colormap(float t) {
    Color c;
    if(t < 1e-2){
        c.r = 1.0f;
        c.g = 1.0f;
        c.b = 1.0f;
        return c;  
    }
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

float random_float(float min, float max) {
    return min + (max - min) * ((float)rand() / (float)RAND_MAX);
}

void set_pixel(unsigned char *fb, int width, int x, int y,
               float ***velocity, float **density, float **pressure) {
    // float r = sqrt(
    //     velocity[x][y][0]*velocity[x][y][0] + 
    //     velocity[x][y][1]*velocity[x][y][1]
    // );
    // r*= 0.1;
    // r = powf(r, 0.5f);
    // float b = density[x][y];
    // if (b < 0.0f) b = 0.0f;
    // if (b > 1.0f) b = 1.0f;
    // float g = r*b;
    
    float mag = sqrt(
        velocity[x][y][0] * velocity[x][y][0] +
        velocity[x][y][1] * velocity[x][y][1]
    );
    // float factor = density[x][y];
    // if(factor < 0.3f) factor += mag*density[x][y]/3.0f;
    // if(factor > 1.0f) factor = 1.0f;
    // if(factor < 0.0f) factor = 0.0f;
    float factor = mag/10.0f;
    // Color c = colormap(mag);
    int index = (y * width + x) * 3;
    fb[index + 0] = (unsigned char)(factor* 255.0f);
    fb[index + 1] = (unsigned char)(factor* 255.0f);
    fb[index + 2] = (unsigned char)(factor* 255.0f);
}
float rand_float(float min, float max) {
    return min + ((float)rand() / RAND_MAX) * (max - min);
}


int render(
    GLFWwindow* window,
    unsigned char *framebuffer,
    int emission_area,
    float emission_rate,
    Data data
){

    int counter = 0;
    int source_y = emission_area;
    int source_x = WIDTH/2;
    int inc_dec = 5; // 0-> increas, 1->decrease
    glfwMakeContextCurrent(window);
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        

        if(counter < 30){
                    
            for(int a_i = -emission_area; a_i< emission_area; a_i++){
                for(int a_j = -emission_area; a_j< emission_area; a_j++){
                    density[source_x + a_i][source_y + a_j]+=emission_rate;
                    velocity[source_x + a_i][source_y + a_j][1] += 1.0f;
                    // if(density[source_x + a_i][source_y + a_j] > 1.0f) density[source_x][source_y] = 1.0f;
                }
            }

        }
        for(int i = 0; i<HEIGHT; i++){
            for(int j = 0; j<WIDTH; j++){
                velocity[i][j][0] += random_float(-1.0f, 1.0f);
                velocity[i][j][1] += random_float(-1.0f, 1.0f);
            }
        }    
        simulation_step(pressure, pressure_buffer, density, density_buffer, b, forces, velocity, velocity_buffer, data);


        for(int i = 0; i<HEIGHT; i++){
            for(int j = 0; j<WIDTH; j++){
                set_pixel(framebuffer, WIDTH, i, j, velocity, density, pressure);
            }
        }        
        
        glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, framebuffer);
        glfwSwapBuffers(window);
        glfwPollEvents();
        counter++;
    }

    glfwTerminate();
    return 0;
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
    data.h = 0.01f;
    data.dt = 0.05f;
    data.jacobi_iter = 40;
    data.viscosity = 0.5f;
    data.scalar_diffusion = 0.01f;
    data.buoyancy_coeff = 0.1f;
    data.conf_strenght = 1.0f;
    int emission_area = 10;
    
    float emission_rate = 0.05f;
    stream = freopen("output.txt", "w", stdout);


    render(window, framebuffer, emission_area, emission_rate, data);

    // for(int t = 0; t<4; t++){
    //     simulation_step(pressure, pressure_buffer, density, density_buffer, b, forces, velocity, velocity_buffer, data);
    // }
    
    free_arrays();
    // system( "type output.txt" );

    return 0;

}
