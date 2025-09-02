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
#include <omp.h>
#include <xmmintrin.h>

// Display
#define DISPLAY_WIDTH 900
#define DISPLAY_HEIGHT 600

// Simulation
#define GRID_WIDTH 200
#define GRID_HEIGHT 200



// Performance monitoring
typedef struct {
    double total_time;
    double min_fps;
    double max_fps;
    double avg_fps;
    int frame_count;
    clock_t last_time;
} PerformanceMetrics;

// Initial conditions

FILE *stream;

typedef struct {
    int mode;
    int shader;
} VisData;

typedef struct {
    int dis_width;
    int dis_height;
    int grid_width;
    int grid_height;
} DisplayParams;

typedef struct {
    // 0 -> square
    // 1-> circle
    // 2 -> random??
    int shape;
    // Bounds
    int min_x;
    int max_x;
    int min_y;
    int max_y;

    // Shape generation
    int amount_shapes;
    int emission_area;
    float emission_rate;
    float emmision_velocity_factor;
} InitialCondition;

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
// Optimized performance monitoring
static PerformanceMetrics perf = {0};

void init_performance_monitor() {
    perf.total_time = 0.0;
    perf.min_fps = 1000.0;
    perf.max_fps = 0.0;
    perf.avg_fps = 0.0;
    perf.frame_count = 0;
    perf.last_time = clock();
}
void update_performance_metrics(double current_fps) {
    if (current_fps < 1e-3) current_fps = 1e-3; // evita infinitos/negativos
    perf.frame_count++;
    perf.total_time += 1.0 / current_fps;
    
    if (current_fps < perf.min_fps) perf.min_fps = current_fps;
    if (current_fps > perf.max_fps) perf.max_fps = current_fps;
    perf.avg_fps = perf.frame_count / perf.total_time;
}

void print_performance_summary() {
    printf("\n=== Performance Summary ===\n");
    printf("Total frames: %d\n", perf.frame_count);
    printf("Average FPS: %.2f\n", perf.avg_fps);
    printf("Min FPS: %.2f\n", perf.min_fps);
    printf("Max FPS: %.2f\n", perf.max_fps);
    printf("Total simulation time: %.2f seconds\n", perf.total_time);
    printf("===========================\n");
}
float rand_float(float min, float max) {
    return min + ((float)rand() / RAND_MAX) * (max - min);
}


void write_to_window(
    float **density,
    unsigned char *framebuffer,
    VisData *vis_data,
    DisplayParams *dis_par
){
    int grid_width = dis_par->grid_width;
    int grid_height = dis_par->grid_height;
    int display_width = dis_par->dis_width;
    int display_height = dis_par->dis_height;


    for(int i = 0; i < display_width; i++){
        for(int j = 0; j < display_height; j++){
            // Round index
            int x_orig = (int)round( i * (grid_width-1)) / (display_width-1);
            int y_orig = (int)round( j * (grid_width-1)) / (display_width-1);

            // Get density
            float density_val = density[x_orig][y_orig];
            if(density_val > 1.0f) density_val = 1.0f;
            if(density_val < 0.0f) density_val = 0.0f;

            // Shader
            int fb_index = (j * display_width + i) * 3;
            if(vis_data->shader == 0){ // Grayscale
                framebuffer[fb_index + 0] = (unsigned char)(density_val* 255.0f);
                framebuffer[fb_index + 1] = (unsigned char)(density_val* 255.0f);
                framebuffer[fb_index + 2] = (unsigned char)(density_val* 255.0f);
            } else if (vis_data->shader == 1){
                framebuffer[fb_index + 0] = (unsigned char)((1.0f-density_val)* 255.0f);
                framebuffer[fb_index + 1] = (unsigned char)((1.0f-density_val)* 255.0f);
                framebuffer[fb_index + 2] = (unsigned char)((1.0f-density_val)* 255.0f);
            } else {
                Color c = colormap(density_val);
                framebuffer[fb_index + 0] = (unsigned char)(c.r* 255.0f);
                framebuffer[fb_index + 1] = (unsigned char)(c.g* 255.0f);
                framebuffer[fb_index + 2] = (unsigned char)(c.b* 255.0f);
            }
            // Write to buffer
            

        }
    }
}

void apply_condition(
    InitialCondition *init_cond
){
    int amount = init_cond->amount_shapes;
    for( int q = 0; q <amount; q++){
        // Random source between bounds
        int x_s = init_cond->min_x + rand() % ( init_cond->max_x - init_cond->min_x);
        int y_s = init_cond->min_y + rand() % ( init_cond->max_y - init_cond->min_y);
        int emission_area = init_cond->emission_area;
        int shape = init_cond->shape;
        float emission_rate = init_cond->emission_rate;
        float vel_factor = init_cond->emmision_velocity_factor;

        double theta = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
        float dir_x = (float)cos(theta);
        float dir_y = (float)sin(theta);

        for(int a_i = -emission_area; a_i< emission_area; a_i++){
            for(int a_j = -emission_area; a_j< emission_area; a_j++){
                if(shape == 0){ // Square
                    density[x_s + a_i][y_s + a_j]+= emission_rate;
                    velocity[x_s + a_i][y_s + a_j][0] += vel_factor * dir_x;
                    velocity[x_s + a_i][y_s + a_j][1] += vel_factor* dir_y;
                } else if (shape == 1){ // Circle
                    if (a_i*a_i + a_j*a_j <= emission_area*emission_area) {
                        density[x_s + a_i][y_s + a_j]+= emission_rate;
                        velocity[x_s + a_i][y_s + a_j][0] += vel_factor * dir_x;
                        velocity[x_s + a_i][y_s + a_j][1] += vel_factor * dir_y;
                    }
                }
            }   
        }
    }
}

int render(
    GLFWwindow* window,
    unsigned char *framebuffer,
    DisplayParams *dis_par,
    VisData *vis_data,
    InitialCondition *init_cond,
    Data *data
){
        // OpenMP: equipo estable y afinidad
        omp_set_dynamic(0);
        omp_set_max_active_levels(1);
        #if defined(_OPENMP) && (_OPENMP>=200805)
        omp_set_schedule(omp_sched_static, 0);
        #endif
        #ifdef _WIN32
        _putenv_s("OMP_PROC_BIND","close");
        _putenv_s("OMP_PLACES","cores");
        #endif
        omp_set_num_threads( choose_threads(data) ); 
    init_performance_monitor();
    int counter = 0;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(0);
        while (!glfwWindowShouldClose(window)) {
            double t0 = glfwGetTime();
            // Clear buffer
                glClear(GL_COLOR_BUFFER_BIT);
                // Add source of fluid
                if (counter == 0){
                    apply_condition(init_cond);
                }
            // Simulate next step
            simulation_step(pressure, pressure_buffer, density, density_buffer, b, forces, velocity, velocity_buffer, data);
            // Write to framebuffer
            write_to_window(density, framebuffer, vis_data, dis_par);
            
            // Framebuffer to window
            glDrawPixels(dis_par->dis_width, dis_par->dis_height, GL_RGB, GL_UNSIGNED_BYTE, framebuffer);
            glfwSwapBuffers(window);
            glfwPollEvents();

            // fps
                double t1 = glfwGetTime();
                double dt = t1 - t0;
                if (dt <= 0.0) dt = 1e-9;          // ← blinda delta
                double current_fps = 1.0 / dt;
                update_performance_metrics(current_fps);
                counter++;
            if (counter % 5 == 0) {
                printf("Frame %d - FPS: %.1f (Avg: %.1f)\n",
                    counter, current_fps, perf.avg_fps);
            }
            if (counter >= 200) {
                printf("Completed 200 frames for performance testing.\n");
                break;
            }
        }
    print_performance_summary();
    glfwTerminate();
    return 0;
}


void parse_config(
    Data *data,
    DisplayParams *dis_par,
    InitialCondition *init_cond,
    VisData *vis_data
){
    // Open config file
    FILE *f = fopen("config.txt", "r");

    // Temporal buffers
    char line[512];
    float ratio = 1.0f;
    char shape[64];
    char mode[64];
    char shader[64];
    float area_factor;
    float vel_factor;

    while (fgets(line, sizeof(line), f)) {
        // Sim data
        if (sscanf(line, "h = %f", &data->h) == 1) continue;
        if (sscanf(line, "dt = %f", &data->dt) == 1) continue;
        if (sscanf(line, "jacobi_iter = %d", &data->jacobi_iter) == 1) continue;
        if (sscanf(line, "viscosity = %f", &data->viscosity) == 1) continue;
        if (sscanf(line, "scalar_diffusion = %f", &data->scalar_diffusion) == 1) continue;
        if (sscanf(line, "buoyancy_coeff = %f", &data->buoyancy_coeff) == 1) continue;
        if (sscanf(line, "conf_strenght = %f", &data->conf_strenght) == 1) continue;
        
        // Display data
        if (sscanf(line, "display_width = %d", &dis_par->dis_width) == 1) continue;
        if (sscanf(line, "display_height = %d", &dis_par->dis_height) == 1) continue;
        if (sscanf(line, "grid_ratio = %f", &ratio) == 1) continue;

        // // Initial condition
        if (sscanf(line, "shape = %63s", &shape) == 1) continue;
        if (sscanf(line, "amount_shapes = %d", &init_cond->amount_shapes) == 1) continue;
        if (sscanf(line, "emission_area = %f", &area_factor) == 1) continue;
        if (sscanf(line, "emission_rate = %f", &init_cond->emission_rate) == 1) continue;
        if (sscanf(line, "emmision_velocity_factor = %f", &vel_factor) == 1) continue;

        // Visualization
        if (sscanf(line, "mode = %63s", &mode) == 1) continue;
        if (sscanf(line, "shader = %63s", &shader) == 1) continue;
    }
    int x =  (int)roundf(dis_par->dis_width * ratio);
    int y =  (int)roundf(dis_par->dis_height * ratio);
    dis_par->grid_width = x;
    dis_par->grid_height= y;
    data->x = x;
    data->y = y;

    // Mode
    if (strcmp(mode, "2d") == 0) {
        vis_data->mode = 0;
    } else if (strcmp(mode, "3d") == 0) {
        vis_data->mode = 1;
    }
    // Shader
    if (strcmp(shader, "grayscale") == 0) {
        vis_data->shader = 0;
    } else if (strcmp(shader, "gray_inverted") == 0) {
        vis_data->shader= 1;
    } else if (strcmp(shader, "lerp") == 0) {
        vis_data->shader= 2;
    }

    // Init cond
    if (strcmp(shape, "square") == 0) {
        init_cond->shape = 0;
    } else if (strcmp(shape, "circle") == 0) {
        init_cond->shape = 1;
    } else if (strcmp(shape, "rand") == 0) {
        init_cond->shape = 2;
    }
    init_cond->emission_area = (int)roundf(max(x, y)*area_factor);
    init_cond->min_x = init_cond->emission_area;
    init_cond->max_x = x - init_cond->emission_area;
    init_cond->min_y = init_cond->emission_area;
    init_cond->max_y = y - init_cond->emission_area;
    init_cond->emmision_velocity_factor = vel_factor * max(x, y);


    // // Data check
    printf("# Data #\n");
    printf("\t X: %f \n",data->x);
    printf("\t Y: %f \n",data->y);
    printf("\t h: %f \n",data->h);
    printf("\t dt: %f \n",data->dt);
    printf("\t jacobi_iter: %d \n",data->jacobi_iter);
    printf("\t viscosity: %f \n",data->viscosity);
    printf("\t scalar_diffusion: %f \n",data->scalar_diffusion);
    printf("\t buoyancy_coeff: %f \n",data->buoyancy_coeff);
    printf("\t conf_strenght: %f \n\n",data->conf_strenght);

    // // Display
    printf("# Display #\n");
    printf("\t Display Width: %d \n",dis_par->dis_width);
    printf("\t Display Height: %d \n",dis_par->dis_height);
    printf("\t Grid Width: %d \n",dis_par->grid_width);
    printf("\t Grid Height: %d \n",dis_par->grid_height);
    
    // // Initial Condition
    printf("# Initial Condition #\n");
    printf("\t Shape: %d \n",init_cond->shape);
    printf("\t amount shapes: %d \n",init_cond->amount_shapes);
    printf("\t emission area: %d \n",init_cond->emission_area);
    printf("\t emission rate: %f \n",init_cond->emission_rate);
    printf("\t velocity factor: %f \n",init_cond->emmision_velocity_factor);
    printf("\t X Bounds: (%d, %d) \n",init_cond->min_x, init_cond->max_x);
    printf("\t Y Bounds: (%d, %d) \n",init_cond->min_y, init_cond->max_y);
    

    // printf("# Visualization #\n");
    printf("\t mode: %d \n",vis_data->mode);
    printf("\t shader: %d \n",vis_data->shader);
    
    fclose(f);
}

int main() {
    // Init data
    Data data;
    DisplayParams dis_par;
    InitialCondition init_cond;
    VisData vis_data;

    // Fetch config data
    parse_config(&data, &dis_par, &init_cond, &vis_data);

    // Display Init
    if (!glfwInit())
        return -1;

    GLFWwindow* window = glfwCreateWindow(
        dis_par.dis_width, dis_par.dis_height, "Fluid SImulation", NULL, NULL
    );
    if (!window) {
        glfwTerminate();
        return -1;
    }
    unsigned char *framebuffer = malloc(dis_par.dis_width * dis_par.dis_height * 3);

    // Simulation Init
    allocate_arrays(dis_par.grid_width, dis_par.grid_height);

    // Rendering
    render(window, framebuffer, &dis_par, &vis_data, &init_cond, &data);
    free_arrays(dis_par.grid_width, dis_par.grid_height);
    free(framebuffer);

    return 0;

}
