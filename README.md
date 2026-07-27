# Screensaver - Smoke Simulator

![Demo](demo.gif)

## Introduction
This project is a screensaver that simulates smoke dissipating through a 2D medium, built first as a sequential implementation and then parallelized with OpenMP to compare performance. It uses an Eulerian grid method, discretizing the domain and applying per-cell operations — advection, diffusion, and pressure projection — that are naturally parallelizable across threads.

Beyond producing an appealing visual effect, the goal is to apply parallel computing concepts in practice, using a classic fluid/gas dispersion problem to highlight the performance difference between a sequential and an optimized, multi-threaded approach.

## Technologies used
- GLFW (graphics)
- OpenMP (parallelization)
- C / gcc (baseline)
- Make (automation)

## Setting up
### GLFW
To get started, first download the library we'll use to open windows on Windows.
1. Go to this link and download the pre-compiled 64-bit Windows binary.
```
https://www.glfw.org/download.html
```
2. This should have downloaded a zip folder. Extract it to the path ```C:dev/```
3. Done! The headers and "a" files are already in this repository, and a window should open when you run it.
### Make
To make compiling and running everything easier, it'll be simpler if we automate it all with a makefile. So, in case you don't have it downloaded, follow these steps:
1. Download msys
```
https://www.msys2.org/
```
Click ok on everything and let it download
2. Run the MSYS2 MinGW program, which will open a terminal. In this terminal type:
```
pacman -Sy
pacman -S mingw-w64-x86_64-make
```
To verify:
```
mingw32-make --version
```
3. In your system's environment variables, under path, add
```
C:\msys64\mingw64\bin
```
4. You may need to restart your computer first, but as a last step, set up an alias for the make executable to make things easier:
```
Set-Alias make "C:\msys64\mingw64\bin\mingw32-make.exe" 
```
5. That's it!!

## Running
### Using Make
**Run Sequential**
```
make seq
```
and to run it
```
make run-seq
```
**Run Parallel**
```
make omp
```
and to run it
```
make run-omp
```

### Manual Sequential
If make gives you trouble, you can obviously do it manually, you just need to run the following.
```
gcc main.c -Iinclude -Llib -lglfw3 -lopengl32 -lgdi32 -o bin/main.exe
```

run:
```
bin/main.exe
```

## Function Reference

### Rendering and color

| Function | Inputs | Outputs | Description |
|---|---|---|---|
| `float lerp(float a, float b, float t)` | `a` (float): start value.<br>`b` (float): end value.<br>`t` (float): factor in [0,1]. | (float): interpolated value. | Linear interpolation; used for gradients and colors. |
| `Color colormap(float t)` | `t` (float): normalized intensity in [0,1]. | `Color` (struct with r,g,b floats). | Maps a scalar value to a color (blue → green → yellow → red). |
| `void write_to_window(float **density, unsigned char *framebuffer, VisData *vis_data, DisplayParams *dis_par)` | `density` (float\*\*): density grid.<br>`framebuffer` (uchar\*): RGB buffer.<br>`vis_data` (VisData\*): mode and shader.<br>`dis_par` (DisplayParams\*): display/grid dimensions. | None (side effect). | Converts density into colored pixels and writes them to the framebuffer for display. |

### Simulation setup and initial conditions

| Function | Inputs | Outputs | Description |
|---|---|---|---|
| `void parse_config(Data *data, DisplayParams *dis_par, InitialCondition *init_cond, VisData *vis_data)` | `data`: physical parameters (h, dt, jacobi_iter, viscosity, etc.).<br>`dis_par`: window and grid dimensions.<br>`init_cond`: emission parameters.<br>`vis_data`: mode and shader. | None (side effect). | Reads `config.txt`, validates values, and loads the simulation and visualization structs. |
| `void apply_condition(InitialCondition *init_cond)` | `init_cond`: defines emission shape, amount, area, and velocity. | None (side effect). | Generates smoke sources, injecting density and initial velocity into the grid. |

### Main render loop

| Function | Inputs | Outputs | Description |
|---|---|---|---|
| `int render(GLFWwindow *window, unsigned char *framebuffer, DisplayParams *dis_par, VisData *vis_data, InitialCondition *init_cond, Data *data)` | `window`: OpenGL context.<br>`framebuffer`: RGB buffer.<br>`dis_par`, `vis_data`, `init_cond`, `data`: runtime structs. | (int): 0 on clean exit. | Runs the simulation loop: applies initial conditions, calls `simulation_step`, generates and draws pixels, measures FPS, and shows stats. |

### State and memory management

| Function | Inputs | Outputs | Description |
|---|---|---|---|
| `void allocate_arrays(int X, int Y)` | `X`, `Y`: grid dimensions. | None (side effect). | Allocates contiguous memory for all fields (velocity, forces, pressure, density, b) and zero-initializes them. |
| `void free_arrays(int X, int Y)` | `X`, `Y`: dimensions used. | None (side effect). | Safely frees all associated arrays and buffers. |

### Copy and debug utilities

| Function | Inputs | Outputs | Description |
|---|---|---|---|
| `void copy_vec(float ***target, float ***origin, Data *data)` | `target`: destination.<br>`origin`: source.<br>`data`: dimensions. | None. | Copies a vector field cell by cell. Uses OpenMP with static scheduling in the parallel version. |
| `void copy_scalar(float **target, float **origin, Data *data)` | `target`: destination.<br>`origin`: source.<br>`data`: dimensions. | None. | Copies a full scalar field; conditionally parallelized on large grids. |
| `void print_all(int stage, float ***force, float ***velocity, float **density, float **pressure, Data *data)` | `stage`: stage number.<br>`force`, `velocity`: vector fields.<br>`density`, `pressure`: scalar fields.<br>`data`: dimensions. | None. | Prints full matrices to the console for debugging. |

### Simulation core (smoke flow)

| Function | Inputs | Outputs | Description |
|---|---|---|---|
| `void update_forces(float ***velocity, float ***forces, float ***buffer_vec, float **buffer_scalar, Data *data)` | `velocity`, `forces`, `buffer_vec`, `buffer_scalar`, `data` with physical parameters. | None. | Computes vorticity and buoyancy/confinement forces, updating velocity. |
| `void density_steps(float ***velocity, float **density, float **density_buffer, Data *data)` | `velocity`, `density`, `density_buffer`, `data`. | None. | Applies semi-Lagrangian advection and diffusion (Jacobi) to density. |
| `void velocity_steps(float ***velocity, float ***velocity_buffer, Data *data)` | `velocity`, `velocity_buffer`, `data`. | None. | Advects velocity by itself and applies viscous diffusion with Jacobi. |
| `void pressure_projection(float **pressure, float **pressure_buffer, float **b, float **density, float ***velocity, Data *data)` | Pressure, density, velocity, and `b` fields.<br>`data`: parameters. | None. | Solves the Poisson equation (Jacobi) and projects velocity onto an incompressible field. |
| `void simulation_step(float **pressure, float **pressure_buffer, float **density, float **density_buffer, float **b, float ***forces, float ***velocity, float ***velocity_buffer, Data *data)` | All fields and buffers. | None. | Runs one full step: `update_forces` → `velocity_steps` → `density_steps` → `pressure_projection`. |

### Performance metrics

| Function | Inputs | Outputs | Description |
|---|---|---|---|
| `void init_performance_monitor(void)` | None. | None. | Initializes FPS counters and accumulated time. |
| `void update_performance_metrics(double current_fps)` | `current_fps`: FPS of the current frame. | None. | Updates min/avg/max FPS and total time metrics. |
| `void print_performance_summary(void)` | None. | None. | Prints a performance summary to the console at the end of the simulation. |

### Miscellaneous utilities

| Function | Inputs | Outputs | Description |
|---|---|---|---|
| `float rand_float(float min, float max)` | `min`, `max`: range bounds. | (float) random value in [min,max]. | Helper generator for randomized initial conditions (e.g., emission direction). |

## Conclusions
The OpenMP version consistently cuts execution time and, for most grid sizes, raises FPS compared to the sequential baseline. This tracks with the per-cell nature of the computation (advection, diffusion, pressure projection), which splits cleanly across threads. At small grid sizes the gain shrinks due to thread creation/sync overhead; as the grid grows, useful work per thread dominates that overhead and speedup improves — until memory bandwidth and cache locality become the limiting factor, particularly in the Jacobi solver used for diffusion and pressure.

Key implementation choices support this: `schedule(static, chunk)` keeps load balanced and locality-friendly, `omp simd` vectorizes the inner loop on top of thread-level parallelism, and thread affinity (`OMP_PROC_BIND`, `OMP_PLACES`) reduces migrations and stabilizes memory latency. Buffer swaps and boundary updates use a `single` + `barrier` pattern to minimize synchronization while avoiding data races. Numerically, the parallel version preserves the same method (semi-Lagrangian advection, Jacobi for diffusion/Poisson, pressure projection) — floating-point summation order can differ slightly from the sequential run, but this has no visible effect on stability or the physical behavior of the smoke.

**Takeaways**
- Parallelization is effective from medium grid sizes onward, improving both runtime and FPS while keeping visual fidelity.
- The per-cell workload is highly parallelizable; `omp for (static)` + `omp simd` + `single`/`barrier` in the Jacobi step gives a solid performance/simplicity tradeoff.
- The practical ceiling is memory (bandwidth/locality) and the unavoidable sequential fraction (Amdahl's law), so speedup grows then plateaus.
- Structural optimizations — contiguous memory blocks, a single thread-count policy, core affinity — improved consistency and overall performance.

**Future work**
- Explore `collapse(2)` or tiling in the 2D loops for better cache use.
- Try `guided` scheduling for cases with uneven active regions.
- Use multigrid or a preconditioner to accelerate the Poisson solve instead of plain Jacobi.
- Move vector fields to a Structure-of-Arrays layout to further favor vectorization.
- Consider MPI + OpenMP if scaling beyond a single node is needed.

## Credits
- [MarchMol](https://github.com/MarchMol) — math and visualization
- [JPS4321](https://github.com/JPS4321) — sequential and parallel approaches
- [Sofiamishel2003](https://github.com/Sofiamishel2003) — optimization

## References
- Admin. (2022, September 2). *Jacobian Method*. BYJUS. https://byjus.com/maths/jacobian-method/
- Colin, B., & Sandu, A. (n.d.). *Fluid Simulation For Computer Graphics: A Tutorial in Grid Based and Particle Based Methods* (Virginia Tech). https://cg.informatik.uni-freiburg.de/intern/seminar/gridFluids_fluid-EulerParticle.pdf
- Kuhn, M., Latu, G., Crouseilles, N., & Genaud, S. (2015). Parallelization of an Advection-Diffusion problem arising in edge plasma physics using hybrid MPI/OpenMP programming. In *Lecture Notes in Computer Science* (pp. 545–557). https://doi.org/10.1007/978-3-662-48096-0_42
