# Screensaver - Smoke Simulator

![Demo](demo.gif)

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

## Credits
- [MarchMol](https://github.com/MarchMol) — math and visualization
- [JPS4321](https://github.com/JPS4321) — sequential and parallel approaches
- [Sofiamishel2003](https://github.com/Sofiamishel2003) — optimization
