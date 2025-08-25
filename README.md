# Screensaver - Smoke Simulator
## Tecnologias utilizadas
- GLFW (graficas)
- OpenMP (paralelizacion)
- C / gcc (baseline)
- Make (automatizar)

## Setting up
### GLFW
Para comenzar, primero descarguen la libreria que vamos a usar para abrir pantallas desde windows.
1. Vayan a este link y decarguen el binario pre-compilado de windows de 64-bits.
```
https://www.glfw.org/download.html
```
2. Esto les tuvo que haber descargado una carpeta zip. Entonces ponganle extraer en la direccion ```C:dev/```
3. Listo! ya los headers y arcivos "a" estan en este repositorio y ya les deberia de abrir una pestaña cuando lo corran.
### Make
Para facilitar la compilacion y corrida de todo, va a ser mas facil si lo automatizamos todo con un makefile. Entonces, por si no lo tienen descargado sigan estos pasos:
1. Descargar el msys
```
https://www.msys2.org/
```
darle ok a todas las cosas y que se descarguen
2. Corran el programa MYSYS2 MinGW que les va a abrir una terminal. En esta terminal pongan:
```
pacman -Sy
pacman -S mingw-w64-x86_64-make
```
Para verificar:
```
mingw32-make --version
```
3. En variables de sistema, dentro de path, pongan la direccion
```
C:\msys64\mingw64\bin
```
4. Puede que necesiten reiniciar su compu antes, pero como ultimo paso ponganle un alias al ejecutable de make para hacerlo mas facil:
```
Set-Alias make "C:\msys64\mingw64\bin\mingw32-make.exe" 
```
5. Ya estuvo!!

## Running
### Usando Make
**Correr Secuencial**
```
make seq
```
y para correr
```
make run-seq
```
**Correr Paralelo**
```
make omp
```
y para correr
```
make run-omp
```

### Manual Secuencial
Si les pelo make, pueden hacerlo manual obvio, solo que necesitan correr lo siguiente.
```
gcc main.c -Iinclude -Llib -lglfw3 -lopengl32 -lgdi32 -o bin/main.exe
```

correr:
```
bin/main.exe
```