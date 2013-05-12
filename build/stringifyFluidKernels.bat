
rem @echo off

premake4 --file=stringifyKernel.lua --kernelfile="../src/BulletFluidsOpenCL/fluidSph.cl" --headerfile="../src/BulletFluidsOpenCL/fluidSphCL.h" --stringname="fluidSphCL" stringify

pause