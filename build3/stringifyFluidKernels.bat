
rem @echo off

premake4 --file=stringifyKernel.lua --kernelfile="../src/Bullet3FluidsOpenCL/fluidSph.cl" --headerfile="../src/Bullet3FluidsOpenCL/fluidSphCL.h" --stringname="fluidSphCL" stringify
premake4 --file=stringifyKernel.lua --kernelfile="../src/Bullet3FluidsOpenCL/fluidSphRigid.cl" --headerfile="../src/Bullet3FluidsOpenCL/fluidSphRigidCL.h" --stringname="fluidSphRigidCL" stringify

pause