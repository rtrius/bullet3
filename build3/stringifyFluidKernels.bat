
rem @echo off

premake4 --file=stringifyKernel.lua --kernelfile="../src/BulletFluidsOpenCL/fluidSph.cl" --headerfile="../src/BulletFluidsOpenCL/fluidSphCL.h" --stringname="fluidSphCL" stringify
premake4 --file=stringifyKernel.lua --kernelfile="../src/BulletFluidsOpenCL/fluidSphRigid.cl" --headerfile="../src/BulletFluidsOpenCL/fluidSphRigidCL.h" --stringname="fluidSphRigidCL" stringify

pause