/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef B3_FLUID_SPH_OPENCL
#define B3_FLUID_SPH_OPENCL

#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"

class b3Vector3;
struct b3FluidSphParametersLocal;
struct b3FluidParticles;

///@brief Manages OpenCL buffers corresponding to b3FluidParticles and b3FluidSphParametersLocal.
class b3FluidSphOpenCL
{
public:
	bool m_initialized;

	b3OpenCLArray<b3FluidSphParametersLocal> m_localParameters;
	
	b3OpenCLArray<b3Vector3> m_pos;
	b3OpenCLArray<b3Vector3> m_vel;
	b3OpenCLArray<b3Vector3> m_vel_eval;
	b3OpenCLArray<b3Vector3> m_accumulatedForce;
	
	b3OpenCLArray<b3Vector3> m_sph_force;
	b3OpenCLArray<b3Scalar> m_density;
	b3OpenCLArray<int> m_cellIndex;
	
	b3FluidSphOpenCL(cl_context context, cl_command_queue queue) :
		m_initialized(0),
	
		m_localParameters(context, queue),
		m_pos(context, queue),
		m_vel(context, queue),
		m_vel_eval(context, queue),
		m_accumulatedForce(context, queue),
		m_sph_force(context, queue),
		m_density(context, queue),
		m_cellIndex(context, queue) {}
	
	void writeToOpenCL(cl_command_queue queue, const b3FluidSphParametersLocal& FL, b3FluidParticles& particles);
	void readFromOpenCL(cl_command_queue queue, b3FluidParticles& particles, b3AlignedObjectArray<b3Vector3>& sphForce);
};

#endif
