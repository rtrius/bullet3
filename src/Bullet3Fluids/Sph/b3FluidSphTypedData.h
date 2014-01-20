/*
BulletFluids 
Copyright (c) 2013 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef B3_FLUID_SPH_TYPED_DATA
#define B3_FLUID_SPH_TYPED_DATA

enum b3FluidSphDataType
{
	FSDT_Invalid,
	FSDT_b3FluidSphSolverDefault_SphParticles,		///<b3FluidSphSolverDefault::SphParticles
	
	FSDT_b3FluidSphOpenCL,							///<b3FluidSphOpenCL
	FSDT_b3FluidSortingGridOpenCL,					///<b3FluidSortingGridOpenCL
	
	FSDT_b3FluidHashGridOpenCL,						///<b3FluidHashGridOpenCL
	
	
	FSDT_Custom1,
	FSDT_Custom2,
	FSDT_Custom3,
	FSDT_Custom4
};

///Used for run time type identification of solver-specific data
class b3FluidSphTypedData
{
public:
	virtual ~b3FluidSphTypedData() {}

	virtual b3FluidSphDataType getType() const { return FSDT_Invalid; }
};

#endif
