/*
BulletFluids 
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
//Portions of this file based on FLUIDS v.2 - SPH Fluid Simulator for CPU and GPU
//Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com
#ifndef B3_FLUID_SPH_PARAMETERS_H
#define B3_FLUID_SPH_PARAMETERS_H
	
#include "Bullet3Common/b3Vector3.h"

///@brief Contains parameters for fluids and fluid-rigid body interaction.
struct b3FluidSphParameters
{
	b3Scalar m_timeStep;				///<Seconds; simulation becomes unstable at > ~0.004s timestep( with setDefaultParameters() ).
	
	///@brief N*m_simulationScale converts N into simulation scale; N/m_simulationScale converts N into world scale.
	///@remarks As SPH fluid simulations are scale sensitive, the simulation
	///(fluid-fluid interaction) is performed at a physically-correct
	///'simulation scale', which is typically much smaller than the 
	///'world scale' at which the particles are rendered.
	b3Scalar m_simulationScale;
	b3Scalar m_sphSmoothRadius;			///<SPH particle interaction radius; use setSphInteractionRadius() to set this; simulation scale; meters.
	
	///@name Kernel function coefficients; dependent on m_sphSmoothRadius; use setSphInteractionRadius() to set these.
	///@{
	b3Scalar m_sphRadiusSquared;		///<m_sphSmoothRadius^2.
	b3Scalar m_poly6KernCoeff;			///<Coefficient of the poly6 kernel; for density calculation.
	b3Scalar m_spikyKernGradCoeff;		///<Coefficient of the gradient of the spiky kernel; for pressure force calculation.
	b3Scalar m_viscosityKernLapCoeff;	///<Coefficient of the Laplacian of the viscosity kernel; for viscosity force calculation.
	b3Scalar m_initialSum; 				///<Self-contributed particle density; should generally be within [0.0, m_sphRadiusSquared^3] (for Wpoly6).
	///@}

	b3Vector3 m_aabbBoundaryMin;		///<Particles cannot move below this boundary; world scale; meters.
	b3Vector3 m_aabbBoundaryMax;		///<Particles cannot move above this boundary; world scale; meters.
	int m_enableAabbBoundary;			///<If nonzero, the particles are confined to m_aabbBoundaryMin and m_aabbBoundaryMax.
	
	b3Vector3 m_gravity;				///<Simulation scale; meters / seconds^2.
	b3Scalar m_sphAccelLimit;			///<Acceleration caused by SPH forces is clamped to this value; world scale; meters / second.
	b3Scalar m_speedLimit;				///<If nonzero, particle speeds are clamped to this value; world scale; meters / second.
	
	b3Scalar m_viscosity;				///<Higher values increase the fluid's resistance to flow; force calculation; pascal*seconds(Pa*s).
	b3Scalar m_restDensity;				///<Used for pressure calculation; kilograms/meters^3
	b3Scalar m_sphParticleMass;			///<Mass of a single particle when calculating SPH density and force; kilograms.
	b3Scalar m_stiffness;				///<Gas constant; higher values make a less compressible, more unstable fluid; pressure calculation; joules.
	
	b3Scalar m_particleDist;			///<Used to determine particle spacing for b3FluidEmitter; simulation scale; meters. 
	b3Scalar m_particleRadius;			///<For collision detection and collision response; world scale; meters.
	b3Scalar m_particleMargin;			///<World scale meters of allowed penetration when colliding with rigids; [0.0, m_particleRadius).
	b3Scalar m_particleMass;			///<Mass of a single particle when colliding with rigid bodies and applying forces; kilograms.
	
	b3Scalar m_boundaryStiff;			///<Spring coefficient; controls the magnitude of the boundary repulsion force.
	b3Scalar m_boundaryDamp;			///<Damping coefficient; controls the influence of relative velocity on the boundary repulsion force.
	b3Scalar m_boundaryFriction;		///<Fraction of tangential velocity removed per frame; [0.0, 1.0]; higher values more unstable.
	b3Scalar m_boundaryRestitution;		///<Fraction of reflected velocity(bounciness); [0.0, 1.0]; higher values more unstable.
	b3Scalar m_boundaryErp;				///<Controls how quickly penetration is removed(per frame impulse: penetration_depth*m_boundaryErp).
	
	b3FluidSphParameters() { setDefaultParameters(); }
	void setDefaultParameters()
	{
		m_timeStep = b3Scalar(0.003);
		m_simulationScale = b3Scalar(0.004);
		setSphInteractionRadius( b3Scalar(0.01) );
	
		m_aabbBoundaryMin.setValue(-B3_LARGE_FLOAT, -B3_LARGE_FLOAT, -B3_LARGE_FLOAT);
		m_aabbBoundaryMax.setValue(B3_LARGE_FLOAT, B3_LARGE_FLOAT, B3_LARGE_FLOAT);
		m_enableAabbBoundary = 0;
	
		m_gravity.setValue(0, b3Scalar(-9.8), 0);
		m_sphAccelLimit = b3Scalar(50000.0);	
		m_speedLimit 	= b3Scalar(0.0);	
		
		m_viscosity 	= b3Scalar(0.2);
		m_restDensity 	= b3Scalar(600.0);
		m_sphParticleMass  = b3Scalar(0.00020543);
		m_stiffness 	= b3Scalar(1.5);
		
		m_particleDist = b3Pow( m_sphParticleMass/m_restDensity, b3Scalar(1.0/3.0) );
		m_particleRadius = b3Scalar(1.0);
		m_particleMargin = b3Scalar(0.05);
		m_particleMass = b3Scalar(0.00020543);
		
		m_boundaryStiff	= b3Scalar(20000.0);
		m_boundaryDamp 	= b3Scalar(256.0);
		m_boundaryFriction 	= b3Scalar(0.0);
		m_boundaryRestitution = b3Scalar(0.0);
		m_boundaryErp = b3Scalar(0.25);
	}
	
	///The grid cell size is dependent on this radius, so b3FluidSph::setGridCellSize() should also be called after this.
	void setSphInteractionRadius(b3Scalar radius)
	{
		m_sphSmoothRadius = radius;
		m_sphRadiusSquared = m_sphSmoothRadius * m_sphSmoothRadius;
		
		m_poly6KernCoeff = b3Scalar(315.0) / ( b3Scalar(64.0) * B3_PI * b3Pow(m_sphSmoothRadius, 9) );
		m_spikyKernGradCoeff = b3Scalar(-45.0) / ( B3_PI * b3Pow(m_sphSmoothRadius, 6) );
		m_viscosityKernLapCoeff = b3Scalar(45.0) / ( B3_PI * b3Pow(m_sphSmoothRadius, 6) );
		
		//
		//m_initialSum = b3Scalar(0.0);
		//m_initialSum = m_sphRadiusSquared*m_sphRadiusSquared*m_sphRadiusSquared;	//poly6 kernel partial result
		m_initialSum = m_sphRadiusSquared*m_sphRadiusSquared*m_sphRadiusSquared * b3Scalar(0.25);
	}
};


#endif


