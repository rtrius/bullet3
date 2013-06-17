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
//Portions of this file based on FLUIDS v.2 - SPH Fluid Simulator for CPU and GPU
//Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com
#include "b3FluidSphSolver.h"

#include "b3FluidSortingGrid.h"

void b3FluidSphSolver::applyForcesSingleFluid(const b3FluidSphParametersGlobal& FG, b3FluidSph* fluid)
{
	B3_PROFILE("b3FluidSphSolver::applyForcesSingleFluid()");

	const b3FluidSphParametersLocal& FL = fluid->getLocalParameters();
	b3FluidParticles& particles = fluid->internalGetParticles();
	
	const b3Scalar invParticleMass = b3Scalar(1.0) / FL.m_particleMass;
	
	for(int i = 0; i < particles.size(); ++i)
	{
		b3Vector3& vel = particles.m_vel[i];
	
		b3Vector3 acceleration = FL.m_gravity + (particles.m_accumulatedForce[i] * invParticleMass);

		//Leapfrog integration
		b3Vector3 vnext = vel + acceleration * FG.m_timeStep;	//v(t+1/2) = v(t-1/2) + a(t) dt
		vel = vnext;
	}
	
	for(int i = 0; i < particles.size(); ++i) particles.m_accumulatedForce[i].setValue(0, 0, 0);
}
void b3FluidSphSolver::integratePositionsSingleFluid(const b3FluidSphParametersGlobal& FG, const b3FluidSphParametersLocal& FL, 
														b3FluidParticles& particles)
{
	B3_PROFILE("b3FluidSphSolver::integratePositionsSingleFluid()");
		
	//Velocity is at simulation scale; divide by simulation scale to convert to world scale
	b3Scalar timeStepDivSimScale = FG.m_timeStep / FG.m_simulationScale;
	
	b3Scalar simulationScaleSpeedLimit = FL.m_speedLimit * FG.m_simulationScale;
	
	if( simulationScaleSpeedLimit != b3Scalar(0.0) )
	{
		for(int i = 0; i < particles.size(); ++i)
		{
			b3Vector3 prevVelocity = particles.m_vel_eval[i];	//Velocity at (t-1/2)
			b3Vector3 nextVelocity = particles.m_vel[i];		//Velccity at (t+1/2)
			
			b3Scalar speed = nextVelocity.length();
			if(speed > simulationScaleSpeedLimit) 
			{
				nextVelocity *= simulationScaleSpeedLimit / speed;
					
				particles.m_vel_eval[i] = (prevVelocity + nextVelocity) * b3Scalar(0.5);	//v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute (sph)forces later
				particles.m_vel[i] = nextVelocity;
			}
		}	
	}
		
	//Leapfrog integration
	//p(t+1) = p(t) + v(t+1/2)*dt
	for(int i = 0; i < particles.size(); ++i) particles.m_pos[i] += particles.m_vel[i] * timeStepDivSimScale;
}

void b3FluidSphSolverDefault::sphComputePressure(const b3FluidSphParametersGlobal& FG, b3FluidSph* fluid, b3FluidSphSolverDefault::SphParticles& sphData)
{
	B3_PROFILE("b3FluidSphSolverDefault::sphComputePressure()");
	
	const int numParticles = fluid->numParticles();
	
	const b3FluidSphParametersLocal& FL = fluid->getLocalParameters();
	const b3FluidSortingGrid& grid = fluid->getGrid();
	b3FluidParticles& particles = fluid->internalGetParticles();
	
	{
		B3_PROFILE("sphComputePressure() - reset sums, clear table");
		
		for(int i = 0; i < numParticles; ++i) sphData.m_invDensity[i] = FG.m_initialSum;
		for(int i = 0; i < numParticles; ++i) sphData.m_neighborTable[i].clear();
	}
	
	{
		B3_PROFILE("sphComputePressure() - compute sums");
		
		for(int group = 0; group < b3FluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
		{
			const b3AlignedObjectArray<int>& currentGroup = grid.internalGetMultithreadingGroup(group);
			if( !currentGroup.size() ) continue;
			
			computeSumsInMultithreadingGroup(FG, currentGroup, grid, particles, sphData);
		}
	}
	
	{
		B3_PROFILE("sphComputePressure() - compute pressure/density");
		
		for(int i = 0; i < numParticles; ++i)
		{
			b3Scalar density = sphData.m_invDensity[i] * FL.m_sphParticleMass * FG.m_poly6KernCoeff;
			sphData.m_pressure[i] = (density - FL.m_restDensity) * FL.m_stiffness;
			sphData.m_invDensity[i] = b3Scalar(1.0) / density;
		}
	}
}

void b3FluidSphSolverDefault::sphComputeForce(const b3FluidSphParametersGlobal& FG, b3FluidSph* fluid, b3FluidSphSolverDefault::SphParticles& sphData)
{
	B3_PROFILE("b3FluidSphSolverDefault::sphComputeForce()");
	
	const b3FluidSphParametersLocal& FL = fluid->getLocalParameters();
	const b3FluidSortingGrid& grid = fluid->getGrid();
	b3FluidParticles& particles = fluid->internalGetParticles();
	
	b3Scalar vterm = FG.m_viscosityKernLapCoeff * FL.m_viscosity;
	
	for(int i = 0; i < particles.size(); ++i)sphData.m_sphForce[i].setValue(0, 0, 0);
	
	for(int group = 0; group < b3FluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
	{
		const b3AlignedObjectArray<int>& currentGroup = grid.internalGetMultithreadingGroup(group);
		if( !currentGroup.size() ) continue;
		
		computeForcesInMultithreadingGroup(FG, vterm, currentGroup, grid, particles, sphData);
	}
	
	for(int i = 0; i < particles.size(); ++i)sphData.m_sphForce[i] *= FL.m_sphParticleMass;
}

void b3FluidSphSolverDefault::calculateSumsInCellSymmetric(const b3FluidSphParametersGlobal& FG, int gridCellIndex, 
															const b3FluidSortingGrid& grid, b3FluidParticles& particles,
															b3FluidSphSolverDefault::SphParticles& sphData)
{
	
	b3FluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	if(currentCell.m_firstIndex <= currentCell.m_lastIndex)	//if cell is not empty
	{
		b3FluidSortingGrid::FoundCells foundCells;
		grid.findCellsSymmetric(particles.m_pos[currentCell.m_firstIndex], foundCells);
		
		for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
		{
			//Remove particle, with index i, from grid cell to prevent self-particle interactions
			++foundCells.m_iterators[0].m_firstIndex;	//Local cell; currentCell == foundCells.m_iterators[0]
			
			for(int cell = 0; cell < b3FluidSortingGrid::NUM_FOUND_CELLS_SYMMETRIC; cell++) 
			{
				b3FluidGridIterator& FI = foundCells.m_iterators[cell];
				
				for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
				{
					//Simulation-scale distance
					b3Vector3 difference = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		
					b3Scalar distanceSquared = difference.length2();
					
					if(FG.m_sphRadiusSquared > distanceSquared)
					{
						b3Scalar c = FG.m_sphRadiusSquared - distanceSquared;
						b3Scalar poly6KernPartialResult = c * c * c;
						sphData.m_invDensity[i] += poly6KernPartialResult;
						sphData.m_invDensity[n] += poly6KernPartialResult;
						
						b3Scalar distance = b3Sqrt(distanceSquared);
						if( !sphData.m_neighborTable[i].isFilled() ) sphData.m_neighborTable[i].addNeighbor(n, distance);
						else if( !sphData.m_neighborTable[n].isFilled() ) sphData.m_neighborTable[n].addNeighbor(i, distance);
						else 
						{
							cell = b3FluidSortingGrid::NUM_FOUND_CELLS_SYMMETRIC;	//Break out of outer loop 
							break;
						}
					}
				}
			}
		}
	}
}

void computeForceNeighborTableSymmetric(const b3FluidSphParametersGlobal& FG, const b3Scalar vterm, int particleIndex, 
										b3FluidParticles& particles, b3FluidSphSolverDefault::SphParticles& sphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < sphData.m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = sphData.m_neighborTable[i].getNeighborIndex(j);
		
		b3Vector3 difference = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		b3Scalar distance = sphData.m_neighborTable[i].getDistance(j);
		
		b3Scalar c = FG.m_sphSmoothRadius - distance;
		b3Scalar pterm = b3Scalar(-0.5) * c * FG.m_spikyKernGradCoeff * (sphData.m_pressure[i] + sphData.m_pressure[n]);
		pterm /= (distance < B3_EPSILON) ? B3_EPSILON : distance;
		
		b3Scalar dterm = c * sphData.m_invDensity[i] * sphData.m_invDensity[n];

		b3Vector3 force(  (pterm * difference.getX() + vterm * (particles.m_vel_eval[n].getX() - particles.m_vel_eval[i].getX())) * dterm,
						  (pterm * difference.getY() + vterm * (particles.m_vel_eval[n].getY() - particles.m_vel_eval[i].getY())) * dterm,
						  (pterm * difference.getZ() + vterm * (particles.m_vel_eval[n].getZ() - particles.m_vel_eval[i].getZ())) * dterm );
		
		sphData.m_sphForce[i] += force;
		sphData.m_sphForce[n] += -force;
	}
}
void b3FluidSphSolverDefault::calculateForcesInCellSymmetric(const b3FluidSphParametersGlobal& FG, const b3Scalar vterm,
															int gridCellIndex, const b3FluidSortingGrid& grid, b3FluidParticles& particles,
															b3FluidSphSolverDefault::SphParticles& sphData)
{
	b3FluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
	{
		computeForceNeighborTableSymmetric(FG, vterm, i, particles, sphData);
	}
}
