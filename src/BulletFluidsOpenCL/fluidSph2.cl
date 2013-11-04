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

#ifdef cl_amd_printf
	#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

typedef float b3Scalar;
typedef float4 b3Vector3;
#define b3Max max
#define b3Min min


//Note that these are vector3 functions -- OpenCL functions are vector4 functions
inline b3Scalar b3Vector3_length2(b3Vector3 v) { return v.x*v.x + v.y*v.y + v.z*v.z; }
inline b3Scalar b3Vector3_dot(b3Vector3 a, b3Vector3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline b3Vector3 b3Vector3_normalize(b3Vector3 v)
{
	b3Scalar length2 = b3Vector3_length2(v);
	if( length2 != (b3Scalar)0.0f ) v /= sqrt(length2);
	
	return v;
}
//Defined in b3FluidSortingGrid.h
#define INVALID_FIRST_INDEX -1
#define INVALID_LAST_INDEX -2


//Syncronize with 'struct b3FluidSphParametersGlobal' in b3FluidSphParameters.h
typedef struct
{
	b3Scalar m_timeStep;
	b3Scalar m_simulationScale;
	b3Scalar m_sphSmoothRadius;
	b3Scalar m_sphRadiusSquared;
	b3Scalar m_poly6KernCoeff;
	b3Scalar m_spikyKernGradCoeff;
	b3Scalar m_viscosityKernLapCoeff;
	b3Scalar m_initialSum;
} b3FluidSphParametersGlobal;

//Syncronize with 'struct b3FluidSphParametersLocal' in b3FluidSphParameters.h
typedef struct
{
	b3Vector3 m_aabbBoundaryMin;
	b3Vector3 m_aabbBoundaryMax;
	int m_enableAabbBoundary;
	b3Vector3 m_gravity;
	b3Scalar m_sphAccelLimit;
	b3Scalar m_speedLimit;
	b3Scalar m_viscosity;
	b3Scalar m_restDensity;
	b3Scalar m_sphParticleMass;
	b3Scalar m_stiffness;
	b3Scalar m_particleDist;
	b3Scalar m_particleRadius;
	b3Scalar m_particleMargin;
	b3Scalar m_particleMass;
	b3Scalar m_boundaryStiff;
	b3Scalar m_boundaryDamp;
	b3Scalar m_boundaryFriction;
	b3Scalar m_boundaryRestitution;
	b3Scalar m_boundaryErp;
} b3FluidSphParametersLocal;


typedef unsigned int b3FluidGridCombinedPos;
typedef int b3FluidGridCoordinate;

typedef struct
{
	int m_firstIndex;
	int m_lastIndex;
	
} b3FluidGridIterator;

typedef struct 
{
	b3FluidGridCombinedPos m_value;
	int m_index;
	
} b3FluidGridValueIndexPair;

typedef struct
{
	b3FluidGridCoordinate x;		
	b3FluidGridCoordinate y;
	b3FluidGridCoordinate z;
	b3FluidGridCoordinate padding;
	
} b3FluidGridPosition;

b3FluidGridPosition getDiscretePosition(b3Scalar cellSize, b3Vector3 position)	//b3FluidSortingGrid::getDiscretePosition()
{
	b3Vector3 discretePosition = position / cellSize;
	
	b3FluidGridPosition result;
	result.x = (b3FluidGridCoordinate)( (position.x >= 0.0f) ? discretePosition.x : floor(discretePosition.x) );
	result.y = (b3FluidGridCoordinate)( (position.y >= 0.0f) ? discretePosition.y : floor(discretePosition.y) );
	result.z = (b3FluidGridCoordinate)( (position.z >= 0.0f) ? discretePosition.z : floor(discretePosition.z) );
	
	return result;
}

#define B3_FLUID_HASH_GRID_COORD_RANGE 64
b3FluidGridCombinedPos getCombinedPositionModulo(b3FluidGridPosition quantizedPosition)
{
	//as_uint() requires that sizeof(b3FluidGridCombinedPos) == sizeof(b3FluidGridCoordinate)
	//This presents an issue if B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT is #defined
	b3FluidGridCombinedPos unsignedX = as_uint(quantizedPosition.x) % B3_FLUID_HASH_GRID_COORD_RANGE;
	b3FluidGridCombinedPos unsignedY = as_uint(quantizedPosition.y) % B3_FLUID_HASH_GRID_COORD_RANGE;
	b3FluidGridCombinedPos unsignedZ = as_uint(quantizedPosition.z) % B3_FLUID_HASH_GRID_COORD_RANGE;
	
	return unsignedX 
		+ unsignedY * B3_FLUID_HASH_GRID_COORD_RANGE
		+ unsignedZ * B3_FLUID_HASH_GRID_COORD_RANGE* B3_FLUID_HASH_GRID_COORD_RANGE;
}

__kernel void generateValueIndexPairs(__global b3Vector3* fluidPositions, __global b3FluidGridValueIndexPair* out_pairs, 
										b3Scalar cellSize, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	b3FluidGridValueIndexPair result;
	result.m_index = index;
	result.m_value = getCombinedPositionModulo( getDiscretePosition(cellSize, fluidPositions[index]) );
	
	out_pairs[index] = result;
}

__kernel void rearrangeParticleArrays(__global b3FluidGridValueIndexPair* sortedPairs, __global b3Vector3* rearrange, 
										__global b3Vector3* temporary, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	//
	int oldIndex = sortedPairs[index].m_index;
	int newIndex = index;
	
	temporary[newIndex] = rearrange[oldIndex];
}


__kernel void resetGridCells(__global b3FluidGridIterator* out_iterators, int numGridCells)
{
	int index = get_global_id(0);
	if(index >= numGridCells) return;
	
	//Crashes on compiling with Catalyst 13.1 if
	//(b3FluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX} is used directly
	int invalidLowerIndex = INVALID_FIRST_INDEX;
	int invalidUpperIndex = INVALID_LAST_INDEX;
	out_iterators[index] = (b3FluidGridIterator){invalidLowerIndex, invalidUpperIndex};
	//out_iterators[index] = (b3FluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX};
}
__kernel void detectIndexRanges(__global b3Vector3* fluidPosition, __global b3FluidGridValueIndexPair* valueIndexPairs, 
								__global b3FluidGridIterator* out_iterators, b3Scalar gridCellSize, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	b3FluidGridCombinedPos gridCellValue = valueIndexPairs[index].m_value;
	
	int lastValidIndex = numFluidParticles - 1;
	
	//if the next particle has a different b3FluidGridCombinedPos(is in another cell),
	//then this particle has the highest index in its cell
	int isLastParticleInCell = (index < lastValidIndex) ? (gridCellValue != valueIndexPairs[index+1].m_value) : 1;
	
	if(isLastParticleInCell)
	{
		int lowerParticleIndex = index;
		int upperParticleIndex = index;
		while( lowerParticleIndex > 0 && valueIndexPairs[lowerParticleIndex - 1].m_value == gridCellValue ) --lowerParticleIndex;
		
		int gridCellIndex = getCombinedPositionModulo( getDiscretePosition(gridCellSize, fluidPosition[index]) );
		out_iterators[gridCellIndex] = (b3FluidGridIterator){ lowerParticleIndex, upperParticleIndex };
	}
}


//
#define B3_EPSILON FLT_EPSILON
__kernel void sphComputePressure(__constant b3FluidSphParametersGlobal* FG,  __constant b3FluidSphParametersLocal* FL,
								__global b3Vector3* fluidPosition, __global b3Scalar* fluidDensity,
								__global b3FluidGridIterator* cellContents, b3Scalar gridCellSize, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar sum = FG->m_initialSum;
	
	b3FluidGridPosition centerCell = getDiscretePosition(gridCellSize, fluidPosition[i]);
	centerCell.x--;
	centerCell.y--;
	centerCell.z--;
	
	for(b3FluidGridCoordinate offsetZ = 0; offsetZ < 3; ++offsetZ)
		for(b3FluidGridCoordinate offsetY = 0; offsetY < 3; ++offsetY)
			for(b3FluidGridCoordinate offsetX = 0; offsetX < 3; ++offsetX)
			{
				b3FluidGridPosition currentCell = centerCell;
				currentCell.z += offsetZ;
				currentCell.y += offsetY;
				currentCell.x += offsetX;
			
				int gridCellIndex = getCombinedPositionModulo(currentCell);
				b3FluidGridIterator gridCell = cellContents[gridCellIndex];
				
				for(int n = gridCell.m_firstIndex; n <= gridCell.m_lastIndex; ++n)
				{
					b3Vector3 delta = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
					b3Scalar distanceSquared = b3Vector3_length2(delta);
					
					b3Scalar c = FG->m_sphRadiusSquared - distanceSquared;
					sum += (c > 0.0f && i != n) ? c*c*c : 0.0f;		//If c is positive, the particle is within interaction radius(poly6 kernel radius)
				}
			}
	
	fluidDensity[i] = sum * FL->m_sphParticleMass * FG->m_poly6KernCoeff;
}


__kernel void sphComputeForce(__constant b3FluidSphParametersGlobal* FG, __constant b3FluidSphParametersLocal* FL,
							__global b3Vector3* fluidPosition, __global b3Vector3* fluidVelEval, 
							__global b3Vector3* fluidSphForce, __global b3Scalar* fluidDensity,
							__global b3FluidGridIterator* cellContents, b3Scalar gridCellSize, int numFluidParticles)
{
	b3Scalar vterm = FG->m_viscosityKernLapCoeff * FL->m_viscosity;
	
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar density_i = fluidDensity[i];
	b3Scalar invDensity_i = 1.0f / density_i;
	b3Scalar pressure_i = (density_i - FL->m_restDensity) * FL->m_stiffness;
	
	b3Vector3 force = {0.0f, 0.0f, 0.0f, 0.0f};
	
	b3FluidGridPosition centerCell = getDiscretePosition(gridCellSize, fluidPosition[i]);
	centerCell.x--;
	centerCell.y--;
	centerCell.z--;
	
	for(b3FluidGridCoordinate offsetZ = 0; offsetZ < 3; ++offsetZ)
		for(b3FluidGridCoordinate offsetY = 0; offsetY < 3; ++offsetY)
			for(b3FluidGridCoordinate offsetX = 0; offsetX < 3; ++offsetX)
			{
				b3FluidGridPosition currentCell = centerCell;
				currentCell.z += offsetZ;
				currentCell.y += offsetY;
				currentCell.x += offsetX;
				
				int gridCellIndex = getCombinedPositionModulo(currentCell);
				b3FluidGridIterator gridCell = cellContents[gridCellIndex];
				
				for(int n = gridCell.m_firstIndex; n <= gridCell.m_lastIndex; ++n)
				{	
					b3Vector3 delta = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
					b3Scalar distanceSquared = b3Vector3_length2(delta);
					
					if(FG->m_sphRadiusSquared > distanceSquared && i != n)
					{
						b3Scalar density_n = fluidDensity[n];
						b3Scalar invDensity_n = 1.0f / density_n;
						b3Scalar pressure_n = (density_n - FL->m_restDensity) * FL->m_stiffness;
					
						b3Scalar distance = sqrt(distanceSquared);
						b3Scalar c = FG->m_sphSmoothRadius - distance;
						b3Scalar pterm = -0.5f * c * FG->m_spikyKernGradCoeff * (pressure_i + pressure_n);
						pterm /= (distance < B3_EPSILON) ? B3_EPSILON : distance;
						
						b3Scalar dterm = c * invDensity_i * invDensity_n;
						
						force += (delta * pterm + (fluidVelEval[n] - fluidVelEval[i]) * vterm) * dterm;
					}
				}
			}
	
	fluidSphForce[i] = force * FL->m_sphParticleMass;
}

__kernel void applyForces(__constant b3FluidSphParametersGlobal* FG,  __constant b3FluidSphParametersLocal* FL, 
						__global b3Vector3* fluidExternalForce, __global b3Vector3* fluidSphAcceleration,
						__global b3Vector3* fluidVel, __global b3Vector3* fluidVelEval, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Vector3 sphAcceleration = fluidSphAcceleration[i];
	{
		b3Scalar accelMagnitude = sqrt( b3Vector3_length2(sphAcceleration) );
		
		b3Scalar simulationScaleAccelLimit = FL->m_sphAccelLimit * FG->m_simulationScale;
		if(accelMagnitude > simulationScaleAccelLimit) sphAcceleration *= simulationScaleAccelLimit / accelMagnitude;
	}

	b3Vector3 acceleration = FL->m_gravity + sphAcceleration + fluidExternalForce[i] / FL->m_particleMass;
	
	b3Vector3 vel = fluidVel[i];
	
	b3Vector3 vnext = vel + acceleration * FG->m_timeStep;		//v(t+1/2) = v(t-1/2) + a(t) dt	
	fluidVel[i] = vnext;
	
	fluidExternalForce[i] = (b3Vector3){0.0f, 0.0f, 0.0f, 0.0f};
}

inline void resolveAabbCollision_impulse(__constant b3FluidSphParametersGlobal* FG,  __constant b3FluidSphParametersLocal* FL, 
										b3Vector3 velocity, b3Vector3 normal, b3Scalar distance, b3Vector3* out_impulse)
{
	if( distance < 0.0f )	//Negative distance indicates penetration
	{
		b3Scalar penetratingMagnitude = b3Vector3_dot(velocity, -normal);
		if( penetratingMagnitude < 0.0f ) penetratingMagnitude = 0.0f;
		
		b3Vector3 penetratingVelocity = -normal * penetratingMagnitude;
		b3Vector3 tangentialVelocity = velocity - penetratingVelocity;
		
		penetratingVelocity *= 1.0f + FL->m_boundaryRestitution;
		
		b3Scalar positionError = (-distance) * (FG->m_simulationScale/FG->m_timeStep) * FL->m_boundaryErp;
		
		*out_impulse += -( penetratingVelocity + (-normal*positionError) + tangentialVelocity * FL->m_boundaryFriction );
	}
}
inline void accumulateBoundaryImpulse(__constant b3FluidSphParametersGlobal* FG,  __constant b3FluidSphParametersLocal* FL, 
								b3Scalar simScaleParticleRadius, b3Vector3 pos, b3Vector3 vel, b3Vector3* out_impulse)
{
	b3Scalar radius = simScaleParticleRadius;
	b3Scalar simScale = FG->m_simulationScale;
	
	b3Vector3 boundaryMin = FL->m_aabbBoundaryMin;
	b3Vector3 boundaryMax = FL->m_aabbBoundaryMax;
	
	resolveAabbCollision_impulse( FG, FL, vel, (b3Vector3){ 1.0f, 0.0f, 0.0f, 0.0f}, ( pos.x - boundaryMin.x )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FG, FL, vel, (b3Vector3){-1.0f, 0.0f, 0.0f, 0.0f}, ( boundaryMax.x - pos.x )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FG, FL, vel, (b3Vector3){0.0f,  1.0f, 0.0f, 0.0f}, ( pos.y - boundaryMin.y )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FG, FL, vel, (b3Vector3){0.0f, -1.0f, 0.0f, 0.0f}, ( boundaryMax.y - pos.y )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FG, FL, vel, (b3Vector3){0.0f, 0.0f,  1.0f, 0.0f}, ( pos.z - boundaryMin.z )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FG, FL, vel, (b3Vector3){0.0f, 0.0f, -1.0f, 0.0f}, ( boundaryMax.z - pos.z )*simScale - radius, out_impulse );
}
__kernel void collideAabbImpulse(__constant b3FluidSphParametersGlobal* FG,  __constant b3FluidSphParametersLocal* FL, 
								__global b3Vector3* fluidPosition, __global b3Vector3* fluidVel, __global b3Vector3* fluidVelEval, 
								int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Vector3 pos = fluidPosition[i];
	b3Vector3 vel = fluidVel[i];
	
	b3Scalar simScaleParticleRadius = FL->m_particleRadius * FG->m_simulationScale;
	
	b3Vector3 aabbImpulse = (b3Vector3){0.0f, 0.0f, 0.0f, 0.0f};
	accumulateBoundaryImpulse(FG, FL, simScaleParticleRadius, pos, vel, &aabbImpulse);
	
	//Leapfrog integration
	b3Vector3 vnext = vel + aabbImpulse;
	fluidVel[i] = vnext;
}

__kernel void integratePositions(__constant b3FluidSphParametersGlobal* FG, __constant b3FluidSphParametersLocal* FL, 
								__global b3Vector3* fluidPosition, __global b3Vector3* fluidVel, __global b3Vector3* fluidVelEval, 
								int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar timeStepDivSimScale = FG->m_timeStep / FG->m_simulationScale;
	
	b3Vector3 prevVelocity = fluidVelEval[i];	//Velocity at (t-1/2)
	b3Vector3 nextVelocity = fluidVel[i];		//Velocity at (t+1/2)
	
	if(FL->m_speedLimit != 0.0f)
	{
		b3Scalar simulationScaleSpeedLimit = FL->m_speedLimit * FG->m_simulationScale;
	
		b3Scalar speed = sqrt( b3Vector3_length2(nextVelocity) );
		if(speed > simulationScaleSpeedLimit) 
		{
			nextVelocity *= simulationScaleSpeedLimit / speed;
			
			fluidVelEval[i] = (prevVelocity + nextVelocity) * 0.5f;		//v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute (sph)forces later
			fluidVel[i] = nextVelocity;
		}
	}
	
	//Leapfrog integration
	//p(t+1) = p(t) + v(t+1/2)*dt
	fluidPosition[i] += nextVelocity * timeStepDivSimScale;
}

