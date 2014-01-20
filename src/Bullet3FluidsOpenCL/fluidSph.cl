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

//Syncronize with 'struct b3FluidSphParameters' in b3FluidSphParameters.h
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
} b3FluidSphParameters;


//#define B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT	//Ensure that this is also #defined in b3FluidSortingGrid.h
#ifdef B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT	
	typedef unsigned long b3FluidGridUint64;
	typedef b3FluidGridUint64 b3FluidGridCombinedPos;	//Range must contain B3_FLUID_GRID_COORD_RANGE^3
	#define B3_FLUID_GRID_COORD_RANGE 2097152		//2^21
	
	inline void splitCombinedPosition(b3FluidGridUint64 resolutionX, b3FluidGridUint64 resolutionY, 
										b3FluidGridUint64 value, int* out_x, int* out_y, int* out_z)
	{
		b3FluidGridUint64 cellsPerLine = resolutionX;
		b3FluidGridUint64 cellsPerPlane = resolutionX * resolutionY;
		
		b3FluidGridUint64 x = value % cellsPerLine;
		b3FluidGridUint64 z = value / cellsPerPlane;
		b3FluidGridUint64 y = (value - z*cellsPerPlane) / cellsPerLine;
		
		*out_x = (int)x;
		*out_z = (int)z;
		*out_y = (int)y;
	}
#else
	typedef unsigned int b3FluidGridCombinedPos;		//Range must contain B3_FLUID_GRID_COORD_RANGE^3
	#define B3_FLUID_GRID_COORD_RANGE 1024			//2^10	
	
	inline void splitCombinedPosition(int resolutionX, int resolutionY, int value, int* out_x, int* out_y, int* out_z)
	{
		int x = value % resolutionX;
		int z = value / (resolutionX*resolutionY);
		int y = (value - z*resolutionX*resolutionY) / resolutionX;
		
		*out_x = (int)x;
		*out_z = (int)z;
		*out_y = (int)y;
	}
#endif

typedef int b3FluidGridCoordinate;
#define B3_FLUID_GRID_COORD_RANGE_HALVED B3_FLUID_GRID_COORD_RANGE/2



typedef struct
{
	int m_firstIndex;
	int m_lastIndex;
	
} b3FluidGridIterator;


//Since the hash function used to determine the 'value' of particles is simply 
//(x + y*CELLS_PER_ROW + z*CELLS_PER_PLANE), adjacent cells have a value 
//that is 1 greater and lesser than the current cell. 
//This makes it possible to query 3 cells simultaneously(as a 3 cell bar extended along the x-axis) 
//by using a 'binary range search' in the range [current_cell_value-1, current_cell_value+1]. 
//Furthermore, as the 3 particle index ranges returned are also adjacent, it is also possible to 
//stitch them together to form a single index range.
#define b3FluidSortingGrid_NUM_FOUND_CELLS_GPU 9

typedef struct
{
	b3FluidGridIterator m_iterators[b3FluidSortingGrid_NUM_FOUND_CELLS_GPU];
	
} b3FluidSortingGridFoundCellsGpu;		//b3FluidSortingGrid::FoundCellsGpu in b3FluidSortingGrid.h

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
b3FluidGridCombinedPos getCombinedPosition(b3FluidGridPosition quantizedPosition)	//b3FluidGridPosition::getCombinedPosition()
{
	b3FluidGridCoordinate signedX = quantizedPosition.x + B3_FLUID_GRID_COORD_RANGE_HALVED;
	b3FluidGridCoordinate signedY = quantizedPosition.y + B3_FLUID_GRID_COORD_RANGE_HALVED;
	b3FluidGridCoordinate signedZ = quantizedPosition.z + B3_FLUID_GRID_COORD_RANGE_HALVED;
	
	b3FluidGridCombinedPos unsignedX = (b3FluidGridCombinedPos)signedX;
	b3FluidGridCombinedPos unsignedY = (b3FluidGridCombinedPos)signedY * B3_FLUID_GRID_COORD_RANGE;
	b3FluidGridCombinedPos unsignedZ = (b3FluidGridCombinedPos)signedZ * B3_FLUID_GRID_COORD_RANGE * B3_FLUID_GRID_COORD_RANGE;
	
	return unsignedX + unsignedY + unsignedZ;
}

__kernel void generateValueIndexPairs(__global b3Vector3* fluidPositions, __global b3FluidGridValueIndexPair* out_pairs, 
										b3Scalar cellSize, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	b3FluidGridValueIndexPair result;
	result.m_index = index;
	result.m_value = getCombinedPosition( getDiscretePosition(cellSize, fluidPositions[index]) );
	
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


__kernel void markUniques(__global b3FluidGridValueIndexPair* valueIndexPairs, __global int* out_retainValueAtThisIndex, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	int lastValidIndex = numFluidParticles - 1;
	
	//Retain if the next particle has a different b3FluidGridCombinedPos(is in another cell)
	int isRetained = (index < lastValidIndex) ? (valueIndexPairs[index].m_value != valueIndexPairs[index+1].m_value) : 1;
	
	out_retainValueAtThisIndex[index] = isRetained;
}

__kernel void storeUniquesAndIndexRanges(__global b3FluidGridValueIndexPair* valueIndexPairs, __global int* retainValue, __global int* scanResults, 
							__global b3FluidGridCombinedPos* out_sortGridValues, __global b3FluidGridIterator* out_iterators, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	if(retainValue[index])
	{
		int gridCellIndex = scanResults[index];
		b3FluidGridCombinedPos gridCellValue = valueIndexPairs[index].m_value;
		
		out_sortGridValues[gridCellIndex] = gridCellValue;
		
		//Perform a linear search for the lower range of this grid cell.
		//When r, the SPH interaction radius, is equivalent to the grid cell size, 
		//the average number of particles per cell is frequently 4 to 7.
		int lowerParticleIndex = index;
		int upperParticleIndex = index;
		while( lowerParticleIndex > 0 && valueIndexPairs[lowerParticleIndex - 1].m_value == gridCellValue ) --lowerParticleIndex;
		
		out_iterators[gridCellIndex] = (b3FluidGridIterator){ lowerParticleIndex, upperParticleIndex };
	}
}

b3FluidGridCombinedPos getCombinedPosition_yAxisOriented(b3FluidGridPosition quantizedPosition)
{
	b3FluidGridCoordinate signedX = quantizedPosition.x + B3_FLUID_GRID_COORD_RANGE_HALVED;
	b3FluidGridCoordinate signedY = quantizedPosition.y + B3_FLUID_GRID_COORD_RANGE_HALVED;
	b3FluidGridCoordinate signedZ = quantizedPosition.z + B3_FLUID_GRID_COORD_RANGE_HALVED;
	
	b3FluidGridCombinedPos unsignedX = (b3FluidGridCombinedPos)signedX * B3_FLUID_GRID_COORD_RANGE;
	b3FluidGridCombinedPos unsignedY = (b3FluidGridCombinedPos)signedY;
	b3FluidGridCombinedPos unsignedZ = (b3FluidGridCombinedPos)signedZ * B3_FLUID_GRID_COORD_RANGE * B3_FLUID_GRID_COORD_RANGE;
	
	return unsignedX + unsignedY + unsignedZ;
}
b3FluidGridCombinedPos getCombinedPosition_zAxisOriented(b3FluidGridPosition quantizedPosition)
{
	b3FluidGridCoordinate signedX = quantizedPosition.x + B3_FLUID_GRID_COORD_RANGE_HALVED;
	b3FluidGridCoordinate signedY = quantizedPosition.y + B3_FLUID_GRID_COORD_RANGE_HALVED;
	b3FluidGridCoordinate signedZ = quantizedPosition.z + B3_FLUID_GRID_COORD_RANGE_HALVED;
	
	b3FluidGridCombinedPos unsignedX = (b3FluidGridCombinedPos)signedX * B3_FLUID_GRID_COORD_RANGE * B3_FLUID_GRID_COORD_RANGE;
	b3FluidGridCombinedPos unsignedY = (b3FluidGridCombinedPos)signedY * B3_FLUID_GRID_COORD_RANGE;
	b3FluidGridCombinedPos unsignedZ = (b3FluidGridCombinedPos)signedZ;
	
	return unsignedX + unsignedY + unsignedZ;
}
__kernel void convertCellValuesAndLoadCellIndex(__global b3FluidGridCombinedPos* activeCells, __global b3FluidGridValueIndexPair* yOrientedPairs, 
								__global b3FluidGridValueIndexPair* zOrientedPairs, int numActiveCells)
{
	int index = get_global_id(0);
	if(index >= numActiveCells) return;
	
	b3FluidGridCombinedPos xAxisOrientedValue = activeCells[index];
	
	b3FluidGridPosition splitPosition;
	splitCombinedPosition(B3_FLUID_GRID_COORD_RANGE, B3_FLUID_GRID_COORD_RANGE, xAxisOrientedValue, 
							&splitPosition.x, &splitPosition.y, &splitPosition.z);
	splitPosition.x -= B3_FLUID_GRID_COORD_RANGE_HALVED;
	splitPosition.y -= B3_FLUID_GRID_COORD_RANGE_HALVED;
	splitPosition.z -= B3_FLUID_GRID_COORD_RANGE_HALVED;
	
	yOrientedPairs[index].m_value = getCombinedPosition_yAxisOriented(splitPosition);
	yOrientedPairs[index].m_index = index;
	
	zOrientedPairs[index].m_value = getCombinedPosition_zAxisOriented(splitPosition);
	zOrientedPairs[index].m_index = index;
}
__kernel void writebackReorientedCellIndicies(__global b3FluidGridValueIndexPair* yOrientedPairs, __global b3FluidGridValueIndexPair* zOrientedPairs,
											__global int* yIndex, __global int* zIndex, int numActiveCells)
{
	int index = get_global_id(0);
	if(index >= numActiveCells) return;
	
	yIndex[ yOrientedPairs[index].m_index ] = index;
	zIndex[ zOrientedPairs[index].m_index ] = index;
}

__kernel void generateUniques(__global b3FluidGridValueIndexPair* sortedPairs, 
							  __global b3FluidGridCombinedPos* out_activeCells, __global b3FluidGridIterator* out_cellContents,
							  __global int* out_numActiveCells, int numSortedPairs )
{
	//Assuming that out_activeCells[] is large enough to contain
	//all active cells( out_activeCells.size() >= numSortedPairs ).

	//Iterate from sortedPairs[0] to sortedPairs[numSortedPairs-1],
	//adding unique b3FluidGridCombinedPos(s) and b3FluidGridIterator(s) to 
	//out_activeCells and out_cellContents, respectively.
	
	if( get_global_id(0) == 0 )
	{
		int numActiveCells = 0;
		
		if( numSortedPairs ) 
		{
			//Crashes on compiling with Catalyst 13.1 if
			//(b3FluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX} is used directly
			int invalidLowerIndex = INVALID_FIRST_INDEX;
			int invalidUpperIndex = INVALID_LAST_INDEX;
		
			out_activeCells[numActiveCells] = sortedPairs[0].m_value;
			//out_cellContents[numActiveCells] = (b3FluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX};
			out_cellContents[numActiveCells] = (b3FluidGridIterator){invalidLowerIndex, invalidUpperIndex};
			++numActiveCells;
			
			out_cellContents[0].m_firstIndex = 0;
			out_cellContents[0].m_lastIndex = 0;
			
			for(int i = 1; i < numSortedPairs; ++i)
			{
				if( sortedPairs[i].m_value != sortedPairs[i - 1].m_value )
				{
					out_activeCells[numActiveCells] = sortedPairs[i].m_value;
					//out_cellContents[numActiveCells] = (b3FluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX};
					out_cellContents[numActiveCells] = (b3FluidGridIterator){invalidLowerIndex, invalidUpperIndex};
					++numActiveCells;
			
					int lastIndex = numActiveCells - 1;
					out_cellContents[lastIndex].m_firstIndex = i;
					out_cellContents[lastIndex].m_lastIndex = i;
					
					//
					out_cellContents[lastIndex - 1].m_lastIndex = i - 1;
				}
			}
			
			int valuesLastIndex = numSortedPairs - 1;
			if( sortedPairs[valuesLastIndex].m_value == sortedPairs[valuesLastIndex - 1].m_value )
			{
				int uniqueLastIndex = numActiveCells - 1;
				out_cellContents[uniqueLastIndex].m_lastIndex = valuesLastIndex;
			}
		}
		
		*out_numActiveCells = numActiveCells;
	}
}

inline void binaryRangeSearch(int numActiveCells, __global b3FluidGridCombinedPos* cellValues,
							  b3FluidGridCombinedPos lowerValue, b3FluidGridCombinedPos upperValue, int* out_lowerIndex, int* out_upperIndex)
{
	int first = 0;
	int last = numActiveCells - 1;
	
	while(first <= last)
	{
//#define USE_OPTIMIZED_BINARY_RANGE_SEARCH
#ifdef USE_OPTIMIZED_BINARY_RANGE_SEARCH
		int mid = (first + last) / 2;
		int lowerValueGreaterThanMidValue = (lowerValue > cellValues[mid]);
		int upperValueLesserThanMidValue = (upperValue < cellValues[mid]);
		
		first = (lowerValueGreaterThanMidValue) ? mid + 1 : first;
		last = (!lowerValueGreaterThanMidValue && upperValueLesserThanMidValue) ? last = mid - 1 : last;
		
		if( !lowerValueGreaterThanMidValue && !upperValueLesserThanMidValue )
		{
			int lowerIndex = mid;
			int upperIndex = mid;
			
			//Perform a linear search to find the lower and upper index range
			while(lowerIndex-1 >= 0 && cellValues[lowerIndex-1] >= lowerValue) lowerIndex--;
			while(upperIndex+1 < numActiveCells && cellValues[upperIndex+1] <= upperValue) upperIndex++;
			
			*out_lowerIndex = lowerIndex;
			*out_upperIndex = upperIndex;
			return;
		}
#else
		int mid = (first + last) / 2;
		if( lowerValue > cellValues[mid] )
		{
			first = mid + 1;
		}
		else if( upperValue < cellValues[mid] )
		{
			last = mid - 1;
		}
		else 
		{
			//At this point, (lowerValue <= cellValues[mid] <= upperValue)
			//Perform a linear search to find the lower and upper index range
		
			int lowerIndex = mid;
			int upperIndex = mid;
			while(lowerIndex-1 >= 0 && cellValues[lowerIndex-1] >= lowerValue) lowerIndex--;
			while(upperIndex+1 < numActiveCells && cellValues[upperIndex+1] <= upperValue) upperIndex++;
		
			*out_lowerIndex = lowerIndex;
			*out_upperIndex = upperIndex;
			return;
		}
#endif
	}
	
	*out_lowerIndex = numActiveCells;
	*out_upperIndex = numActiveCells;
}

inline void findCellsFromGridPosition(int numActiveCells, __global b3FluidGridCombinedPos* cellValues, __global b3FluidGridIterator* cellContents, 
										b3FluidGridPosition combinedPosition, b3FluidGridIterator* out_cells)
{
	b3FluidGridPosition cellIndicies[b3FluidSortingGrid_NUM_FOUND_CELLS_GPU];
	
	b3FluidGridPosition indicies = combinedPosition;

	for(int i = 0; i < b3FluidSortingGrid_NUM_FOUND_CELLS_GPU; ++i) cellIndicies[i] = indicies;
	cellIndicies[1].y++;
	cellIndicies[2].z++;
	cellIndicies[3].y++;
	cellIndicies[3].z++;
	
	cellIndicies[4].y--;
	cellIndicies[5].z--;
	cellIndicies[6].y--;
	cellIndicies[6].z--;
	
	cellIndicies[7].y++;
	cellIndicies[7].z--;
	
	cellIndicies[8].y--;
	cellIndicies[8].z++;
	
	for(int i = 0; i < b3FluidSortingGrid_NUM_FOUND_CELLS_GPU; ++i) 
	{
		//Crashes on compiling with Catalyst 13.1 if
		//(b3FluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX} is used directly
		int invalidLowerIndex = INVALID_FIRST_INDEX;
		int invalidUpperIndex = INVALID_LAST_INDEX;
		out_cells[i] = (b3FluidGridIterator){invalidLowerIndex, invalidUpperIndex};
		//out_cells[i] = (b3FluidGridIterator){INVALID_FIRST_INDEX, INVALID_LAST_INDEX};
	}
	for(int i = 0; i < b3FluidSortingGrid_NUM_FOUND_CELLS_GPU; ++i)
	{
	
		b3FluidGridPosition lower = cellIndicies[i];
		lower.x--;
	
		b3FluidGridPosition upper = cellIndicies[i];
		upper.x++;
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(numActiveCells, cellValues, getCombinedPosition(lower), getCombinedPosition(upper), &lowerIndex, &upperIndex);
		
		if(lowerIndex != numActiveCells)
		{
			out_cells[i] = (b3FluidGridIterator){cellContents[lowerIndex].m_firstIndex, cellContents[upperIndex].m_lastIndex};
		}
	
	}
}


__kernel void findNeighborCellsPerCell( __constant int* numActiveCells, __global b3FluidGridCombinedPos* cellValues, 
										__global b3FluidGridIterator* cellContents, __global b3FluidSortingGridFoundCellsGpu* out_foundCells)
{
	int gridCellIndex = get_global_id(0);
	if(gridCellIndex >= *numActiveCells) return;
	
	b3FluidGridCombinedPos combinedPosition = cellValues[gridCellIndex];
	
	b3FluidGridPosition splitPosition;
	splitCombinedPosition(B3_FLUID_GRID_COORD_RANGE, B3_FLUID_GRID_COORD_RANGE, combinedPosition, 
							&splitPosition.x, &splitPosition.y, &splitPosition.z);
	splitPosition.x -= B3_FLUID_GRID_COORD_RANGE_HALVED;
	splitPosition.y -= B3_FLUID_GRID_COORD_RANGE_HALVED;
	splitPosition.z -= B3_FLUID_GRID_COORD_RANGE_HALVED;
	b3FluidGridIterator foundCells[b3FluidSortingGrid_NUM_FOUND_CELLS_GPU];
	findCellsFromGridPosition(*numActiveCells, cellValues, cellContents, splitPosition, foundCells);
	for(int cell = 0; cell < b3FluidSortingGrid_NUM_FOUND_CELLS_GPU; ++cell) out_foundCells[gridCellIndex].m_iterators[cell] = foundCells[cell];
}

__kernel void findGridCellIndexPerParticle(__constant int* numActiveCells, __global b3FluidGridIterator* cellContents, 
											__global int* out_gridCellIndicies)
{
	int gridCellIndex = get_global_id(0);
	if(gridCellIndex >= *numActiveCells) return;
	
	b3FluidGridIterator foundCell = cellContents[gridCellIndex];
	for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n) out_gridCellIndicies[n] = gridCellIndex;
}

//
#define B3_EPSILON FLT_EPSILON
__kernel void sphComputePressure(__constant b3FluidSphParameters* FP,
								  __global b3Vector3* fluidPosition, __global b3Scalar* fluidDensity,
								  __global b3FluidSortingGridFoundCellsGpu* foundCells, __global int* foundCellIndex, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar sum = FP->m_initialSum;
	
	for(int cell = 0; cell < b3FluidSortingGrid_NUM_FOUND_CELLS_GPU; ++cell) 
	{
		b3FluidGridIterator foundCell = foundCells[ foundCellIndex[i] ].m_iterators[cell];
		
		for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n)
		{
			b3Vector3 delta = (fluidPosition[i] - fluidPosition[n]) * FP->m_simulationScale;	//Simulation scale distance
			b3Scalar distanceSquared = b3Vector3_length2(delta);
			
			b3Scalar c = FP->m_sphRadiusSquared - distanceSquared;
			sum += (c > 0.0f && i != n) ? c*c*c : 0.0f;		//If c is positive, the particle is within interaction radius(poly6 kernel radius)
		}
	}
	
	fluidDensity[i] = sum * FP->m_sphParticleMass * FP->m_poly6KernCoeff;
}


__kernel void sphComputeForce(__constant b3FluidSphParameters* FP,
							   __global b3Vector3* fluidPosition, __global b3Vector3* fluidVelEval, 
							   __global b3Vector3* fluidSphForce, __global b3Scalar* fluidDensity,
							   __global b3FluidSortingGridFoundCellsGpu* foundCells, __global int* foundCellIndex, int numFluidParticles)
{
	b3Scalar vterm = FP->m_viscosityKernLapCoeff * FP->m_viscosity;
	
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar density_i = fluidDensity[i];
	b3Scalar invDensity_i = 1.0f / density_i;
	b3Scalar pressure_i = (density_i - FP->m_restDensity) * FP->m_stiffness;
	
	b3Vector3 force = {0.0f, 0.0f, 0.0f, 0.0f};
	
	for(int cell = 0; cell < b3FluidSortingGrid_NUM_FOUND_CELLS_GPU; ++cell) 
	{
		b3FluidGridIterator foundCell = foundCells[ foundCellIndex[i] ].m_iterators[cell];
		
		for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n)
		{	
			b3Vector3 delta = (fluidPosition[i] - fluidPosition[n]) * FP->m_simulationScale;	//Simulation scale distance
			b3Scalar distanceSquared = b3Vector3_length2(delta);
			
			if(FP->m_sphRadiusSquared > distanceSquared && i != n)
			{
				b3Scalar density_n = fluidDensity[n];
				b3Scalar invDensity_n = 1.0f / density_n;
				b3Scalar pressure_n = (density_n - FP->m_restDensity) * FP->m_stiffness;
			
				b3Scalar distance = sqrt(distanceSquared);
				b3Scalar c = FP->m_sphSmoothRadius - distance;
				b3Scalar pterm = -0.5f * c * FP->m_spikyKernGradCoeff * (pressure_i + pressure_n);
				pterm /= (distance < B3_EPSILON) ? B3_EPSILON : distance;
				
				b3Scalar dterm = c * invDensity_i * invDensity_n;
				
				force += (delta * pterm + (fluidVelEval[n] - fluidVelEval[i]) * vterm) * dterm;
			}
		}
	}
	
	fluidSphForce[i] = force * FP->m_sphParticleMass;
}

__kernel void applyForces(__constant b3FluidSphParameters* FP, 
						__global b3Vector3* fluidExternalForce, __global b3Vector3* fluidSphAcceleration,
						__global b3Vector3* fluidVel, __global b3Vector3* fluidVelEval, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Vector3 sphAcceleration = fluidSphAcceleration[i];
	{
		b3Scalar accelMagnitude = sqrt( b3Vector3_length2(sphAcceleration) );
		
		b3Scalar simulationScaleAccelLimit = FP->m_sphAccelLimit * FP->m_simulationScale;
		if(accelMagnitude > simulationScaleAccelLimit) sphAcceleration *= simulationScaleAccelLimit / accelMagnitude;
	}

	b3Vector3 acceleration = FP->m_gravity + sphAcceleration + fluidExternalForce[i] / FP->m_particleMass;
	
	b3Vector3 vel = fluidVel[i];
	
	b3Vector3 vnext = vel + acceleration * FP->m_timeStep;		//v(t+1/2) = v(t-1/2) + a(t) dt	
	fluidVel[i] = vnext;
	
	fluidExternalForce[i] = (b3Vector3){0.0f, 0.0f, 0.0f, 0.0f};
}

inline void resolveAabbCollision_impulse(__constant b3FluidSphParameters* FP, 
										b3Vector3 velocity, b3Vector3 normal, b3Scalar distance, b3Vector3* out_impulse)
{
	if( distance < 0.0f )	//Negative distance indicates penetration
	{
		b3Scalar penetratingMagnitude = b3Vector3_dot(velocity, -normal);
		if( penetratingMagnitude < 0.0f ) penetratingMagnitude = 0.0f;
		
		b3Vector3 penetratingVelocity = -normal * penetratingMagnitude;
		b3Vector3 tangentialVelocity = velocity - penetratingVelocity;
		
		penetratingVelocity *= 1.0f + FP->m_boundaryRestitution;
		
		b3Scalar positionError = (-distance) * (FP->m_simulationScale/FP->m_timeStep) * FP->m_boundaryErp;
		
		*out_impulse += -( penetratingVelocity + (-normal*positionError) + tangentialVelocity * FP->m_boundaryFriction );
	}
}
inline void accumulateBoundaryImpulse(__constant b3FluidSphParameters* FP, 
								b3Scalar simScaleParticleRadius, b3Vector3 pos, b3Vector3 vel, b3Vector3* out_impulse)
{
	b3Scalar radius = simScaleParticleRadius;
	b3Scalar simScale = FP->m_simulationScale;
	
	b3Vector3 boundaryMin = FP->m_aabbBoundaryMin;
	b3Vector3 boundaryMax = FP->m_aabbBoundaryMax;
	
	resolveAabbCollision_impulse( FP, vel, (b3Vector3){ 1.0f, 0.0f, 0.0f, 0.0f}, ( pos.x - boundaryMin.x )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FP, vel, (b3Vector3){-1.0f, 0.0f, 0.0f, 0.0f}, ( boundaryMax.x - pos.x )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FP, vel, (b3Vector3){0.0f,  1.0f, 0.0f, 0.0f}, ( pos.y - boundaryMin.y )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FP, vel, (b3Vector3){0.0f, -1.0f, 0.0f, 0.0f}, ( boundaryMax.y - pos.y )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FP, vel, (b3Vector3){0.0f, 0.0f,  1.0f, 0.0f}, ( pos.z - boundaryMin.z )*simScale - radius, out_impulse );
	resolveAabbCollision_impulse( FP, vel, (b3Vector3){0.0f, 0.0f, -1.0f, 0.0f}, ( boundaryMax.z - pos.z )*simScale - radius, out_impulse );
}
__kernel void collideAabbImpulse(__constant b3FluidSphParameters* FP, 
								__global b3Vector3* fluidPosition, __global b3Vector3* fluidVel, __global b3Vector3* fluidVelEval, 
								int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Vector3 pos = fluidPosition[i];
	b3Vector3 vel = fluidVel[i];
	
	b3Scalar simScaleParticleRadius = FP->m_particleRadius * FP->m_simulationScale;
	
	b3Vector3 aabbImpulse = (b3Vector3){0.0f, 0.0f, 0.0f, 0.0f};
	accumulateBoundaryImpulse(FP, simScaleParticleRadius, pos, vel, &aabbImpulse);
	
	//Leapfrog integration
	b3Vector3 vnext = vel + aabbImpulse;
	fluidVel[i] = vnext;
}

__kernel void integratePositions(__constant b3FluidSphParameters* FP, 
								__global b3Vector3* fluidPosition, __global b3Vector3* fluidVel, __global b3Vector3* fluidVelEval, 
								int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar timeStepDivSimScale = FP->m_timeStep / FP->m_simulationScale;
	
	b3Vector3 prevVelocity = fluidVelEval[i];	//Velocity at (t-1/2)
	b3Vector3 nextVelocity = fluidVel[i];		//Velocity at (t+1/2)
	
	if(FP->m_speedLimit != 0.0f)
	{
		b3Scalar simulationScaleSpeedLimit = FP->m_speedLimit * FP->m_simulationScale;
	
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


// ////////////////////////////////////////////////////////////////////////////
// Modulo Hash Grid
// ////////////////////////////////////////////////////////////////////////////
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

__kernel void generateValueIndexPairsModulo(__global b3Vector3* fluidPositions, __global b3FluidGridValueIndexPair* out_pairs, 
										b3Scalar cellSize, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	b3FluidGridValueIndexPair result;
	result.m_index = index;
	result.m_value = getCombinedPositionModulo( getDiscretePosition(cellSize, fluidPositions[index]) );
	
	out_pairs[index] = result;
}

__kernel void resetGridCellsModulo(__global b3FluidGridIterator* out_iterators, int numGridCells)
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
__kernel void detectIndexRangesModulo(__global b3Vector3* fluidPosition, __global b3FluidGridValueIndexPair* valueIndexPairs, 
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

__kernel void sphComputePressureModulo(__constant b3FluidSphParameters* FP,
								__global b3Vector3* fluidPosition, __global b3Scalar* fluidDensity,
								__global b3FluidGridIterator* cellContents, b3Scalar gridCellSize, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar sum = FP->m_initialSum;
	
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
					b3Vector3 delta = (fluidPosition[i] - fluidPosition[n]) * FP->m_simulationScale;	//Simulation scale distance
					b3Scalar distanceSquared = b3Vector3_length2(delta);
					
					b3Scalar c = FP->m_sphRadiusSquared - distanceSquared;
					sum += (c > 0.0f && i != n) ? c*c*c : 0.0f;		//If c is positive, the particle is within interaction radius(poly6 kernel radius)
				}
			}
	
	fluidDensity[i] = sum * FP->m_sphParticleMass * FP->m_poly6KernCoeff;
}
__kernel void sphComputeForceModulo(__constant b3FluidSphParameters* FP,
							__global b3Vector3* fluidPosition, __global b3Vector3* fluidVelEval, 
							__global b3Vector3* fluidSphForce, __global b3Scalar* fluidDensity,
							__global b3FluidGridIterator* cellContents, b3Scalar gridCellSize, int numFluidParticles)
{
	b3Scalar vterm = FP->m_viscosityKernLapCoeff * FP->m_viscosity;
	
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar density_i = fluidDensity[i];
	b3Scalar invDensity_i = 1.0f / density_i;
	b3Scalar pressure_i = (density_i - FP->m_restDensity) * FP->m_stiffness;
	
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
					b3Vector3 delta = (fluidPosition[i] - fluidPosition[n]) * FP->m_simulationScale;	//Simulation scale distance
					b3Scalar distanceSquared = b3Vector3_length2(delta);
					
					if(FP->m_sphRadiusSquared > distanceSquared && i != n)
					{
						b3Scalar density_n = fluidDensity[n];
						b3Scalar invDensity_n = 1.0f / density_n;
						b3Scalar pressure_n = (density_n - FP->m_restDensity) * FP->m_stiffness;
					
						b3Scalar distance = sqrt(distanceSquared);
						b3Scalar c = FP->m_sphSmoothRadius - distance;
						b3Scalar pterm = -0.5f * c * FP->m_spikyKernGradCoeff * (pressure_i + pressure_n);
						pterm /= (distance < B3_EPSILON) ? B3_EPSILON : distance;
						
						b3Scalar dterm = c * invDensity_i * invDensity_n;
						
						force += (delta * pterm + (fluidVelEval[n] - fluidVelEval[i]) * vterm) * dterm;
					}
				}
			}
	
	fluidSphForce[i] = force * FP->m_sphParticleMass;
}
