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
#include "b3FluidSortingGrid.h"

#include "Bullet3Common/b3Vector3.h"
#include "Bullet3Common/b3Quickprof.h"		//B3_PROFILE(name) macro

#include "b3FluidParticles.h"

#define SWAP_REARRANGED_ARRAY
#ifdef SWAP_REARRANGED_ARRAY
#include <cstring>	//memcpy()
#endif

template<typename T>
void rearrangeToMatchSortedValues(const b3AlignedObjectArray<b3FluidGridValueIndexPair>& sortedValues, 
									b3AlignedObjectArray<T>& temp, b3AlignedObjectArray<T>& out_rearranged)
{
	{
		B3_PROFILE("Resize");
		temp.resize( sortedValues.size() );
	}
	
	{
		B3_PROFILE("Rearrange");
		for(int i = 0; i < sortedValues.size(); ++i)
		{
			int oldIndex = sortedValues[i].m_index;
			int newIndex = i;
				
			temp[newIndex] = out_rearranged[oldIndex];
		}
	}
	
	{
		B3_PROFILE("Copy back");
	
#ifndef SWAP_REARRANGED_ARRAY
		out_rearranged = temp;
#else
		const int SIZEOF_ARRAY = sizeof(b3AlignedObjectArray<T>);
	
		char swap[SIZEOF_ARRAY];
		
		memcpy(swap, &temp, SIZEOF_ARRAY);
		memcpy(&temp, &out_rearranged, SIZEOF_ARRAY);
		memcpy(&out_rearranged, swap, SIZEOF_ARRAY);
#endif
	}
}

struct ValueIndexPairSortPredicate 
{
	inline bool operator() (const b3FluidGridValueIndexPair& a, const b3FluidGridValueIndexPair& b) const 
	{
		return (a.m_value < b.m_value);
	}
};
void sortParticlesByValues(b3FluidParticles& particles, b3AlignedObjectArray<b3FluidGridValueIndexPair>& values,
							 b3AlignedObjectArray<b3Vector3>& tempVector, b3AlignedObjectArray<void*>& tempVoid)
{
	{
		B3_PROFILE("sortParticlesByValues() - quickSort");
		
		values.quickSort( ValueIndexPairSortPredicate() );
	}
	
	{
		B3_PROFILE("sortParticlesByValues() - move data");
		
		rearrangeToMatchSortedValues(values, tempVector, particles.m_pos);
		rearrangeToMatchSortedValues(values, tempVector, particles.m_vel);
		rearrangeToMatchSortedValues(values, tempVector, particles.m_vel_eval);
		rearrangeToMatchSortedValues(values, tempVector, particles.m_accumulatedForce);
		rearrangeToMatchSortedValues(values, tempVoid, particles.m_userPointer);
	}
}


void b3FluidSortingGrid::insertParticles(b3FluidParticles& particles)
{
	{
		B3_PROFILE("b3FluidSortingGrid() - generate");
		m_tempPairs.resize( particles.size() );
		
		if( particles.size() )
		{
			m_pointMin.setValue(B3_LARGE_FLOAT, B3_LARGE_FLOAT, B3_LARGE_FLOAT);
			m_pointMax.setValue(-B3_LARGE_FLOAT, -B3_LARGE_FLOAT, -B3_LARGE_FLOAT);
			for(int i = 0; i < particles.size(); ++i) 
			{
				const b3Vector3& position = particles.m_pos[i];
			
				m_pointMin.setMin(position);
				m_pointMax.setMax(position);
			
				b3FluidGridPosition indicies = getDiscretePosition(position);
				m_tempPairs[i] = b3FluidGridValueIndexPair( indicies.getCombinedPosition(), i );
			}
		}
		else
		{
			m_pointMin.setValue(0, 0, 0);
			m_pointMax.setValue(0, 0, 0);
		}
	}
	
	
	//Sort fluidSystem and values by m_value(s) in m_tempPairs
	{
		B3_PROFILE("b3FluidSortingGrid() - sort");
		sortParticlesByValues(particles, m_tempPairs, m_tempBufferVector, m_tempBufferVoid);
	}
	
	m_activeCells.resize(0);
	m_cellContents.resize(0);
	{
		B3_PROFILE("b3FluidSortingGrid() - find unique");
		
		//Iterate through m_tempPairs to find the unique b3FluidGridCombinedPos(s),
		//and the index ranges(fluids[] index) at which each value appears
		if( m_tempPairs.size() ) 
		{
			m_activeCells.push_back( m_tempPairs[0].m_value );
			m_cellContents.push_back( b3FluidGridIterator() );
			m_cellContents[0].m_firstIndex = 0;
			m_cellContents[0].m_lastIndex = 0;
			
			for(int i = 1; i < m_tempPairs.size(); ++i)
			{
				if( m_tempPairs[i].m_value != m_tempPairs[i - 1].m_value )
				{
					m_activeCells.push_back( m_tempPairs[i].m_value );
					m_cellContents.push_back( b3FluidGridIterator() );
					
					int lastIndex = m_cellContents.size() - 1;
					m_cellContents[lastIndex].m_firstIndex = i;
					m_cellContents[lastIndex].m_lastIndex = i;
					
					//
					m_cellContents[lastIndex - 1].m_lastIndex = i - 1;
				}
			}
			
			int valuesLastIndex = m_tempPairs.size() - 1;
			if( m_tempPairs[valuesLastIndex].m_value == m_tempPairs[valuesLastIndex - 1].m_value )
			{
				int uniqueLastIndex = m_cellContents.size() - 1;
				m_cellContents[uniqueLastIndex].m_lastIndex = valuesLastIndex;
			}
		}
	}
	
	generateMultithreadingGroups();
}

//Based on b3AlignedObjectArray::findBinarySearch()
//instead of finding a single value, it finds a range of values
//and returns the index range containing that value range.
//Assumes that cellValues is sorted in ascending order.
//Returns cellValues.size() on failure.
void binaryRangeSearch(const b3AlignedObjectArray<b3FluidGridCombinedPos>& cellValues,
					   b3FluidGridCombinedPos lowerValue, b3FluidGridCombinedPos upperValue, int& out_lowerIndex, int& out_upperIndex)
{
	int first = 0;
	int last = cellValues.size() - 1;
	
	while(first <= last)
	{
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
			while(upperIndex+1 < cellValues.size() && cellValues[upperIndex+1] <= upperValue) upperIndex++;
		
			out_lowerIndex = lowerIndex;
			out_upperIndex = upperIndex;
			return;
		}
	}

	out_lowerIndex = cellValues.size();
	out_upperIndex = cellValues.size();
}
void b3FluidSortingGrid::forEachGridCell(const b3Vector3& aabbMin, const b3Vector3& aabbMax, b3FluidSortingGrid::AabbCallback& callback) const
{
	b3FluidGridPosition minIndicies = getDiscretePosition(aabbMin);
	b3FluidGridPosition maxIndicies = getDiscretePosition(aabbMax);

	for(b3FluidGridCoordinate z = minIndicies.z; z <= maxIndicies.z; ++z)
		for(b3FluidGridCoordinate y = minIndicies.y; y <= maxIndicies.y; ++y)
		{
			const bool USE_BINARY_RANGE_SEARCH = true;
			if(USE_BINARY_RANGE_SEARCH)
			{
				b3FluidGridPosition lower;
				lower.x = minIndicies.x;
				lower.y = y;
				lower.z = z;
					
				b3FluidGridPosition upper;
				upper.x = maxIndicies.x;
				upper.y = y;
				upper.z = z;
				
				int lowerIndex, upperIndex;
				binaryRangeSearch( m_activeCells, lower.getCombinedPosition(), upper.getCombinedPosition(), lowerIndex, upperIndex );
			
				if( lowerIndex != m_activeCells.size() && upperIndex != m_activeCells.size() )
				{
					b3FluidGridIterator FI(m_cellContents[lowerIndex].m_firstIndex, m_cellContents[upperIndex].m_lastIndex);
					if( !callback.processParticles(FI, aabbMin, aabbMax) ) return;
				}
			}
			else
			{
				for(b3FluidGridCoordinate x = minIndicies.x; x <= maxIndicies.x; ++x)
				{
					b3FluidGridPosition current;
					current.x = x;
					current.y = y;
					current.z = z;
				
					//findBinarySearch() returns m_activeCells.size() on failure
					int gridCellIndex = m_activeCells.findBinarySearch( current.getCombinedPosition() );
					if( gridCellIndex != m_activeCells.size() )
					{
						if( !callback.processParticles(m_cellContents[gridCellIndex], aabbMin, aabbMax) ) return;
					}
				}
			}
		}
}

b3FluidGridPosition b3FluidSortingGrid::getDiscretePosition(const b3Vector3& position) const
{
	//Using only 'result.x = static_cast<b3FluidGridCoordinate>( position.getX() / m_gridCellSize )'
	//would cause positions in (-m_gridCellSize, m_gridCellSize) to convert to 0;
	//that is, cell (0,0,0) would be twice as large as desired.
	//
	//To resolve this, define the indicies such that:
	//[0, m_gridCellSize) converts to 0
	//[-m_gridCellSize, 0) converts to -1
	
	b3Vector3 discretePosition = position / m_gridCellSize;
	
	//Worlds larger than 2^21^3 cells are unsupported
	const b3Scalar MIN = static_cast<b3Scalar>(-B3_FLUID_GRID_COORD_RANGE_HALVED);
	const b3Scalar MAX = static_cast<b3Scalar>(B3_FLUID_GRID_COORD_RANGE_HALVED - 1);
	b3Assert( MIN <= discretePosition.getX() && discretePosition.getX() <= MAX );
	b3Assert( MIN <= discretePosition.getY() && discretePosition.getY() <= MAX );
	b3Assert( MIN <= discretePosition.getZ() && discretePosition.getZ() <= MAX );
	
	b3FluidGridPosition result;
	result.x = static_cast<b3FluidGridCoordinate>( (position.getX() >= b3Scalar(0.0)) ? discretePosition.getX() : floor(discretePosition.getX()) );
	result.y = static_cast<b3FluidGridCoordinate>( (position.getY() >= b3Scalar(0.0)) ? discretePosition.getY() : floor(discretePosition.getY()) );
	result.z = static_cast<b3FluidGridCoordinate>( (position.getZ() >= b3Scalar(0.0)) ? discretePosition.getZ() : floor(discretePosition.getZ()) );
	
	return result;
}

void b3FluidSortingGrid::findAdjacentGridCells(b3FluidGridPosition indicies, b3FluidSortingGrid::FoundCells& out_gridCells) const
{	
	const b3FluidGridIterator INVALID_ITERATOR(INVALID_FIRST_INDEX, INVALID_LAST_INDEX);
	
	b3FluidGridPosition cellIndicies[b3FluidSortingGrid::NUM_FOUND_CELLS];
	
	for(int i = 0; i < 9; ++i) cellIndicies[i] = indicies;
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

	for(int i = 0; i < b3FluidSortingGrid::NUM_FOUND_CELLS; ++i) out_gridCells.m_iterators[i] = INVALID_ITERATOR;
	for(int i = 0; i < 9; ++i)
	{
		b3FluidGridPosition lower = cellIndicies[i];
		lower.x--;
		
		b3FluidGridPosition upper = cellIndicies[i];
		upper.x++;
			
		int lowerIndex, upperIndex;
		binaryRangeSearch( m_activeCells, lower.getCombinedPosition(), upper.getCombinedPosition(), lowerIndex, upperIndex );
		
		if( lowerIndex != m_activeCells.size() )
		{
			int range = upperIndex - lowerIndex;
			
			for(int n = 0; n <= range; ++n) out_gridCells.m_iterators[i*3 + n] = m_cellContents[lowerIndex + n];
		}
	}
}

void b3FluidSortingGrid::findAdjacentGridCellsSymmetric(b3FluidGridPosition indicies, b3FluidSortingGrid::FoundCells& out_gridCells) const
{
	const b3FluidGridIterator INVALID_ITERATOR(INVALID_FIRST_INDEX, INVALID_LAST_INDEX);
	
	//Only 14 of 27 cells need to be checked due to symmetry
	//out_gridCells.m_iterators[0] must correspond to the center grid cell
	//
	//6 binary ranges(extended along the x-axis)
	//Upper: 5 cells, 2 binary ranges(1 3-cell bar; 1 2-cell bar)
	//Center: 5 cells, 2 binary ranges(1 3-cell bar; 1 2-cell bar)
	//Lower: 4 cells, 2 binary ranges(1 3-cell bar; 1 cell)
	//
	//	C = Checked
	//	N = Not checked
	//
	//	The cell marked as '[C]' corresponds to the parameter 'b3FluidGridPosition indicies',
	//	which is the center cell from which the other cells are calculated.
	//
	//	Upper(Y+)           Center              Lower(Y-)
	//	C    C    C			C    C    C			C    C    C			<--- 3 3-cell bars on this line
	//	C    C    N			C   [C]   N			C    N    N			<--- 2 2-cell bars + 1 cell on this line
	//	N    N    N			N    N    N			N    N    N
	//
	//		       Z-
	// 		       |
	//		X- <---|---> X+
	//		       |
	//		       Z+
	//
	for(int i = 0; i < b3FluidSortingGrid::NUM_FOUND_CELLS; ++i) out_gridCells.m_iterators[i] = INVALID_ITERATOR;

	b3FluidGridPosition centers[6];	//Center cells of the 6 bars
	for(int i = 0; i < 6; ++i) centers[i] = indicies;
	
	//centers[0] and centers[1] contain 2-cell bars
	centers[1].y++;
	
	//centers[2] to centers[4] contain 3-cell bars
	centers[2].z--;
	
	centers[3].z--;
	centers[3].y++;
	
	centers[4].z--;
	centers[4].y--;
	
	//centers[5] contains 1 cell
	centers[5].y--;
	centers[5].x--;
	
	//centers[0]
	{
		b3FluidGridPosition lower = centers[0];
		lower.x--;
	
		b3FluidGridPosition upper = centers[0];
		b3FluidGridCombinedPos centerValue = centers[0].getCombinedPosition();
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(m_activeCells, lower.getCombinedPosition(), upper.getCombinedPosition(), lowerIndex, upperIndex);
		if( lowerIndex != m_activeCells.size() )
		{
			//out_gridCells.m_iterators[0] must be the center grid cell if it exists, and INVALID_ITERATOR otherwise
			if(lowerIndex != upperIndex) 
			{
				out_gridCells.m_iterators[0] = m_cellContents[upperIndex];
				out_gridCells.m_iterators[1] = m_cellContents[lowerIndex];
			}
			else
			{
				if(centerValue == m_activeCells[lowerIndex])out_gridCells.m_iterators[0] = m_cellContents[lowerIndex];
				else out_gridCells.m_iterators[1] = m_cellContents[lowerIndex];
			}
		}
	}
	
	//centers[1]
	{
		b3FluidGridPosition lower = centers[1];
		lower.x--;
	
		b3FluidGridPosition upper = centers[1];
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(m_activeCells, lower.getCombinedPosition(), upper.getCombinedPosition(), lowerIndex, upperIndex);
		if( lowerIndex != m_activeCells.size() )
		{
			out_gridCells.m_iterators[2] = m_cellContents[lowerIndex];
			if(lowerIndex != upperIndex) out_gridCells.m_iterators[3] = m_cellContents[upperIndex];
		}
	}
	
	//centers[2-4]
	for(int i = 0; i < 3; ++i)
	{
		b3FluidGridPosition lower = centers[i+2];
		lower.x--;
	
		b3FluidGridPosition upper = centers[i+2];
		upper.x++;
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(m_activeCells, lower.getCombinedPosition(), upper.getCombinedPosition(), lowerIndex, upperIndex);
		
		if( lowerIndex != m_activeCells.size() )
		{
			int range = upperIndex - lowerIndex;
			
			for(int n = 0; n <= range; ++n) out_gridCells.m_iterators[4 + i*3 + n] = m_cellContents[lowerIndex + n];
		}
	}
	
	//centers[5]
	{
		b3FluidGridPosition lowerAndUpper = centers[5];
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(m_activeCells, lowerAndUpper.getCombinedPosition(), lowerAndUpper.getCombinedPosition(), lowerIndex, upperIndex);
		if( lowerIndex != m_activeCells.size() ) out_gridCells.m_iterators[13] = m_cellContents[lowerIndex];
	}
}


#ifndef B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT	//#defined in "b3FluidSortingGrid.h"
	//Extracts the (x,y,z) indicies from value, where
	//value == x + y*resolutionX + z*resolutionX*resolutionY
	static void splitCombinedPosition(int resolutionX, int resolutionY, int value, int& out_x, int& out_y, int& out_z)
	{
		int x = value % resolutionX;
		int z = value / (resolutionX*resolutionY);
		int y = (value - z*resolutionX*resolutionY) / resolutionX;
				
		out_x = x;
		out_z = z;
		out_y = y;
	}
#else
	static void splitCombinedPosition(unsigned long long int resolutionX, unsigned long long int resolutionY, 
									unsigned long long int value, int& out_x, int& out_y, int& out_z)
	{
		unsigned long long int cellsPerLine = resolutionX;
		unsigned long long int cellsPerPlane = resolutionX * resolutionY;
											 
		unsigned long long int x = value % cellsPerLine;
		unsigned long long int z = value / cellsPerPlane;
		unsigned long long int y = (value - z*cellsPerPlane) / cellsPerLine;
				
		out_x = static_cast<int>(x);
		out_z = static_cast<int>(z);
		out_y = static_cast<int>(y);
	}
#endif
void b3FluidSortingGrid::generateMultithreadingGroups()
{	
	//Processing particle-particle interactions in a single grid cell may access 
	//up to 3^3 grid cells if the interaction is symmetric(that is, if both the 
	//particle in the center cell and the particle in a neighbor cell are written to).
	
	//In order to prevent a thread from accessing a cell that may be written to by 
	//another thread, the cells are grouped into cell blocks of 27 cells each.
	//Multiple threads cannot collide if all of the cells at the same position in
	//all blocks are operated on simultaneously. Overall, this method results in
	//27 sequentially processed groups of cells, where the cells in each group 
	//may be split among multiple threads.
	
	B3_PROFILE("generateMultithreadingGroups()");
	
	for(int i = 0; i < b3FluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++i) m_multithreadingGroups[i].resize(0);
		
	for(int cell = 0; cell < getNumGridCells(); ++cell)
	{
		b3FluidGridIterator FI = getGridCell(cell);
		if( !(FI.m_firstIndex <= FI.m_lastIndex) ) continue;
	
		int index_x, index_y, index_z;
		splitCombinedPosition(B3_FLUID_GRID_COORD_RANGE, B3_FLUID_GRID_COORD_RANGE, m_activeCells[cell], index_x, index_y, index_z);

		//For b3FluidSortingGrid: Convert range from [0, 1023] to [1, 1024]
		//(b3FluidGridPosition::getCombinedPosition() already converts from [-512, 511], to [0, 1023])
		index_x += 1;
		index_y += 1;
		index_z += 1;
		
		//For each dimension, place indicies into one of 3 categories such that
		//indicies (1, 2, 3, 4, 5, 6, ...) correspond to categories (1, 2, 3, 1, 2, 3, ...)
		int group = 0;
		
		if(index_x % 3 == 0) group += 0;
		else if(index_x % 2 == 0) group += 1;
		else group += 2;
		
		if(index_y % 3 == 0) group += 0;
		else if(index_y % 2 == 0) group += 3;
		else group += 6;
		
		if(index_z % 3 == 0) group += 0;
		else if(index_z % 2 == 0) group += 9;
		else group += 18;
		
		m_multithreadingGroups[group].push_back(cell);
	}
}
