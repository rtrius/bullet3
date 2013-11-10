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
#ifndef B3_FLUID_SORTING_GRID_H
#define B3_FLUID_SORTING_GRID_H

#include "Bullet3Common/b3Vector3.h"
#include "Bullet3Common/b3AlignedObjectArray.h"
#include "Bullet3Common/b3Logging.h"

#include "b3FluidParticles.h"

class b3Vector3;
struct b3FluidParticles;


///@brief Used to iterate through all particles in a single b3FluidSortingGrid cell.
///@remarks
///The standard method for iterating through a grid cell is:
///@code
///		b3FluidGridIterator FI;
///		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
///@endcode
struct b3FluidGridIterator
{
	int m_firstIndex;
	int m_lastIndex;
	
	b3FluidGridIterator() {}
	b3FluidGridIterator(int firstIndex, int lastIndex) : m_firstIndex(firstIndex), m_lastIndex(lastIndex) {}
};

//B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT may be disabled in order to use the faster OpenCL grid update
//in order to do this, '#define B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT' must be commented out in 3 places:
//	b3FluidSortingGrid.h,
//	fluidSph.cl,
//	fluidSphCL.h (or re-stringify from fluidSph.cl)
//
//#define B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT	//Ensure that this is also #defined in "fluidSph.cl"
#ifdef B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT
	typedef unsigned long long int b3FluidGridUint64;
	typedef b3FluidGridUint64 b3FluidGridCombinedPos;					//Range must contain B3_FLUID_GRID_COORD_RANGE^3
	const b3FluidGridCombinedPos B3_FLUID_GRID_COORD_RANGE = 2097152;	//2^21
#else
	typedef unsigned int b3FluidGridCombinedPos;						//Range must contain B3_FLUID_GRID_COORD_RANGE^3
	const b3FluidGridCombinedPos B3_FLUID_GRID_COORD_RANGE = 1024;		//2^10
#endif

typedef int b3FluidGridCoordinate;
const b3FluidGridCoordinate B3_FLUID_GRID_COORD_RANGE_HALVED = B3_FLUID_GRID_COORD_RANGE/2;

///For sorting; contains a b3FluidSortingGrid grid cell id and fluid particle index.
struct b3FluidGridValueIndexPair
{
	b3FluidGridCombinedPos m_value;		///<Grid cell id
	int m_index;						///<Fluid particle index
	
	b3FluidGridValueIndexPair() {}
	b3FluidGridValueIndexPair(b3FluidGridCombinedPos value, int index) : m_value(value), m_index(index) {}
};

///@brief Contains a world scale position quantized to units of b3FluidSortingGrid.m_gridCellSize.
struct b3FluidGridPosition
{
	b3FluidGridCoordinate x;		
	b3FluidGridCoordinate y;
	b3FluidGridCoordinate z;
private:
	b3FluidGridCoordinate padding;
	
public:
	bool operator==(const b3FluidGridPosition& GI) const { return (x == GI.x && y == GI.y && z == GI.z); }
	bool operator!=(const b3FluidGridPosition& GI) const { return (x != GI.x || y != GI.y || z != GI.z); }
	
	bool operator>(const b3FluidGridPosition& GI) const
	{
		if(z != GI.z) return (z > GI.z);
		if(y != GI.y) return (y > GI.y);
		return (x > GI.x);
	}	
	bool operator<(const b3FluidGridPosition& GI) const
	{
		if(z != GI.z) return (z < GI.z);
		if(y != GI.y) return (y < GI.y);
		return (x < GI.x);
	}
	
	b3FluidGridCombinedPos getCombinedPosition() const
	{
		//Convert range 
		//from [-B3_FLUID_GRID_COORD_RANGE_HALVED, B3_FLUID_GRID_COORD_RANGE_HALVED - 1] 
		//  to [0, B3_FLUID_GRID_COORD_RANGE - 1] before combining
		//e.g. from [-512, 511] to [0, 1023] before combining	
		b3FluidGridCoordinate signedX = x + B3_FLUID_GRID_COORD_RANGE_HALVED;
		b3FluidGridCoordinate signedY = y + B3_FLUID_GRID_COORD_RANGE_HALVED;
		b3FluidGridCoordinate signedZ = z + B3_FLUID_GRID_COORD_RANGE_HALVED;
		
		b3FluidGridCombinedPos unsignedX = static_cast<b3FluidGridCombinedPos>(signedX);
		b3FluidGridCombinedPos unsignedY = static_cast<b3FluidGridCombinedPos>(signedY) * B3_FLUID_GRID_COORD_RANGE;
		b3FluidGridCombinedPos unsignedZ = static_cast<b3FluidGridCombinedPos>(signedZ) * B3_FLUID_GRID_COORD_RANGE * B3_FLUID_GRID_COORD_RANGE;
		
		return unsignedX + unsignedY + unsignedZ;
	}
};



///@brief Uniform grid broadphase for b3FluidSph particles.
///@remarks
///A fundamental operation in SPH fluid simulations is the detection
///of collisions between fluid particles.
///@par
///Since testing each fluid pair would require O(n^2) operations, a grid based 
///broadphase is implemented to accelerate the search to O(kn), where k
///is the average number of particles in each cell. Each grid cell has a size 
///of r, where r is the SPH interaction radius at world scale. When a particle 
///is queried, it searches a 3x3x3 grid cell volume surrounding its position.
///Note that, for each particle, a 'collision' is detected if the center of
///other particles is within r; that is, the effective radius of collision,
///were all particles to be treated as spheres(and not as points), is r/2.
///@par
///Particles are stored as a set of index ranges. First, the world scale positions
///of particles are quantized into grid cell coordinates(b3FluidGridPosition). Those
///integer coordinates are then converted into single values that are used for sorting
///(b3FluidGridCombinedPos). After sorting the particles by the grid cell they are contained in, 
///the lower and upper indicies of each nonempty cell is detected and stored(b3FluidGridIterator).
///@par
///Effective size: B3_FLUID_GRID_COORD_RANGE^3, which is currently 1024^3 
///or 2^21^3(with #define B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT) grid cells.
///Worlds larger than 2^21^3 are unsupported.
class b3FluidSortingGrid
{
	//INVALID_LAST_INDEX must be lower than INVALID_FIRST_INDEX,
	//such that the below loop will not execute.
	//		for(int i = FI.m_firstIndex; i <= FI.m_lastIndex; ++i)
	static const int INVALID_FIRST_INDEX = -1;
	static const int INVALID_LAST_INDEX = INVALID_FIRST_INDEX - 1;
public:
	static const int NUM_MULTITHREADING_GROUPS = 27; 	///<Number of grid cells that may be accessed when iterating through a single grid cell
	static const int NUM_FOUND_CELLS = 27;				///<Number of grid cells returned from b3FluidSortingGrid::findCells()
	static const int NUM_FOUND_CELLS_SYMMETRIC = 14;	///<Number of grid cells returned from b3FluidSortingGrid::findCellsSymmetric()
	static const int NUM_FOUND_CELLS_GPU = 9;			///<OpenCL solver represents 27 cells as 9 3-cell bars
	
	struct FoundCells { b3FluidGridIterator m_iterators[b3FluidSortingGrid::NUM_FOUND_CELLS]; }; ///<Contains results of b3FluidSortingGrid::findCells()
	struct FoundCellsGpu { b3FluidGridIterator m_iterators[b3FluidSortingGrid::NUM_FOUND_CELLS_GPU]; };
	
private:
	b3Vector3 m_pointMin;	//AABB calculated from the center of fluid particles, without considering particle radius
	b3Vector3 m_pointMax;
	
	///Each array contains a set of grid cell indicies that may be simultaneously processed
	b3AlignedObjectArray<int> m_multithreadingGroups[b3FluidSortingGrid::NUM_MULTITHREADING_GROUPS];
	
	b3Scalar m_gridCellSize;

	b3AlignedObjectArray<b3FluidGridCombinedPos> m_activeCells;		//Stores the value of each nonempty grid cell
	b3AlignedObjectArray<b3FluidGridIterator> m_cellContents;	//Stores the range of indicies that correspond to the values in m_activeCells
	
	b3AlignedObjectArray<b3FluidGridValueIndexPair> m_tempPairs;
	b3AlignedObjectArray<b3Vector3> m_tempBufferVector;
	b3AlignedObjectArray<void*> m_tempBufferVoid;
	
public:
	b3FluidSortingGrid()
	{
		m_pointMin.setValue(0,0,0);
		m_pointMax.setValue(0,0,0); 
		m_gridCellSize = b3Scalar(1.0);
	}

	void insertParticles(b3FluidParticles& fluids);
	
	void clear() 
	{ 
		m_activeCells.resize(0);
		m_cellContents.resize(0);
		for(int i = 0; i < b3FluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++i) m_multithreadingGroups[i].resize(0);
	}
	
	///Returns a 3x3x3 group of b3FluidGridIterator, which is the maximum extent of cells
	///that may interact with an AABB defined by (position - radius, position + radius). 
	///Where radius is the SPH smoothing radius, in b3FluidSphParametersGlobal, converted to world scale.
	///@param position Center of the AABB defined by (position - radius, position + radius).
	void findCells(const b3Vector3& position, b3FluidSortingGrid::FoundCells& out_gridCells) const
	{
		findAdjacentGridCells( getDiscretePosition(position), out_gridCells );
	}
	
	///Returns 14 grid cells, with out_gridCells->m_iterator[0] as the center cell corresponding to position.
	void findCellsSymmetric(const b3Vector3& position, b3FluidSortingGrid::FoundCells& out_gridCells) const
	{
		findAdjacentGridCellsSymmetric( getDiscretePosition(position), out_gridCells );
	}
	
	int getNumGridCells() const { return m_activeCells.size(); }	///<Returns the number of nonempty grid cells.
	b3FluidGridIterator getGridCell(int gridCellIndex) const { return m_cellContents[gridCellIndex]; }
	
	struct AabbCallback
	{
		AabbCallback() {}
		virtual ~AabbCallback() {}
		
		///b3FluidSortingGrid::forEachGridCell() will continue calling processParticles() if this returns true
		virtual bool processParticles(const b3FluidGridIterator FI, const b3Vector3& aabbMin, const b3Vector3& aabbMax) = 0;
	};
	void forEachGridCell(const b3Vector3& aabbMin, const b3Vector3& aabbMax, b3FluidSortingGrid::AabbCallback& callback) const;

	b3Scalar getCellSize() const { return m_gridCellSize; }
	void setCellSize(b3Scalar simulationScale, b3Scalar sphSmoothRadius) 
	{
		m_gridCellSize = sphSmoothRadius / simulationScale; 	//Divide by simulationScale to convert to world scale
	}
	
	///Returns the AABB calculated from the center of each fluid particle in the grid, without considering particle radius
	void getPointAabb(b3Vector3& out_pointMin, b3Vector3& out_pointMax) const
	{
		out_pointMin = m_pointMin;
		out_pointMax = m_pointMax;
	}
	
	b3AlignedObjectArray<b3FluidGridCombinedPos>& internalGetActiveCells() { return m_activeCells; }
	b3AlignedObjectArray<b3FluidGridIterator>& internalGetCellContents() { return m_cellContents; }
	
	const b3AlignedObjectArray<int>& internalGetMultithreadingGroup(int index) const { return m_multithreadingGroups[index]; }
	
	b3Vector3& internalGetPointAabbMin() { return m_pointMin; }
	b3Vector3& internalGetPointAabbMax() { return m_pointMax; }
	
private:
	b3FluidGridPosition getDiscretePosition(const b3Vector3& position) const;

	void findAdjacentGridCells(b3FluidGridPosition indicies, b3FluidSortingGrid::FoundCells& out_gridCells) const;
	void findAdjacentGridCellsSymmetric(b3FluidGridPosition indicies, b3FluidSortingGrid::FoundCells& out_gridCells) const;
	
	void generateMultithreadingGroups();
};

#endif



