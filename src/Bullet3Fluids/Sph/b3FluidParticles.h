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
#ifndef B3_FLUID_PARTICLES_H
#define B3_FLUID_PARTICLES_H

#include "Bullet3Common/b3AlignedObjectArray.h"

class b3Vector3;


///@brief Coordinates the parallel arrays used to store fluid particles.
///@remarks
///Members of this struct should not be accessed directly, except for calling:
/// - b3AlignedObjectArray::operator[]() 
/// - b3AlignedObjectArray::at()
///@par
///If using a GPU solver, make sure that the b3FluidSph::getGpuSyncFlags() are set correctly
///or these arrays will not be updated.
struct b3FluidParticles
{
	int m_maxParticles;
	
	//Parallel arrays
	b3AlignedObjectArray<b3Vector3> m_position;				///<World scale.
	b3AlignedObjectArray<b3Vector3> m_velocity;				///<Simulation scale; actually the velocity at (t+1/2t) for leapfrog integration.
	b3AlignedObjectArray<b3Vector3> m_velocityEval;			///<Simulation scale; average of the current(t+1/2t) and last frame(t-1/2t) m_velocity. 
	b3AlignedObjectArray<void*> m_userPointer;
	
	b3AlignedObjectArray<b3Vector3> m_accumulatedForce;		///<Applied during stepSimulation(), then set to 0; simulation scale.
	
	b3FluidParticles() : m_maxParticles(0) {}
	
	int	size() const	{ return m_position.size(); }

	int addParticle(const b3Vector3& position);		///<Returns size() if size() == getMaxParticles().
	void removeParticle(int index);					///<Swaps indicies if index does not correspond to the last index; invalidates grid.
	
	void resize(int newSize);						///<Does not initialize particles if( newSize > size() ).
	
	void setMaxParticles(int maxNumParticles);
	int getMaxParticles() const { return m_maxParticles; }
};


#endif


