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
#include "b3FluidParticles.h"

#include "Bullet3Common/b3Vector3.h"

int b3FluidParticles::addParticle(const b3Vector3& position)
{
	if( size() < m_maxParticles )
	{
		m_position.push_back( b3Vector3() );
		m_velocity.push_back( b3Vector3() );
		m_velocityEval.push_back( b3Vector3() );
		m_accumulatedForce.push_back( b3Vector3() );
		m_userPointer.push_back(0);
		
		int index = size() - 1;
		
		m_position[index] = position;
		m_velocity[index].setValue(0,0,0);
		m_velocityEval[index].setValue(0,0,0);
		m_accumulatedForce[index].setValue(0,0,0);
		
		return index;
	}
	
	return size();
}
void b3FluidParticles::removeParticle(int index)
{
	b3Assert(0 <= index);
	b3Assert( index < size() );
	
	int lastIndex = size() - 1;
	
	if(index < lastIndex) 
	{
		m_position[index] = m_position[lastIndex];
		m_velocity[index] = m_velocity[lastIndex];
		m_velocityEval[index] = m_velocityEval[lastIndex];
		m_accumulatedForce[index] = m_accumulatedForce[lastIndex];
		m_userPointer[index] = m_userPointer[lastIndex];
	}
	m_position.pop_back();
	m_velocity.pop_back();
	m_velocityEval.pop_back();
	m_accumulatedForce.pop_back();
	m_userPointer.pop_back();
}

void b3FluidParticles::resize(int newSize)
{
	if(newSize > m_maxParticles) m_maxParticles = newSize;

	m_position.resize(newSize);
	m_velocity.resize(newSize);
	m_velocityEval.resize(newSize);
	m_accumulatedForce.resize(newSize);
	m_userPointer.resize(newSize);
}

void b3FluidParticles::setMaxParticles(int maxNumParticles)
{
	m_maxParticles = maxNumParticles;
	
	m_position.reserve(maxNumParticles);
	m_velocity.reserve(maxNumParticles);
	m_velocityEval.reserve(maxNumParticles);
	m_accumulatedForce.reserve(maxNumParticles);
	m_userPointer.reserve(maxNumParticles);
}
