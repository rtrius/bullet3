/*
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef B3_STATE_RIGID_CONTACTS
#define B3_STATE_RIGID_CONTACTS

#include "Bullet3Collision/NarrowPhaseCollision/b3Contact4.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLInclude.h"
#include "Bullet3OpenCL/ParallelPrimitives/b3OpenCLArray.h"

struct b3StateRigidContacts
{
	int m_maxContacts;

	b3OpenCLArray<b3Contact4>* m_pContactBuffersGPU[2];
	int	m_currentContactBuffer;

	b3AlignedObjectArray<b3Contact4> m_contactOutCPU;


	b3StateRigidContacts(cl_context context, cl_command_queue queue, int maxContacts);
	virtual ~b3StateRigidContacts();


	const struct b3Contact4* getContactsCPU();

	cl_mem	getContactsGpu();
	int	getNumContactsGpu() const;

};

#endif

