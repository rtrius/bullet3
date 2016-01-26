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
#include "b3StateRigidContacts.h"

b3StateRigidContacts::b3StateRigidContacts(cl_context context, cl_command_queue queue, int maxContacts) :
	m_maxContacts(maxContacts),
	m_currentContactBuffer(0)
{
	m_contactOutCPU.resize(maxContacts);

	m_pContactBuffersGPU[0] = new b3OpenCLArray<b3Contact4>(context, queue, maxContacts, true);
	m_pContactBuffersGPU[1] = new b3OpenCLArray<b3Contact4>(context, queue, maxContacts, true);
}
b3StateRigidContacts::~b3StateRigidContacts()
{
	delete m_pContactBuffersGPU[0];
	delete m_pContactBuffersGPU[1];
}

int	b3StateRigidContacts::getNumContactsGpu() const
{
	return m_pContactBuffersGPU[m_currentContactBuffer]->size();
}
cl_mem b3StateRigidContacts::getContactsGpu()
{
	return m_pContactBuffersGPU[m_currentContactBuffer]->getBufferCL();
}

const b3Contact4* b3StateRigidContacts::getContactsCPU()
{
	m_pContactBuffersGPU[m_currentContactBuffer]->copyToHost(m_contactOutCPU);
	return &m_contactOutCPU[0];
}