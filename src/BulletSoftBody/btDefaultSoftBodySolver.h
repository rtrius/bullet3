/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  http://continuousphysics.com/Bullet/

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef BT_SOFT_BODY_DEFAULT_SOLVER_H
#define BT_SOFT_BODY_DEFAULT_SOLVER_H


#include "BulletSoftBody/btSoftBodySolvers.h"
#include "btSoftBodySolverVertexBuffer.h"
struct btCollisionObjectWrapper;

class btDefaultSoftBodySolver : public btSoftBodySolver
{
protected:
	btAlignedObjectArray<btSoftBody*> m_softBodySet;


public:
	btDefaultSoftBodySolver();
	virtual ~btDefaultSoftBodySolver();

	virtual void updateSoftBodies( );

	virtual void optimize( btAlignedObjectArray< btSoftBody * > &softBodies,bool forceUpdate=false );

	virtual void solveConstraints( float solverdt );

	virtual void predictMotion( float solverdt );

	static void VSolve_Links(btSoftBody* psb,btScalar kst);
	static void PSolve_Anchors(btSoftBody* psb,btScalar kst,btScalar ti);
	static void PSolve_RigidContacts(btSoftBody* psb,btScalar kst,btScalar ti);
	static void PSolve_SoftContacts(btSoftBody* psb,btScalar,btScalar ti);
	static void PSolve_Links(btSoftBody* psb,btScalar kst,btScalar ti);
};

#endif // #ifndef BT_ACCELERATED_SOFT_BODY_CPU_SOLVER_H
