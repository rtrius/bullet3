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

#include "btDefaultSoftBodySolver.h"
#include "BulletSoftBody/btSoftBody.h"
#include "BulletSoftBody/btSoftBodyInternals.h"


btDefaultSoftBodySolver::btDefaultSoftBodySolver()
{
}

btDefaultSoftBodySolver::~btDefaultSoftBodySolver()
{
}

void btDefaultSoftBodySolver::optimize( btAlignedObjectArray< btSoftBody * > &softBodies , bool forceUpdate)
{
	m_softBodySet.copyFromArray( softBodies );
}

void btDefaultSoftBodySolver::updateSoftBodies( )
{
	for ( int i=0; i < m_softBodySet.size(); i++)
	{
		btSoftBody*	psb = m_softBodySet[i];
		if (psb->isActive())
		{
			psb->updateNormals();	
		}
	}
}


void solveConstraintsSingleSoftBody(btSoftBody* softBody)
{
	//Prepare anchors
	for(int i = 0; i < softBody->m_anchors.size(); ++i)
	{
		btSoftBodyAnchor& a = softBody->m_anchors[i];
		btSoftBodyNode& node = softBody->m_nodes[a.m_nodeIndex];
		
		const btVector3	ra = a.m_body->getWorldTransform().getBasis() * a.m_local;
		a.m_impulseMatrix = ImpulseMatrix(	softBody->m_timeStep, node.m_invMass, a.m_body->getInvMass(),
									a.m_body->getInvInertiaTensorWorld(), ra);
		a.m_rotatedPosition = ra;
		a.m_invMassDt = softBody->m_timeStep * node.m_invMass;
		a.m_body->activate();
	}
	
	//Solve positions
	if(softBody->m_cfg.m_positionIterations > 0)
	{
		for(int iteration = 0; iteration < softBody->m_cfg.m_positionIterations; ++iteration)
		{
			const btScalar ti = iteration / (btScalar)softBody->m_cfg.m_positionIterations;
			btDefaultSoftBodySolver::PSolve_Anchors(softBody,1,ti);
			btDefaultSoftBodySolver::PSolve_RigidContacts(softBody,1,ti);
			btDefaultSoftBodySolver::PSolve_SoftContacts(softBody,1,ti);
			btDefaultSoftBodySolver::PSolve_Links(softBody,1,ti);
		}
		const btScalar	vc = (1 - softBody->m_cfg.m_damping) /  softBody->m_timeStep;
		for(int i = 0; i < softBody->m_nodes.size(); ++i)
		{
			btSoftBodyNode& n = softBody->m_nodes[i];
			n.m_velocity = (n.m_position - n.m_prevPosition)*vc;
			n.m_accumulatedForce = btVector3(0,0,0);		
		}
	}
}

void btDefaultSoftBodySolver::solveConstraints( float solverdt )
{
	for(int i = 0; i < m_softBodySet.size(); ++i)
	{
		btSoftBody*	psb = m_softBodySet[i];
		if ( psb->isActive() ) solveConstraintsSingleSoftBody(psb);
	}	
}

void btDefaultSoftBodySolver::predictMotion( float timeStep )
{
	for ( int i=0; i < m_softBodySet.size(); ++i)
	{
		btSoftBody*	psb = m_softBodySet[i];

		if (psb->isActive())
		{
			psb->predictMotion(timeStep);		
		}
	}
}

void btDefaultSoftBodySolver::PSolve_Anchors(btSoftBody* psb,btScalar kst,btScalar ti)
{
	const btScalar anchorHardness = psb->m_cfg.m_anchorHardness * kst;
	const btScalar	dt=psb->m_timeStep;
	for(int i=0,ni=psb->m_anchors.size();i<ni;++i)
	{
		const btSoftBodyAnchor& a = psb->m_anchors[i];
		const btTransform& t = a.m_body->getWorldTransform();
		btSoftBodyNode& n = psb->m_nodes[a.m_nodeIndex];
		const btVector3 wa = t*a.m_local;
		const btVector3 va = a.m_body->getVelocityInLocalPoint(a.m_rotatedPosition)*dt;
		const btVector3 vb = n.m_position - n.m_prevPosition;
		const btVector3 vr = (va - vb) + (wa - n.m_position) * anchorHardness;
		const btVector3 impulse = a.m_impulseMatrix * vr * a.m_influence;
		n.m_position += impulse * a.m_invMassDt;
		a.m_body->applyImpulse(-impulse,a.m_rotatedPosition);
	}
}

void btDefaultSoftBodySolver::PSolve_RigidContacts(btSoftBody* psb, btScalar kst, btScalar ti)
{
	const btScalar	dt = psb->m_timeStep;
	const btScalar	mrg = psb->getCollisionShape()->getMargin();
	for(int i=0,ni=psb->m_rigidContacts.size();i<ni;++i)
	{
		const btSoftRigidContact& c = psb->m_rigidContacts[i];
		if (c.m_colObj->hasContactResponse()) 
		{
			btSoftBodyNode& node = psb->m_nodes[c.m_nodeIndex];
		
			btRigidBody* tmpRigid = (btRigidBody*)btRigidBody::upcast(c.m_colObj);
			const btVector3		va = tmpRigid ? tmpRigid->getVelocityInLocalPoint(c.m_relativeNodePosition) * dt : btVector3(0,0,0);
			const btVector3		vb = node.m_position - node.m_prevPosition;	
			const btVector3		vr = vb-va;
			const btScalar		dn = btDot(vr, c.m_normal);		
			if(dn<=SIMD_EPSILON)
			{
				const btScalar		dp = btMin( (btDot(node.m_position, c.m_normal) + c.m_offset), mrg );
				const btVector3		fv = vr - (c.m_normal * dn);
				
				const btVector3		impulse = c.m_impulseMatrix * ( (vr - (fv * c.m_combinedFriction) + (c.m_normal * (dp * c.m_hardness))) * kst );
				node.m_position -= impulse * c.m_invMassDt;
				if (tmpRigid)
					tmpRigid->applyImpulse(impulse,c.m_relativeNodePosition);
			}
		}
	}
}

void btDefaultSoftBodySolver::PSolve_SoftContacts(btSoftBody* psb,btScalar,btScalar ti)
{

	for(int i=0,ni=psb->m_softContacts.size();i<ni;++i)
	{
		const btSoftSoftContact& c = psb->m_softContacts[i];
		const btVector3& nr = c.m_normal;
		btSoftBodyNode& n = psb->m_nodes[c.m_nodeIndex];
		
		btAlignedObjectArray<btSoftBodyNode>& faceNodes = c.m_faceSoftBody->m_nodes;
		const btSoftBodyFace& face = c.m_faceSoftBody->m_faces[c.m_faceIndex];
		btSoftBodyNode& faceNode0 = faceNodes[ face.m_indicies[0] ];
		btSoftBodyNode& faceNode1 = faceNodes[ face.m_indicies[1] ];
		btSoftBodyNode& faceNode2 = faceNodes[ face.m_indicies[2] ];
	
		btVector3 p = BaryEval(faceNode0.m_position, faceNode1.m_position, faceNode2.m_position, c.m_weights);
		btVector3 q = BaryEval(faceNode0.m_prevPosition, faceNode1.m_prevPosition, faceNode2.m_prevPosition, c.m_weights);											
		const btVector3		vr=(n.m_position-n.m_prevPosition)-(p-q);
		btVector3			corr(0,0,0);
		btScalar dot = btDot(vr,nr);
		if(dot<0)
		{
			const btScalar	j = c.m_margin - ( btDot(nr, n.m_position) - btDot(nr, p) );
			corr+=c.m_normal*j;
		}
		corr			-=	ProjectOnPlane(vr,nr)*c.m_friction;
		n.m_position += corr*c.m_nodeCfm;
		faceNode0.m_position -= corr*(c.m_faceCfm*c.m_weights.x());
		faceNode1.m_position -= corr*(c.m_faceCfm*c.m_weights.y());
		faceNode2.m_position -= corr*(c.m_faceCfm*c.m_weights.z());
	}
}

void btDefaultSoftBodySolver::PSolve_Links(btSoftBody* psb,btScalar kst,btScalar ti)
{
	for(int i=0,ni=psb->m_links.size();i<ni;++i)
	{			
		btSoftBodyLink& l = psb->m_links[i];
		if(l.m_scaledCombinedInvMass > 0)
		{
			btSoftBodyNode& a = psb->m_nodes[ l.m_linkIndicies[0] ];
			btSoftBodyNode& b = psb->m_nodes[ l.m_linkIndicies[1] ];
			const btVector3	del = b.m_position - a.m_position;
			const btScalar	len=del.length2();
			if (l.m_restLengthSquared+len > SIMD_EPSILON)
			{
				const btScalar k = ( (l.m_restLengthSquared-len)/(l.m_scaledCombinedInvMass*(l.m_restLengthSquared+len)) ) * kst;
				a.m_position -= del * (k * a.m_invMass);
				b.m_position += del * (k * b.m_invMass);
			}
		}
	}
}
