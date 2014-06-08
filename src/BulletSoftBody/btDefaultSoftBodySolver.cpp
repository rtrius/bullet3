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
	//Prepare links
	for(int i = 0; i < softBody->m_links.size();++i)
	{
		btSoftBody::Link& l = softBody->m_links[i];
		l.m_gradient0to1 = l.m_n[1]->m_q - l.m_n[0]->m_q;
		l.m_impulseScaling = 1 / (l.m_gradient0to1.length2() * l.m_scaledCombinedInvMass);
	}
	
	//Prepare anchors
	for(int i = 0; i < softBody->m_anchors.size(); ++i)
	{
		btSoftBody::Anchor& a = softBody->m_anchors[i];
		const btVector3	ra = a.m_body->getWorldTransform().getBasis() * a.m_local;
		a.m_impulseMatrix = ImpulseMatrix(	softBody->m_sst.sdt, a.m_node->m_invMass, a.m_body->getInvMass(),
									a.m_body->getInvInertiaTensorWorld(), ra);
		a.m_rotatedPosition = ra;
		a.m_invMassDt = softBody->m_sst.sdt * a.m_node->m_invMass;
		a.m_body->activate();
	}
	
	//Solve velocities
	if(softBody->m_cfg.m_velocityIterations > 0)
	{
		// Solve
		for(int iteration = 0; iteration < softBody->m_cfg.m_velocityIterations; ++iteration)
		{
			btDefaultSoftBodySolver::VSolve_Links(softBody,1);
		}
		
		// Update
		for(int i = 0; i < softBody->m_nodes.size();++i)
		{
			btSoftBody::Node& n = softBody->m_nodes[i];
			n.m_x =	n.m_q + n.m_v * softBody->m_sst.sdt;
		}
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
		const btScalar	vc = (1 - softBody->m_cfg.m_damping) /  softBody->m_sst.sdt;
		for(int i = 0; i < softBody->m_nodes.size(); ++i)
		{
			btSoftBody::Node& n = softBody->m_nodes[i];
			n.m_v = (n.m_x - n.m_q)*vc;
			n.m_f = btVector3(0,0,0);		
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




void btDefaultSoftBodySolver::VSolve_Links(btSoftBody* psb,btScalar kst)
{
	for(int i=0,ni=psb->m_links.size();i<ni;++i)
	{			
		btSoftBody::Link& l = psb->m_links[i];
		btSoftBody::Node** n = l.m_n;
		const btScalar j = -btDot(l.m_gradient0to1, n[0]->m_v - n[1]->m_v) * l.m_impulseScaling * kst;
		n[0]->m_v += l.m_gradient0to1 * (j * n[0]->m_invMass);
		n[1]->m_v -= l.m_gradient0to1 * (j * n[1]->m_invMass);
	}
}

void btDefaultSoftBodySolver::PSolve_Anchors(btSoftBody* psb,btScalar kst,btScalar ti)
{
	const btScalar anchorHardness = psb->m_cfg.m_anchorHardness * kst;
	const btScalar	dt=psb->m_sst.sdt;
	for(int i=0,ni=psb->m_anchors.size();i<ni;++i)
	{
		const btSoftBody::Anchor& a = psb->m_anchors[i];
		const btTransform& t = a.m_body->getWorldTransform();
		btSoftBody::Node& n = *a.m_node;
		const btVector3 wa = t*a.m_local;
		const btVector3 va = a.m_body->getVelocityInLocalPoint(a.m_rotatedPosition)*dt;
		const btVector3 vb = n.m_x-n.m_q;
		const btVector3 vr = (va - vb) + (wa - n.m_x) * anchorHardness;
		const btVector3 impulse = a.m_impulseMatrix * vr * a.m_influence;
		n.m_x += impulse * a.m_invMassDt;
		a.m_body->applyImpulse(-impulse,a.m_rotatedPosition);
	}
}

void btDefaultSoftBodySolver::PSolve_RigidContacts(btSoftBody* psb, btScalar kst, btScalar ti)
{
	const btScalar	dt = psb->m_sst.sdt;
	const btScalar	mrg = psb->getCollisionShape()->getMargin();
	for(int i=0,ni=psb->m_rigidContacts.size();i<ni;++i)
	{
		const btSoftBody::RigidContact& c = psb->m_rigidContacts[i];
		const btSoftBody::sCti&			cti = c.m_cti;	
		if (cti.m_colObj->hasContactResponse()) 
		{
			btRigidBody* tmpRigid = (btRigidBody*)btRigidBody::upcast(cti.m_colObj);
			const btVector3		va = tmpRigid ? tmpRigid->getVelocityInLocalPoint(c.m_relativeNodePosition) * dt : btVector3(0,0,0);
			const btVector3		vb = c.m_node->m_x-c.m_node->m_q;	
			const btVector3		vr = vb-va;
			const btScalar		dn = btDot(vr, cti.m_normal);		
			if(dn<=SIMD_EPSILON)
			{
				const btScalar		dp = btMin( (btDot(c.m_node->m_x, cti.m_normal) + cti.m_offset), mrg );
				const btVector3		fv = vr - (cti.m_normal * dn);
				// c0 is the impulse matrix, c3 is 1 - the friction coefficient or 0, c4 is the contact hardness coefficient
				const btVector3		impulse = c.m_impulseMatrix * ( (vr - (fv * c.m_combinedFriction) + (cti.m_normal * (dp * c.m_hardness))) * kst );
				c.m_node->m_x -= impulse * c.m_invMassDt;
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
		const btSoftBody::SoftContact& c = psb->m_softContacts[i];
		const btVector3&	nr=c.m_normal;
		btSoftBody::Node&				n=*c.m_node;
		btSoftBody::Face&				f=*c.m_face;
		const btVector3		p=BaryEval(	f.m_n[0]->m_x,
			f.m_n[1]->m_x,
			f.m_n[2]->m_x,
			c.m_weights);
		const btVector3		q=BaryEval(	f.m_n[0]->m_q,
			f.m_n[1]->m_q,
			f.m_n[2]->m_q,
			c.m_weights);											
		const btVector3		vr=(n.m_x-n.m_q)-(p-q);
		btVector3			corr(0,0,0);
		btScalar dot = btDot(vr,nr);
		if(dot<0)
		{
			const btScalar	j=c.m_margin-(btDot(nr,n.m_x)-btDot(nr,p));
			corr+=c.m_normal*j;
		}
		corr			-=	ProjectOnPlane(vr,nr)*c.m_friction;
		n.m_x			+=	corr*c.m_cfm[0];
		f.m_n[0]->m_x	-=	corr*(c.m_cfm[1]*c.m_weights.x());
		f.m_n[1]->m_x	-=	corr*(c.m_cfm[1]*c.m_weights.y());
		f.m_n[2]->m_x	-=	corr*(c.m_cfm[1]*c.m_weights.z());
	}
}

void btDefaultSoftBodySolver::PSolve_Links(btSoftBody* psb,btScalar kst,btScalar ti)
{
	for(int i=0,ni=psb->m_links.size();i<ni;++i)
	{			
		btSoftBody::Link&	l=psb->m_links[i];
		if(l.m_scaledCombinedInvMass > 0)
		{
			btSoftBody::Node&			a=*l.m_n[0];
			btSoftBody::Node&			b=*l.m_n[1];
			const btVector3	del=b.m_x-a.m_x;
			const btScalar	len=del.length2();
			if (l.m_restLengthSquared+len > SIMD_EPSILON)
			{
				const btScalar k = ( (l.m_restLengthSquared-len)/(l.m_scaledCombinedInvMass*(l.m_restLengthSquared+len)) ) * kst;
				a.m_x -= del * (k * a.m_invMass);
				b.m_x += del * (k * b.m_invMass);
			}
		}
	}
}
