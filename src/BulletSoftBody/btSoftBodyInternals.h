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
///btSoftBody implementation by Nathanael Presson

#ifndef _BT_SOFT_BODY_INTERNALS_H
#define _BT_SOFT_BODY_INTERNALS_H

#include "btSoftBody.h"


#include "LinearMath/btQuickprof.h"
#include "LinearMath/btPolarDecomposition.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseInterface.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btConvexInternalShape.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpa2.h"
#include <string.h> //for memset
//
// btSymMatrix
//
template <typename T>
struct btSymMatrix
{
	btSymMatrix() : dim(0)					{}
	btSymMatrix(int n,const T& init=T())	{ resize(n,init); }
	void					resize(int n,const T& init=T())			{ dim=n;store.resize((n*(n+1))/2,init); }
	int						index(int c,int r) const				{ if(c>r) btSwap(c,r);btAssert(r<dim);return((r*(r+1))/2+c); }
	T&						operator()(int c,int r)					{ return(store[index(c,r)]); }
	const T&				operator()(int c,int r) const			{ return(store[index(c,r)]); }
	btAlignedObjectArray<T>	store;
	int						dim;
};	

//
// btSoftBodyCollisionShape
//
class btSoftBodyCollisionShape : public btConcaveShape
{
public:
	btSoftBody*						m_body;

	btSoftBodyCollisionShape(btSoftBody* backptr)
	{
		m_shapeType = SOFTBODY_SHAPE_PROXYTYPE;
		m_body=backptr;
	}

	virtual ~btSoftBodyCollisionShape()
	{

	}

	void	processAllTriangles(btTriangleCallback* /*callback*/,const btVector3& /*aabbMin*/,const btVector3& /*aabbMax*/) const
	{
		//not yet
		btAssert(0);
	}

	///getAabb returns the axis aligned bounding box in the coordinate frame of the given transform t.
	virtual void getAabb(const btTransform& t,btVector3& aabbMin,btVector3& aabbMax) const
	{
		// t is usually identity, except when colliding against btCompoundShape. See Issue 512 
		const btVector3	mins=m_body->m_aabbMin;
		const btVector3	maxs=m_body->m_aabbMax;
		const btVector3	crns[]={t*btVector3(mins.x(),mins.y(),mins.z()),
			t*btVector3(maxs.x(),mins.y(),mins.z()),
			t*btVector3(maxs.x(),maxs.y(),mins.z()),
			t*btVector3(mins.x(),maxs.y(),mins.z()),
			t*btVector3(mins.x(),mins.y(),maxs.z()),
			t*btVector3(maxs.x(),mins.y(),maxs.z()),
			t*btVector3(maxs.x(),maxs.y(),maxs.z()),
			t*btVector3(mins.x(),maxs.y(),maxs.z())};
		aabbMin=aabbMax=crns[0];
		for(int i=1;i<8;++i)
		{
			aabbMin.setMin(crns[i]);
			aabbMax.setMax(crns[i]);
		}
	}


	virtual void	setLocalScaling(const btVector3&)
	{		
		///na
	}
	virtual const btVector3& getLocalScaling() const
	{
		static const btVector3 dummy(1,1,1);
		return dummy;
	}
	virtual void	calculateLocalInertia(btScalar, btVector3&) const
	{
		///not yet
		btAssert(0);
	}
	virtual const char*	getName()const
	{
		return "SoftBody";
	}

};


//
// Inline's
//

//
template <typename T>
static inline void			ZeroInitialize(T& value)
{
	memset(&value,0,sizeof(T));
}
//
template <typename T>
static inline T				Lerp(const T& a,const T& b,btScalar t)
{ return(a+(b-a)*t); }

//
static inline btMatrix3x3	Lerp(	const btMatrix3x3& a,
								 const btMatrix3x3& b,
								 btScalar t)
{
	btMatrix3x3	r;
	r[0]=Lerp(a[0],b[0],t);
	r[1]=Lerp(a[1],b[1],t);
	r[2]=Lerp(a[2],b[2],t);
	return(r);
}

//
template <typename T>
static inline bool			SameSign(const T& x,const T& y)
{ return((x*y)>0); }

//
static inline btMatrix3x3	Cross(const btVector3& v)
{
	btMatrix3x3	m;
	m[0]=btVector3(0,-v.z(),+v.y());
	m[1]=btVector3(+v.z(),0,-v.x());
	m[2]=btVector3(-v.y(),+v.x(),0);
	return(m);
}
//
static inline btMatrix3x3	Diagonal(btScalar x)
{
	btMatrix3x3	m;
	m[0]=btVector3(x,0,0);
	m[1]=btVector3(0,x,0);
	m[2]=btVector3(0,0,x);
	return(m);
}
//
static inline btMatrix3x3	Add(const btMatrix3x3& a,
								const btMatrix3x3& b)
{
	btMatrix3x3	r;
	for(int i=0;i<3;++i) r[i]=a[i]+b[i];
	return(r);
}
//
static inline btMatrix3x3	Sub(const btMatrix3x3& a,
								const btMatrix3x3& b)
{
	btMatrix3x3	r;
	for(int i=0;i<3;++i) r[i]=a[i]-b[i];
	return(r);
}
//
static inline btMatrix3x3	Mul(const btMatrix3x3& a,
								btScalar b)
{
	btMatrix3x3	r;
	for(int i=0;i<3;++i) r[i]=a[i]*b;
	return(r);
}

//
static inline btMatrix3x3	MassMatrix(btScalar im,const btMatrix3x3& iwi,const btVector3& r)
{
	const btMatrix3x3	cr=Cross(r);
	return(Sub(Diagonal(im),cr*iwi*cr));
}

//
static inline btMatrix3x3	ImpulseMatrix(	btScalar dt,
										  btScalar ima,
										  btScalar imb,
										  const btMatrix3x3& iwi,
										  const btVector3& r)
{
	return(Diagonal(1/dt)*Add(Diagonal(ima),MassMatrix(imb,iwi,r)).inverse());
}

//
static inline btMatrix3x3	ImpulseMatrix(	btScalar ima,const btMatrix3x3& iia,const btVector3& ra,
										  btScalar imb,const btMatrix3x3& iib,const btVector3& rb)	
{
	return(Add(MassMatrix(ima,iia,ra),MassMatrix(imb,iib,rb)).inverse());
}

//
static inline btVector3		ProjectOnAxis(	const btVector3& v,
										  const btVector3& a)
{
	return(a*btDot(v,a));
}
//
static inline btVector3		ProjectOnPlane(	const btVector3& v,
										   const btVector3& a)
{
	return(v-ProjectOnAxis(v,a));
}

//
static inline void			ProjectOrigin(	const btVector3& a,
										  const btVector3& b,
										  btVector3& prj,
										  btScalar& sqd)
{
	const btVector3	d=b-a;
	const btScalar	m2=d.length2();
	if(m2>SIMD_EPSILON)
	{	
		const btScalar	t = btMax( btScalar(0.0), btMin(-btDot(a,d)/m2, btScalar(1.0)) );
		const btVector3	p=a+d*t;
		const btScalar	l2=p.length2();
		if(l2<sqd)
		{
			prj=p;
			sqd=l2;
		}
	}
}
//
static inline void			ProjectOrigin(	const btVector3& a,
										  const btVector3& b,
										  const btVector3& c,
										  btVector3& prj,
										  btScalar& sqd)
{
	const btVector3&	q=btCross(b-a,c-a);
	const btScalar		m2=q.length2();
	if(m2>SIMD_EPSILON)
	{
		const btVector3	n=q/btSqrt(m2);
		const btScalar	k=btDot(a,n);
		const btScalar	k2=k*k;
		if(k2<sqd)
		{
			const btVector3	p=n*k;
			if(	(btDot(btCross(a-p,b-p),q)>0)&&
				(btDot(btCross(b-p,c-p),q)>0)&&
				(btDot(btCross(c-p,a-p),q)>0))
			{			
				prj=p;
				sqd=k2;
			}
			else
			{
				ProjectOrigin(a,b,prj,sqd);
				ProjectOrigin(b,c,prj,sqd);
				ProjectOrigin(c,a,prj,sqd);
			}
		}
	}
}

//
template <typename T>
static inline T				BaryEval(		const T& a,
									 const T& b,
									 const T& c,
									 const btVector3& coord)
{
	return(a*coord.x()+b*coord.y()+c*coord.z());
}
//
static inline btVector3		BaryCoord(	const btVector3& a,
									  const btVector3& b,
									  const btVector3& c,
									  const btVector3& p)
{
	const btScalar	w[]={	btCross(a-p,b-p).length(),
		btCross(b-p,c-p).length(),
		btCross(c-p,a-p).length()};
	const btScalar	isum=1/(w[0]+w[1]+w[2]);
	return(btVector3(w[1]*isum,w[2]*isum,w[0]*isum));
}

///Calculates the fraction, in range [0, 1], at which the line segment a-b intersects shape; 
///returns -1 if both vertices are on the same side of the shape. 
static btScalar				ImplicitSolve(btSoftImplicitShape* shape,
										  const btVector3& a,
										  const btVector3& b,
										  const btScalar accuracy,
										  const int maxiterations=256)
{
	btScalar	span[2]={0,1};
	btScalar	values[2]={shape->signedDistance(a),shape->signedDistance(b)};
	if(values[0]>values[1])
	{
		btSwap(span[0],span[1]);
		btSwap(values[0],values[1]);
	}
	if(values[0]>-accuracy) return(-1);
	if(values[1]<+accuracy) return(-1);
	for(int i=0;i<maxiterations;++i)
	{
		const btScalar	t=Lerp(span[0],span[1],values[0]/(values[0]-values[1]));
		const btScalar	v=shape->signedDistance(Lerp(a,b,t));
		if((t<=0)||(t>=1))		break;
		if(btFabs(v)<accuracy)	return(t);
		if(v<0)
		{ span[0]=t;values[0]=v; }
		else
		{ span[1]=t;values[1]=v; }
	}
	return(-1);
}

///Returns the area of a parallelogram defined by (x1 - x0) and (x2 - x0);
///if x0, x1, and x2 form a triangle, then the result should be halved.
static inline btScalar AreaOfParallelogram(const btVector3& x0,  const btVector3& x1, const btVector3& x2)
{
	return btCross(x1 - x0, x2 - x0).length();
}

//
// btSoftColliders
//
struct btSoftColliders
{
	///Detects collision between rigid body and soft body using SDF(Signed Distance Field)
	struct	CollideSDF_RS : btDbvt::ICollide
	{
		void Process(const btDbvtNode* leaf)
		{
			int nodeIndex = reinterpret_cast<int>(leaf->data);
			DoNode(nodeIndex);
		}
		void DoNode(int nodeIndex) const
		{
			 btSoftBodyNode& n = psb->m_nodes[nodeIndex];
		
			const btScalar m = n.m_invMass > 0 ? dynmargin : stamargin;
			btSoftRigidContact c;

			if( !n.m_isAttachedToAnchor && psb->checkContact(m_colObj1Wrap,n.m_position,m,c) )
			{
				const btScalar	ima = n.m_invMass;
				const btScalar	imb= m_rigidBody? m_rigidBody->getInvMass() : 0.f;
				const btScalar	ms=ima+imb;
				if(ms>0)
				{
					const btTransform&	wtr=m_rigidBody?m_rigidBody->getWorldTransform() : m_colObj1Wrap->getCollisionObject()->getWorldTransform();
					static const btMatrix3x3	iwiStatic(0,0,0,0,0,0,0,0,0);
					const btMatrix3x3&	iwi=m_rigidBody?m_rigidBody->getInvInertiaTensorWorld() : iwiStatic;
					const btVector3		ra=n.m_position-wtr.getOrigin();
					const btVector3		va=m_rigidBody ? m_rigidBody->getVelocityInLocalPoint(ra)*psb->m_timeStep : btVector3(0,0,0);
					const btVector3		vb = n.m_position - n.m_prevPosition;	
					const btVector3		vr=vb-va;
					const btScalar		dn=btDot(vr,c.m_normal);
					const btVector3		fv=vr-c.m_normal*dn;
					const btScalar		fc=psb->m_cfg.m_dynamicFriction*m_colObj1Wrap->getCollisionObject()->getFriction();
					c.m_nodeIndex	=	nodeIndex;
					c.m_impulseMatrix = ImpulseMatrix(psb->m_timeStep,ima,imb,iwi,ra);
					c.m_relativeNodePosition = ra;
					c.m_invMassDt = ima * psb->m_timeStep;
			        c.m_combinedFriction =	fv.length2() < (dn*fc*dn*fc) ? 0 : 1 - fc;
					c.m_hardness =	m_colObj1Wrap->getCollisionObject()->isStaticOrKinematicObject()
									? psb->m_cfg.m_kinematicContactHardness : psb->m_cfg.m_rigidContactHardness;
					psb->m_rigidContacts.push_back(c);
					if (m_rigidBody)
						m_rigidBody->activate();
				}
			}
		}
		btSoftBody*		psb;
		const btCollisionObjectWrapper*	m_colObj1Wrap;
		btRigidBody*	m_rigidBody;
		btScalar		dynmargin;
		btScalar		stamargin;
	};
	
	///Detects collision between 2 soft bodies by testing vertices(nodes) against triangles(faces)
	struct SoftSoftVertexFaceCollider : btDbvt::ICollide
	{
		void Process(const btDbvtNode* lnode, const btDbvtNode* lface)
		{
			int nodeIndex = reinterpret_cast<int>(lnode->data);
			btSoftBodyNode& node = m_nodeSoftBody->m_nodes[nodeIndex];
			
			int faceIndex = reinterpret_cast<int>(lface->data);
			btSoftBodyFace&	face = m_faceSoftBody->m_faces[faceIndex];
			btVector3			o = node.m_position;
			btVector3			p;
			btScalar			d=SIMD_INFINITY;
			
			const btSoftBodyNode&	faceNode0 = m_faceSoftBody->m_nodes[ face.m_indicies[0] ];
			const btSoftBodyNode&	faceNode1 = m_faceSoftBody->m_nodes[ face.m_indicies[1] ];
			const btSoftBodyNode&	faceNode2 = m_faceSoftBody->m_nodes[ face.m_indicies[2] ];
			
			ProjectOrigin(faceNode0.m_position - o, faceNode1.m_position - o, faceNode2.m_position - o, p, d);
			const btScalar	m = mrg + (o - node.m_prevPosition).length() * 2;
			if(d<(m*m))
			{
				btVector3 w = BaryCoord(faceNode0.m_position, faceNode1.m_position, faceNode2.m_position, p + o);
				btScalar ma = node.m_invMass;
				btScalar mb = BaryEval(faceNode0.m_invMass, faceNode1.m_invMass, faceNode2.m_invMass, w);
				if(	(faceNode0.m_invMass <= 0) || (faceNode1.m_invMass <= 0)|| (faceNode2.m_invMass <= 0) )
				{
					mb=0;
				}
				const btScalar	ms=ma+mb;
				if(ms>0)
				{
					btSoftSoftContact	c;
					c.m_nodeSoftBody = m_nodeSoftBody;
					c.m_faceSoftBody = m_faceSoftBody;
					c.m_normal		=	p/-btSqrt(d);
					c.m_margin		=	m;
					c.m_nodeIndex	=	nodeIndex;
					c.m_faceIndex	=	faceIndex;
					c.m_weights		=	w;
					c.m_friction	=	btMax(m_nodeSoftBody->m_cfg.m_dynamicFriction, m_faceSoftBody->m_cfg.m_dynamicFriction);
					c.m_nodeCfm		=	ma / ms * m_nodeSoftBody->m_cfg.m_softContactHardness;
					c.m_faceCfm		=	mb / ms * m_faceSoftBody->m_cfg.m_softContactHardness;
					m_nodeSoftBody->m_softContacts.push_back(c);
				}
			}	
		}
		btSoftBody* m_nodeSoftBody;
		btSoftBody* m_faceSoftBody;
		btScalar		mrg;
	};
};

#endif //_BT_SOFT_BODY_INTERNALS_H
