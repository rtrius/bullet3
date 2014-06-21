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
static btScalar				ImplicitSolve(	btSoftBody::ImplicitFn* shape,
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

//
static inline btDbvtVolume	VolumeOf(	const btSoftBody::Face& f,
									 btScalar margin)
{
	const btVector3*	pts[]={	&f.m_n[0]->m_x,
		&f.m_n[1]->m_x,
		&f.m_n[2]->m_x};
	btDbvtVolume		vol=btDbvtVolume::FromPoints(pts,3);
	vol.Expand(btVector3(margin,margin,margin));
	return(vol);
}

//
static inline btScalar			AreaOf(		const btVector3& x0,
									   const btVector3& x1,
									   const btVector3& x2)
{
	const btVector3	a=x1-x0;
	const btVector3	b=x2-x0;
	const btVector3	cr=btCross(a,b);
	const btScalar	area=cr.length();
	return(area);
}

//
static inline btScalar		VolumeOf(	const btVector3& x0,
									 const btVector3& x1,
									 const btVector3& x2,
									 const btVector3& x3)
{
	const btVector3	a=x1-x0;
	const btVector3	b=x2-x0;
	const btVector3	c=x3-x0;
	return(btDot(a,btCross(b,c)));
}



//
static inline int		MatchEdge(	const btSoftBody::Node* a,
								  const btSoftBody::Node* b,
								  const btSoftBody::Node* ma,
								  const btSoftBody::Node* mb)
{
	if((a==ma)&&(b==mb)) return(0);
	if((a==mb)&&(b==ma)) return(1);
	return(-1);
}

//
// btSoftColliders
//
struct btSoftColliders
{
	///Detects collision between rigid body and soft body using SDF(Signed Distance Field)
	struct	CollideSDF_RS : btDbvt::ICollide
	{
		void		Process(const btDbvtNode* leaf)
		{
			btSoftBody::Node*	node=(btSoftBody::Node*)leaf->data;
			DoNode(*node);
		}
		void		DoNode(btSoftBody::Node& n) const
		{
			const btScalar			m = n.m_invMass > 0 ? dynmargin : stamargin;
			btSoftBody::RigidContact c;

			if(	(!n.m_battach)&&
				psb->checkContact(m_colObj1Wrap,n.m_x,m,c))
			{
				const btScalar	ima = n.m_invMass;
				const btScalar	imb= m_rigidBody? m_rigidBody->getInvMass() : 0.f;
				const btScalar	ms=ima+imb;
				if(ms>0)
				{
					const btTransform&	wtr=m_rigidBody?m_rigidBody->getWorldTransform() : m_colObj1Wrap->getCollisionObject()->getWorldTransform();
					static const btMatrix3x3	iwiStatic(0,0,0,0,0,0,0,0,0);
					const btMatrix3x3&	iwi=m_rigidBody?m_rigidBody->getInvInertiaTensorWorld() : iwiStatic;
					const btVector3		ra=n.m_x-wtr.getOrigin();
					const btVector3		va=m_rigidBody ? m_rigidBody->getVelocityInLocalPoint(ra)*psb->m_timeStep : btVector3(0,0,0);
					const btVector3		vb=n.m_x-n.m_q;	
					const btVector3		vr=vb-va;
					const btScalar		dn=btDot(vr,c.m_normal);
					const btVector3		fv=vr-c.m_normal*dn;
					const btScalar		fc=psb->m_cfg.m_dynamicFriction*m_colObj1Wrap->getCollisionObject()->getFriction();
					c.m_node	=	&n;
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
	struct	CollideVF_SS : btDbvt::ICollide
	{
		void		Process(const btDbvtNode* lnode,
			const btDbvtNode* lface)
		{
			btSoftBody::Node*	node=(btSoftBody::Node*)lnode->data;
			btSoftBody::Face*	face=(btSoftBody::Face*)lface->data;
			btVector3			o=node->m_x;
			btVector3			p;
			btScalar			d=SIMD_INFINITY;
			ProjectOrigin(	face->m_n[0]->m_x-o,
				face->m_n[1]->m_x-o,
				face->m_n[2]->m_x-o,
				p,d);
			const btScalar	m=mrg+(o-node->m_q).length()*2;
			if(d<(m*m))
			{
				const btSoftBody::Node*	n[]={face->m_n[0],face->m_n[1],face->m_n[2]};
				const btVector3			w=BaryCoord(n[0]->m_x,n[1]->m_x,n[2]->m_x,p+o);
				const btScalar			ma=node->m_invMass;
				btScalar				mb=BaryEval(n[0]->m_invMass, n[1]->m_invMass, n[2]->m_invMass, w);
				if(	(n[0]->m_invMass <= 0) || (n[1]->m_invMass <= 0)|| (n[2]->m_invMass <= 0) )
				{
					mb=0;
				}
				const btScalar	ms=ma+mb;
				if(ms>0)
				{
					btSoftBody::SoftContact	c;
					c.m_normal		=	p/-btSqrt(d);
					c.m_margin		=	m;
					c.m_node		=	node;
					c.m_face		=	face;
					c.m_weights		=	w;
					c.m_friction	=	btMax(psb[0]->m_cfg.m_dynamicFriction,psb[1]->m_cfg.m_dynamicFriction);
					c.m_cfm[0]		=	ma / ms * psb[0]->m_cfg.m_softContactHardness;
					c.m_cfm[1]		=	mb / ms * psb[1]->m_cfg.m_softContactHardness;
					psb[0]->m_softContacts.push_back(c);
				}
			}	
		}
		btSoftBody*		psb[2];
		btScalar		mrg;
	};
};

#endif //_BT_SOFT_BODY_INTERNALS_H
