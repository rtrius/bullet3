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

#include "btSoftBodyInternals.h"
#include "BulletSoftBody/btSoftBodySolvers.h"
#include "btSoftBodyData.h"
#include "LinearMath/btSerializer.h"


//
btSoftBody::btSoftBody(btSoftBodyWorldInfo*	worldInfo, int numNodes, const btVector3* nodePositions, const btScalar* nodeMasses)
:m_softBodySolver(0),m_worldInfo(worldInfo)
{	
	//
	initDefaults();

	// Default material 
	Material*	pm=appendMaterial();
	pm->m_linearStiffness	=	1;
	pm->m_debugDraw = true;

	// Nodes
	const btScalar		margin=getCollisionShape()->getMargin();
	m_nodes.resize(numNodes);
	for(int i = 0; i < numNodes; ++i)
	{	
		Node&	n=m_nodes[i];
		ZeroInitialize(n);
		n.m_x		=	nodePositions ? *nodePositions++ : btVector3(0,0,0);
		n.m_q		=	n.m_x;
		n.m_invMass		= nodeMasses ? *nodeMasses++ : 1;
		n.m_invMass		= (n.m_invMass > 0) ? 1 / n.m_invMass : 0;
		n.m_leaf	=	m_nodeBvh.insert(btDbvtVolume::FromCR(n.m_x,margin),&n);
		n.m_material=	pm;
	}
	updateBounds();	

}

btSoftBody::btSoftBody(btSoftBodyWorldInfo*	worldInfo)
:m_worldInfo(worldInfo)
{
	initDefaults();
}


void	btSoftBody::initDefaults()
{
	m_internalType		=	CO_SOFT_BODY;
	m_bUpdateRtCst		=	true;
	m_aabbMin = btVector3(0,0,0);
	m_aabbMax = btVector3(0,0,0);
	m_worldTransform.setIdentity();
	
	///for now, create a collision shape internally
	m_collisionShape = new btSoftBodyCollisionShape(this);
	m_collisionShape->setMargin(0.25f);

	m_restLengthScale = btScalar(1.0);
}

//
btSoftBody::~btSoftBody()
{
	//for now, delete the internal shape
	delete m_collisionShape;	
	int i;

	for(i=0;i<m_materials.size();++i) 
		btAlignedFree(m_materials[i]);
}

//
bool			btSoftBody::checkLink(int node0,int node1) const
{
	return(checkLink(&m_nodes[node0],&m_nodes[node1]));
}

//
bool			btSoftBody::checkLink(const Node* node0,const Node* node1) const
{
	const Node*	n[]={node0,node1};
	for(int i=0,ni=m_links.size();i<ni;++i)
	{
		const Link&	l=m_links[i];
		if(	(l.m_n[0]==n[0]&&l.m_n[1]==n[1])||
			(l.m_n[0]==n[1]&&l.m_n[1]==n[0]))
		{
			return(true);
		}
	}
	return(false);
}

//
bool			btSoftBody::checkFace(int node0,int node1,int node2) const
{
	const Node*	n[]={	&m_nodes[node0],
		&m_nodes[node1],
		&m_nodes[node2]};
	for(int i=0,ni=m_faces.size();i<ni;++i)
	{
		const Face&	f=m_faces[i];
		int			c=0;
		for(int j=0;j<3;++j)
		{
			if(	(f.m_n[j]==n[0])||
				(f.m_n[j]==n[1])||
				(f.m_n[j]==n[2])) c|=1<<j; else break;
		}
		if(c==7) return(true);
	}
	return(false);
}

//
btSoftBody::Material*		btSoftBody::appendMaterial()
{
	Material*	pm=new(btAlignedAlloc(sizeof(Material),16)) Material();
	if(m_materials.size()>0)
		*pm=*m_materials[0];
	else
		ZeroInitialize(*pm);
	m_materials.push_back(pm);
	return(pm);
}

//
void btSoftBody::appendNode(const btVector3& position, btScalar mass)
{
	if(m_nodes.capacity()==m_nodes.size())
	{
		pointersToIndices();
		m_nodes.reserve(m_nodes.size()*2+1);
		indicesToPointers();
	}
	const btScalar	margin=getCollisionShape()->getMargin();
	m_nodes.push_back(Node());
	Node&			n=m_nodes[m_nodes.size()-1];
	ZeroInitialize(n);
	n.m_x			=	position;
	n.m_q			=	n.m_x;
	n.m_invMass 	= (mass > 0) ? 1 / mass : 0;
	n.m_material	=	m_materials[0];
	n.m_leaf		=	m_nodeBvh.insert(btDbvtVolume::FromCR(n.m_x,margin),&n);
}

//
void			btSoftBody::appendLink(int model,Material* mat)
{
	Link	l;
	if(model>=0)
		l=m_links[model];
	else
	{ ZeroInitialize(l);l.m_material=mat?mat:m_materials[0]; }
	m_links.push_back(l);
}

//
void			btSoftBody::appendLink(	int node0,
									   int node1,
									   Material* mat,
									   bool bcheckexist)
{
	appendLink(&m_nodes[node0],&m_nodes[node1],mat,bcheckexist);
}

//
void			btSoftBody::appendLink(	Node* node0,
									   Node* node1,
									   Material* mat,
									   bool bcheckexist)
{
	if((!bcheckexist)||(!checkLink(node0,node1)))
	{
		appendLink(-1,mat);
		Link&	l=m_links[m_links.size()-1];
		l.m_n[0]		=	node0;
		l.m_n[1]		=	node1;
		l.m_restLength = (l.m_n[0]->m_x - l.m_n[1]->m_x).length();
		m_bUpdateRtCst=true;
	}
}

//
void			btSoftBody::appendFace(int model,Material* mat)
{
	Face	f;
	if(model>=0)
	{ f=m_faces[model]; }
	else
	{ ZeroInitialize(f);f.m_material=mat?mat:m_materials[0]; }
	m_faces.push_back(f);
}

//
void			btSoftBody::appendFace(int node0,int node1,int node2,Material* mat)
{
	if (node0==node1)
		return;
	if (node1==node2)
		return;
	if (node2==node0)
		return;

	appendFace(-1,mat);
	Face&	f=m_faces[m_faces.size()-1];
	btAssert(node0!=node1);
	btAssert(node1!=node2);
	btAssert(node2!=node0);
	f.m_n[0]	=	&m_nodes[node0];
	f.m_n[1]	=	&m_nodes[node1];
	f.m_n[2]	=	&m_nodes[node2];
	f.m_area = AreaOf(f.m_n[0]->m_x, f.m_n[1]->m_x, f.m_n[2]->m_x);	
	m_bUpdateRtCst=true;
}

//

void			btSoftBody::appendAnchor(int node,btRigidBody* body, bool disableCollisionBetweenLinkedBodies,btScalar influence)
{
	btVector3 local = body->getWorldTransform().inverse()*m_nodes[node].m_x;
	appendAnchor(node,body,local,disableCollisionBetweenLinkedBodies,influence);
}

//
void			btSoftBody::appendAnchor(int node,btRigidBody* body, const btVector3& localPivot,bool disableCollisionBetweenLinkedBodies,btScalar influence)
{
	if (disableCollisionBetweenLinkedBodies)
	{
		if (m_collisionDisabledObjects.findLinearSearch(body)==m_collisionDisabledObjects.size())
		{
			m_collisionDisabledObjects.push_back(body);
		}
	}

	Anchor	a;
	a.m_node			=	&m_nodes[node];
	a.m_body			=	body;
	a.m_local			=	localPivot;
	a.m_node->m_battach	=	1;
	a.m_influence = influence;
	m_anchors.push_back(a);
}

//
void btSoftBody::addForceAllNodes(const btVector3& force)
{
	for(int i = 0; i < m_nodes.size(); ++i)  addForce(force,i);
}

//
void			btSoftBody::addForce(const btVector3& force,int node)
{
	Node&	n=m_nodes[node];
	if(n.m_invMass > 0) n.m_f += force;
}


//
void btSoftBody::addVelocityAllNodes(const btVector3& velocity)
{
	for(int i = 0; i < m_nodes.size(); ++i) addVelocity(velocity,i);
}

void btSoftBody::setVelocityAllNodes(const btVector3& velocity)
{
	for(int i = 0; i < m_nodes.size(); ++i) 
	{
		Node& n = m_nodes[i];
		if(n.m_invMass > 0) n.m_v = velocity;
	}
}


//
void btSoftBody::addVelocity(const btVector3& velocity, int node)
{
	Node& n = m_nodes[node];
	if(n.m_invMass > 0) n.m_v += velocity;
}

//
void			btSoftBody::setMass(int node,btScalar mass)
{
	m_nodes[node].m_invMass = (mass > 0) ? 1 / mass : 0;
	m_bUpdateRtCst=true;
}

//
btScalar		btSoftBody::getMass(int node) const
{
	return (m_nodes[node].m_invMass > 0) ? 1 / m_nodes[node].m_invMass : 0;
}

//
btScalar		btSoftBody::getTotalMass() const
{
	btScalar	mass=0;
	for(int i=0;i<m_nodes.size();++i)
	{
		mass+=getMass(i);
	}
	return(mass);
}

//
void			btSoftBody::setTotalMass(btScalar mass,bool fromfaces)
{
	int i;

	if(fromfaces)
	{

		for(i=0;i<m_nodes.size();++i)
		{
			m_nodes[i].m_invMass = 0;
		}
		for(i=0;i<m_faces.size();++i)
		{
			const Face&		f=m_faces[i];
			const btScalar	twicearea=AreaOf(	f.m_n[0]->m_x,
				f.m_n[1]->m_x,
				f.m_n[2]->m_x);
			for(int j=0;j<3;++j)
			{
				f.m_n[j]->m_invMass += twicearea;
			}
		}
		for( i=0;i<m_nodes.size();++i)
		{
			m_nodes[i].m_invMass = 1 / m_nodes[i].m_invMass;
		}
	}
	const btScalar	tm=getTotalMass();
	const btScalar	itm=1/tm;
	for( i=0;i<m_nodes.size();++i)
	{
		m_nodes[i].m_invMass /= itm * mass;
	}
	m_bUpdateRtCst=true;
}

//
void			btSoftBody::setTotalDensity(btScalar density)
{
	setTotalMass( getClosedTrimeshVolume() * density, true);
}

//
void			btSoftBody::transform(const btTransform& trs)
{
	const btScalar	margin=getCollisionShape()->getMargin();
	ATTRIBUTE_ALIGNED16(btDbvtVolume)	vol;
	
	for(int i=0,ni=m_nodes.size();i<ni;++i)
	{
		Node&	n=m_nodes[i];
		n.m_x=trs*n.m_x;
		n.m_q=trs*n.m_q;
		n.m_normal = trs.getBasis() * n.m_normal;
		vol = btDbvtVolume::FromCR(n.m_x,margin);
		
		m_nodeBvh.update(n.m_leaf,vol);
	}
	updateNormals();
	updateBounds();
	updateConstants();
}

//
void			btSoftBody::translate(const btVector3& trs)
{
	btTransform	t;
	t.setIdentity();
	t.setOrigin(trs);
	transform(t);
}

//
void			btSoftBody::rotate(	const btQuaternion& rot)
{
	btTransform	t;
	t.setIdentity();
	t.setRotation(rot);
	transform(t);
}

//
void			btSoftBody::scale(const btVector3& scl)
{

	const btScalar	margin=getCollisionShape()->getMargin();
	ATTRIBUTE_ALIGNED16(btDbvtVolume)	vol;
	
	for(int i=0,ni=m_nodes.size();i<ni;++i)
	{
		Node&	n=m_nodes[i];
		n.m_x*=scl;
		n.m_q*=scl;
		vol = btDbvtVolume::FromCR(n.m_x,margin);
		m_nodeBvh.update(n.m_leaf, vol);
	}
	updateNormals();
	updateBounds();
	updateConstants();
}

//
btScalar btSoftBody::getRestLengthScale()
{
	return m_restLengthScale;
}

//
void btSoftBody::setRestLengthScale(btScalar restLengthScale)
{
	for(int i=0, ni=m_links.size(); i<ni; ++i)
	{
		Link&		l=m_links[i];
		l.m_restLength = l.m_restLength / m_restLengthScale * restLengthScale;
		l.m_restLengthSquared = l.m_restLength * l.m_restLength;
	}
	m_restLengthScale = restLengthScale;
	
	if (getActivationState() == ISLAND_SLEEPING)
		activate();
}

//
void btSoftBody::setPose(btScalar poseMatching)
{
	m_pose.m_poseMatching = poseMatching;

	//Set per-node weights
	{
		int numKinematicNodes = 0;
		for(int i = 0; i < m_nodes.size(); ++i)
			if(m_nodes[i].m_invMass <= 0) ++numKinematicNodes;
		
		btScalar dynamicMass = getTotalMass();
		btScalar kinematicMass = dynamicMass * m_nodes.size() * btScalar(1000.0);
		btScalar totalMass = dynamicMass + btScalar(numKinematicNodes) * kinematicMass;
		
		m_pose.m_weights.resize( m_nodes.size() );
		for(int i = 0; i < m_nodes.size(); ++i)
		{
			Node& n = m_nodes[i];
			m_pose.m_weights[i] = (n.m_invMass > 0) ? 1 / (m_nodes[i].m_invMass * totalMass) : kinematicMass / totalMass;
		}
	}
	
	//Initialize initial positions, center of mass, rotation
	{
		btVector3 centerOfMass = evaluateCenterOfMass();
		
		m_pose.m_referencePositions.resize( m_nodes.size() );
		for(int i = 0; i < m_nodes.size(); ++i) m_pose.m_referencePositions[i] = m_nodes[i].m_x - centerOfMass;
		
		m_pose.m_centerOfMass = centerOfMass;
		m_pose.m_rotation.setIdentity();
	}
	
	updateConstants();
}

//
void btSoftBody::updateAndApplyPose()
{
	if( m_pose.m_poseMatching != btScalar(0.0) )
	{
		//Update center of mass
		btVector3 centerOfMass = evaluateCenterOfMass();
		
		//Update rotation
		btMatrix3x3 rotation;
		{
			btMatrix3x3 Apq, scaling;
			
			Apq[0] = Apq[1] = Apq[2] = btVector3(0,0,0);
			Apq[0].setX(SIMD_EPSILON);
			Apq[1].setY(SIMD_EPSILON*2);
			Apq[2].setZ(SIMD_EPSILON*3);
			
			for(int i = 0; i < m_nodes.size(); ++i)
			{
				const btVector3 a = m_pose.m_weights[i] * (m_nodes[i].m_x - centerOfMass);
				const btVector3& b = m_pose.m_referencePositions[i];
				Apq[0] += a.x() * b;
				Apq[1] += a.y() * b;
				Apq[2] += a.z() * b;
			}
			
			const btPolarDecomposition polar;  
			polar.decompose(Apq, rotation, scaling);
		}
		
		//
		m_pose.m_centerOfMass = centerOfMass;
		m_pose.m_rotation = rotation;
		
		//Apply pose matching constraint
		{
			for(int i = 0; i < m_nodes.size(); ++i)
			{
				Node& n = m_nodes[i];
				if(n.m_invMass > 0)
				{
					btVector3 x = rotation * m_pose.m_referencePositions[i] + centerOfMass;
					n.m_x = Lerp(n.m_x, x, m_pose.m_poseMatching);
				}
			}
		}
	}
}

//
btVector3 btSoftBody::evaluateCenterOfMass() const
{
	btVector3 centerOfMass(0,0,0);
	
	if( m_pose.m_poseMatching != btScalar(0.0) )
	{
		for(int i = 0; i < m_nodes.size(); ++i) centerOfMass += m_nodes[i].m_x * m_pose.m_weights[i];
	}
	
	return centerOfMass;
}

//
btScalar btSoftBody::getClosedTrimeshVolume() const
{
	btScalar volume(0.0);
	
	if(m_nodes.size() > 0)
	{
		//Could also use (0,0,0) for origin; here, we assume that m_nodes[0] is closer so that
		//there will be less loss of floating point precision if the soft body is far from (0,0,0).
		const btVector3& origin = m_nodes[0].m_x;	
		for(int i = 0; i < m_faces.size(); ++i)
		{
			const Face&	f = m_faces[i];
			
			//Each triangle is combined with the origin to form a tetrahedron.
			//By summing up the signed volume of all tetrahedra, we get the volume of the (closed) triangle mesh.
			//Signed volume of tetrahedron is positive if the triangle's normal points away from the origin, negative otherwise.
			btScalar signedTetrahedronVolume = btDot( f.m_n[0]->m_x - origin, btCross(f.m_n[1]->m_x - origin, f.m_n[2]->m_x - origin) );
			
			volume += signedTetrahedronVolume;
		}
		volume /= btScalar(6.0);
	}
	
	return volume;
}


struct NodeLinks
{
    btAlignedObjectArray<int> m_links;
};



//
int				btSoftBody::generateBendingConstraints(int distance,Material* mat)
{
	int i,j;

	if(distance>1)
	{
		/* Build graph	*/ 
		const int		n=m_nodes.size();
		const unsigned	inf=(~(unsigned)0)>>1;
		unsigned*		adj=new unsigned[n*n];
		

#define IDX(_x_,_y_)	((_y_)*n+(_x_))
		for(j=0;j<n;++j)
		{
			for(i=0;i<n;++i)
			{
				if(i!=j)
				{
					adj[IDX(i,j)]=adj[IDX(j,i)]=inf;
				}
				else
				{
					adj[IDX(i,j)]=adj[IDX(j,i)]=0;
				}
			}
		}
		for( i=0;i<m_links.size();++i)
		{
			const int	ia=(int)(m_links[i].m_n[0]-&m_nodes[0]);
			const int	ib=(int)(m_links[i].m_n[1]-&m_nodes[0]);
			adj[IDX(ia,ib)]=1;
			adj[IDX(ib,ia)]=1;
		}


		//special optimized case for distance == 2
		if (distance == 2)
		{

			btAlignedObjectArray<NodeLinks> nodeLinks;


			/* Build node links */
			nodeLinks.resize(m_nodes.size());

			for( i=0;i<m_links.size();++i)
			{
				const int	ia=(int)(m_links[i].m_n[0]-&m_nodes[0]);
				const int	ib=(int)(m_links[i].m_n[1]-&m_nodes[0]);
				if (nodeLinks[ia].m_links.findLinearSearch(ib)==nodeLinks[ia].m_links.size())
					nodeLinks[ia].m_links.push_back(ib);

				if (nodeLinks[ib].m_links.findLinearSearch(ia)==nodeLinks[ib].m_links.size())
					nodeLinks[ib].m_links.push_back(ia);
			}
			for (int ii=0;ii<nodeLinks.size();ii++)
			{
				int i=ii;

				for (int jj=0;jj<nodeLinks[ii].m_links.size();jj++)
				{
					int k = nodeLinks[ii].m_links[jj];
					for (int kk=0;kk<nodeLinks[k].m_links.size();kk++)
					{
						int j = nodeLinks[k].m_links[kk];
						if (i!=j)
						{
							const unsigned	sum=adj[IDX(i,k)]+adj[IDX(k,j)];
							btAssert(sum==2);
							if(adj[IDX(i,j)]>sum)
							{
								adj[IDX(i,j)]=adj[IDX(j,i)]=sum;
							}
						}

					}
				}
			}
		} else
		{
			///generic Floyd's algorithm
			for(int k=0;k<n;++k)
			{
				for(j=0;j<n;++j)
				{
					for(i=j+1;i<n;++i)
					{
						const unsigned	sum=adj[IDX(i,k)]+adj[IDX(k,j)];
						if(adj[IDX(i,j)]>sum)
						{
							adj[IDX(i,j)]=adj[IDX(j,i)]=sum;
						}
					}
				}
			}
		}


		/* Build links	*/ 
		int	nlinks=0;
		for(j=0;j<n;++j)
		{
			for(i=j+1;i<n;++i)
			{
				if(adj[IDX(i,j)]==(unsigned)distance)
				{
					appendLink(i,j,mat);
					m_links[m_links.size()-1].m_bbending=1;
					++nlinks;
				}
			}
		}
		delete[] adj;		
		return(nlinks);
	}
	return(0);
}

//
void			btSoftBody::randomizeConstraints()
{
	unsigned long	seed=243703;
#define NEXTRAND (seed=(1664525L*seed+1013904223L)&0xffffffff)
	int i,ni;

	for(i=0,ni=m_links.size();i<ni;++i)
	{
		btSwap(m_links[i],m_links[NEXTRAND%ni]);
	}
	for(i=0,ni=m_faces.size();i<ni;++i)
	{
		btSwap(m_faces[i],m_faces[NEXTRAND%ni]);
	}
#undef NEXTRAND
}

//
void			btSoftBody::refine(ImplicitFn* ifn,btScalar accurary,bool cut)
{
	const Node*			nbase = &m_nodes[0];
	int					ncount = m_nodes.size();
	btSymMatrix<int>	edges(ncount,-2);
	int					newnodes=0;
	int i,j,k,ni;

	/* Filter out		*/ 
	for(i=0;i<m_links.size();++i)
	{
		Link&	l=m_links[i];
		if(l.m_bbending)
		{
			if(!SameSign(ifn->Eval(l.m_n[0]->m_x),ifn->Eval(l.m_n[1]->m_x)))
			{
				btSwap(m_links[i],m_links[m_links.size()-1]);
				m_links.pop_back();--i;
			}
		}	
	}
	/* Fill edges		*/ 
	for(i=0;i<m_links.size();++i)
	{
		Link&	l=m_links[i];
		edges(int(l.m_n[0]-nbase),int(l.m_n[1]-nbase))=-1;
	}
	for(i=0;i<m_faces.size();++i)
	{	
		Face&	f=m_faces[i];
		edges(int(f.m_n[0]-nbase),int(f.m_n[1]-nbase))=-1;
		edges(int(f.m_n[1]-nbase),int(f.m_n[2]-nbase))=-1;
		edges(int(f.m_n[2]-nbase),int(f.m_n[0]-nbase))=-1;
	}
	/* Intersect		*/ 
	for(i=0;i<ncount;++i)
	{
		for(j=i+1;j<ncount;++j)
		{
			if(edges(i,j)==-1)
			{
				Node&			a=m_nodes[i];
				Node&			b=m_nodes[j];
				const btScalar	t=ImplicitSolve(ifn,a.m_x,b.m_x,accurary);
				if(t>0)
				{
					const btVector3	x=Lerp(a.m_x,b.m_x,t);
					const btVector3	v=Lerp(a.m_v,b.m_v,t);
					btScalar		m=0;
					if(a.m_invMass > 0)
					{
						if(b.m_invMass > 0)
						{
							const btScalar	ma = 1 / a.m_invMass;
							const btScalar	mb = 1 / b.m_invMass;
							const btScalar	mc=Lerp(ma,mb,t);
							const btScalar	f=(ma+mb)/(ma+mb+mc);
							a.m_invMass = 1 / (ma * f);
							b.m_invMass = 1 / (mb * f);
							m=mc*f;
						}
						else
						{
							a.m_invMass /= 0.5f;
							m = 1 / a.m_invMass;
						}
					}
					else
					{
						if(b.m_invMass > 0)
						{
							b.m_invMass /= 0.5f;
							m = 1 / b.m_invMass;
						}
						else
							m=0;
					}
					appendNode(x,m);
					edges(i,j)=m_nodes.size()-1;
					m_nodes[edges(i,j)].m_v=v;
					++newnodes;
				}
			}
		}
	}
	nbase=&m_nodes[0];
	/* Refine links		*/ 
	for(i=0,ni=m_links.size();i<ni;++i)
	{
		Link&		feat=m_links[i];
		const int	idx[]={	int(feat.m_n[0]-nbase),
			int(feat.m_n[1]-nbase)};
		if((idx[0]<ncount)&&(idx[1]<ncount))
		{
			const int ni=edges(idx[0],idx[1]);
			if(ni>0)
			{
				appendLink(i);
				Link*		pft[]={	&m_links[i],
					&m_links[m_links.size()-1]};			
				pft[0]->m_n[0]=&m_nodes[idx[0]];
				pft[0]->m_n[1]=&m_nodes[ni];
				pft[1]->m_n[0]=&m_nodes[ni];
				pft[1]->m_n[1]=&m_nodes[idx[1]];
			}
		}
	}
	/* Refine faces		*/ 
	for(i=0;i<m_faces.size();++i)
	{
		const Face&	feat=m_faces[i];
		const int	idx[]={	int(feat.m_n[0]-nbase),
			int(feat.m_n[1]-nbase),
			int(feat.m_n[2]-nbase)};
		for(j=2,k=0;k<3;j=k++)
		{
			if((idx[j]<ncount)&&(idx[k]<ncount))
			{
				const int ni=edges(idx[j],idx[k]);
				if(ni>0)
				{
					appendFace(i);
					const int	l=(k+1)%3;
					Face*		pft[]={	&m_faces[i],
						&m_faces[m_faces.size()-1]};
					pft[0]->m_n[0]=&m_nodes[idx[l]];
					pft[0]->m_n[1]=&m_nodes[idx[j]];
					pft[0]->m_n[2]=&m_nodes[ni];
					pft[1]->m_n[0]=&m_nodes[ni];
					pft[1]->m_n[1]=&m_nodes[idx[k]];
					pft[1]->m_n[2]=&m_nodes[idx[l]];
					appendLink(ni,idx[l],pft[0]->m_material);
					--i;break;
				}
			}
		}
	}
	/* Cut				*/ 
	if(cut)
	{	
		btAlignedObjectArray<int>	cnodes;
		const int					pcount=ncount;
		int							i;
		ncount=m_nodes.size();
		cnodes.resize(ncount,0);
		/* Nodes		*/ 
		for(i=0;i<ncount;++i)
		{
			const btVector3	x=m_nodes[i].m_x;
			if((i>=pcount)||(btFabs(ifn->Eval(x))<accurary))
			{
				const btVector3	v=m_nodes[i].m_v;
				btScalar		m=getMass(i);
				if(m > 0)
				{
					m*=0.5f;
					m_nodes[i].m_invMass /= 0.5f;
				}
				appendNode(x,m);
				cnodes[i]=m_nodes.size()-1;
				m_nodes[cnodes[i]].m_v=v;
			}
		}
		nbase=&m_nodes[0];
		/* Links		*/ 
		for(i=0,ni=m_links.size();i<ni;++i)
		{
			const int		id[]={	int(m_links[i].m_n[0]-nbase),
				int(m_links[i].m_n[1]-nbase)};
			int				todetach=0;
			if(cnodes[id[0]]&&cnodes[id[1]])
			{
				appendLink(i);
				todetach=m_links.size()-1;
			}
			else
			{
				if((	(ifn->Eval(m_nodes[id[0]].m_x)<accurary)&&
					(ifn->Eval(m_nodes[id[1]].m_x)<accurary)))
					todetach=i;
			}
			if(todetach)
			{
				Link&	l=m_links[todetach];
				for(int j=0;j<2;++j)
				{
					int cn=cnodes[int(l.m_n[j]-nbase)];
					if(cn) l.m_n[j]=&m_nodes[cn];
				}			
			}
		}
		/* Faces		*/ 
		for(i=0,ni=m_faces.size();i<ni;++i)
		{
			Node**			n=	m_faces[i].m_n;
			if(	(ifn->Eval(n[0]->m_x)<accurary)&&
				(ifn->Eval(n[1]->m_x)<accurary)&&
				(ifn->Eval(n[2]->m_x)<accurary))
			{
				for(int j=0;j<3;++j)
				{
					int cn=cnodes[int(n[j]-nbase)];
					if(cn) n[j]=&m_nodes[cn];
				}
			}
		}
		/* Clean orphans	*/ 
		int							nnodes=m_nodes.size();
		btAlignedObjectArray<int>	ranks;
		btAlignedObjectArray<int>	todelete;
		ranks.resize(nnodes,0);
		for(i=0,ni=m_links.size();i<ni;++i)
		{
			for(int j=0;j<2;++j) ranks[int(m_links[i].m_n[j]-nbase)]++;
		}
		for(i=0,ni=m_faces.size();i<ni;++i)
		{
			for(int j=0;j<3;++j) ranks[int(m_faces[i].m_n[j]-nbase)]++;
		}
		for(i=0;i<m_links.size();++i)
		{
			const int	id[]={	int(m_links[i].m_n[0]-nbase),
				int(m_links[i].m_n[1]-nbase)};
			const bool	sg[]={	ranks[id[0]]==1,
				ranks[id[1]]==1};
			if(sg[0]||sg[1])
			{
				--ranks[id[0]];
				--ranks[id[1]];
				btSwap(m_links[i],m_links[m_links.size()-1]);
				m_links.pop_back();--i;
			}
		}
#if 0	
		for(i=nnodes-1;i>=0;--i)
		{
			if(!ranks[i]) todelete.push_back(i);
		}	
		if(todelete.size())
		{		
			btAlignedObjectArray<int>&	map=ranks;
			for(int i=0;i<nnodes;++i) map[i]=i;
			PointersToIndices(this);
			for(int i=0,ni=todelete.size();i<ni;++i)
			{
				int		j=todelete[i];
				int&	a=map[j];
				int&	b=map[--nnodes];
				m_nodeBvh.remove(m_nodes[a].m_leaf);m_nodes[a].m_leaf=0;
				btSwap(m_nodes[a],m_nodes[b]);
				j=a;a=b;b=j;			
			}
			IndicesToPointers(this,&map[0]);
			m_nodes.resize(nnodes);
		}
#endif
	}
	m_bUpdateRtCst=true;
}

//
bool			btSoftBody::cutLink(const Node* node0,const Node* node1,btScalar position)
{
	return(cutLink(int(node0-&m_nodes[0]),int(node1-&m_nodes[0]),position));
}

//
bool			btSoftBody::cutLink(int node0,int node1,btScalar position)
{
	bool			done=false;
	int i,ni;
//	const btVector3	d=m_nodes[node0].m_x-m_nodes[node1].m_x;
	const btVector3	x=Lerp(m_nodes[node0].m_x,m_nodes[node1].m_x,position);
	const btVector3	v=Lerp(m_nodes[node0].m_v,m_nodes[node1].m_v,position);
	const btScalar	m=1;
	appendNode(x,m);
	appendNode(x,m);
	Node*			pa=&m_nodes[node0];
	Node*			pb=&m_nodes[node1];
	Node*			pn[2]={	&m_nodes[m_nodes.size()-2],
		&m_nodes[m_nodes.size()-1]};
	pn[0]->m_v=v;
	pn[1]->m_v=v;
	for(i=0,ni=m_links.size();i<ni;++i)
	{
		const int mtch=MatchEdge(m_links[i].m_n[0],m_links[i].m_n[1],pa,pb);
		if(mtch!=-1)
		{
			appendLink(i);
			Link*	pft[]={&m_links[i],&m_links[m_links.size()-1]};
			pft[0]->m_n[1]=pn[mtch];
			pft[1]->m_n[0]=pn[1-mtch];
			done=true;
		}
	}
	for(i=0,ni=m_faces.size();i<ni;++i)
	{
		for(int k=2,l=0;l<3;k=l++)
		{
			const int mtch=MatchEdge(m_faces[i].m_n[k],m_faces[i].m_n[l],pa,pb);
			if(mtch!=-1)
			{
				appendFace(i);
				Face*	pft[]={&m_faces[i],&m_faces[m_faces.size()-1]};
				pft[0]->m_n[l]=pn[mtch];
				pft[1]->m_n[k]=pn[1-mtch];
				appendLink(pn[0],pft[0]->m_n[(l+1)%3],pft[0]->m_material,true);
				appendLink(pn[1],pft[0]->m_n[(l+1)%3],pft[0]->m_material,true);
			}
		}
	}
	if(!done)
	{
		m_nodeBvh.remove(pn[0]->m_leaf);
		m_nodeBvh.remove(pn[1]->m_leaf);
		m_nodes.pop_back();
		m_nodes.pop_back();
	}
	return(done);
}

//
bool btSoftBody::rayTest(const btVector3& rayFrom, const btVector3& rayTo, sRayCast& results)
{
	if( m_faces.size() && m_faceBvh.empty() ) initializeFaceTree();

	results.body	=	this;
	results.fraction = 1.f;
	results.index	=	-1;

	return (rayTest(rayFrom, rayTo, results.fraction, results.index, false) != 0);
}

//
void btSoftBody::predictMotion(btScalar timeStep)
{
	if(m_bUpdateRtCst)
	{
		m_bUpdateRtCst=false;
		updateConstants();
		m_faceBvh.clear();
		initializeFaceTree();
	}

	// Prepare
	m_timeStep = timeStep;
	btScalar velocityMargin	= timeStep * 3;
	btScalar radialMargin = getCollisionShape()->getMargin();
	btScalar updateMargin = radialMargin*(btScalar)0.25;
	// Forces
	addVelocityAllNodes(m_worldInfo->m_gravity * timeStep);
	applyForces();
	// Integrate
	for(int i = 0; i < m_nodes.size(); ++i)
	{
		Node&	n=m_nodes[i];
		n.m_q	=	n.m_x;
		btVector3 deltaV = n.m_f * n.m_invMass * timeStep;
		{
			btScalar maxDisplacement = m_worldInfo->m_maxDisplacement;
			btScalar clampDeltaV = maxDisplacement / timeStep;
			for (int c = 0; c < 3; c++)
			{
				if (deltaV[c]>clampDeltaV)
				{
					deltaV[c] = clampDeltaV;
				}
				if (deltaV[c]<-clampDeltaV)
				{
					deltaV[c]=-clampDeltaV;
				}
			}
		}
		n.m_v	+=	deltaV;
		n.m_x	+=	n.m_v * timeStep;
		n.m_f	=	btVector3(0,0,0);
	}
	
	// Clusters
	//updateClusters();
	
	// Bounds
	updateBounds();	
	// Nodes
	ATTRIBUTE_ALIGNED16(btDbvtVolume)	vol;
	for(int i = 0; i < m_nodes.size(); ++i)
	{
		Node&	n=m_nodes[i];
		vol = btDbvtVolume::FromCR(n.m_x, radialMargin);
		m_nodeBvh.update(n.m_leaf, vol, n.m_v * velocityMargin, updateMargin);
	}
	// Faces
	if( !m_faceBvh.empty() )
	{
		for(int i = 0; i < m_faces.size(); ++i)
		{
			Face& f = m_faces[i];
			const btVector3	v = (f.m_n[0]->m_v + f.m_n[1]->m_v + f.m_n[2]->m_v) / 3;
			vol = VolumeOf(f, radialMargin);
			m_faceBvh.update(f.m_leaf, vol, v * velocityMargin, updateMargin);
		}
	}
	
	// Pose
	updateAndApplyPose();
	
	// Clear contacts
	m_rigidContacts.resize(0);
	m_softContacts.resize(0);
	// Optimize dbvt's
	m_nodeBvh.optimizeIncremental(1);
	m_faceBvh.optimizeIncremental(1);
}


/// RayFromToCaster takes a ray from, ray to (instead of direction!)
struct	RayFromToCaster : btDbvt::ICollide
{
	btVector3			m_rayFrom;
	btVector3			m_rayTo;
	btVector3			m_rayNormalizedDirection;
	btScalar			m_mint;
	btSoftBody::Face*	m_face;
	int					m_tests;
	RayFromToCaster(const btVector3& rayFrom,const btVector3& rayTo,btScalar mxt);
	void Process(const btDbvtNode* leaf);

	static inline btScalar	rayFromToTriangle(const btVector3& rayFrom, const btVector3& rayTo,
		const btVector3& rayNormalizedDirection,
		const btVector3& a, const btVector3& b, const btVector3& c,
		btScalar maxt=SIMD_INFINITY);
};
RayFromToCaster::RayFromToCaster(const btVector3& rayFrom,const btVector3& rayTo,btScalar mxt)
{
	m_rayFrom = rayFrom;
	m_rayNormalizedDirection = (rayTo-rayFrom);
	m_rayTo = rayTo;
	m_mint	=	mxt;
	m_face	=	0;
	m_tests	=	0;
}
void RayFromToCaster::Process(const btDbvtNode* leaf)
{
	btSoftBody::Face&	f=*(btSoftBody::Face*)leaf->data;
	const btScalar		t=rayFromToTriangle(	m_rayFrom,m_rayTo,m_rayNormalizedDirection,
		f.m_n[0]->m_x,
		f.m_n[1]->m_x,
		f.m_n[2]->m_x,
		m_mint);
	if((t>0)&&(t<m_mint)) 
	{ 
		m_mint=t;m_face=&f; 
	}
	++m_tests;
}
btScalar RayFromToCaster::rayFromToTriangle(const btVector3& rayFrom, const btVector3& rayTo,
											const btVector3& rayNormalizedDirection,
											const btVector3& a, const btVector3& b, const btVector3& c,
											btScalar maxt)
{
	static const btScalar	ceps=-SIMD_EPSILON*10;
	static const btScalar	teps=SIMD_EPSILON*10;

	const btVector3			n=btCross(b-a,c-a);
	const btScalar			d=btDot(a,n);
	const btScalar			den=btDot(rayNormalizedDirection,n);
	if(!btFuzzyZero(den))
	{
		const btScalar		num=btDot(rayFrom,n)-d;
		const btScalar		t=-num/den;
		if((t>teps)&&(t<maxt))
		{
			const btVector3	hit=rayFrom+rayNormalizedDirection*t;
			if(	(btDot(n,btCross(a-hit,b-hit))>ceps)	&&			
				(btDot(n,btCross(b-hit,c-hit))>ceps)	&&
				(btDot(n,btCross(c-hit,a-hit))>ceps))
			{
				return(t);
			}
		}
	}
	return(-1);
}

//
void				btSoftBody::pointersToIndices()
{
#define	PTR2IDX(_p_,_b_)	reinterpret_cast<btSoftBody::Node*>((_p_)-(_b_))
	btSoftBody::Node*	base=m_nodes.size() ? &m_nodes[0] : 0;
	int i,ni;

	for(i=0,ni=m_nodes.size();i<ni;++i)
	{
		if(m_nodes[i].m_leaf)
		{
			m_nodes[i].m_leaf->data=*(void**)&i;
		}
	}
	for(i=0,ni=m_links.size();i<ni;++i)
	{
		m_links[i].m_n[0]=PTR2IDX(m_links[i].m_n[0],base);
		m_links[i].m_n[1]=PTR2IDX(m_links[i].m_n[1],base);
	}
	for(i=0,ni=m_faces.size();i<ni;++i)
	{
		m_faces[i].m_n[0]=PTR2IDX(m_faces[i].m_n[0],base);
		m_faces[i].m_n[1]=PTR2IDX(m_faces[i].m_n[1],base);
		m_faces[i].m_n[2]=PTR2IDX(m_faces[i].m_n[2],base);
		if(m_faces[i].m_leaf)
		{
			m_faces[i].m_leaf->data=*(void**)&i;
		}
	}
	for(i=0,ni=m_anchors.size();i<ni;++i)
	{
		m_anchors[i].m_node=PTR2IDX(m_anchors[i].m_node,base);
	}
#undef	PTR2IDX
}

//
void				btSoftBody::indicesToPointers(const int* map)
{
#define	IDX2PTR(_p_,_b_)	map?(&(_b_)[map[(((char*)_p_)-(char*)0)]]):	\
	(&(_b_)[(((char*)_p_)-(char*)0)])
	btSoftBody::Node*	base=m_nodes.size() ? &m_nodes[0]:0;
	int i,ni;

	for(i=0,ni=m_nodes.size();i<ni;++i)
	{
		if(m_nodes[i].m_leaf)
		{
			m_nodes[i].m_leaf->data=&m_nodes[i];
		}
	}
	for(i=0,ni=m_links.size();i<ni;++i)
	{
		m_links[i].m_n[0]=IDX2PTR(m_links[i].m_n[0],base);
		m_links[i].m_n[1]=IDX2PTR(m_links[i].m_n[1],base);
	}
	for(i=0,ni=m_faces.size();i<ni;++i)
	{
		m_faces[i].m_n[0]=IDX2PTR(m_faces[i].m_n[0],base);
		m_faces[i].m_n[1]=IDX2PTR(m_faces[i].m_n[1],base);
		m_faces[i].m_n[2]=IDX2PTR(m_faces[i].m_n[2],base);
		if(m_faces[i].m_leaf)
		{
			m_faces[i].m_leaf->data=&m_faces[i];
		}
	}
	for(i=0,ni=m_anchors.size();i<ni;++i)
	{
		m_anchors[i].m_node=IDX2PTR(m_anchors[i].m_node,base);
	}
#undef	IDX2PTR
}

//
int					btSoftBody::rayTest(const btVector3& rayFrom,const btVector3& rayTo,
										btScalar& mint,int& index,bool bcountonly) const
{
	int	cnt=0;
	btVector3 dir = rayTo-rayFrom;
	

	if( bcountonly || m_faceBvh.empty() )
	{
		
		for(int i=0,ni=m_faces.size();i<ni;++i)
		{
			const btSoftBody::Face&	f=m_faces[i];

			const btScalar			t=RayFromToCaster::rayFromToTriangle(	rayFrom,rayTo,dir,
				f.m_n[0]->m_x,
				f.m_n[1]->m_x,
				f.m_n[2]->m_x,
				mint);
			if(t>0)
			{
				++cnt;
				if(!bcountonly)
				{
					index=i;
					mint=t;
				}
			}
		}
	}
	else	//Use dbvt
	{
		RayFromToCaster	collider(rayFrom,rayTo,mint);

		btDbvt::rayTest(m_faceBvh.m_root, rayFrom, rayTo, collider);
		if(collider.m_face)
		{
			mint=collider.m_mint;
			index=(int)(collider.m_face-&m_faces[0]);
			cnt=1;
		}
	}

	return(cnt);
}

//
void btSoftBody::initializeFaceTree()
{
	m_faceBvh.clear();
	for(int i = 0; i < m_faces.size(); ++i)
	{
		Face& f = m_faces[i];
		f.m_leaf = m_faceBvh.insert( VolumeOf(f,0), &f );
	}
}



//
bool btSoftBody::checkContact(const btCollisionObjectWrapper* colObjWrap,
							const btVector3& worldSpaceNodePosition,
							btScalar margin,
							btSoftBody::RigidContact& contact) const
{
	btVector3 normal;
	const btCollisionShape *shp = colObjWrap->getCollisionShape();
	const btTransform &wtr = colObjWrap->getWorldTransform();

	btScalar dst = m_worldInfo->m_sparsesdf.Evaluate(wtr.invXform(worldSpaceNodePosition), shp, normal, margin);
	if(dst < 0)
	{
		contact.m_colObj = colObjWrap->getCollisionObject();
		contact.m_normal = wtr.getBasis() * normal;
		contact.m_offset = -btDot( contact.m_normal, worldSpaceNodePosition - contact.m_normal * dst );
		return true;
	}
	return false;
}

//
void btSoftBody::updateNormals()
{
	for(int i = 0; i < m_nodes.size(); ++i)
	{
		m_nodes[i].m_normal = btVector3(0,0,0);
	}
	for(int i = 0; i < m_faces.size(); ++i)
	{
		btSoftBody::Face& f = m_faces[i];
		const btVector3 n = btCross(f.m_n[1]->m_x - f.m_n[0]->m_x, f.m_n[2]->m_x - f.m_n[0]->m_x);
		f.m_normal = n.normalized();
		f.m_n[0]->m_normal += n;
		f.m_n[1]->m_normal += n;
		f.m_n[2]->m_normal += n;
	}
	for(int i = 0; i < m_nodes.size(); ++i)
	{
		btScalar len = m_nodes[i].m_normal.length();
		if (len > SIMD_EPSILON) m_nodes[i].m_normal /= len;
	}
}

//
void btSoftBody::updateBounds()
{
	if(m_nodeBvh.m_root)
	{
		const btVector3& mins = m_nodeBvh.m_root->volume.Mins();
		const btVector3& maxs = m_nodeBvh.m_root->volume.Maxs();
		btScalar margin = getCollisionShape()->getMargin();
		btVector3 margin3 = btVector3(margin, margin, margin);
		
		m_aabbMin = mins - margin3;
		m_aabbMax = maxs + margin3;
		if( 0 != getBroadphaseHandle() )
			m_worldInfo->m_broadphase->setAabb(	getBroadphaseHandle(), m_aabbMin, m_aabbMax, m_worldInfo->m_dispatcher);
	}
	else
	{
		m_aabbMin = m_aabbMax = btVector3(0,0,0);
	}
}


//
void				btSoftBody::updateArea(bool averageArea)
{
	int i,ni;

	/* Face area		*/ 
	for(i=0,ni=m_faces.size();i<ni;++i)
	{
		Face&		f=m_faces[i];
		f.m_area = AreaOf(f.m_n[0]->m_x, f.m_n[1]->m_x, f.m_n[2]->m_x);
	}
	
	/* Node area		*/ 

	if (averageArea)
	{
		btAlignedObjectArray<int>	counts;
		counts.resize(m_nodes.size(),0);
		for(i=0,ni=m_nodes.size();i<ni;++i)
		{
			m_nodes[i].m_area	=	0;
		}
		for(i=0,ni=m_faces.size();i<ni;++i)
		{
			btSoftBody::Face&	f=m_faces[i];
			for(int j=0;j<3;++j)
			{
				const int index=(int)(f.m_n[j]-&m_nodes[0]);
				counts[index]++;
				f.m_n[j]->m_area += btFabs(f.m_area);
			}
		}
		for(i=0,ni=m_nodes.size();i<ni;++i)
		{
			if(counts[i]>0)
				m_nodes[i].m_area/=(btScalar)counts[i];
			else
				m_nodes[i].m_area=0;
		}
	}
	else
	{
		// initialize node area as zero
		for(i=0,ni=m_nodes.size();i<ni;++i)
		{
			m_nodes[i].m_area=0;	
		}

		for(i=0,ni=m_faces.size();i<ni;++i)
		{
			btSoftBody::Face&	f=m_faces[i];

			for(int j=0;j<3;++j)
			{
				f.m_n[j]->m_area += f.m_area;
			}
		}

		for(i=0,ni=m_nodes.size();i<ni;++i)
		{
			m_nodes[i].m_area *= 0.3333333f;
		}
	}
}


void btSoftBody::updateConstants()
{
	//Set current link lengths as resting lengths
	for(int i = 0; i < m_links.size(); ++i)
	{
		Link& l = m_links[i];
		l.m_restLength = (l.m_n[0]->m_x - l.m_n[1]->m_x).length();
		l.m_restLengthSquared = l.m_restLength * l.m_restLength;
	}
	
	//Update link constants
	for(int i = 0; i < m_links.size(); ++i)
	{
		Link& l= m_links[i];
		Material& m = *l.m_material;
		l.m_scaledCombinedInvMass = (l.m_n[0]->m_invMass + l.m_n[1]->m_invMass) / m.m_linearStiffness;
	}
	
	updateArea();
}


//
void btSoftBody::applyForces()
{
	BT_PROFILE("SoftBody applyForces");
	
	m_closedTrimeshForce.applyForcesToNodes( getClosedTrimeshVolume(), m_nodes );

	addAeroForces(m_aeroForce, m_timeStep, m_nodes, m_faces);
}

//
void btSoftBody::defaultCollisionHandler(const btCollisionObjectWrapper* pcoWrap)
{
	btSoftColliders::CollideSDF_RS	docollide;
	btRigidBody* prb1=(btRigidBody*) btRigidBody::upcast(pcoWrap->getCollisionObject());
	btTransform	wtr=pcoWrap->getWorldTransform();

	const btTransform ctr = pcoWrap->getWorldTransform();
	const btScalar timemargin = (wtr.getOrigin() - ctr.getOrigin()).length();
	const btScalar basemargin = getCollisionShape()->getMargin();
	btVector3 mins;
	btVector3 maxs;
	ATTRIBUTE_ALIGNED16(btDbvtVolume) volume;
	pcoWrap->getCollisionShape()->getAabb(pcoWrap->getWorldTransform(), mins, maxs);
	volume=btDbvtVolume::FromMM(mins,maxs);
	volume.Expand(btVector3(basemargin,basemargin,basemargin));
	docollide.psb =	this;
	docollide.m_colObj1Wrap = pcoWrap;
	docollide.m_rigidBody = prb1;

	docollide.dynmargin	= basemargin + timemargin;
	docollide.stamargin	= basemargin;
	m_nodeBvh.collideTV(m_nodeBvh.m_root, volume, docollide);
}

//
void btSoftBody::defaultCollisionHandler(btSoftBody* psb)
{
	//only self-collision for Cluster, not Vertex-Face yet
	if (this != psb)
	{
		btSoftColliders::CollideVF_SS docollide;
		
		// common
		docollide.mrg =	getCollisionShape()->getMargin()+ psb->getCollisionShape()->getMargin();
		
		// psb0 nodes vs psb1 faces
		docollide.psb[0] = this;
		docollide.psb[1] = psb;
		docollide.psb[0]->m_nodeBvh.collideTT(docollide.psb[0]->m_nodeBvh.m_root, docollide.psb[1]->m_faceBvh.m_root, docollide);
		
		// psb1 nodes vs psb0 faces
		docollide.psb[0] = psb;
		docollide.psb[1] = this;
		docollide.psb[0]->m_nodeBvh.collideTT(docollide.psb[0]->m_nodeBvh.m_root, docollide.psb[1]->m_faceBvh.m_root, docollide);
	}
}

/************************************************************************************
///Aero force
************************************************************************************/
static inline void ApplyClampedForce(btSoftBody::Node& n, const btVector3& f, btScalar dt)
{
	const btScalar	dtim=dt * n.m_invMass;
	if((f*dtim).length2()>n.m_v.length2()) n.m_f-=ProjectOnAxis(n.m_v,f.normalized())/dtim;	//Apply clamped force
	else n.m_f+=f;	//Unclamped force 
}

void btSoftBody::addAeroForces(const AeroForce& aeroForce, btScalar timeStep, btAlignedObjectArray<Node>& nodes, btAlignedObjectArray<Face>& faces)
{
	const bool use_medium =	(aeroForce.m_liftCoeff > 0) || (aeroForce.m_dragCoeff > 0);
	if(use_medium)
	{
		for(int i = 0; i < nodes.size(); ++i)
			if(nodes[i].m_invMass > 0) addAeroForceToNode(aeroForce, timeStep, nodes, i);
		
		for(int i = 0; i < faces.size(); ++i) addAeroForceToFace(aeroForce, timeStep, faces, i);	 
	}
}

void btSoftBody::addAeroForceToNode(const AeroForce& aeroForce, btScalar timeStep, btAlignedObjectArray<Node>& nodes, int nodeIndex)
{
	const btScalar kLF = aeroForce.m_liftCoeff;
	const btScalar kDG = aeroForce.m_dragCoeff;
	
	const bool as_lift = kLF>0;
	const bool as_drag = kDG>0;
	const bool as_aero = as_lift || as_drag;
	const bool as_vaero = as_aero && (aeroForce.m_model < btSoftBody::eAeroModel::F_TwoSided);

	Node& n = nodes[nodeIndex];

	if(n.m_invMass > 0)
	{
		//Aerodynamics
		if(as_vaero)
		{				
			const btVector3	rel_v = n.m_v - aeroForce.m_windVelocity;
			const btScalar rel_v_len = rel_v.length();
			const btScalar	rel_v2 = rel_v.length2();

			if(rel_v2>SIMD_EPSILON)
			{
				const btVector3 rel_v_nrm = rel_v.normalized();
				btVector3	nrm = n.m_normal;

				if (aeroForce.m_model == btSoftBody::eAeroModel::V_TwoSidedLiftDrag)
				{
					nrm *= (btScalar)( (btDot(nrm,rel_v) < 0) ? -1 : +1);
					btVector3 fDrag(0, 0, 0);
					btVector3 fLift(0, 0, 0);

					btScalar n_dot_v = nrm.dot(rel_v_nrm);
					btScalar tri_area = 0.5f * n.m_area;
							
					fDrag = 0.5f * kDG * aeroForce.m_airDensity * rel_v2 * tri_area * n_dot_v * (-rel_v_nrm);
							
					// Check angle of attack
					// cos(10º) = 0.98480
					if ( 0 < n_dot_v && n_dot_v < 0.98480f)
						fLift = 0.5f * kLF * aeroForce.m_airDensity * rel_v_len * tri_area * btSqrt(1.0f-n_dot_v*n_dot_v) * (nrm.cross(rel_v_nrm).cross(rel_v_nrm));

					// Check if the velocity change resulted by aero drag force exceeds the current velocity of the node.
					btVector3 del_v_by_fDrag = fDrag * n.m_invMass * timeStep;
					btScalar del_v_by_fDrag_len2 = del_v_by_fDrag.length2();
					btScalar v_len2 = n.m_v.length2();

					if (del_v_by_fDrag_len2 >= v_len2 && del_v_by_fDrag_len2 > 0)
					{
						btScalar del_v_by_fDrag_len = del_v_by_fDrag.length();
						btScalar v_len = n.m_v.length();
						fDrag *= btScalar(0.8)*(v_len / del_v_by_fDrag_len);
					}

					n.m_f += fDrag;
					n.m_f += fLift;
				}
				else if (aeroForce.m_model == btSoftBody::eAeroModel::V_Point || aeroForce.m_model == btSoftBody::eAeroModel::V_TwoSided)
				{
					if (btSoftBody::eAeroModel::V_TwoSided)
						nrm *= (btScalar)( (btDot(nrm,rel_v) < 0) ? -1 : +1);

					const btScalar dvn = btDot(rel_v,nrm);
					//Compute forces
					if(dvn>0)
					{
						btVector3		force(0,0,0);
						const btScalar	c0	=	n.m_area * dvn * rel_v2/2;
						const btScalar	c1	=	c0 * aeroForce.m_airDensity;
						force	+=	nrm*(-c1*kLF);
						force	+=	rel_v.normalized() * (-c1 * kDG);
						ApplyClampedForce(n, force, timeStep);
					}
				}	
			}
		}
	}
}

void btSoftBody::addAeroForceToFace(const AeroForce& aeroForce, btScalar timeStep, btAlignedObjectArray<Face>& faces, int faceIndex)
{
	const btScalar kLF = aeroForce.m_liftCoeff;
	const btScalar kDG = aeroForce.m_dragCoeff;
	
	const bool as_lift = kLF>0;
	const bool as_drag = kDG>0;
	const bool as_aero = as_lift || as_drag;
	const bool as_faero = as_aero && (aeroForce.m_model >= btSoftBody::eAeroModel::F_TwoSided);

	if(as_faero)
	{
		btSoftBody::Face&	f=faces[faceIndex];

		const btVector3	v=(f.m_n[0]->m_v+f.m_n[1]->m_v+f.m_n[2]->m_v)/3;
		
		const btVector3	rel_v = v - aeroForce.m_windVelocity;
		const btScalar rel_v_len = rel_v.length();
		const btScalar	rel_v2=rel_v.length2();

		if(rel_v2>SIMD_EPSILON)
		{
			const btVector3 rel_v_nrm = rel_v.normalized();
			btVector3	nrm = f.m_normal;

			if (aeroForce.m_model == btSoftBody::eAeroModel::F_TwoSidedLiftDrag)
			{
				nrm *= (btScalar)( (btDot(nrm,rel_v) < 0) ? -1 : +1);

				btVector3 fDrag(0, 0, 0);
				btVector3 fLift(0, 0, 0);

				btScalar n_dot_v = nrm.dot(rel_v_nrm);
				btScalar tri_area = 0.5f * f.m_area;
					
				fDrag = 0.5f * kDG * aeroForce.m_airDensity * rel_v2 * tri_area * n_dot_v * (-rel_v_nrm);

				// Check angle of attack
				// cos(10º) = 0.98480
				if ( 0 < n_dot_v && n_dot_v < 0.98480f)
					fLift = 0.5f * kLF * aeroForce.m_airDensity * rel_v_len * tri_area * btSqrt(1.0f-n_dot_v*n_dot_v) * (nrm.cross(rel_v_nrm).cross(rel_v_nrm));

				fDrag /= 3;
				fLift /= 3;

				for(int j=0;j<3;++j) 
				{
					if (f.m_n[j]->m_invMass > 0)
					{
						// Check if the velocity change resulted by aero drag force exceeds the current velocity of the node.
						btVector3 del_v_by_fDrag = fDrag * f.m_n[j]->m_invMass * timeStep;
						btScalar del_v_by_fDrag_len2 = del_v_by_fDrag.length2();
						btScalar v_len2 = f.m_n[j]->m_v.length2();

						if (del_v_by_fDrag_len2 >= v_len2 && del_v_by_fDrag_len2 > 0)
						{
							btScalar del_v_by_fDrag_len = del_v_by_fDrag.length();
							btScalar v_len = f.m_n[j]->m_v.length();
							fDrag *= btScalar(0.8)*(v_len / del_v_by_fDrag_len);
						}

						f.m_n[j]->m_f += fDrag; 
						f.m_n[j]->m_f += fLift;
					}
				}
			}
			else if (aeroForce.m_model == btSoftBody::eAeroModel::F_OneSided || aeroForce.m_model == btSoftBody::eAeroModel::F_TwoSided)
			{
				if (btSoftBody::eAeroModel::F_TwoSided)
					nrm *= (btScalar)( (btDot(nrm,rel_v) < 0) ? -1 : +1);

				const btScalar	dvn=btDot(rel_v,nrm);
				/* Compute forces	*/ 
				if(dvn>0)
				{
					btVector3		force(0,0,0);
					const btScalar	c0	=	f.m_area * dvn * rel_v2;
					const btScalar	c1	=	c0*aeroForce.m_airDensity;
					force	+=	nrm*(-c1*kLF);
					force	+=	rel_v.normalized()*(-c1*kDG);
					force	/=	3;
					for(int j=0;j<3;++j) ApplyClampedForce(*f.m_n[j],force,timeStep);
				}
			}
		}
	}

}
/************************************************************************************
///Aero
************************************************************************************/