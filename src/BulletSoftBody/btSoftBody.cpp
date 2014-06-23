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
		n.m_leaf	=	m_nodeBvh.insert( btDbvtVolume::FromCR(n.m_x,margin), reinterpret_cast<void*>(i) );
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
	for(int i = 0; i < m_links.size(); ++i)
	{
		const Link&	l = m_links[i];
		if(	(l.m_linkIndicies[0] == node0 && l.m_linkIndicies[1] == node1) ||
			(l.m_linkIndicies[0] == node1 && l.m_linkIndicies[1] == node0)) return true;
	}
	
	return false;
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
	if( m_nodes.capacity() == m_nodes.size() ) m_nodes.reserve( m_nodes.size() * 2 + 1 );
	
	const btScalar	margin=getCollisionShape()->getMargin();
	m_nodes.push_back(Node());
	
	int nodeIndex = m_nodes.size() - 1;
	Node&			n = m_nodes[nodeIndex];
	ZeroInitialize(n);
	n.m_x			=	position;
	n.m_q			=	n.m_x;
	n.m_invMass 	= (mass > 0) ? 1 / mass : 0;
	n.m_material	=	m_materials[0];
	n.m_leaf		=	m_nodeBvh.insert( btDbvtVolume::FromCR(n.m_x,margin), reinterpret_cast<void*>(nodeIndex) );
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

	if((!bcheckexist)||(!checkLink(node0,node1)))
	{
		appendLink(-1,mat);
		Link& l = m_links[ m_links.size()-1 ];
		l.m_linkIndicies[0] = node0;
		l.m_linkIndicies[1] = node1;
		l.m_restLength = (m_nodes[node0].m_x - m_nodes[node1].m_x).length();
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
	btAssert(node0!=node1);
	btAssert(node1!=node2);
	btAssert(node2!=node0);
	
	if (node0==node1)
		return;
	if (node1==node2)
		return;
	if (node2==node0)
		return;

	appendFace(-1,mat);
	Face&	f=m_faces[m_faces.size()-1];
	f.m_indicies[0] = node0;
	f.m_indicies[1] = node1;
	f.m_indicies[2] = node2;
	
	f.m_area = AreaOfParallelogram(m_nodes[node0].m_x, m_nodes[node1].m_x, m_nodes[node2].m_x) * btScalar(0.5f);
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
	a.m_nodeIndex		=	node;
	a.m_body			=	body;
	a.m_local			=	localPivot;
	a.m_influence = influence;
	m_anchors.push_back(a);

	m_nodes[node].m_battach = 1;
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
			const Face& face = m_faces[i];
			btSoftBody::Node* faceNodes[3] = { &m_nodes[ face.m_indicies[0] ], &m_nodes[ face.m_indicies[1] ], &m_nodes[ face.m_indicies[2] ] };
			
			btScalar twicearea = AreaOfParallelogram(faceNodes[0]->m_x, faceNodes[1]->m_x, faceNodes[2]->m_x);
			for(int j = 0; j < 3; ++j) faceNodes[j]->m_invMass += twicearea;
			
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
			
			btVector3 a = m_nodes[ f.m_indicies[0] ].m_x - origin;
			btVector3 b = m_nodes[ f.m_indicies[1] ].m_x - origin;
			btVector3 c = m_nodes[ f.m_indicies[2] ].m_x - origin;
			
			//Each triangle is combined with the origin to form a tetrahedron.
			//By summing up the signed volume of all tetrahedra, we get the volume of the (closed) triangle mesh.
			//Signed volume of tetrahedron is positive if the triangle's normal points away from the origin, negative otherwise.
			btScalar signedTetrahedronVolume = btDot( a, btCross(b, c) );
			
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
		// Build graph	 
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
			const int ia = m_links[i].m_linkIndicies[0];
			const int ib = m_links[i].m_linkIndicies[1];
			adj[IDX(ia,ib)]=1;
			adj[IDX(ib,ia)]=1;
		}


		//special optimized case for distance == 2
		if (distance == 2)
		{

			btAlignedObjectArray<NodeLinks> nodeLinks;


			// Build node links 
			nodeLinks.resize(m_nodes.size());

			for( i=0;i<m_links.size();++i)
			{
				const int ia = m_links[i].m_linkIndicies[0];
				const int ib = m_links[i].m_linkIndicies[1];
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


		// Build links	 
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
void btSoftBody::refine(ImplicitFn* shape, btScalar accuracy, bool cut)
{
	int numNodesInitial = m_nodes.size();
	
	//Symmetric matrix of size (numNodesInitial x numNodesInitial), initialized with all values = -2
	btSymMatrix<int> edges(numNodesInitial, -2);
	
	//Remove all bending links that intersect the shape; links are not removed if entirely inside or outside the shape
	for(int i = 0; i < m_links.size(); ++i)
	{
		Link& l = m_links[i];
		
		const btVector3& nodePosition0 = m_nodes[ l.m_linkIndicies[0] ].m_x;
		const btVector3& nodePosition1 = m_nodes[ l.m_linkIndicies[1] ].m_x;
		
		if(l.m_bbending)
		{
			if( !SameSign(shape->signedDistance(nodePosition0), shape->signedDistance(nodePosition1)) )
			{
				btSwap( m_links[i], m_links[m_links.size() - 1] );
				m_links.pop_back();
				--i;
			}
		}	
	}
	
	//If 2 nodes share a link, set their entries in the matrix to -1
	for(int i = 0; i < m_links.size(); ++i)
	{
		Link& l = m_links[i];
		edges(l.m_linkIndicies[0], l.m_linkIndicies[1]) = -1;
	}
	
	//If 3 nodes share a face, set their entries in the matrix to -1
	for(int i = 0; i < m_faces.size(); ++i)
	{	
		Face& f = m_faces[i];
		edges(f.m_indicies[0], f.m_indicies[1]) = -1;
		edges(f.m_indicies[1], f.m_indicies[2]) = -1;
		edges(f.m_indicies[2], f.m_indicies[0]) = -1;
	}
	
	//Intersect
	//If 2 nodes a and b share a link or face, and the line from a to b intersects the implicit shape's surface,
	//create a new node at the point of intersection.
	for(int i = 0; i < numNodesInitial; ++i)
	{
		for(int j = i + 1; j < numNodesInitial; ++j)
		{
			if( edges(i,j) == -1 )		//If nodes i and j share a link or face
			{
				Node& a = m_nodes[i];
				Node& b = m_nodes[j];
				const btScalar	t=ImplicitSolve(shape,a.m_x,b.m_x,accuracy);
				
				if(t>0)
				{
					btScalar mass=0;
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
							mass=mc*f;
						}
						else
						{
							a.m_invMass /= 0.5f;
							mass = 1 / a.m_invMass;
						}
					}
					else
					{
						if(b.m_invMass > 0)
						{
							b.m_invMass /= 0.5f;
							mass = 1 / b.m_invMass;
						}
						else mass = 0;
					}
					
					btVector3 position = Lerp(a.m_x,b.m_x,t);
					btVector3 velocity = Lerp(a.m_v,b.m_v,t);
					
					appendNode(position, mass);
					
					int newNodeIndex = m_nodes.size() - 1;
					edges(i,j) = newNodeIndex;
					m_nodes[newNodeIndex].m_v = velocity;
				}
			}
		}
	}
	
	//Refine links - split each link that intersects the implicit shape
	for(int i = 0, numLinks = m_links.size(); i < numLinks; ++i)
	{
		Link& link = m_links[i];
		int linkIndex0 = link.m_linkIndicies[0];
		int linkIndex1 = link.m_linkIndicies[1];
		
		if(linkIndex0 < numNodesInitial && linkIndex1 < numNodesInitial)		//If both nodes are not newly created in the intersect stage above
		{
			const int ni = edges(linkIndex0, linkIndex1);
			if(ni > 0)		//If a new node was created between these nodes
			{
				appendLink(i);
				Link* pft[] = {	&m_links[i], &m_links[m_links.size()-1] };			
				pft[0]->m_linkIndicies[0] = linkIndex0;
				pft[0]->m_linkIndicies[1] = ni;
				pft[1]->m_linkIndicies[0] = ni;
				pft[1]->m_linkIndicies[1] = linkIndex1;
			}
		}
	}
	
	//Refine faces
	for(int i = 0; i < m_faces.size(); ++i)
	{
		const Face&	face = m_faces[i];
		const int idx[] = {	face.m_indicies[0], face.m_indicies[1], face.m_indicies[2] };
		for(int j = 2, k = 0; k < 3; j = k++)
		{
			if(idx[j] < numNodesInitial && idx[k] < numNodesInitial)
			{
				const int ni=edges(idx[j],idx[k]);
				if(ni>0)
				{
					appendFace(i);
					const int	l=(k+1)%3;
					Face* pft[] = { &m_faces[i], &m_faces[m_faces.size()-1] };
					pft[0]->m_indicies[0] = idx[l];
					pft[0]->m_indicies[1] = idx[j];
					pft[0]->m_indicies[2] = ni;
					pft[1]->m_indicies[0] = ni;
					pft[1]->m_indicies[1] = idx[k];
					pft[1]->m_indicies[2] = idx[l];
					appendLink(ni,idx[l],pft[0]->m_material);
					--i;
					break;
				}
			}
		}
	}
	
	// Cut
	if(cut)
	{	
		int currNodeCount = m_nodes.size();
		
		btAlignedObjectArray<int> cnodes;
		cnodes.resize(currNodeCount,0);
		
		// Nodes		 
		for(int i = 0; i < currNodeCount; ++i)
		{
			const btVector3	nodePosition = m_nodes[i].m_x;
			if( (i >= numNodesInitial) || (btFabs(shape->signedDistance(nodePosition)) < accuracy) )
			{
				btScalar m = getMass(i);
				if(m > 0)
				{
					m *= 0.5f;
					m_nodes[i].m_invMass /= 0.5f;
				}
				appendNode(nodePosition, m);
				cnodes[i] = m_nodes.size() - 1;
				m_nodes[cnodes[i]].m_v = m_nodes[i].m_v;
			}
		}
		
		// Links
		for(int i = 0, numLinks = m_links.size(); i < numLinks; ++i)
		{
			int linkIndex0 = m_links[i].m_linkIndicies[0];
			int linkIndex1 = m_links[i].m_linkIndicies[1];
			
			int	todetach = 0;
			if(cnodes[linkIndex0] && cnodes[linkIndex1])
			{
				appendLink(i);
				todetach=m_links.size()-1;
			}
			else
			{
				if(	(shape->signedDistance(m_nodes[linkIndex0].m_x) < accuracy) &&
					(shape->signedDistance(m_nodes[linkIndex1].m_x) < accuracy) ) todetach = i;
			}
			
			if(todetach)
			{
				Link&	l=m_links[todetach];
				for(int j = 0; j < 2; ++j)
				{
					int nodeIndex = cnodes[ l.m_linkIndicies[j] ];
					if(nodeIndex) l.m_linkIndicies[j] = nodeIndex;
				}			
			}
		}
		
		// Faces		 
		for(int i = 0; i < m_faces.size(); ++i)
		{	
			btSoftBody::Face& face = m_faces[i];
			btVector3 nodePosition0 = m_nodes[ face.m_indicies[0] ].m_x;
			btVector3 nodePosition1 = m_nodes[ face.m_indicies[1] ].m_x;
			btVector3 nodePosition2 = m_nodes[ face.m_indicies[2] ].m_x;
			
			if(	(shape->signedDistance(nodePosition0) < accuracy) &&
				(shape->signedDistance(nodePosition1) < accuracy) &&
				(shape->signedDistance(nodePosition2) < accuracy) )
			{
				for(int j = 0; j < 3; ++j)
				{
					int nodeIndex = cnodes[ face.m_indicies[j] ];
					if(nodeIndex) face.m_indicies[j] = nodeIndex;
				}
			}
		}
		
		// Clean orphans
		btAlignedObjectArray<int> referencesPerNode;
		referencesPerNode.resize( m_nodes.size(), 0 );
		for(int i = 0; i < m_links.size(); ++i)
		{
			for(int j = 0; j < 2; ++j) referencesPerNode[ m_links[i].m_linkIndicies[j] ]++;
		}
		for(int i = 0; i < m_faces.size(); ++i)
		{
			for(int j = 0; j < 3; ++j) referencesPerNode[ m_faces[i].m_indicies[j] ]++;
		}
		
		for(int i = 0; i < m_links.size(); ++i)
		{
			int linkIndex0 = m_links[i].m_linkIndicies[0];
			int linkIndex1 = m_links[i].m_linkIndicies[1];
			
			if(referencesPerNode[linkIndex0] == 1 || referencesPerNode[linkIndex1] == 1)
			{
				--referencesPerNode[linkIndex0];
				--referencesPerNode[linkIndex1];
				btSwap(m_links[i],m_links[m_links.size()-1]);
				m_links.pop_back();
				--i;
			}
		}
	}
	
	m_bUpdateRtCst=true;
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
			const btSoftBody::Node& faceNode0 = m_nodes[ f.m_indicies[0] ];
			const btSoftBody::Node& faceNode1 = m_nodes[ f.m_indicies[1] ];
			const btSoftBody::Node& faceNode2 = m_nodes[ f.m_indicies[2] ];
			
			const btVector3* points[3] = { &faceNode0.m_x, &faceNode1.m_x, &faceNode2.m_x} ;
			btDbvtVolume aabb = btDbvtVolume::FromPoints(points, 3);
			aabb.Expand( btVector3(radialMargin, radialMargin, radialMargin) );
			
			btVector3 averageVelocity = (faceNode0.m_v + faceNode1.m_v + faceNode2.m_v) / 3;
			
			m_faceBvh.update(f.m_leaf, aabb, averageVelocity * velocityMargin, updateMargin);
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
	const btAlignedObjectArray<btSoftBody::Node>& m_nodes;

	btVector3			m_rayFrom;
	btVector3			m_rayTo;
	btVector3			m_rayNormalizedDirection;
	btScalar			m_mint;
	btSoftBody::Face*	m_face;
	int					m_tests;
	
	
	RayFromToCaster(const btAlignedObjectArray<btSoftBody::Node>& nodes, const btVector3& rayFrom,const btVector3& rayTo,btScalar mxt);
	void Process(const btDbvtNode* leaf);

	static inline btScalar	rayFromToTriangle(const btVector3& rayFrom, const btVector3& rayTo,
		const btVector3& rayNormalizedDirection,
		const btVector3& a, const btVector3& b, const btVector3& c,
		btScalar maxt=SIMD_INFINITY);
};
RayFromToCaster::RayFromToCaster(const btAlignedObjectArray<btSoftBody::Node>& nodes, 
								const btVector3& rayFrom, const btVector3& rayTo, btScalar mxt) : m_nodes(nodes)
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
	btSoftBody::Face& f = *(btSoftBody::Face*)leaf->data;
	btVector3 node0 = m_nodes[ f.m_indicies[0] ].m_x;
	btVector3 node1 = m_nodes[ f.m_indicies[1] ].m_x;
	btVector3 node2 = m_nodes[ f.m_indicies[2] ].m_x;
	
	btScalar t = rayFromToTriangle(	m_rayFrom, m_rayTo, m_rayNormalizedDirection, node0, node1, node2, m_mint);
	if((t>0)&&(t<m_mint)) 
	{ 
		m_mint = t;
		m_face = &f; 
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
int btSoftBody::rayTest(const btVector3& rayFrom, const btVector3& rayTo, btScalar& hitFraction, int& faceIndex, bool bcountonly) const
{
	int numHits = 0;
	btVector3 dir = rayTo-rayFrom;
	

	if( bcountonly || m_faceBvh.empty() )
	{
		
		for(int i = 0; i < m_faces.size(); ++i)
		{
			const btSoftBody::Face&	f=m_faces[i];

			const btSoftBody::Node& faceNode0 = m_nodes[ f.m_indicies[0] ];
			const btSoftBody::Node& faceNode1 = m_nodes[ f.m_indicies[1] ];
			const btSoftBody::Node& faceNode2 = m_nodes[ f.m_indicies[2] ];

			btScalar t = RayFromToCaster::rayFromToTriangle(rayFrom, rayTo, dir, faceNode0.m_x, faceNode1.m_x, faceNode2.m_x, hitFraction);
			if(t > 0)
			{
				++numHits;
				if(!bcountonly)
				{
					faceIndex = i;
					hitFraction = t;
				}
			}
		}
	}
	else	//Use dbvt
	{
		RayFromToCaster	collider(m_nodes, rayFrom, rayTo, hitFraction);

		btDbvt::rayTest(m_faceBvh.m_root, rayFrom, rayTo, collider);
		if(collider.m_face)
		{
			hitFraction = collider.m_mint;
			faceIndex = (int)(collider.m_face-&m_faces[0]);
			numHits = 1;
		}
	}

	return numHits;
}

//
void btSoftBody::initializeFaceTree()
{
	m_faceBvh.clear();
	for(int i = 0; i < m_faces.size(); ++i)
	{
		Face& f = m_faces[i];
		const btSoftBody::Node& faceNode0 = m_nodes[ f.m_indicies[0] ];
		const btSoftBody::Node& faceNode1 = m_nodes[ f.m_indicies[1] ];
		const btSoftBody::Node& faceNode2 = m_nodes[ f.m_indicies[2] ];
		
		const btVector3* points[3] = { &faceNode0.m_x, &faceNode1.m_x, &faceNode2.m_x} ;
		btDbvtVolume faceAabb = btDbvtVolume::FromPoints(points, 3);
			
		f.m_leaf = m_faceBvh.insert( faceAabb, &f );
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
		btSoftBody::Face& face = m_faces[i];
		btSoftBody::Node& faceNode0 = m_nodes[ face.m_indicies[0] ];
		btSoftBody::Node& faceNode1 = m_nodes[ face.m_indicies[1] ];
		btSoftBody::Node& faceNode2 = m_nodes[ face.m_indicies[2] ];
		
		btVector3 normal = btCross(faceNode1.m_x - faceNode0.m_x, faceNode2.m_x - faceNode0.m_x).normalized();
		face.m_normal = normal;
		faceNode0.m_normal += normal;
		faceNode1.m_normal += normal;
		faceNode2.m_normal += normal;
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


void btSoftBody::updateArea()
{
	// Face area
	for(int i = 0; i < m_faces.size(); ++i)
	{
		Face& f = m_faces[i];
		const btVector3& position0 = m_nodes[ f.m_indicies[0] ].m_x;
		const btVector3& position1 = m_nodes[ f.m_indicies[1] ].m_x;
		const btVector3& position2 = m_nodes[ f.m_indicies[2] ].m_x;
		
		f.m_area = AreaOfParallelogram(position0, position1, position2) * btScalar(0.5f);
	}
	
	// Node area
	{
		for(int i = 0; i < m_nodes.size(); ++i) m_nodes[i].m_area = btScalar(0.0);

		for(int i = 0; i < m_faces.size(); ++i)
		{
			btSoftBody::Face& f = m_faces[i];
			
			btScalar area = btFabs(f.m_area);
			m_nodes[ f.m_indicies[0] ].m_area += area;
			m_nodes[ f.m_indicies[1] ].m_area += area;
			m_nodes[ f.m_indicies[2] ].m_area += area;
		}

		for(int i = 0; i < m_nodes.size(); ++i) m_nodes[i].m_area *= btScalar(0.3333333);
	}
}


void btSoftBody::updateConstants()
{
	//Set current link lengths as resting lengths
	for(int i = 0; i < m_links.size(); ++i)
	{
		Link& l = m_links[i];
		
		const btVector3& nodePosition0 = m_nodes[ l.m_linkIndicies[0] ].m_x;
		const btVector3& nodePosition1 = m_nodes[ l.m_linkIndicies[1] ].m_x;
		
		l.m_restLength = (nodePosition0 - nodePosition1).length();
		l.m_restLengthSquared = l.m_restLength * l.m_restLength;
	}
	
	//Update link constants
	for(int i = 0; i < m_links.size(); ++i)
	{
		Link& l= m_links[i];
		Material& m = *l.m_material;
		
		btScalar invMass0 = m_nodes[ l.m_linkIndicies[0] ].m_invMass;
		btScalar invMass1 = m_nodes[ l.m_linkIndicies[1] ].m_invMass;
		l.m_scaledCombinedInvMass = (invMass0 + invMass1) / m.m_linearStiffness;
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
		btSoftColliders::SoftSoftVertexFaceCollider callback;
		
		// common
		callback.mrg =	getCollisionShape()->getMargin()+ psb->getCollisionShape()->getMargin();
		
		//
		callback.m_nodeSoftBody = this;
		callback.m_faceSoftBody = psb;
		callback.m_nodeSoftBody->m_nodeBvh.collideTT(callback.m_nodeSoftBody->m_nodeBvh.m_root, callback.m_faceSoftBody->m_faceBvh.m_root, callback);
		
		//
		callback.m_nodeSoftBody = psb;
		callback.m_faceSoftBody = this;
		callback.m_nodeSoftBody->m_nodeBvh.collideTT(callback.m_nodeSoftBody->m_nodeBvh.m_root, callback.m_faceSoftBody->m_faceBvh.m_root, callback);
	}
}

/************************************************************************************
///Aero force
************************************************************************************/
static inline void ApplyClampedForce(btSoftBody::Node* n, const btVector3& f, btScalar dt)
{
	btScalar dtim = dt * n->m_invMass;
	if( (f * dtim).length2() > n->m_v.length2() ) n->m_f -= ProjectOnAxis( n->m_v, f.normalized() ) / dtim;	//Apply clamped force
	else n->m_f += f;	//Unclamped force 
}

void btSoftBody::addAeroForces(const AeroForce& aeroForce, btScalar timeStep, btAlignedObjectArray<Node>& nodes, btAlignedObjectArray<Face>& faces)
{
	const bool use_medium =	(aeroForce.m_liftCoeff > 0) || (aeroForce.m_dragCoeff > 0);
	if(use_medium)
	{
		for(int i = 0; i < nodes.size(); ++i)
			if(nodes[i].m_invMass > 0) addAeroForceToNode(aeroForce, timeStep, nodes, i);
		
		for(int i = 0; i < faces.size(); ++i) addAeroForceToFace(aeroForce, timeStep, nodes, faces, i);	 
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
					btScalar tri_area = n.m_area;
							
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
						ApplyClampedForce(&n, force, timeStep);
					}
				}	
			}
		}
	}
}

void btSoftBody::addAeroForceToFace(const AeroForce& aeroForce, btScalar timeStep, btAlignedObjectArray<Node>& nodes, btAlignedObjectArray<Face>& faces, int faceIndex)
{
	const btScalar kLF = aeroForce.m_liftCoeff;
	const btScalar kDG = aeroForce.m_dragCoeff;
	
	const bool as_lift = kLF>0;
	const bool as_drag = kDG>0;
	const bool as_aero = as_lift || as_drag;
	const bool as_faero = as_aero && (aeroForce.m_model >= btSoftBody::eAeroModel::F_TwoSided);

	if(as_faero)
	{
		btSoftBody::Face& f = faces[faceIndex];
		btSoftBody::Node* faceNodes[3] = { &nodes[ f.m_indicies[0] ], &nodes[ f.m_indicies[1] ], &nodes[ f.m_indicies[2] ] };

		const btVector3	v=(faceNodes[0]->m_v+faceNodes[1]->m_v+faceNodes[2]->m_v)/3;
		
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
				btScalar tri_area = f.m_area;
					
				fDrag = 0.5f * kDG * aeroForce.m_airDensity * rel_v2 * tri_area * n_dot_v * (-rel_v_nrm);

				// Check angle of attack
				// cos(10º) = 0.98480
				if ( 0 < n_dot_v && n_dot_v < 0.98480f)
					fLift = 0.5f * kLF * aeroForce.m_airDensity * rel_v_len * tri_area * btSqrt(1.0f-n_dot_v*n_dot_v) * (nrm.cross(rel_v_nrm).cross(rel_v_nrm));

				fDrag /= 3;
				fLift /= 3;

				for(int j=0;j<3;++j) 
				{
					if (faceNodes[j]->m_invMass > 0)
					{
						// Check if the velocity change resulted by aero drag force exceeds the current velocity of the node.
						btVector3 del_v_by_fDrag = fDrag * faceNodes[j]->m_invMass * timeStep;
						btScalar del_v_by_fDrag_len2 = del_v_by_fDrag.length2();
						btScalar v_len2 = faceNodes[j]->m_v.length2();

						if (del_v_by_fDrag_len2 >= v_len2 && del_v_by_fDrag_len2 > 0)
						{
							btScalar del_v_by_fDrag_len = del_v_by_fDrag.length();
							btScalar v_len = faceNodes[j]->m_v.length();
							fDrag *= btScalar(0.8)*(v_len / del_v_by_fDrag_len);
						}

						faceNodes[j]->m_f += fDrag; 
						faceNodes[j]->m_f += fLift;
					}
				}
			}
			else if (aeroForce.m_model == btSoftBody::eAeroModel::F_OneSided || aeroForce.m_model == btSoftBody::eAeroModel::F_TwoSided)
			{
				if (btSoftBody::eAeroModel::F_TwoSided)
					nrm *= (btScalar)( (btDot(nrm,rel_v) < 0) ? -1 : +1);

				const btScalar	dvn=btDot(rel_v,nrm);
				// Compute forces	 
				if(dvn>0)
				{
					btVector3		force(0,0,0);
					const btScalar	c0	=	f.m_area * dvn * rel_v2;
					const btScalar	c1	=	c0*aeroForce.m_airDensity;
					force	+=	nrm*(-c1*kLF);
					force	+=	rel_v.normalized()*(-c1*kDG);
					force	/=	3;
					for(int j = 0; j < 3; ++j) ApplyClampedForce(faceNodes[j], force, timeStep);
				}
			}
		}
	}

}
/************************************************************************************
///Aero
************************************************************************************/