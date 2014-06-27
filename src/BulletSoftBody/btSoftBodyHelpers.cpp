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
///btSoftBodyHelpers.cpp by Nathanael Presson

#include "btSoftBodyInternals.h"
#include <stdio.h>
#include <string.h>
#include "btSoftBodyHelpers.h"
#include "LinearMath/btConvexHull.h"
#include "LinearMath/btConvexHullComputer.h"


//
static void				drawVertex(	btIDebugDraw* idraw,
								   const btVector3& x,btScalar s,const btVector3& c)
{
	idraw->drawLine(x-btVector3(s,0,0),x+btVector3(s,0,0),c);
	idraw->drawLine(x-btVector3(0,s,0),x+btVector3(0,s,0),c);
	idraw->drawLine(x-btVector3(0,0,s),x+btVector3(0,0,s),c);
}

//
static void				drawBox(	btIDebugDraw* idraw,
								const btVector3& mins,
								const btVector3& maxs,
								const btVector3& color)
{
	const btVector3	c[]={	btVector3(mins.x(),mins.y(),mins.z()),
		btVector3(maxs.x(),mins.y(),mins.z()),
		btVector3(maxs.x(),maxs.y(),mins.z()),
		btVector3(mins.x(),maxs.y(),mins.z()),
		btVector3(mins.x(),mins.y(),maxs.z()),
		btVector3(maxs.x(),mins.y(),maxs.z()),
		btVector3(maxs.x(),maxs.y(),maxs.z()),
		btVector3(mins.x(),maxs.y(),maxs.z())};
	idraw->drawLine(c[0],c[1],color);idraw->drawLine(c[1],c[2],color);
	idraw->drawLine(c[2],c[3],color);idraw->drawLine(c[3],c[0],color);
	idraw->drawLine(c[4],c[5],color);idraw->drawLine(c[5],c[6],color);
	idraw->drawLine(c[6],c[7],color);idraw->drawLine(c[7],c[4],color);
	idraw->drawLine(c[0],c[4],color);idraw->drawLine(c[1],c[5],color);
	idraw->drawLine(c[2],c[6],color);idraw->drawLine(c[3],c[7],color);
}

//
static void				drawTree(	btIDebugDraw* idraw,
								 const btDbvtNode* node,
								 int depth,
								 const btVector3& ncolor,
								 const btVector3& lcolor,
								 int mindepth,
								 int maxdepth)
{
	if(node)
	{
		if(node->isinternal()&&((depth<maxdepth)||(maxdepth<0)))
		{
			drawTree(idraw,node->childs[0],depth+1,ncolor,lcolor,mindepth,maxdepth);
			drawTree(idraw,node->childs[1],depth+1,ncolor,lcolor,mindepth,maxdepth);
		}
		if(depth>=mindepth)
		{
			const btScalar	scl=(btScalar)(node->isinternal()?1:1);
			const btVector3	mi=node->volume.Center()-node->volume.Extents()*scl;
			const btVector3	mx=node->volume.Center()+node->volume.Extents()*scl;
			drawBox(idraw,mi,mx,node->isleaf()?lcolor:ncolor);
		}
	}
}

//
template <typename T>
static inline T				sum(const btAlignedObjectArray<T>& items)
{
	T	v;
	if(items.size())
	{
		v=items[0];
		for(int i=1,ni=items.size();i<ni;++i)
		{
			v+=items[i];
		}
	}
	return(v);
}

//
template <typename T,typename Q>
static inline void			add(btAlignedObjectArray<T>& items,const Q& value)
{
	for(int i=0,ni=items.size();i<ni;++i)
	{
		items[i]+=value;
	}
}

//
template <typename T,typename Q>
static inline void			mul(btAlignedObjectArray<T>& items,const Q& value)
{
	for(int i=0,ni=items.size();i<ni;++i)
	{
		items[i]*=value;
	}
}

//
template <typename T>
static inline T				average(const btAlignedObjectArray<T>& items)
{
	const btScalar	n=(btScalar)(items.size()>0?items.size():1);
	return(sum(items)/n);
}

//
static inline btScalar		tetravolume(const btVector3& x0,
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
#if 0
static btVector3		stresscolor(btScalar stress)
{
	static const btVector3	spectrum[]=	{	btVector3(1,0,1),
		btVector3(0,0,1),
		btVector3(0,1,1),
		btVector3(0,1,0),
		btVector3(1,1,0),
		btVector3(1,0,0),
		btVector3(1,0,0)};
	static const int		ncolors=sizeof(spectrum)/sizeof(spectrum[0])-1;
	static const btScalar	one=1;
	stress=btMax<btScalar>(0,btMin<btScalar>(1,stress))*ncolors;
	const int				sel=(int)stress;
	const btScalar			frc=stress-sel;
	return(spectrum[sel]+(spectrum[sel+1]-spectrum[sel])*frc);
}
#endif

//
void			btSoftBodyHelpers::Draw(btSoftBody* psb, btIDebugDraw* idraw, int drawflags)
{
	const btScalar		scl=(btScalar)0.1;
	const btScalar		nscl=scl*5;
	const btVector3		lcolor=btVector3(0,0,0);
	const btVector3		ncolor=btVector3(1,1,1);
	const btVector3		ccolor=btVector3(1,0,0);
	int i,j,nj;

	{
		// Nodes	 
		if(0!=(drawflags&fDrawFlags::Nodes))
		{
			for(i=0;i<psb->m_nodes.size();++i)
			{
				const btSoftBodyNode& n = psb->m_nodes[i];
				idraw->drawLine(n.m_position-btVector3(scl,0,0),n.m_position+btVector3(scl,0,0),btVector3(1,0,0));
				idraw->drawLine(n.m_position-btVector3(0,scl,0),n.m_position+btVector3(0,scl,0),btVector3(0,1,0));
				idraw->drawLine(n.m_position-btVector3(0,0,scl),n.m_position+btVector3(0,0,scl),btVector3(0,0,1));
			}
		}
		// Links	 
		if(0!=(drawflags&fDrawFlags::Links))
		{
			for(i=0;i<psb->m_links.size();++i)
			{
				const btSoftBodyLink& l = psb->m_links[i];
				if(!l.m_material->m_debugDraw) continue;
				idraw->drawLine( psb->m_nodes[ l.m_linkIndicies[0] ].m_position, psb->m_nodes[ l.m_linkIndicies[1] ].m_position, lcolor);
			}
		}
		// Normals	 
		if(0!=(drawflags&fDrawFlags::Normals))
		{
			for(i=0;i<psb->m_nodes.size();++i)
			{
				const btSoftBodyNode& n = psb->m_nodes[i];
				const btVector3 d = n.m_normal * nscl;
				idraw->drawLine(n.m_position,n.m_position+d,ncolor);
				idraw->drawLine(n.m_position,n.m_position-d,ncolor*0.5);
			}
		}
		// Contacts	 
		if(0!=(drawflags&fDrawFlags::Contacts))
		{
			static const btVector3		axis[]={btVector3(1,0,0),
				btVector3(0,1,0),
				btVector3(0,0,1)};
			for(i=0;i<psb->m_rigidContacts.size();++i)
			{		
				const btSoftRigidContact& c = psb->m_rigidContacts[i];
				const btVector3& nodePosition = psb->m_nodes[c.m_nodeIndex].m_position;
				
				btVector3 o = nodePosition - c.m_normal * ( btDot(nodePosition, c.m_normal) + c.m_offset );
				btVector3 x = btCross(c.m_normal,axis[c.m_normal.minAxis()]).normalized();
				btVector3 y = btCross(x,c.m_normal).normalized();
				idraw->drawLine(o-x*nscl,o+x*nscl,ccolor);
				idraw->drawLine(o-y*nscl,o+y*nscl,ccolor);
				idraw->drawLine(o,o+c.m_normal*nscl*3,btVector3(1,1,0));
			}
		}
		// Faces	 
		if(0!=(drawflags&fDrawFlags::Faces))
		{
			const btAlignedObjectArray<btSoftBodyNode>& nodes = psb->m_nodes;
		
			const btScalar	scl=(btScalar)0.8;
			const btScalar	alp=(btScalar)1;
			const btVector3	col(0,(btScalar)0.7,0);
			for(i=0;i<psb->m_faces.size();++i)
			{
				const btSoftBodyFace&	face = psb->m_faces[i];
				if(!face.m_material->m_debugDraw) continue;
				
				const btVector3 x[] = { nodes[ face.m_indicies[0] ].m_position, nodes[ face.m_indicies[1] ].m_position, nodes[ face.m_indicies[2] ].m_position };
				const btVector3 c = (x[0] + x[1] + x[2]) / 3;
				idraw->drawTriangle((x[0]-c)*scl+c, (x[1]-c)*scl+c, (x[2]-c)*scl+c, col,alp);
			}	
		}
	}
	// Anchors	 
	if(0!=(drawflags&fDrawFlags::Anchors))
	{
		for(i=0;i<psb->m_anchors.size();++i)
		{
			const btSoftBodyAnchor& a = psb->m_anchors[i];
			const btSoftBodyNode& node = psb->m_nodes[a.m_nodeIndex];
			
			const btVector3 q = a.m_body->getWorldTransform() * a.m_local;
			drawVertex( idraw, node.m_position, 0.25, btVector3(1,0,0) );
			drawVertex( idraw, q, 0.25, btVector3(0,1,0) );
			idraw->drawLine( node.m_position, q, btVector3(1,1,1) );
		}
		for(i=0;i<psb->m_nodes.size();++i)
		{
			const btSoftBodyNode&	n=psb->m_nodes[i];
			if(n.m_invMass <= 0)
			{
				drawVertex(idraw,n.m_position,0.25,btVector3(1,0,0));
			}
		}
	}
	
	// Node tree	 
	if(0!=(drawflags&fDrawFlags::NodeTree))		DrawNodeTree(psb,idraw);
	// Face tree	 
	if(0!=(drawflags&fDrawFlags::FaceTree))		DrawFaceTree(psb,idraw);
}

//
void			btSoftBodyHelpers::DrawInfos(		btSoftBody* psb,
											 btIDebugDraw* idraw,
											 bool masses,
											 bool areas,
											 bool)
{
	for(int i=0;i<psb->m_nodes.size();++i)
	{
		const btSoftBodyNode&	n=psb->m_nodes[i];
		char					text[2048]={0};
		char					buff[1024];
		if(masses)
		{
			sprintf(buff," M(%.2f)",1/n.m_invMass);
			strcat(text,buff);
		}
		if(areas)
		{
			sprintf(buff," A(%.2f)",n.m_area);
			strcat(text,buff);
		}
		if(text[0]) idraw->draw3dText(n.m_position,text);
	}
}

//
void			btSoftBodyHelpers::DrawNodeTree(	btSoftBody* psb,
												btIDebugDraw* idraw,
												int mindepth,
												int maxdepth)
{
	drawTree(idraw, psb->m_nodeBvh.m_root, 0, btVector3(1,0,1), btVector3(1,1,1), mindepth, maxdepth);
}

//
void			btSoftBodyHelpers::DrawFaceTree(	btSoftBody* psb,
												btIDebugDraw* idraw,
												int mindepth,
												int maxdepth)
{
	drawTree(idraw,psb->m_faceBvh.m_root,0,btVector3(0,1,0),btVector3(1,0,0),mindepth,maxdepth);
}

//The btSoftBody object from the BulletSDK includes an array of Nodes and Links. These links appear
// to be first set up to connect a node to between 5 and 6 of its neighbors [480 links], 
//and then to the rest of the nodes after the execution of the Floyd-Warshall graph algorithm 
//[another 930 links]. 
//The way the links are stored by default, we have a number of cases where adjacent links share a node in common
// - this leads to the creation of a data dependency through memory. 
//The PSolve_Links() function reads and writes nodes as it iterates over each link. 
//So, we now have the possibility of a data dependency between iteration X 
//that processes link L with iteration X+1 that processes link L+1 
//because L and L+1 have one node in common, and iteration X updates the positions of that node, 
//and iteration X+1 reads in the position of that shared node.
//
//Such a memory dependency limits the ability of a modern CPU to speculate beyond 
//a certain point because it has to respect a possible dependency 
//- this prevents the CPU from making full use of its out-of-order resources. 
//If we re-order the links such that we minimize the cases where a link L and L+1 share a common node, 
//we create a temporal gap between when the node position is written, 
//and when it is subsequently read. This in turn allows the CPU to continue execution without 
//risking a dependency violation. Such a reordering would result in significant speedups on 
//modern CPUs with lots of execution resources. 
//In our testing, we see it have a tremendous impact not only on the A7, 
//but also on all x86 cores that ship with modern Macs. 
//The attached source file includes a single function (ReoptimizeLinkOrder) which can be called on a 
//btSoftBody object in the solveConstraints() function before the actual solver is invoked, 
//or right after generateBendingConstraints() once we have all 1410 links.


//===================================================================
//
//
// This function takes in a list of interdependent Links and tries 
// to maximize the distance between calculation
// of dependent links.  This increases the amount of parallelism that can
// be exploited by out-of-order instruction processors with large but
// (inevitably) finite instruction windows.
//
//===================================================================

// A small structure to track lists of dependent link calculations
class LinkDeps_t {
	public:
	int value;			// A link calculation that is dependent on this one
		// Positive values = "input A" while negative values = "input B"
	LinkDeps_t *next;	// Next dependence in the list
};
typedef LinkDeps_t *LinkDepsPtr_t;

// Dependency list constants
#define REOP_NOT_DEPENDENT	-1
#define REOP_NODE_COMPLETE	-2	// Must be less than REOP_NOT_DEPENDENT


void btSoftBodyHelpers::ReoptimizeLinkOrder(btSoftBody *psb)
{
	int i, nLinks=psb->m_links.size(), nNodes=psb->m_nodes.size();
	btSoftBodyLink *lr;
	int ar, br;
	btSoftBodyNode *node0 = &(psb->m_nodes[0]);
	btSoftBodyNode *node1 = &(psb->m_nodes[1]);
	LinkDepsPtr_t linkDep;
	int readyListHead, readyListTail, linkNum, linkDepFrees, depLink;
	
	// Allocate temporary buffers
	int *nodeWrittenAt = new int[nNodes+1];	// What link calculation produced this node's current values?
	int *linkDepA = new int[nLinks];			// Link calculation input is dependent upon prior calculation #N
	int *linkDepB = new int[nLinks];
	int *readyList = new int[nLinks];		// List of ready-to-process link calculations (# of links, maximum)
	LinkDeps_t *linkDepFreeList = new LinkDeps_t[2*nLinks];		// Dependent-on-me list elements (2x# of links, maximum)
	LinkDepsPtr_t *linkDepListStarts = new LinkDepsPtr_t[nLinks];	// Start nodes of dependent-on-me lists, one for each link
		
	// Copy the original, unsorted links to a side buffer
	btSoftBodyLink *linkBuffer = new btSoftBodyLink[nLinks];
	memcpy(linkBuffer, &(psb->m_links[0]), sizeof(btSoftBodyLink)*nLinks);

	// Clear out the node setup and ready list
	for (i=0; i < nNodes+1; i++) {
		nodeWrittenAt[i] = REOP_NOT_DEPENDENT;
	}
	for (i=0; i < nLinks; i++) {
		linkDepListStarts[i] = NULL;
	}
	readyListHead = readyListTail = linkDepFrees = 0;

	// Initial link analysis to set up data structures
	for (i=0; i < nLinks; i++) {
	
		// Note which prior link calculations we are dependent upon & build up dependence lists
		lr = &(psb->m_links[i]);
		ar = ( (&psb->m_nodes[ lr->m_linkIndicies[0] ]) - node0 ) / (node1 - node0);
		br = ( (&psb->m_nodes[ lr->m_linkIndicies[1] ]) - node0 ) / (node1 - node0);
		if (nodeWrittenAt[ar] > REOP_NOT_DEPENDENT) {
			linkDepA[i] = nodeWrittenAt[ar];
			linkDep = &linkDepFreeList[linkDepFrees++];
			linkDep->value = i;
			linkDep->next = linkDepListStarts[nodeWrittenAt[ar]];
			linkDepListStarts[nodeWrittenAt[ar]] = linkDep;
		} else {
			linkDepA[i] = REOP_NOT_DEPENDENT;
		}
		if (nodeWrittenAt[br] > REOP_NOT_DEPENDENT) {
			linkDepB[i] = nodeWrittenAt[br];
			linkDep = &linkDepFreeList[linkDepFrees++];
			linkDep->value = -(i+1);
			linkDep->next = linkDepListStarts[nodeWrittenAt[br]];
			linkDepListStarts[nodeWrittenAt[br]] = linkDep;
		} else {
			linkDepB[i] = REOP_NOT_DEPENDENT;
		}
		
		// Add this link to the initial ready list, if it is not dependent on any other links
		if ((linkDepA[i] == REOP_NOT_DEPENDENT) && (linkDepB[i] == REOP_NOT_DEPENDENT)) {
			readyList[readyListTail++] = i;
			linkDepA[i] = linkDepB[i] = REOP_NODE_COMPLETE;	// Probably not needed now
		}
		
		// Update the nodes to mark which ones are calculated by this link
		nodeWrittenAt[ar] = nodeWrittenAt[br] = i;
	}
	
	// Process the ready list and create the sorted list of links
	// -- By treating the ready list as a queue, we maximize the distance between any
	//    inter-dependent node calculations
	// -- All other (non-related) nodes in the ready list will automatically be inserted
	//    in between each set of inter-dependent link calculations by this loop
	i = 0;
	while (readyListHead != readyListTail) {
		// Use ready list to select the next link to process
		linkNum = readyList[readyListHead++];
		// Copy the next-to-calculate link back into the original link array
		psb->m_links[i++] = linkBuffer[linkNum];

		// Free up any link inputs that are dependent on this one
		linkDep = linkDepListStarts[linkNum];
		while (linkDep) {
			depLink = linkDep->value;
			if (depLink >= 0) {
				linkDepA[depLink] = REOP_NOT_DEPENDENT;
			} else {
				depLink = -depLink - 1;
				linkDepB[depLink] = REOP_NOT_DEPENDENT;
			}
			// Add this dependent link calculation to the ready list if *both* inputs are clear
			if ((linkDepA[depLink] == REOP_NOT_DEPENDENT) && (linkDepB[depLink] == REOP_NOT_DEPENDENT)) {
				readyList[readyListTail++] = depLink;
				linkDepA[depLink] = linkDepB[depLink] = REOP_NODE_COMPLETE;	// Probably not needed now
			}
			linkDep = linkDep->next;
		}
	}

	// Delete the temporary buffers
	delete [] nodeWrittenAt;
	delete [] linkDepA;
	delete [] linkDepB;
	delete [] readyList;
	delete [] linkDepFreeList;
	delete [] linkDepListStarts;
	delete [] linkBuffer;
}


//
void			btSoftBodyHelpers::DrawFrame(		btSoftBody* psb,
											 btIDebugDraw* idraw)
{
	if( psb->m_pose.m_poseMatching != btScalar(0.0) )
	{
		static const btScalar	axisScaling = 10;
		static const btScalar	nodeScaling = (btScalar)0.1;
		const btVector3			centerOfMass = psb->m_pose.m_centerOfMass;
		const btMatrix3x3		rotation = psb->m_pose.m_rotation;
		const btVector3			Xaxis = ( rotation * btVector3(1,0,0) ).normalized();
		const btVector3			Yaxis = ( rotation * btVector3(0,1,0) ).normalized();
		const btVector3			Zaxis = ( rotation * btVector3(0,0,1) ).normalized();
		idraw->drawLine( centerOfMass, centerOfMass + Xaxis * axisScaling, btVector3(1,0,0) );
		idraw->drawLine( centerOfMass, centerOfMass + Yaxis * axisScaling, btVector3(0,1,0) );
		idraw->drawLine( centerOfMass, centerOfMass + Zaxis * axisScaling, btVector3(0,0,1) );
		for(int i=0;i<psb->m_pose.m_referencePositions.size();++i)
		{
			const btVector3	x = centerOfMass + rotation * psb->m_pose.m_referencePositions[i];
			drawVertex( idraw, x, nodeScaling, btVector3(1,0,1) );
		}
	}
}

//
btSoftBody*		btSoftBodyHelpers::CreateRope(	btSoftBodyWorldInfo& worldInfo, const btVector3& from,
											  const btVector3& to,
											  int res,
											  int fixeds)
{
	// Create nodes	 
	const int		r=res+2;
	btVector3*		x=new btVector3[r];
	btScalar*		m=new btScalar[r];
	int i;

	for(i=0;i<r;++i)
	{
		const btScalar	t=i/(btScalar)(r-1);
		x[i]=lerp(from,to,t);
		m[i]=1;
	}
	btSoftBody*		psb= new btSoftBody(&worldInfo,r,x,m);
	if(fixeds&1) psb->setMass(0,0);
	if(fixeds&2) psb->setMass(r-1,0);
	delete[] x;
	delete[] m;
	// Create links	 
	for(i=1;i<r;++i)
	{
		psb->appendLink(i-1,i);
	}
	// Finished		 
	return(psb);
}

//
btSoftBody*		btSoftBodyHelpers::CreatePatch(btSoftBodyWorldInfo& worldInfo,const btVector3& corner00,
											   const btVector3& corner10,
											   const btVector3& corner01,
											   const btVector3& corner11,
											   int resx,
											   int resy,
											   int fixeds,
											   bool gendiags)
{
#define IDX(_x_,_y_)	((_y_)*rx+(_x_))
	// Create nodes	 
	if((resx<2)||(resy<2)) return(0);
	const int	rx=resx;
	const int	ry=resy;
	const int	tot=rx*ry;
	btVector3*	x=new btVector3[tot];
	btScalar*	m=new btScalar[tot];
	int iy;

	for(iy=0;iy<ry;++iy)
	{
		const btScalar	ty=iy/(btScalar)(ry-1);
		const btVector3	py0=lerp(corner00,corner01,ty);
		const btVector3	py1=lerp(corner10,corner11,ty);
		for(int ix=0;ix<rx;++ix)
		{
			const btScalar	tx=ix/(btScalar)(rx-1);
			x[IDX(ix,iy)]=lerp(py0,py1,tx);
			m[IDX(ix,iy)]=1;
		}
	}
	btSoftBody*		psb=new btSoftBody(&worldInfo,tot,x,m);
	if(fixeds&1)	psb->setMass(IDX(0,0),0);
	if(fixeds&2)	psb->setMass(IDX(rx-1,0),0);
	if(fixeds&4)	psb->setMass(IDX(0,ry-1),0);
	if(fixeds&8)	psb->setMass(IDX(rx-1,ry-1),0);
	delete[] x;
	delete[] m;
	// Create links	and faces  
	for(iy=0;iy<ry;++iy)
	{
		for(int ix=0;ix<rx;++ix)
		{
			const int	idx=IDX(ix,iy);
			const bool	mdx=(ix+1)<rx;
			const bool	mdy=(iy+1)<ry;
			if(mdx) psb->appendLink(idx,IDX(ix+1,iy));
			if(mdy) psb->appendLink(idx,IDX(ix,iy+1));
			if(mdx&&mdy)
			{
				if((ix+iy)&1)
				{
					psb->appendFace(IDX(ix,iy),IDX(ix+1,iy),IDX(ix+1,iy+1));
					psb->appendFace(IDX(ix,iy),IDX(ix+1,iy+1),IDX(ix,iy+1));
					if(gendiags)
					{
						psb->appendLink(IDX(ix,iy),IDX(ix+1,iy+1));
					}
				}
				else
				{
					psb->appendFace(IDX(ix,iy+1),IDX(ix,iy),IDX(ix+1,iy));
					psb->appendFace(IDX(ix,iy+1),IDX(ix+1,iy),IDX(ix+1,iy+1));
					if(gendiags)
					{
						psb->appendLink(IDX(ix+1,iy),IDX(ix,iy+1));
					}
				}
			}
		}
	}
	// Finished		 
#undef IDX
	return(psb);
}

//
btSoftBody*		btSoftBodyHelpers::CreatePatchUV(btSoftBodyWorldInfo& worldInfo,
												 const btVector3& corner00,
												 const btVector3& corner10,
												 const btVector3& corner01,
												 const btVector3& corner11,
												 int resx,
												 int resy,
												 int fixeds,
												 bool gendiags,
												 float* tex_coords)
{

	/*
	*
	*  corners:
	*
	*  [0][0]     corner00 ------- corner01   [resx][0]
	*                |                |
	*                |                |
	*  [0][resy]  corner10 -------- corner11  [resx][resy]
	*
	*
	*
	*
	*
	*
	*   "fixedgs" map:
	*
	*  corner00     -->   +1
	*  corner01     -->   +2
	*  corner10     -->   +4
	*  corner11     -->   +8
	*  upper middle -->  +16
	*  left middle  -->  +32
	*  right middle -->  +64
	*  lower middle --> +128
	*  center       --> +256
	*
	*
	*   tex_coords size   (resx-1)*(resy-1)*12
	*
	*
	*
	*     SINGLE QUAD INTERNALS
	*
	*  1) btSoftBody's nodes and links,
	*     diagonal link is optional ("gendiags")
	*
	*
	*    node00 ------ node01
	*      | .              
	*      |   .            
	*      |     .          
	*      |       .        
	*      |         .      
	*    node10        node11
	*
	*
	*
	*   2) Faces:
	*      two triangles,
	*      UV Coordinates (hier example for single quad)
	*      
	*     (0,1)          (0,1)  (1,1)
	*     1 |\            3 \-----| 2
	*       | \              \    |
	*       |  \              \   |
	*       |   \              \  |
	*       |    \              \ |
	*     2 |-----\ 3            \| 1
	*     (0,0)    (1,0)       (1,0)
	*
	*
	*
	*
	*
	*
	*/

#define IDX(_x_,_y_)	((_y_)*rx+(_x_))
	// Create nodes		 
	if((resx<2)||(resy<2)) return(0);
	const int	rx=resx;
	const int	ry=resy;
	const int	tot=rx*ry;
	btVector3*	x=new btVector3[tot];
	btScalar*	m=new btScalar[tot];

	int iy;

	for(iy=0;iy<ry;++iy)
	{
		const btScalar	ty=iy/(btScalar)(ry-1);
		const btVector3	py0=lerp(corner00,corner01,ty);
		const btVector3	py1=lerp(corner10,corner11,ty);
		for(int ix=0;ix<rx;++ix)
		{
			const btScalar	tx=ix/(btScalar)(rx-1);
			x[IDX(ix,iy)]=lerp(py0,py1,tx);
			m[IDX(ix,iy)]=1;
		}
	}
	btSoftBody*	psb=new btSoftBody(&worldInfo,tot,x,m);
	if(fixeds&1)		psb->setMass(IDX(0,0),0);
	if(fixeds&2)		psb->setMass(IDX(rx-1,0),0);
	if(fixeds&4)		psb->setMass(IDX(0,ry-1),0);
	if(fixeds&8)		psb->setMass(IDX(rx-1,ry-1),0);
	if(fixeds&16)		psb->setMass(IDX((rx-1)/2,0),0);
	if(fixeds&32)		psb->setMass(IDX(0,(ry-1)/2),0);
	if(fixeds&64)		psb->setMass(IDX(rx-1,(ry-1)/2),0);
	if(fixeds&128)		psb->setMass(IDX((rx-1)/2,ry-1),0);
	if(fixeds&256)		psb->setMass(IDX((rx-1)/2,(ry-1)/2),0);
	delete[] x;
	delete[] m;


	int z = 0;
	// Create links	and faces	 
	for(iy=0;iy<ry;++iy)
	{
		for(int ix=0;ix<rx;++ix)
		{
			const bool	mdx=(ix+1)<rx;
			const bool	mdy=(iy+1)<ry;

			int node00=IDX(ix,iy);
			int node01=IDX(ix+1,iy);
			int node10=IDX(ix,iy+1);
			int node11=IDX(ix+1,iy+1);

			if(mdx) psb->appendLink(node00,node01);
			if(mdy) psb->appendLink(node00,node10);
			if(mdx&&mdy)
			{
				psb->appendFace(node00,node10,node11);
				if (tex_coords) {
					tex_coords[z+0]=CalculateUV(resx,resy,ix,iy,0);
					tex_coords[z+1]=CalculateUV(resx,resy,ix,iy,1);
					tex_coords[z+2]=CalculateUV(resx,resy,ix,iy,0);
					tex_coords[z+3]=CalculateUV(resx,resy,ix,iy,2);
					tex_coords[z+4]=CalculateUV(resx,resy,ix,iy,3);
					tex_coords[z+5]=CalculateUV(resx,resy,ix,iy,2);
				}
				psb->appendFace(node11,node01,node00);
				if (tex_coords) {
					tex_coords[z+6 ]=CalculateUV(resx,resy,ix,iy,3);
					tex_coords[z+7 ]=CalculateUV(resx,resy,ix,iy,2);
					tex_coords[z+8 ]=CalculateUV(resx,resy,ix,iy,3);
					tex_coords[z+9 ]=CalculateUV(resx,resy,ix,iy,1);
					tex_coords[z+10]=CalculateUV(resx,resy,ix,iy,0);
					tex_coords[z+11]=CalculateUV(resx,resy,ix,iy,1);
				}
				if (gendiags) psb->appendLink(node00,node11);
				z += 12;
			}
		}
	}
	// Finished	 
#undef IDX
	return(psb);
}

float   btSoftBodyHelpers::CalculateUV(int resx,int resy,int ix,int iy,int id)
{

	/*
	*
	*
	*    node00 --- node01
	*      |          |
	*    node10 --- node11
	*
	*
	*   ID map:
	*
	*   node00 s --> 0
	*   node00 t --> 1
	*
	*   node01 s --> 3
	*   node01 t --> 1
	*
	*   node10 s --> 0
	*   node10 t --> 2
	*
	*   node11 s --> 3
	*   node11 t --> 2
	*
	*
	*/

	float tc=0.0f;
	if (id == 0) {
		tc = (1.0f/((resx-1))*ix);
	}
	else if (id==1) {
		tc = (1.0f/((resy-1))*(resy-1-iy));
	}
	else if (id==2) {
		tc = (1.0f/((resy-1))*(resy-1-iy-1));
	}
	else if (id==3) {
		tc = (1.0f/((resx-1))*(ix+1));
	}
	return tc;
}
//
btSoftBody*		btSoftBodyHelpers::CreateEllipsoid(btSoftBodyWorldInfo& worldInfo,const btVector3& center,
												   const btVector3& radius,
												   int res)
{
	struct	Hammersley
	{
		static void	Generate(btVector3* x,int n)
		{
			for(int i=0;i<n;i++)
			{
				btScalar	p=0.5,t=0;
				for(int j=i;j;p*=0.5,j>>=1) if(j&1) t+=p;
				btScalar	w=2*t-1;
				btScalar	a=(SIMD_PI+2*i*SIMD_PI)/n;
				btScalar	s=btSqrt(1-w*w);
				*x++=btVector3(s*btCos(a),s*btSin(a),w);
			}
		}
	};
	btAlignedObjectArray<btVector3>	vtx;
	vtx.resize(3+res);
	Hammersley::Generate(&vtx[0],vtx.size());
	for(int i=0;i<vtx.size();++i)
	{
		vtx[i]=vtx[i]*radius+center;
	}
	return(CreateFromConvexHull(worldInfo,&vtx[0],vtx.size()));
}



//
btSoftBody*		btSoftBodyHelpers::CreateFromTriMesh(btSoftBodyWorldInfo& worldInfo,const btScalar*	vertices,
													 const int* triangles,
													 int ntriangles, bool randomizeConstraints)
{
	int		maxidx=0;
	int i,j,ni;

	for(i=0,ni=ntriangles*3;i<ni;++i)
	{
		maxidx=btMax(triangles[i],maxidx);
	}
	++maxidx;
	btAlignedObjectArray<bool>		chks;
	btAlignedObjectArray<btVector3>	vtx;
	chks.resize(maxidx*maxidx,false);
	vtx.resize(maxidx);
	for(i=0,j=0,ni=maxidx*3;i<ni;++j,i+=3)
	{
		vtx[j]=btVector3(vertices[i],vertices[i+1],vertices[i+2]);
	}
	btSoftBody*		psb=new btSoftBody(&worldInfo,vtx.size(),&vtx[0],0);
	for( i=0,ni=ntriangles*3;i<ni;i+=3)
	{
		const int idx[]={triangles[i],triangles[i+1],triangles[i+2]};
#define IDX(_x_,_y_) ((_y_)*maxidx+(_x_))
		for(int j=2,k=0;k<3;j=k++)
		{
			if(!chks[IDX(idx[j],idx[k])])
			{
				chks[IDX(idx[j],idx[k])]=true;
				chks[IDX(idx[k],idx[j])]=true;
				psb->appendLink(idx[j],idx[k]);
			}
		}
#undef IDX
		psb->appendFace(idx[0],idx[1],idx[2]);
	}

	if (randomizeConstraints)
	{
		psb->randomizeConstraints();
	}

	return(psb);
}

//
btSoftBody*		btSoftBodyHelpers::CreateFromConvexHull(btSoftBodyWorldInfo& worldInfo,	const btVector3* vertices,
														int nvertices, bool randomizeConstraints)
{
	HullDesc		hdsc(QF_TRIANGLES,nvertices,vertices);
	HullResult		hres;
	HullLibrary		hlib;//?? 
	hdsc.mMaxVertices=nvertices;
	hlib.CreateConvexHull(hdsc,hres);
	btSoftBody*		psb=new btSoftBody(&worldInfo,(int)hres.mNumOutputVertices,
		&hres.m_OutputVertices[0],0);
	for(int i=0;i<(int)hres.mNumFaces;++i)
	{
		const int idx[]={	static_cast<int>(hres.m_Indices[i*3+0]),
							static_cast<int>(hres.m_Indices[i*3+1]),
							static_cast<int>(hres.m_Indices[i*3+2])};
		if(idx[0]<idx[1]) psb->appendLink(	idx[0],idx[1]);
		if(idx[1]<idx[2]) psb->appendLink(	idx[1],idx[2]);
		if(idx[2]<idx[0]) psb->appendLink(	idx[2],idx[0]);
		psb->appendFace(idx[0],idx[1],idx[2]);
	}
	hlib.ReleaseResult(hres);
	if (randomizeConstraints)
	{
		psb->randomizeConstraints();
	}
	return(psb);
}

