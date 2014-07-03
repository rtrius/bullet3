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

#ifndef _BT_SOFT_BODY_H
#define _BT_SOFT_BODY_H

#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btTransform.h"
#include "LinearMath/btIDebugDraw.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"

#include "BulletCollision/CollisionShapes/btConcaveShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
#include "btSparseSDF.h"
#include "BulletCollision/BroadphaseCollision/btDbvt.h"

//#ifdef BT_USE_DOUBLE_PRECISION
//#define btRigidBodyData	btRigidBodyDoubleData
//#define btRigidBodyDataName	"btRigidBodyDoubleData"
//#else
#define btSoftBodyData	btSoftBodyFloatData
#define btSoftBodyDataName	"btSoftBodyFloatData"
//#endif //BT_USE_DOUBLE_PRECISION

class btBroadphaseInterface;
class btDispatcher;
class btSoftBodySolver;

struct	btSoftBodyWorldInfo
{
	btScalar				m_maxDisplacement;
	btBroadphaseInterface*	m_broadphase;
	btDispatcher*	m_dispatcher;
	btVector3				m_gravity;
	btSparseSdf<3>			m_sparsesdf;

	btSoftBodyWorldInfo() :
		m_maxDisplacement(1000.f),//avoid soft body from 'exploding' so use some upper threshold of maximum motion that a node can travel per frame
		m_broadphase(0),
		m_dispatcher(0),
		m_gravity(0,-10,0)
	{
	}
};	

struct btSoftBodyNode
{
	btVector3 m_position;
	btVector3 m_prevPosition;		///<Previous step position; if m_position is overwritten this should be overwritten with the same value.
	btVector3 m_velocity;			///<Rate at which position changes.
	btVector3 m_accumulatedForce;	///<Applied, then set to 0 every frame; if(m_invMass == 0.0) this is set to 0 before forces are applied.
	btVector3 m_normal;				///<Averaged normal over all btSoftBodyFace(s) this node is assigned to.
	btScalar m_invMass;				///<1.0 / mass; if set to 0 the node is fixed and cannot move.
	btScalar m_area;				///<For each btSoftBodyFace that a node is a part of, it receives 1/3 of the area of the face.
	btDbvtNode* m_leaf;				///<BVH Leaf data; m_leaf->data contains the index of this node
	int m_isAttachedToAnchor:1;		///<If nonzero this node is attached to a btSoftBodyAnchor(disables rigid body collision)
};

///Distance constraint between 2 nodes; used for stretch, shear, and bending constraints
struct btSoftBodyLink
{
	int m_linkIndicies[2];		///<Indicies of m_nodes[]
	btScalar m_restLength;		///<Constraint tries to keep nodes at this distance
	btScalar m_linkStiffness;	///<[0,1]; lower stiffness means that Links are easier to stretch.
	
	int m_isBendingLink:1;		///If nonzero this is a bending link, and m_linkStiffness is scaled by m_bendingStiffness instead of m_stretchStiffness.
};

///Triangle; used for rendering, collision and various forces(aerodynamic, pressure)
struct btSoftBodyFace
{
	int m_indicies[3];				///<Indicies of m_nodes[]
};
	


///Applies forces to btSoftBody; requires that the soft body is a single closed triangle mesh.
///The volume calculation will be incorrect if multiple closed triangle meshes are in a single soft body,
///of if the mesh is open.
struct btSoftBodyClosedTrimeshForce
{
	///@name Pressure force; a force that either shrinks or expands the mesh
	///@{
	btScalar m_pressure;				///<[-inf, +inf]; greater than 0 expands, less than 0 shrinks, and 0 == disabled.
	///@}
	
	///@name Volume conservation force; tries to maintain a constant volume
	///@{
	btScalar m_volumeConservation; 		///<[0, +inf]; force magnitude
	btScalar m_restVolume;				///<0.0 == force is disabled; set using btSoftBody::getClosedTrimeshVolume() to enable
	///@}
	
	btSoftBodyClosedTrimeshForce()
	{
		m_pressure = btScalar(0.0);
	
		m_volumeConservation = btScalar(0.0);
		m_restVolume = btScalar(0.0);
	}
	
	///@param closedTrimeshVolume Current volume of the soft body; use btSoftBody::getClosedTrimeshVolume() to obtain.
	void applyForcesToNodes(btScalar closedTrimeshVolume, btAlignedObjectArray<btSoftBodyNode>& nodes) const
	{
		bool applyVolumeForce = m_pressure != btScalar(0.0) || m_volumeConservation > btScalar(0.0);
		if(applyVolumeForce)
		{
			btScalar pressureForce = ( btScalar(1.0) / btFabs(closedTrimeshVolume) ) * m_pressure;
			btScalar volumeForce = (m_restVolume - closedTrimeshVolume) * m_volumeConservation;
			
			btScalar closedTrimeshForce = pressureForce + volumeForce;
			
			for(int i = 0; i < nodes.size(); ++i)
			{
				btSoftBodyNode& n = nodes[i];
				if(n.m_invMass > 0) n.m_accumulatedForce += n.m_normal * (n.m_area * closedTrimeshForce);
			}
		}
	}
};


///Generates and applies an aerodynamic force to a btSoftBody.
///Useful for simulating objects such as flags or falling papers.
struct btSoftBodyAeroForce
{
	enum AeroModel
	{
		V_TwoSided,			///Vertex normals are flipped to match velocity
		V_TwoSidedLiftDrag, ///Vertex normals are flipped to match velocity and lift and drag forces are applied
	};
	
	btSoftBodyAeroForce::AeroModel m_model; 		///<Aerodynamic model (default: V_TwoSided)
	btVector3 m_windVelocity;
	btScalar m_dragCoeff;			///<Range [0, +inf]
	btScalar m_liftCoeff;			///<Range [0, +inf]
	btScalar m_airDensity;
	
	btSoftBodyAeroForce()
	{
		m_model = btSoftBodyAeroForce::V_TwoSided;
		m_windVelocity = btVector3(0,0,0);
		m_dragCoeff = btScalar(0.0);
		m_liftCoeff = btScalar(0.0);
		m_airDensity = btScalar(1.2);
	}
	
	void addAeroForces(btScalar timeStep, btAlignedObjectArray<btSoftBodyNode>& nodes);
	void addAeroForceToNode(btScalar timeStep, btAlignedObjectArray<btSoftBodyNode>& nodes, int nodeIndex);
};

///Moves the soft body nodes to match a certain shape
struct btSoftBodyPose
{
	btScalar m_poseMatching;								///<[0, 1]; 0.0 == disabled; use setPose() to set this
	btAlignedObjectArray<btVector3> m_referencePositions;	///<Pose constraint moves each node to these positions
	btAlignedObjectArray<btScalar> m_weights;				///<Per-node weights; used for computing the center of mass
	btVector3 m_centerOfMass;
	btMatrix3x3 m_rotation;
	
	btSoftBodyPose()
	{
		m_poseMatching = btScalar(0.0);
		m_centerOfMass = btVector3(0,0,0);
		m_rotation.setIdentity();
	}
	
	void initializePose(const btAlignedObjectArray<btSoftBodyNode>& nodes, btScalar totalMass, btScalar poseMatching);
	void updateAndApplyPose(btAlignedObjectArray<btSoftBodyNode>& nodes);
	btVector3 evaluateCenterOfMass(const btAlignedObjectArray<btSoftBodyNode>& nodes) const;
};

struct btSoftBodyConfig
{
	btScalar m_damping;						///<[0, 1]; fraction of velocity removed per timestep(0.05 means that velocity is scaled by 0.95).
	btScalar m_dynamicFriction;				///<[0, 1]
	
	btScalar m_rigidContactHardness;		///<[0, 1]; contact hardness for dynamic rigid bodies.
	btScalar m_kinematicContactHardness;	///<[0, 1]; contact hardness for static and kinematic(inverse_mass == 0) rigid bodies.
	btScalar m_softContactHardness;			///<[0, 1]; hardness(ERP) determines how quickly penetration is resolved.
	btScalar m_anchorHardness;				///<[0, 1]
	
	btScalar m_stretchStiffness;			///<[0, 1]; btSoftBodyLink.m_linkStiffness is multiplied by this to get the actual stiffness.
	btScalar m_bendingStiffness;			///<[0, 1]; used in place of m_stretchStiffness if (btSoftBodyLink.m_isBendingLink != 0).
	
	int m_positionIterations;				///<More iterations improves the quality, but also makes the constraints more stiff
	
	btSoftBodyConfig()
	{
		m_damping = btScalar(0.0);
		m_dynamicFriction = btScalar(0.2);
		
		m_rigidContactHardness = btScalar(1.0);
		m_kinematicContactHardness = btScalar(0.1);
		m_softContactHardness = btScalar(1.0);
		m_anchorHardness = btScalar(0.7);
		
		m_stretchStiffness = btScalar(1.0);
		m_bendingStiffness = btScalar(1.0);
		
		m_positionIterations = 1;
	}
};

class btSoftBody;

struct btSoftBodyRaycastResult
{
	btSoftBody*	body;
	int index;		///<Triangle index
	btScalar fraction;	///<Time of impact fraction (rayorg+(rayto-rayfrom)*fraction); 1.f == no hit
};
	


///Contact between a btRigidBody and node
struct btSoftRigidContact
{
	const btCollisionObject*	m_colObj;		// Rigid body
	btVector3					m_normal;	// Outward normal
	btScalar					m_offset;	// Offset from origin
	int							m_nodeIndex;
	btMatrix3x3					m_impulseMatrix;
	btVector3					m_relativeNodePosition;		///<Position of m_node relative to the btCollisionObject
	btScalar					m_invMassDt;				///<inverse_mass * timeStep
	btScalar					m_combinedFriction;
	btScalar					m_hardness;
};
	
///Node-Face contact between 2 soft bodies
struct btSoftSoftContact
{
	btSoftBody* m_nodeSoftBody;
	btSoftBody* m_faceSoftBody;
	int			m_nodeIndex;	///<Index to m_nodes[] of m_nodeSoftBody
	int			m_faceIndex;	///<Index to m_faces[] of m_faceSoftBody
	btVector3	m_weights;
	btVector3	m_normalOnFace;
	btScalar	m_margin;
	btScalar	m_friction;	
	btScalar	m_nodeCfm;		// Constraint force mixing
	btScalar	m_faceCfm;		// Constraint force mixing
};

///Attaches a btSoftBodyNode to a btRigidBody
struct btSoftBodyAnchor
{
	int						m_nodeIndex;
	btVector3				m_local;		// Anchor position in body space
	btRigidBody*			m_body;			// Body
	btScalar				m_influence;
	btMatrix3x3				m_impulseMatrix;
	btVector3				m_rotatedPosition;
	btScalar				m_invMassDt;		///<inverse_mass * timeStep
};

struct btSoftImplicitShape
{
	virtual btScalar signedDistance(const btVector3& x) = 0;		///<Returns negative if penetrating
};
class btSoftBodyMeshModifier
{
public:
	///Each soft body should have its own shape when calling this
	static void refine(btSoftBody* softBody, btSoftImplicitShape* shape, btScalar accuracy, bool cut);
	
	///Generate bending constraints based on distance in the adjency graph
	static int generateBendingConstraints( btSoftBody* softBody, int distance, btScalar stiffness = btScalar(1.0) );

	///Randomize constraints to reduce solver bias
	static void randomizeConstraints(btAlignedObjectArray<btSoftBodyLink>& links);	
};

///Defines how the nodes of a soft body are connected(constraints, triangles).
///Multiple btSoftBody(s) can share a single shape, if the shape is not modified.
struct btSoftBodyShape
{
	///Default positions and masses of vertices(nodes); m_referencePositions.size() == m_referenceMasses.size()
	btAlignedObjectArray<btVector3> m_referencePositions;
	btAlignedObjectArray<btScalar> m_referenceMasses;
	
	btAlignedObjectArray<btSoftBodyLink> m_links;	///<Distance constraints
	btAlignedObjectArray<btSoftBodyFace> m_faces;	///<Triangles
	
	int numNodes() const { return m_referencePositions.size(); }
	
	bool checkLink(int node0, int node1) const;
	void appendLink(int node0, int node1, btScalar stiffness = btScalar(1.0), bool bcheckexist = false);
	void appendFace(int node0, int node1, int node2);
	void scale(const btVector3& scaling);
};

///The btSoftBody is an class to simulate cloth and volumetric soft bodies. 
///There is two-way interaction between btSoftBody and btRigidBody/btCollisionObject.
class btSoftBody : public btCollisionObject
{
public:

	btAlignedObjectArray<const class btCollisionObject*> m_collisionDisabledObjects;	//btRigidBody or btCollisionObject

	btSoftBodySolver* m_softBodySolver;
	
	btSoftBodyClosedTrimeshForce m_closedTrimeshForce;
	btSoftBodyAeroForce m_aeroForce;
	btSoftBodyPose m_pose;		///<Pose matching constraint; can be used to ensure that the soft body maintains a certain shape
	
	btSoftBodyConfig m_cfg;					// Configuration
	btScalar				m_timeStep;
	btSoftBodyWorldInfo*	m_worldInfo;	// World info
	
	btAlignedObjectArray<btSoftBodyNode> m_nodes;
	btAlignedObjectArray<btSoftBodyAnchor> m_anchors;
	btAlignedObjectArray<btSoftRigidContact> m_rigidContacts;
	btAlignedObjectArray<btSoftSoftContact> m_softContacts;
	
	btAlignedObjectArray<btDbvtNode*> m_faceLeaves;		///<m_faceLeaves[]->data contains the index of the face in m_softShape
	
	btSoftBodyShape* m_softShape;
	
	btVector3 m_aabbMin;
	btVector3 m_aabbMax;
	
	bool					m_bUpdateRtCst;	// Update runtime constants
	btDbvt					m_nodeBvh;		// Nodes tree
	btDbvt					m_faceBvh;		// Faces tree

public:
	btSoftBody(btSoftBodyWorldInfo* worldInfo);
	virtual ~btSoftBody();

	void setupNodesAndShape(btSoftBodyShape* softShape);
	
	btSoftBodyWorldInfo* getWorldInfo() { return m_worldInfo; }

	void setPose(btScalar poseMatching);		///<Enables pose matching constraint if (poseMatching != 0.0)
	
	///@todo: avoid internal softbody shape hack and move collision code to collision library
	virtual void setCollisionShape(btCollisionShape* collisionShape) {}

	
	void appendNode(const btVector3& position, btScalar mass);

	void appendAnchor(int node, btRigidBody* body, bool disableCollisionBetweenLinkedBodies=false,btScalar influence = 1);
	void appendAnchor(int node,btRigidBody* body, const btVector3& localPivot,bool disableCollisionBetweenLinkedBodies=false,btScalar influence = 1);

	void addForceAllNodes(const btVector3& force);
	void addVelocityAllNodes(const btVector3& velocity);
	void setVelocityAllNodes(const btVector3& velocity);
	
	void addForce(const btVector3& force, int node);		//Add force to a node of the body
	void addVelocity(const btVector3& velocity, int node);	//Add velocity to a node of the body
	
	void setMass(int node, btScalar mass);
	btScalar getMass(int node) const;
	btScalar getTotalMass() const;
	void setTotalMass(btScalar mass, bool fromfaces=false);	//Set total mass (weighted by previous masses)
	void setTotalDensity(btScalar density);
	
	void transform(const btTransform& trs);
	void translate(const btVector3& trs);
	void rotate(const btQuaternion& rot);
	
	btScalar getClosedTrimeshVolume() const;		///<Returns the volume, assuming that m_faces represents a closed triangle mesh
	
	///Ray casting using rayFrom and rayTo in worldspace, (not direction!)
	bool rayTest(const btVector3& rayFrom, const btVector3& rayTo, btSoftBodyRaycastResult& results);
	void predictMotion(btScalar timeStep);
	void defaultCollisionHandler(const btCollisionObjectWrapper* pcoWrap);
	void defaultCollisionHandler(btSoftBody* psb);
	

	// Set the solver that handles this soft body
	// Should not be allowed to get out of sync with reality
	// Currently called internally on addition to the world
	void setSoftBodySolver(btSoftBodySolver *softBodySolver) { m_softBodySolver = softBodySolver; }
	
	btSoftBodySolver *getSoftBodySolver() { return m_softBodySolver; }
	btSoftBodySolver *getSoftBodySolver() const { return m_softBodySolver; }

	static const btSoftBody* upcast(const btCollisionObject* colObj)
	{
		if (colObj->getInternalType()==CO_SOFT_BODY)
			return (const btSoftBody*)colObj;
		return 0;
	}
	static btSoftBody* upcast(btCollisionObject* colObj)
	{
		if (colObj->getInternalType()==CO_SOFT_BODY)
			return (btSoftBody*)colObj;
		return 0;
	}

	// ::btCollisionObject
	virtual void getAabb(btVector3& aabbMin,btVector3& aabbMax) const
	{
		aabbMin = m_aabbMin;
		aabbMax = m_aabbMax;
	}
	
	// Private

	int rayTest(const btVector3& rayFrom,const btVector3& rayTo, btScalar& hitFraction, int& faceIndex, bool bcountonly) const;
	void initializeFaceTree();
	bool checkContact(const btCollisionObjectWrapper* colObjWrap,const btVector3& worldSpaceNodePosition,btScalar margin,btSoftRigidContact& contact) const;
	void updateNormals();
	void updateBounds();
	void updateArea();
	void applyForces();	
};


#endif //_BT_SOFT_BODY_H
