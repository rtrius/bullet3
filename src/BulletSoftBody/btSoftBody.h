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


///The btSoftBody is an class to simulate cloth and volumetric soft bodies. 
///There is two-way interaction between btSoftBody and btRigidBody/btCollisionObject.
class	btSoftBody : public btCollisionObject
{
public:
	btAlignedObjectArray<const class btCollisionObject*> m_collisionDisabledObjects;

	// The solver object that handles this soft body
	btSoftBodySolver *m_softBodySolver;

	// Enumerations
	struct	eFeature { enum _ {
		None,
		Node,
		Link,
		Face,
		Tetra
	};};
	
	// API Types
	struct sRayCast
	{
		btSoftBody*	body;		/// soft body
		eFeature::_	feature;	/// feature type
		int			index;		/// feature index
		btScalar	fraction;		/// time of impact fraction (rayorg+(rayto-rayfrom)*fraction)
	};

	struct	ImplicitFn
	{
		virtual btScalar	Eval(const btVector3& x)=0;
	};

	// Internal types

	/// sCti is Softbody contact info
	struct	sCti
	{
		const btCollisionObject*	m_colObj;		// Rigid body
		btVector3		m_normal;	// Outward normal
		btScalar		m_offset;	// Offset from origin
	};	
	
	struct	Element
	{
		void*			m_tag;			// User data
		Element() : m_tag(0) {}
	};
	
	struct	Material : Element
	{
		btScalar				m_linearStiffness;	///<[0,1]; lower stiffness means that Links are easier to stretch.
		bool					m_debugDraw;
	};

	
	struct	Feature : Element
	{
		Material*				m_material;
	};
	
	struct	Node : Feature
	{
		btVector3				m_x;			// Position
		btVector3				m_q;			// Previous step position
		btVector3				m_v;			// Velocity
		btVector3				m_f;			// Force accumulator
		btVector3				m_normal;
		btScalar				m_invMass;		///<1 / mass
		btScalar				m_area;			// Area
		btDbvtNode*				m_leaf;			// Leaf data
		int						m_battach:1;	// Attached
	};
	
	struct	Link : Feature
	{
		Node*					m_n[2];						// Node pointers
		btScalar				m_restLength;				///<Constraint tries to keep nodes at this distance
		int						m_bbending:1;				// Bending link
		btScalar				m_scaledCombinedInvMass;	///<(inverse_mass_0 + inverse_mass_1) * linearStiffness
		btScalar				m_restLengthSquared;
		btScalar				m_impulseScaling;			// |gradient|^2/c0
		btVector3				m_gradient0to1;
	};
	
	struct	Face : Feature
	{
		Node*					m_n[3];			// Node pointers
		btVector3				m_normal;		// Normal
		btScalar				m_area;
		btDbvtNode*				m_leaf;			// Leaf data
	};
	
	struct	Tetra : Feature
	{
		Node*					m_n[4];			// Node pointers		
		btScalar				m_restVolume;
		btDbvtNode*				m_leaf;			// Leaf data
	};
	
	struct	RigidContact
	{
		sCti					m_cti;			// Contact infos
		Node*					m_node;			// Owner node
		btMatrix3x3				m_impulseMatrix;
		btVector3				m_relativeNodePosition;		///<Position of m_node relative to the btCollisionObject
		btScalar				m_invMassDt;				///<inverse_mass * timeStep
		btScalar				m_combinedFriction;
		btScalar				m_hardness;
	};
	
	struct	SoftContact
	{
		Node*					m_node;			// Node
		Face*					m_face;			// Face
		btVector3				m_weights;		// Weigths
		btVector3				m_normal;		// Normal
		btScalar				m_margin;		// Margin
		btScalar				m_friction;		// Friction
		btScalar				m_cfm[2];		// Constraint force mixing
	};
	
	struct	Anchor
	{
		Node*					m_node;			// Node pointer
		btVector3				m_local;		// Anchor position in body space
		btRigidBody*			m_body;			// Body
		btScalar				m_influence;
		btMatrix3x3				m_impulseMatrix;
		btVector3				m_rotatedPosition;
		btScalar				m_invMassDt;		///<inverse_mass * timeStep
	};
	
	struct	Pose
	{
		bool					m_bframe;		// Is frame
		btScalar				m_volume;		// Rest volume
		btAlignedObjectArray<btVector3>			m_pos;			// Reference positions
		btAlignedObjectArray<btScalar>			m_wgh;			// Weights
		btVector3				m_com;			// COM
		btMatrix3x3				m_rot;			// Rotation
		btMatrix3x3				m_scl;			// Scale
		btMatrix3x3				m_aqq;			// Base scaling
	};
	
	struct	Config
	{
		btScalar m_damping;						///<[0, 1]; fraction of velocity removed per timestep(0.05 means that velocity is scaled by 0.95).
		btScalar m_pressure;					///<[-inf, +inf]
		btScalar m_volumeConservation; 			///<[0, +inf]
		btScalar m_dynamicFriction;				///<[0, 1]
		btScalar m_poseMatching;				///<[0, 1]
		btScalar m_rigidContactHardness;		///<[0, 1]; contact hardness for dynamic rigid bodies.
		btScalar m_kinematicContactHardness;	///<[0, 1]; contact hardness for static and kinematic(inverse_mass == 0) rigid bodies.
		btScalar m_softContactHardness;			///<[0, 1]
		btScalar m_anchorHardness;				///<[0, 1]
		btScalar				maxvolume;		// Maximum volume ratio for pose
		int m_velocityIterations;
		int m_positionIterations;
	};
	
	struct	SolverState
	{
		btScalar				sdt;			// dt
	};	

	// Fields

	Config					m_cfg;			// Configuration
	SolverState				m_sst;			// Solver state
	Pose					m_pose;			// Pose
	btSoftBodyWorldInfo*	m_worldInfo;	// World info
	btAlignedObjectArray<Node>				m_nodes;
	btAlignedObjectArray<Link>				m_links;
	btAlignedObjectArray<Face>				m_faces;
	btAlignedObjectArray<Tetra>				m_tetras;
	btAlignedObjectArray<Anchor>			m_anchors;
	btAlignedObjectArray<RigidContact>		m_rigidContacts;
	btAlignedObjectArray<SoftContact>		m_softContacts;
	btAlignedObjectArray<Material*>				m_materials;
	btVector3				m_bounds[2];	// Spatial bounds	
	bool					m_bUpdateRtCst;	// Update runtime constants
	btDbvt					m_ndbvt;		// Nodes tree
	btDbvt					m_fdbvt;		// Faces tree

	btScalar        	m_restLengthScale;
	
	// Api
	btSoftBody(	btSoftBodyWorldInfo* worldInfo, int numNodes, const btVector3* nodePositions, const btScalar* nodeMasses);
	btSoftBody(	btSoftBodyWorldInfo* worldInfo);

	virtual ~btSoftBody();
	
	void initDefaults();

	btSoftBodyWorldInfo* getWorldInfo() { return m_worldInfo; }

	///@todo: avoid internal softbody shape hack and move collision code to collision library
	virtual void setCollisionShape(btCollisionShape* collisionShape) {}

	bool checkLink(int node0, int node1) const;
	bool checkLink(const Node* node0, const Node* node1) const;
	bool checkFace(int node0, int node1, int node2) const;	///<Check for existing face
	
	Material* appendMaterial();
	
	void appendNode(const btVector3& position, btScalar mass);
	
	void appendLink(int model=-1,Material* mat=0);
	void appendLink(int node0, int node1, Material* mat=0, bool bcheckexist=false);
	void appendLink(Node* node0, Node* node1, Material* mat=0, bool bcheckexist=false);
	void appendFace(int model=-1,Material* mat=0);
	void appendFace(int node0, int node1, int node2, Material* mat=0);
	void appendTetra(int model, Material* mat);
	void appendTetra(int node0, int node1, int node2, int node3, Material* mat=0);

	void appendAnchor(int node, btRigidBody* body, bool disableCollisionBetweenLinkedBodies=false,btScalar influence = 1);
	void appendAnchor(int node,btRigidBody* body, const btVector3& localPivot,bool disableCollisionBetweenLinkedBodies=false,btScalar influence = 1);

	
	void addForce(const btVector3& force);									//Add force (or gravity) to the entire body
	void addForce(const btVector3& force, int node);						//Add force (or gravity) to a node of the body
	
	void addVelocity(const btVector3& velocity);	//Add velocity to the entire body	
	void setVelocity(const btVector3& velocity);	//Set velocity for the entire body
	void addVelocity(const btVector3& velocity, int node);	//Add velocity to a node of the body
	void setMass(int node, btScalar mass);
	btScalar getMass(int node) const;
	btScalar getTotalMass() const;
	void setTotalMass(btScalar mass, bool fromfaces=false);	//Set total mass (weighted by previous masses)
	void setTotalDensity(btScalar density);
	void setVolumeMass(btScalar mass);			//Set volume mass (using tetrahedrons)
	void setVolumeDensity(btScalar density);	//Set volume density (using tetrahedrons)
	void transform(const btTransform& trs);
	void translate(const btVector3& trs);
	void rotate(const btQuaternion& rot);
	void scale(	const btVector3& scl);
	btScalar getRestLengthScale();					//Get link resting lengths scale
	void setRestLengthScale(btScalar restLength);	//Scale resting length of all springs
	void setPose(bool bvolume, bool bframe);		//Set current state as pose	
	void resetLinkRestLengths();					//Set current link lengths as resting lengths
	btScalar getVolume() const;
	
	int generateBendingConstraints(int distance, Material* mat=0);	//Generate bending constraints based on distance in the adjency graph
	void randomizeConstraints();	//Randomize constraints to reduce solver bias	
	
	void refine(ImplicitFn* ifn,btScalar accurary,bool cut);
	bool cutLink(int node0,int node1,btScalar position);
	bool cutLink(const Node* node0,const Node* node1,btScalar position);

	///Ray casting using rayFrom and rayTo in worldspace, (not direction!)
	bool rayTest(const btVector3& rayFrom, const btVector3& rayTo, sRayCast& results);
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
		aabbMin = m_bounds[0];
		aabbMax = m_bounds[1];
	}
	
	// Private
	void pointersToIndices();
	void indicesToPointers(const int* map=0);

	int rayTest(const btVector3& rayFrom,const btVector3& rayTo, btScalar& mint,eFeature::_& feature,int& index,bool bcountonly) const;
	void initializeFaceTree();
	btVector3 evaluateCom() const;
	bool checkContact(const btCollisionObjectWrapper* colObjWrap,const btVector3& x,btScalar margin,btSoftBody::sCti& cti) const;
	void updateNormals();
	void updateBounds();
	void updatePose();
	void updateConstants();
	void updateLinkConstants();
	void updateArea(bool averageArea = true);
	void applyForces();	
	
	/************************************************************************************
	///Aero force
	************************************************************************************/
	struct eAeroModel 
	{ 
		enum _ 
		{
			V_Point,			///Vertex normals are oriented toward velocity
			V_TwoSided,			///Vertex normals are flipped to match velocity	
			V_TwoSidedLiftDrag, ///Vertex normals are flipped to match velocity and lift and drag forces are applied
			V_OneSided,			///Vertex normals are taken as it is	
			F_TwoSided,			///Face normals are flipped to match velocity
			F_TwoSidedLiftDrag,	///Face normals are flipped to match velocity and lift and drag forces are applied 
			F_OneSided,			///Face normals are taken as it is		
			END
		};
	};
	
	///Generates and applies an aerodynamic force to a btSoftBody.
	///Useful for simulating objects such as flags or falling papers.
	struct AeroForce
	{
		eAeroModel::_ m_model; 		///<Aerodynamic model (default: V_Point)
		btVector3 m_windVelocity;
		btScalar m_dragCoeff;			///<Range [0, +inf]
		btScalar m_liftCoeff;			///<Range [0, +inf]
		btScalar m_airDensity;
		
		AeroForce()
		{
			m_model = eAeroModel::V_Point;
			m_windVelocity = btVector3(0,0,0);
			m_dragCoeff = btScalar(0.0);
			m_liftCoeff = btScalar(0.0);
			m_airDensity = btScalar(1.2);
		}
	};
	static void addAeroForces(const AeroForce& aeroForce, btScalar timeStep, btAlignedObjectArray<Node>& nodes, btAlignedObjectArray<Face>& faces);
	static void addAeroForceToNode(const AeroForce& aeroForce, btScalar timeStep, btAlignedObjectArray<Node>& nodes, int nodeIndex);	//Add aero force to a node of the body
	static void addAeroForceToFace(const AeroForce& aeroForce, btScalar timeStep, btAlignedObjectArray<Face>& faces, int faceIndex);	//Add aero force to a face of the body

	AeroForce m_aeroForce;
	
	/************************************************************************************
	///Aero force
	************************************************************************************/
};


#endif //_BT_SOFT_BODY_H
