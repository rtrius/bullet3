
#ifdef cl_amd_printf
	#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

typedef float b3Scalar;
typedef float4 b3Vector3;
#define b3Sqrt sqrt
#define b3Max max
#define b3Min min


//Note that these are vector3 functions -- OpenCL functions are vector4 functions
inline b3Scalar b3Vector3_length2(b3Vector3 v) { return v.x*v.x + v.y*v.y + v.z*v.z; }
inline b3Scalar b3Vector3_dot(b3Vector3 a, b3Vector3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline b3Vector3 b3Vector3_normalize(b3Vector3 v)
{
	b3Scalar length2 = b3Vector3_length2(v);
	if( length2 != (b3Scalar)0.0f ) v /= sqrt(length2);
	
	return v;
}
//Defined in b3FluidSortingGrid.h
#define INVALID_FIRST_INDEX -1
#define INVALID_LAST_INDEX -2

//Syncronize with 'struct b3FluidSphParameters' in b3FluidSphParameters.h
typedef struct
{
	b3Scalar m_timeStep;
	b3Scalar m_simulationScale;
	b3Scalar m_sphSmoothRadius;
	b3Scalar m_sphRadiusSquared;
	b3Scalar m_poly6KernCoeff;
	b3Scalar m_spikyKernGradCoeff;
	b3Scalar m_viscosityKernLapCoeff;
	b3Scalar m_initialSum;
	b3Vector3 m_aabbBoundaryMin;
	b3Vector3 m_aabbBoundaryMax;
	int m_enableAabbBoundary;
	b3Vector3 m_gravity;
	b3Scalar m_sphAccelLimit;
	b3Scalar m_speedLimit;
	b3Scalar m_viscosity;
	b3Scalar m_restDensity;
	b3Scalar m_sphParticleMass;
	b3Scalar m_stiffness;
	b3Scalar m_particleDist;
	b3Scalar m_particleRadius;
	b3Scalar m_particleMargin;
	b3Scalar m_particleMass;
	b3Scalar m_boundaryStiff;
	b3Scalar m_boundaryDamp;
	b3Scalar m_boundaryFriction;
	b3Scalar m_boundaryRestitution;
	b3Scalar m_boundaryErp;
} b3FluidSphParameters;


//#define B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT	//Ensure that this is also #defined in b3FluidSortingGrid.h
#ifdef B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT	
	typedef unsigned long b3FluidGridUint64;
	typedef b3FluidGridUint64 b3FluidGridCombinedPos;	//Range must contain B3_FLUID_GRID_COORD_RANGE^3
	#define B3_FLUID_GRID_COORD_RANGE 2097152		//2^21
	
	inline void splitCombinedPosition(b3FluidGridUint64 resolutionX, b3FluidGridUint64 resolutionY, 
										b3FluidGridUint64 value, int* out_x, int* out_y, int* out_z)
	{
		b3FluidGridUint64 cellsPerLine = resolutionX;
		b3FluidGridUint64 cellsPerPlane = resolutionX * resolutionY;
		
		b3FluidGridUint64 x = value % cellsPerLine;
		b3FluidGridUint64 z = value / cellsPerPlane;
		b3FluidGridUint64 y = (value - z*cellsPerPlane) / cellsPerLine;
		
		*out_x = (int)x;
		*out_z = (int)z;
		*out_y = (int)y;
	}
#else
	typedef unsigned int b3FluidGridCombinedPos;		//Range must contain B3_FLUID_GRID_COORD_RANGE^3
	#define B3_FLUID_GRID_COORD_RANGE 1024			//2^10	
	
	inline void splitCombinedPosition(int resolutionX, int resolutionY, int value, int* out_x, int* out_y, int* out_z)
	{
		int x = value % resolutionX;
		int z = value / (resolutionX*resolutionY);
		int y = (value - z*resolutionX*resolutionY) / resolutionX;
		
		*out_x = (int)x;
		*out_z = (int)z;
		*out_y = (int)y;
	}
#endif

typedef int b3FluidGridCoordinate;
#define B3_FLUID_GRID_COORD_RANGE_HALVED B3_FLUID_GRID_COORD_RANGE/2



typedef struct
{
	int m_firstIndex;
	int m_lastIndex;
	
} b3FluidGridIterator;


//Since the hash function used to determine the 'value' of particles is simply 
//(x + y*CELLS_PER_ROW + z*CELLS_PER_PLANE), adjacent cells have a value 
//that is 1 greater and lesser than the current cell. 
//This makes it possible to query 3 cells simultaneously(as a 3 cell bar extended along the x-axis) 
//by using a 'binary range search' in the range [current_cell_value-1, current_cell_value+1]. 
//Furthermore, as the 3 particle index ranges returned are also adjacent, it is also possible to 
//stitch them together to form a single index range.
#define b3FluidSortingGrid_NUM_FOUND_CELLS_GPU 9

typedef struct
{
	b3FluidGridIterator m_iterators[b3FluidSortingGrid_NUM_FOUND_CELLS_GPU];
	
} b3FluidSortingGridFoundCellsGpu;		//b3FluidSortingGrid::FoundCellsGpu in b3FluidSortingGrid.h

typedef struct 
{
	b3FluidGridCombinedPos m_value;
	int m_index;
	
} b3FluidGridValueIndexPair;

typedef struct
{
	b3FluidGridCoordinate x;		
	b3FluidGridCoordinate y;
	b3FluidGridCoordinate z;
	b3FluidGridCoordinate padding;
	
} b3FluidGridPosition;

b3FluidGridPosition getDiscretePosition(b3Scalar cellSize, b3Vector3 position)	//b3FluidSortingGrid::getDiscretePosition()
{
	b3Vector3 discretePosition = position / cellSize;
	
	b3FluidGridPosition result;
	result.x = (b3FluidGridCoordinate)( (position.x >= 0.0f) ? discretePosition.x : floor(discretePosition.x) );
	result.y = (b3FluidGridCoordinate)( (position.y >= 0.0f) ? discretePosition.y : floor(discretePosition.y) );
	result.z = (b3FluidGridCoordinate)( (position.z >= 0.0f) ? discretePosition.z : floor(discretePosition.z) );
	
	return result;
}
b3FluidGridCombinedPos getCombinedPosition(b3FluidGridPosition quantizedPosition)	//b3FluidGridPosition::getCombinedPosition()
{
	b3FluidGridCoordinate signedX = quantizedPosition.x + B3_FLUID_GRID_COORD_RANGE_HALVED;
	b3FluidGridCoordinate signedY = quantizedPosition.y + B3_FLUID_GRID_COORD_RANGE_HALVED;
	b3FluidGridCoordinate signedZ = quantizedPosition.z + B3_FLUID_GRID_COORD_RANGE_HALVED;
	
	b3FluidGridCombinedPos unsignedX = (b3FluidGridCombinedPos)signedX;
	b3FluidGridCombinedPos unsignedY = (b3FluidGridCombinedPos)signedY * B3_FLUID_GRID_COORD_RANGE;
	b3FluidGridCombinedPos unsignedZ = (b3FluidGridCombinedPos)signedZ * B3_FLUID_GRID_COORD_RANGE * B3_FLUID_GRID_COORD_RANGE;
	
	return unsignedX + unsignedY + unsignedZ;
}

inline int binarySearch(__global b3FluidGridCombinedPos *sortGridValues, int sortGridValuesSize, b3FluidGridCombinedPos value)
{
	//From b3AlignedObjectArray::findBinarySearch()
	//Assumes sortGridValues[] is sorted
	
	int first = 0;
	int last = sortGridValuesSize - 1;
	
	while(first <= last) 
	{
		int mid = (first + last) / 2;
		if(value > sortGridValues[mid]) first = mid + 1;
		else if(value < sortGridValues[mid]) last = mid - 1;
		else return mid;
	}

	return sortGridValuesSize;
}


#define B3_EPSILON FLT_EPSILON

//Syncronize with defines b3FluidSphRigidInteractorCL.h
#define MAX_FLUID_RIGID_PAIRS 32
#define MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE 16
#define MAX_FLUID_CONTACTS_PER_DYNAMIC_RIGID 256

///Contains the indicies of rigid bodies whose AABB intersects with that of a single fluid particle
typedef struct
{
	int m_numIndicies;
	int m_rigidIndicies[MAX_FLUID_RIGID_PAIRS];
	int m_rigidSubIndicies[MAX_FLUID_RIGID_PAIRS];	//Only used if the rigid has triangle mesh or compound shape
	
} FluidRigidPairs;

typedef struct
{
	int m_numContacts;
	int m_rigidIndicies[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];
	b3Scalar m_distances[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];
	b3Vector3 m_pointsOnRigid[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];		//World space point
	b3Vector3 m_normalsOnRigid[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];		//World space normal
	
} FluidRigidContacts;



///Contains indicies of fluid particles that are colliding with dynamic(invMass != 0.0) rigid bodies
///If the rigid body is static, m_numContacts should be 0; this struct is used to apply fluid-rigid impulses to the rigid bodies
typedef struct
{
	int m_numContacts;
	
	//Fluid particle indicies; actual contact data is stored in the FluidRigidContact struct(1 per particle)
	int m_fluidIndicies[MAX_FLUID_CONTACTS_PER_DYNAMIC_RIGID];
	int m_contactIndicies[MAX_FLUID_CONTACTS_PER_DYNAMIC_RIGID];		//range [0, MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE - 1]
	
} RigidFluidContacts;


// -----------------------------------------------------------------------------
// Bullet3
// -----------------------------------------------------------------------------
#define SHAPE_CONVEX_HULL 3
#define SHAPE_PLANE 4
#define SHAPE_CONCAVE_TRIMESH 5
#define SHAPE_COMPOUND_OF_CONVEX_HULLS 6
#define SHAPE_SPHERE 7

#define mymake_float4 (float4)

typedef unsigned int u32;

typedef struct
{
	float4 m_row[3];
}Matrix3x3;

__inline
float dot3F4(float4 a, float4 b);

__inline
float4 mtMul1(Matrix3x3 a, float4 b)
{
	float4 ans;
	ans.x = dot3F4( a.m_row[0], b );
	ans.y = dot3F4( a.m_row[1], b );
	ans.z = dot3F4( a.m_row[2], b );
	ans.w = 0.f;
	return ans;
}

__inline
float4 mtMul3(float4 a, Matrix3x3 b)
{
	float4 colx = mymake_float4(b.m_row[0].x, b.m_row[1].x, b.m_row[2].x, 0);
	float4 coly = mymake_float4(b.m_row[0].y, b.m_row[1].y, b.m_row[2].y, 0);
	float4 colz = mymake_float4(b.m_row[0].z, b.m_row[1].z, b.m_row[2].z, 0);

	float4 ans;
	ans.x = dot3F4( a, colx );
	ans.y = dot3F4( a, coly );
	ans.z = dot3F4( a, colz );
	return ans;
}

typedef struct
{
	int m_numChildShapes;
	float m_radius;
	int m_shapeType;
	int m_shapeIndex;
	
} btCollidableGpu;

typedef struct
{
	float4	m_childPosition;
	float4	m_childOrientation;
	int m_shapeIndex;
	int m_unused0;
	int m_unused1;
	int m_unused2;
} btGpuChildShape;

typedef struct
{
	float4 m_pos;
	float4 m_quat;
	float4 m_linVel;
	float4 m_angVel;

	u32 m_collidableIdx;
	float m_invMass;
	float m_restituitionCoeff;
	float m_frictionCoeff;
} BodyData;				//b3RigidBodyCL in C++ (Bullet3Collision/NarrowPhaseCollision/b3RigidBodyCL.h)

typedef struct
{
	Matrix3x3 m_invInertia;
	Matrix3x3 m_initInvInertia;
} InertiaTensor;		//b3InertiaCL in C++ (Bullet3Collision/NarrowPhaseCollision/b3RigidBodyCL.h)

typedef struct  
{
	float4		m_localCenter;
	float4		m_extents;
	float4		mC;
	float4		mE;
	
	float			m_radius;
	int	m_faceOffset;
	int m_numFaces;
	int	m_numVertices;

	int m_vertexOffset;
	int	m_uniqueEdgesOffset;
	int	m_numUniqueEdges;
	int m_unused;
} ConvexPolyhedronCL;	//b3ConvexPolyhedronCL in C++ (Bullet3OpenCL/NarrowphaseCollision/b3ConvexPolyhedronCL.h)

typedef struct 
{
	union
	{
		float4	m_min;
		float   m_minElems[4];
		int			m_minIndices[4];
	};
	union
	{
		float4	m_max;
		float   m_maxElems[4];
		int			m_maxIndices[4];
	};
} btAabbCL;		//b3SapAabb in C++ (Bullet3OpenCL/BroadphaseCollision/b3SapAabb.h)

typedef struct
{
	float4 m_plane;
	int m_indexOffset;
	int m_numIndices;
} btGpuFace;

#define make_float4 (float4)

typedef float4 Quaternion;

__inline
Quaternion qtMul(Quaternion a, Quaternion b);

__inline
float4 qtRotate(Quaternion q, float4 vec);

__inline
Quaternion qtInvert(Quaternion q);

__inline
float4 cross3(float4 a, float4 b)
{
	return cross(a,b);
}

__inline
float dot3F4(float4 a, float4 b)
{
	float4 a1 = make_float4(a.xyz,0.f);
	float4 b1 = make_float4(b.xyz,0.f);
	return dot(a1, b1);
}

__inline
Quaternion qtMul(Quaternion a, Quaternion b)
{
	Quaternion ans;
	ans = cross3( a, b );
	ans += a.w*b+b.w*a;
//	ans.w = a.w*b.w - (a.x*b.x+a.y*b.y+a.z*b.z);
	ans.w = a.w*b.w - dot3F4(a, b);
	return ans;
}

__inline
float4 qtRotate(Quaternion q, float4 vec)
{
	Quaternion qInv = qtInvert( q );
	float4 vcpy = vec;
	vcpy.w = 0.f;
	float4 out = qtMul(qtMul(q,vcpy),qInv);
	return out;
}

__inline
Quaternion qtInvert(Quaternion q)
{
	return (Quaternion)(-q.xyz, q.w);
}

__inline
float4 qtInvRotate(const Quaternion q, float4 vec)
{
	return qtRotate( qtInvert( q ), vec );
}

__inline
float4 transform(const float4* p, const float4* translation, const Quaternion* orientation)
{
	return qtRotate( *orientation, *p ) + (*translation);
}

void	trInverse(float4 translationIn, Quaternion orientationIn,
		float4* translationOut, Quaternion* orientationOut)
{
	*orientationOut = qtInvert(orientationIn);
	*translationOut = qtRotate(*orientationOut, -translationIn);
}

void	trMul(float4 translationA, Quaternion orientationA,
						float4 translationB, Quaternion orientationB,
		float4* translationOut, Quaternion* orientationOut)
{
	*orientationOut = qtMul(orientationA,orientationB);
	*translationOut = transform(&translationB,&translationA,&orientationA);
}

float signedDistanceFromPointToPlane(float4 point, float4 planeEqn, float4* closestPointOnFace)
{
	float4 n = (float4)(planeEqn.x, planeEqn.y, planeEqn.z, 0);
	float dist = dot3F4(n, point) + planeEqn.w;
	*closestPointOnFace = point - dist * n;
	return dist;
}

inline bool IsPointInPolygon(float4 p, 
							const btGpuFace* face,
							__global const float4* baseVertex,
							__global const  int* convexIndices,
							float4* out)
{
    float4 a;
    float4 b;
    float4 ab;
    float4 ap;
    float4 v;

	float4 plane = make_float4(face->m_plane.x,face->m_plane.y,face->m_plane.z,0.f);
	
	if (face->m_numIndices<2)
		return false;

	
	float4 v0 = baseVertex[convexIndices[face->m_indexOffset + face->m_numIndices-1]];
	
	b = v0;

    for(unsigned i=0; i != face->m_numIndices; ++i)
    {
		a = b;
		float4 vi = baseVertex[convexIndices[face->m_indexOffset + i]];
		b = vi;
        ab = b-a;
        ap = p-a;
        v = cross3(ab,plane);

        if (dot(ap, v) > 0.f)
        {
            float ab_m2 = dot(ab, ab);
            float rt = ab_m2 != 0.f ? dot(ab, ap) / ab_m2 : 0.f;
            if (rt <= 0.f)
            {
                *out = a;
            }
            else if (rt >= 1.f) 
            {
                *out = b;
            }
            else
            {
            	float s = 1.f - rt;
				out[0].x = s * a.x + rt * b.x;
				out[0].y = s * a.y + rt * b.y;
				out[0].z = s * a.z + rt * b.z;
            }
            return false;
        }
    }
    return true;
}

//Modified computeContactSphereConvex() from primitiveContacts.cl
bool computeContactSphereConvex
(
	int collidableIndexB, 
	__global const btCollidableGpu* collidables,
	__global const ConvexPolyhedronCL* convexShapes,
	__global const float4* convexVertices,
	__global const int* convexIndices,
	__global const btGpuFace* faces,
	float4 spherePos2,
	float radius,
	float4 pos,
	float4 quat,

	float* out_distance,
	float4* out_normalOnRigidWorld,
	float4* out_pointOnRigidWorld
)
{
	float4 invPos;
	float4 invOrn;

	trInverse(pos,quat, &invPos,&invOrn);

	float4 spherePos = transform(&spherePos2,&invPos,&invOrn);

	int shapeIndex = collidables[collidableIndexB].m_shapeIndex;
	int numFaces = convexShapes[shapeIndex].m_numFaces;
	float4 closestPnt = (float4)(0, 0, 0, 0);
	float4 hitNormalWorld = (float4)(0, 0, 0, 0);
	float minDist = -1000000.f;
	bool bCollide = true;

	for ( int f = 0; f < numFaces; f++ )
	{
		btGpuFace face = faces[convexShapes[shapeIndex].m_faceOffset+f];

		// set up a plane equation 
		float4 planeEqn;
		float4 n1 = face.m_plane;
		n1.w = 0.f;
		planeEqn = n1;
		planeEqn.w = face.m_plane.w;
		
	
		// compute a signed distance from the vertex in cloth to the face of rigidbody.
		float4 pntReturn;
		float dist = signedDistanceFromPointToPlane(spherePos, planeEqn, &pntReturn);

		// If the distance is positive, the plane is a separating plane. 
		if ( dist > radius )
		{
			bCollide = false;
			break;
		}


		if (dist>0)
		{
			//might hit an edge or vertex
			float4 out;
			float4 zeroPos = make_float4(0,0,0,0);

			bool isInPoly = IsPointInPolygon(spherePos,
					&face,
					&convexVertices[convexShapes[shapeIndex].m_vertexOffset],
					convexIndices,
           &out);
			if (isInPoly)
			{
				if (dist>minDist)
				{
					minDist = dist;
					closestPnt = pntReturn;
					hitNormalWorld = planeEqn;
					
				}
			} else
			{
				float4 tmp = spherePos-out;
				float l2 = dot(tmp,tmp);
				if (l2<radius*radius)
				{
					dist  = sqrt(l2);
					if (dist>minDist)
					{
						minDist = dist;
						closestPnt = out;
						hitNormalWorld = tmp/dist;
						
					}
					
				} else
				{
					bCollide = false;
					break;
				}
			}
		} else
		{
			if ( dist > minDist )
			{
				minDist = dist;
				closestPnt = pntReturn;
				hitNormalWorld.xyz = planeEqn.xyz;
			}
		}
		
	}

	

	if (bCollide && minDist > -10000)
	{
		float4 normalOnSurfaceB1 = qtRotate(quat,-hitNormalWorld);
		float4 pOnB1 = transform(&closestPnt,&pos,&quat);
		
		float actualDepth = minDist-radius;
		if (actualDepth<=0.f)
		{
			*out_distance = actualDepth;
			*out_normalOnRigidWorld = normalOnSurfaceB1;
			*out_pointOnRigidWorld = pOnB1;
			return true;
		}
	}
	
	return false;
}

bool pointInTriangle(const float4* vertices, const float4* normal, float4 *p )
{

	const float4* p1 = &vertices[0];
	const float4* p2 = &vertices[1];
	const float4* p3 = &vertices[2];

	float4 edge1;	edge1 = (*p2 - *p1);
	float4 edge2;	edge2 = ( *p3 - *p2 );
	float4 edge3;	edge3 = ( *p1 - *p3 );

	
	float4 p1_to_p; p1_to_p = ( *p - *p1 );
	float4 p2_to_p; p2_to_p = ( *p - *p2 );
	float4 p3_to_p; p3_to_p = ( *p - *p3 );

	float4 edge1_normal; edge1_normal = ( cross(edge1,*normal));
	float4 edge2_normal; edge2_normal = ( cross(edge2,*normal));
	float4 edge3_normal; edge3_normal = ( cross(edge3,*normal));

	
	
	float r1, r2, r3;
	r1 = dot(edge1_normal,p1_to_p );
	r2 = dot(edge2_normal,p2_to_p );
	r3 = dot(edge3_normal,p3_to_p );
	
	if ( r1 > 0 && r2 > 0 && r3 > 0 )
		return true;
    if ( r1 <= 0 && r2 <= 0 && r3 <= 0 ) 
		return true;
	return false;

}


float segmentSqrDistance(float4 from, float4 to,float4 p, float4* nearest) 
{
	float4 diff = p - from;
	float4 v = to - from;
	float t = dot(v,diff);
	
	if (t > 0) 
	{
		float dotVV = dot(v,v);
		if (t < dotVV) 
		{
			t /= dotVV;
			diff -= t*v;
		} else 
		{
			t = 1;
			diff -= v;
		}
	} else
	{
		t = 0;
	}
	*nearest = from + t*v;
	return dot(diff,diff);	
}

//Modified computeContactSphereTriangle() from primitiveContacts.cl
bool computeContactSphereTriangle(const float4* triangleVertices, 
								float4 spherePos2, 
								float radius, 
								float4 rigidPos, 
								float4 rigidOrn,
									
								float* out_distance,
								float4* out_normalOnRigidWorld,
								float4* out_pointOnRigidWorld)
{

	float4 invPos;
	float4 invOrn;

	trInverse(rigidPos,rigidOrn, &invPos,&invOrn);
	float4 spherePos = transform(&spherePos2,&invPos,&invOrn);
	float4 closestPnt = (float4)(0, 0, 0, 0);
	float4 hitNormalWorld = (float4)(0, 0, 0, 0);
	float minDist = -1000000.f;
	bool bCollide = true;

	
	//////////////////////////////////////

	float4 sphereCenter;
	sphereCenter = spherePos;

	const float4* vertices = triangleVertices;
	float contactBreakingThreshold = 0.f;//todo?
	float radiusWithThreshold = radius + contactBreakingThreshold;
	float4 edge10;
	edge10 = vertices[1]-vertices[0];
	edge10.w = 0.f;//is this needed?
	float4 edge20;
	edge20 = vertices[2]-vertices[0];
	edge20.w = 0.f;//is this needed?
	float4 normal = cross3(edge10,edge20);
	normal = normalize(normal);
	float4 p1ToCenter;
	p1ToCenter = sphereCenter - vertices[0];
	
	float distanceFromPlane = dot(p1ToCenter,normal);

	if (distanceFromPlane < 0.f)
	{
		//triangle facing the other way
		distanceFromPlane *= -1.f;
		normal *= -1.f;
	}
	hitNormalWorld = normal;

	bool isInsideContactPlane = distanceFromPlane < radiusWithThreshold;
	
	// Check for contact / intersection
	bool hasContact = false;
	float4 contactPoint;
	if (isInsideContactPlane) 
	{
	
		if (pointInTriangle(vertices,&normal, &sphereCenter)) 
		{
			// Inside the contact wedge - touches a point on the shell plane
			hasContact = true;
			contactPoint = sphereCenter - normal*distanceFromPlane;
			
		} else {
			// Could be inside one of the contact capsules
			float contactCapsuleRadiusSqr = radiusWithThreshold*radiusWithThreshold;
			float4 nearestOnEdge;
			int numEdges = 3;
			for (int i = 0; i < numEdges; i++) 
			{
				float4 pa =vertices[i];
				float4 pb = vertices[(i+1)%3];

				float distanceSqr = segmentSqrDistance(pa,pb,sphereCenter, &nearestOnEdge);
				if (distanceSqr < contactCapsuleRadiusSqr) 
				{
					// Yep, we're inside a capsule
					hasContact = true;
					contactPoint = nearestOnEdge;
					
				}
				
			}
		}
	}

	if (hasContact) 
	{

		closestPnt = contactPoint;
		float4 contactToCenter = sphereCenter - contactPoint;
		minDist = length(contactToCenter);
		if (minDist>0.f)
		{
			hitNormalWorld = normalize(contactToCenter);//*(1./minDist);
		}
		bCollide  = true;
	}


	/////////////////////////////////////

	if (bCollide && minDist > -10000)
	{
		
		float4 normalOnSurfaceB1 = qtRotate(rigidOrn,-hitNormalWorld);
		float4 pOnB1 = transform(&closestPnt,&rigidPos,&rigidOrn);
		float actualDepth = minDist-radius;

		
		if (actualDepth<=0.f)
		{
			*out_distance = actualDepth;
			*out_normalOnRigidWorld = normalOnSurfaceB1;
			*out_pointOnRigidWorld = pOnB1;
		}
	}

	return hasContact;
}
// -----------------------------------------------------------------------------
// Bullet3
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Bullet3 - bvhTraversal.cl
// -----------------------------------------------------------------------------
#define MAX_NUM_PARTS_IN_BITS 10

///btQuantizedBvhNode is a compressed aabb node, 16 bytes.
///Node can be used for leafnode or internal node. Leafnodes can point to 32-bit triangle index (non-negative range).
typedef struct
{
	//12 bytes
	unsigned short int	m_quantizedAabbMin[3];
	unsigned short int	m_quantizedAabbMax[3];
	//4 bytes
	int	m_escapeIndexOrTriangleIndex;
} btQuantizedBvhNode;

typedef struct
{
	float4		m_aabbMin;
	float4		m_aabbMax;
	float4		m_quantization;
	int			m_numNodes;
	int			m_numSubTrees;
	int			m_nodeOffset;
	int			m_subTreeOffset;

} b3BvhInfo;


int	getTriangleIndex(const btQuantizedBvhNode* rootNode)
{
	unsigned int x=0;
	unsigned int y = (~(x&0))<<(31-MAX_NUM_PARTS_IN_BITS);
	// Get only the lower bits where the triangle index is stored
	return (rootNode->m_escapeIndexOrTriangleIndex&~(y));
}

int isLeaf(const btQuantizedBvhNode* rootNode)
{
	//skipindex is negative (internal node), triangleindex >=0 (leafnode)
	return (rootNode->m_escapeIndexOrTriangleIndex >= 0)? 1 : 0;
}
	
int getEscapeIndex(const btQuantizedBvhNode* rootNode)
{
	return -rootNode->m_escapeIndexOrTriangleIndex;
}

typedef struct
{
	//12 bytes
	unsigned short int	m_quantizedAabbMin[3];
	unsigned short int	m_quantizedAabbMax[3];
	//4 bytes, points to the root of the subtree
	int			m_rootNodeIndex;
	//4 bytes
	int			m_subtreeSize;
	int			m_padding[3];
} btBvhSubtreeInfo;


int testQuantizedAabbAgainstQuantizedAabb(
								const unsigned short int* aabbMin1,
								const unsigned short int* aabbMax1,
								const unsigned short int* aabbMin2,
								const unsigned short int* aabbMax2)
{
	//int overlap = 1;
	if (aabbMin1[0] > aabbMax2[0])
		return 0;
	if (aabbMax1[0] < aabbMin2[0])
		return 0;
	if (aabbMin1[1] > aabbMax2[1])
		return 0;
	if (aabbMax1[1] < aabbMin2[1])
		return 0;
	if (aabbMin1[2] > aabbMax2[2])
		return 0;
	if (aabbMax1[2] < aabbMin2[2])
		return 0;
	return 1;
	//overlap = ((aabbMin1[0] > aabbMax2[0]) || (aabbMax1[0] < aabbMin2[0])) ? 0 : overlap;
	//overlap = ((aabbMin1[2] > aabbMax2[2]) || (aabbMax1[2] < aabbMin2[2])) ? 0 : overlap;
	//overlap = ((aabbMin1[1] > aabbMax2[1]) || (aabbMax1[1] < aabbMin2[1])) ? 0 : overlap;
	//return overlap;
}


void quantizeWithClamp(unsigned short* out, float4 point2,int isMax, float4 bvhAabbMin, float4 bvhAabbMax, float4 bvhQuantization)
{
	float4 clampedPoint = max(point2,bvhAabbMin);
	clampedPoint = min (clampedPoint, bvhAabbMax);

	float4 v = (clampedPoint - bvhAabbMin) * bvhQuantization;
	if (isMax)
	{
		out[0] = (unsigned short) (((unsigned short)(v.x+1.f) | 1));
		out[1] = (unsigned short) (((unsigned short)(v.y+1.f) | 1));
		out[2] = (unsigned short) (((unsigned short)(v.z+1.f) | 1));
	} else
	{
		out[0] = (unsigned short) (((unsigned short)(v.x) & 0xfffe));
		out[1] = (unsigned short) (((unsigned short)(v.y) & 0xfffe));
		out[2] = (unsigned short) (((unsigned short)(v.z) & 0xfffe));
	}

}

//Modified bvhTraversalKernel() from bvhTraversal.cl
void getIntersectingBvhLeaves(__global const btCollidableGpu* collidables,
							__global const btBvhSubtreeInfo* subtreeHeadersRoot,
							__global const btQuantizedBvhNode* quantizedNodesRoot,
							__global const b3BvhInfo* bvhInfos,
						
							int rigidBodyIndex, 
							int trimeshCollidableIndex,
							int particleIndex, 
							b3Vector3 particleAabbMin, 
							b3Vector3 particleAabbMax,
							
							__global FluidRigidPairs* out_pairs)
{
	btCollidableGpu trimeshCollidable = collidables[trimeshCollidableIndex];
	if(trimeshCollidable.m_shapeType != SHAPE_CONCAVE_TRIMESH) return;


	b3BvhInfo bvhInfo = bvhInfos[trimeshCollidable.m_numChildShapes];

	float4 bvhAabbMin = bvhInfo.m_aabbMin;
	float4 bvhAabbMax = bvhInfo.m_aabbMax;
	float4 bvhQuantization = bvhInfo.m_quantization;
	int numSubtreeHeaders = bvhInfo.m_numSubTrees;
	__global const btBvhSubtreeInfo* subtreeHeaders = &subtreeHeadersRoot[bvhInfo.m_subTreeOffset];
	__global const btQuantizedBvhNode* quantizedNodes = &quantizedNodesRoot[bvhInfo.m_nodeOffset];
	

	unsigned short int quantizedQueryAabbMin[3];
	unsigned short int quantizedQueryAabbMax[3];
	quantizeWithClamp(quantizedQueryAabbMin,particleAabbMin,false,bvhAabbMin, bvhAabbMax,bvhQuantization);
	quantizeWithClamp(quantizedQueryAabbMax,particleAabbMax,true ,bvhAabbMin, bvhAabbMax,bvhQuantization);
	
	for (int i=0;i<numSubtreeHeaders;i++)
	{
		btBvhSubtreeInfo subtree = subtreeHeaders[i];
				
		int overlap = testQuantizedAabbAgainstQuantizedAabb(quantizedQueryAabbMin,quantizedQueryAabbMax,subtree.m_quantizedAabbMin,subtree.m_quantizedAabbMax);
		if (overlap != 0)
		{
			int startNodeIndex = subtree.m_rootNodeIndex;
			int endNodeIndex = subtree.m_rootNodeIndex+subtree.m_subtreeSize;
			int curIndex = startNodeIndex;
			int escapeIndex;
			int isLeafNode;
			int aabbOverlap;
			while (curIndex < endNodeIndex)
			{
				btQuantizedBvhNode rootNode = quantizedNodes[curIndex];
				aabbOverlap = testQuantizedAabbAgainstQuantizedAabb(quantizedQueryAabbMin,quantizedQueryAabbMax,rootNode.m_quantizedAabbMin,rootNode.m_quantizedAabbMax);
				isLeafNode = isLeaf(&rootNode);
				if (aabbOverlap)
				{
					if (isLeafNode)
					{
						int triangleIndex = getTriangleIndex(&rootNode);
						
						int pairIndex = atomic_inc(&out_pairs[particleIndex].m_numIndicies);
						if(pairIndex < MAX_FLUID_RIGID_PAIRS) 
						{
							out_pairs[particleIndex].m_rigidIndicies[pairIndex] = rigidBodyIndex;
							out_pairs[particleIndex].m_rigidSubIndicies[pairIndex] = triangleIndex;
						}
					} 
					curIndex++;
				} else
				{
					if (isLeafNode)
					{
						curIndex++;
					} else
					{
						escapeIndex = getEscapeIndex(&rootNode);
						curIndex += escapeIndex;
					}
				}
			}
		}
	}

}
// -----------------------------------------------------------------------------
// Bullet3 - bvhTraversal.cl
// -----------------------------------------------------------------------------

bool computeContactSpherePlane(float4 spherePos, float sphereRadius,
								float4 rigidPos, float4 rigidOrn, float4 rigidPlaneEquation,
								float* out_distance, float4* out_normalOnRigidWorld, float4* out_pointOnRigidWorld)
{
	//Convert sphere position from world space to rigid space
	float4 invRigidPos, invRigidOrn;
	trInverse(rigidPos, rigidOrn, &invRigidPos, &invRigidOrn);

	float4 spherePosInRigidSpace = transform(&spherePos, &invRigidPos, &invRigidOrn);
	
	//
	float4 pointOnPlaneRigidSpace;
	float distance = signedDistanceFromPointToPlane(spherePosInRigidSpace, rigidPlaneEquation, &pointOnPlaneRigidSpace) - sphereRadius;
	
	if(distance < 0.f)
	{
		float4 normalOnPlaneRigidSpace = make_float4(rigidPlaneEquation.xyz, 0.f);
	
		//Convert contact from rigid space to world space
		*out_distance = distance;
		*out_normalOnRigidWorld = qtRotate(rigidOrn, normalOnPlaneRigidSpace); 
		*out_pointOnRigidWorld = transform(&normalOnPlaneRigidSpace, &rigidPos, &rigidOrn);
		
		return true;
	}
	
	return false;
}

bool computeContactSphereSphere(float4 particlePos, float particleRadius,
								float4 rigidPos, float rigidRadius,
								float* out_distance, float4* out_normalOnRigidWorld, float4* out_pointOnRigidWorld)
{
	float4 rigidToParticle = particlePos - rigidPos;
	float distanceBetweenCenters = length(rigidToParticle);
	float distance = distanceBetweenCenters - (particleRadius + rigidRadius);
	
	if(distance < 0.f)
	{
		*out_distance = distance;
		*out_normalOnRigidWorld = (distanceBetweenCenters > B3_EPSILON) ? rigidToParticle / distanceBetweenCenters : (float4)(0.f, 1.f, 0.f, 0.f);
		*out_pointOnRigidWorld = rigidPos + *out_normalOnRigidWorld*rigidRadius;
		
		return true;
	}
	
	return false;
}


__kernel void clearFluidRigidPairsAndContacts(__global FluidRigidPairs* pairs, __global FluidRigidPairs* midphasePairs,
											__global FluidRigidContacts* contacts, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	pairs[i].m_numIndicies = 0;
	midphasePairs[i].m_numIndicies = 0;
	contacts[i].m_numContacts = 0;
}

__kernel void detectLargeAabbRigids(__constant b3FluidSphParameters* FP,  
									__global btAabbCL* rigidBodyWorldAabbs,
									__global BodyData* rigidBodies, __global btCollidableGpu* collidables,
									__global int* out_numLargeAabbRigids, __global int* out_largeAabbRigidIndicies, 
									int maxLargeRigidAabbs, int numRigidBodies)
{
	int i = get_global_id(0);
	if(i >= numRigidBodies) return;
	
	b3Scalar gridCellSize = FP->m_sphSmoothRadius / FP->m_simulationScale;
	b3Scalar particleRadius = FP->m_particleRadius;
	b3Vector3 radiusAabbExtent = (b3Vector3){ particleRadius, particleRadius, particleRadius, 0.0f };
	
	btAabbCL rigidAabb = rigidBodyWorldAabbs[i];
	
	b3Vector3 expandedRigidAabbMin = rigidAabb.m_min - radiusAabbExtent;
	b3Vector3 expandedRigidAabbMax = rigidAabb.m_max + radiusAabbExtent;
	b3Scalar rigidAabbCornerDistance = b3Sqrt( b3Vector3_length2(expandedRigidAabbMax - expandedRigidAabbMin) );
	
	BodyData rigidBody = rigidBodies[i];
	int collidableIndex = rigidBody.m_collidableIdx;
	
	b3Scalar maxAabbExtent = gridCellSize * (b3Scalar)B3_FLUID_GRID_COORD_RANGE_HALVED;
	if( fabs(expandedRigidAabbMin.x) > maxAabbExtent
	 || fabs(expandedRigidAabbMin.y) > maxAabbExtent
	 || fabs(expandedRigidAabbMin.z) > maxAabbExtent
	 || fabs(expandedRigidAabbMax.x) > maxAabbExtent
	 || fabs(expandedRigidAabbMax.y) > maxAabbExtent
	 || fabs(expandedRigidAabbMax.z) > maxAabbExtent
	 || rigidAabbCornerDistance > maxAabbExtent * 0.05f
	 || collidables[collidableIndex].m_shapeType == SHAPE_CONCAVE_TRIMESH )
	{
		int largeRigidIndex = atomic_inc(out_numLargeAabbRigids);	//out_numLargeAabbRigids is assumed to be 0 before the kernel is executed
		if(largeRigidIndex < maxLargeRigidAabbs) out_largeAabbRigidIndicies[largeRigidIndex] = i;
	}
}
__kernel void fluidLargeRigidBroadphase(__constant b3FluidSphParameters* FP, 
										__global b3Vector3* fluidPosition, __global btAabbCL* rigidBodyWorldAabbs,
										__global BodyData* rigidBodies, __global btCollidableGpu* collidables,
										__global int* numLargeAabbRigids, __global int* largeAabbRigidIndicies,
										__global FluidRigidPairs* out_pairs, __global FluidRigidPairs* out_midphasePairs, 
										int maxLargeRigidAabbs, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar gridCellSize = FP->m_sphSmoothRadius / FP->m_simulationScale;
	b3Scalar particleRadius = FP->m_particleRadius;
	b3Vector3 radiusAabbExtent = (b3Vector3){ particleRadius, particleRadius, particleRadius, 0.0f };
	
	b3Vector3 particlePos = fluidPosition[i];
	
	int numValidLargeAabbRigids = min(*numLargeAabbRigids, maxLargeRigidAabbs);
	if( get_global_id(0) == 0 ) *numLargeAabbRigids = numValidLargeAabbRigids;
	
	for(int n = 0; n < numValidLargeAabbRigids; ++n)
	{
		int rigidIndex = largeAabbRigidIndicies[n];
	
		BodyData rigidBody = rigidBodies[rigidIndex];
		int collidableIndex = rigidBody.m_collidableIdx;
	
		bool needsMidphase = ( collidables[collidableIndex].m_shapeType == SHAPE_CONCAVE_TRIMESH 
							|| collidables[collidableIndex].m_shapeType == SHAPE_COMPOUND_OF_CONVEX_HULLS );
						
		__global FluidRigidPairs* output = (needsMidphase) ? out_midphasePairs : out_pairs;
	
		btAabbCL rigidAabb = rigidBodyWorldAabbs[rigidIndex];
		b3Vector3 expandedRigidAabbMin = rigidAabb.m_min - radiusAabbExtent;
		b3Vector3 expandedRigidAabbMax = rigidAabb.m_max + radiusAabbExtent;
		
		if( expandedRigidAabbMin.x <= particlePos.x && particlePos.x <= expandedRigidAabbMax.x
		 && expandedRigidAabbMin.y <= particlePos.y && particlePos.y <= expandedRigidAabbMax.y
		 && expandedRigidAabbMin.z <= particlePos.z && particlePos.z <= expandedRigidAabbMax.z )
		{
			int pairIndex = output[i].m_numIndicies++;
			if(pairIndex < MAX_FLUID_RIGID_PAIRS) output[i].m_rigidIndicies[pairIndex] = rigidIndex;
			else return;
		}
	}
}

__kernel void fluidSmallRigidBroadphase(__constant b3FluidSphParameters* FP, 
									__global b3Vector3* fluidPosition, __global b3FluidGridCombinedPos* cellValues, 
									__global b3FluidGridIterator* cellContents, __global btAabbCL* rigidBodyWorldAabbs,
									__global BodyData* rigidBodies, __global btCollidableGpu* collidables,
									__global FluidRigidPairs* out_pairs, __global FluidRigidPairs* out_midphasePairs,			
									int numGridCells, int numRigidBodies)
{
	int i = get_global_id(0);
	if(i >= numRigidBodies) return;
	
	b3Scalar gridCellSize = FP->m_sphSmoothRadius / FP->m_simulationScale;
	b3Scalar particleRadius = FP->m_particleRadius;
	b3Vector3 radiusAabbExtent = (b3Vector3){ particleRadius, particleRadius, particleRadius, 0.0f };
	
	btAabbCL rigidAabb = rigidBodyWorldAabbs[i];
	//rigidAabb.m_min.w = 0.0f;	//	check if necessary(if using vector functions for point-AABB test)
	//rigidAabb.m_max.w = 0.0f;
	
	b3Vector3 expandedRigidAabbMin = rigidAabb.m_min - radiusAabbExtent;
	b3Vector3 expandedRigidAabbMax = rigidAabb.m_max + radiusAabbExtent;
	b3Scalar rigidAabbCornerDistance = b3Sqrt( b3Vector3_length2(expandedRigidAabbMax - expandedRigidAabbMin) );
	
	BodyData rigidBody = rigidBodies[i];
	int collidableIndex = rigidBody.m_collidableIdx;
	
	bool needsMidphase = ( collidables[collidableIndex].m_shapeType == SHAPE_CONCAVE_TRIMESH 
						|| collidables[collidableIndex].m_shapeType == SHAPE_COMPOUND_OF_CONVEX_HULLS );
	__global FluidRigidPairs* output = (needsMidphase) ? out_midphasePairs : out_pairs;
	
	b3Scalar maxAabbExtent = gridCellSize * (b3Scalar)B3_FLUID_GRID_COORD_RANGE_HALVED;
	if( fabs(expandedRigidAabbMin.x) > maxAabbExtent
	 || fabs(expandedRigidAabbMin.y) > maxAabbExtent
	 || fabs(expandedRigidAabbMin.z) > maxAabbExtent
	 || fabs(expandedRigidAabbMax.x) > maxAabbExtent
	 || fabs(expandedRigidAabbMax.y) > maxAabbExtent
	 || fabs(expandedRigidAabbMax.z) > maxAabbExtent 
	 || rigidAabbCornerDistance > maxAabbExtent * 0.05f
	 || collidables[collidableIndex].m_shapeType == SHAPE_CONCAVE_TRIMESH ) return;
						
	b3FluidGridPosition quantizedAabbMin = getDiscretePosition( gridCellSize, expandedRigidAabbMin );
	b3FluidGridPosition quantizedAabbMax = getDiscretePosition( gridCellSize, expandedRigidAabbMax );
						
	for(int z = quantizedAabbMin.z; z <= quantizedAabbMax.z; ++z)
		for(int y = quantizedAabbMin.y; y <= quantizedAabbMax.y; ++y)
			for(int x = quantizedAabbMin.x; x <= quantizedAabbMax.x; ++x)
			{
				b3FluidGridPosition currentCell;
				currentCell.x = x;
				currentCell.y = y;
				currentCell.z = z;
				
				b3FluidGridCombinedPos currentCellValue = getCombinedPosition(currentCell);
			
				int cellIndex = binarySearch(cellValues, numGridCells, currentCellValue);
				if(cellIndex != numGridCells)
				{
					b3FluidGridIterator fluidCell = cellContents[cellIndex]; 
					for(int particleIndex = fluidCell.m_firstIndex; particleIndex <= fluidCell.m_lastIndex; ++particleIndex)
					{
						b3Vector3 particlePos = fluidPosition[particleIndex];
						//particlePos.w = 0.0f;	//	check if necessary(if using vector functions for point-AABB test)
						
						if( expandedRigidAabbMin.x <= particlePos.x && particlePos.x <= expandedRigidAabbMax.x
						 && expandedRigidAabbMin.y <= particlePos.y && particlePos.y <= expandedRigidAabbMax.y
						 && expandedRigidAabbMin.z <= particlePos.z && particlePos.z <= expandedRigidAabbMax.z )
						{
							int pairIndex = atomic_inc(&output[particleIndex].m_numIndicies);
							if(pairIndex < MAX_FLUID_RIGID_PAIRS) output[particleIndex].m_rigidIndicies[pairIndex] = i;
						}
					}
				}
			
				
			}
}

// ////////////////////////////////////////////////////////////////////////////
// Modulo Hash Grid Only
// ////////////////////////////////////////////////////////////////////////////
#define B3_FLUID_HASH_GRID_COORD_RANGE 64
b3FluidGridCombinedPos getCombinedPositionModulo(b3FluidGridPosition quantizedPosition)
{
	//as_uint() requires that sizeof(b3FluidGridCombinedPos) == sizeof(b3FluidGridCoordinate)
	//This presents an issue if B3_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT is #defined
	b3FluidGridCombinedPos unsignedX = as_uint(quantizedPosition.x) % B3_FLUID_HASH_GRID_COORD_RANGE;
	b3FluidGridCombinedPos unsignedY = as_uint(quantizedPosition.y) % B3_FLUID_HASH_GRID_COORD_RANGE;
	b3FluidGridCombinedPos unsignedZ = as_uint(quantizedPosition.z) % B3_FLUID_HASH_GRID_COORD_RANGE;
	
	return unsignedX 
		+ unsignedY * B3_FLUID_HASH_GRID_COORD_RANGE
		+ unsignedZ * B3_FLUID_HASH_GRID_COORD_RANGE* B3_FLUID_HASH_GRID_COORD_RANGE;
}
__kernel void fluidSmallRigidBroadphaseModulo(__constant b3FluidSphParameters* FP, 
									__global b3Vector3* fluidPosition, 
									__global b3FluidGridIterator* cellContents, __global btAabbCL* rigidBodyWorldAabbs,
									__global BodyData* rigidBodies, __global btCollidableGpu* collidables,
									__global FluidRigidPairs* out_pairs, __global FluidRigidPairs* out_midphasePairs,			
									int numRigidBodies)
{
	int i = get_global_id(0);
	if(i >= numRigidBodies) return;
	
	b3Scalar gridCellSize = FP->m_sphSmoothRadius / FP->m_simulationScale;
	b3Scalar particleRadius = FP->m_particleRadius;
	b3Vector3 radiusAabbExtent = (b3Vector3){ particleRadius, particleRadius, particleRadius, 0.0f };
	
	btAabbCL rigidAabb = rigidBodyWorldAabbs[i];
	//rigidAabb.m_min.w = 0.0f;	//	check if necessary(if using vector functions for point-AABB test)
	//rigidAabb.m_max.w = 0.0f;
	
	b3Vector3 expandedRigidAabbMin = rigidAabb.m_min - radiusAabbExtent;
	b3Vector3 expandedRigidAabbMax = rigidAabb.m_max + radiusAabbExtent;
	b3Scalar rigidAabbCornerDistance = b3Sqrt( b3Vector3_length2(expandedRigidAabbMax - expandedRigidAabbMin) );
	
	BodyData rigidBody = rigidBodies[i];
	int collidableIndex = rigidBody.m_collidableIdx;
	
	bool needsMidphase = ( collidables[collidableIndex].m_shapeType == SHAPE_CONCAVE_TRIMESH 
						|| collidables[collidableIndex].m_shapeType == SHAPE_COMPOUND_OF_CONVEX_HULLS );
						
	__global FluidRigidPairs* output = (needsMidphase) ? out_midphasePairs : out_pairs;
	
	//	may need to use different criteria for modulo grid(use of B3_FLUID_GRID_COORD_RANGE_HALVED is arbitrary)
	b3Scalar maxAabbExtent = gridCellSize * (b3Scalar)B3_FLUID_GRID_COORD_RANGE_HALVED;	
	if( fabs(expandedRigidAabbMin.x) > maxAabbExtent
	 || fabs(expandedRigidAabbMin.y) > maxAabbExtent
	 || fabs(expandedRigidAabbMin.z) > maxAabbExtent
	 || fabs(expandedRigidAabbMax.x) > maxAabbExtent
	 || fabs(expandedRigidAabbMax.y) > maxAabbExtent
	 || fabs(expandedRigidAabbMax.z) > maxAabbExtent 
	 || rigidAabbCornerDistance > maxAabbExtent * 0.05f
	 || collidables[collidableIndex].m_shapeType == SHAPE_CONCAVE_TRIMESH ) return;
						
	b3FluidGridPosition quantizedAabbMin = getDiscretePosition( gridCellSize, expandedRigidAabbMin );
	b3FluidGridPosition quantizedAabbMax = getDiscretePosition( gridCellSize, expandedRigidAabbMax );
	
	for(int z = quantizedAabbMin.z; z <= quantizedAabbMax.z; ++z)
		for(int y = quantizedAabbMin.y; y <= quantizedAabbMax.y; ++y)
			for(int x = quantizedAabbMin.x; x <= quantizedAabbMax.x; ++x)
			{
				b3FluidGridPosition currentCell;
				currentCell.x = x;
				currentCell.y = y;
				currentCell.z = z;
				
				b3FluidGridCombinedPos gridCellIndex = getCombinedPositionModulo(currentCell);
				
				{
					b3FluidGridIterator fluidCell = cellContents[gridCellIndex]; 
					for(int particleIndex = fluidCell.m_firstIndex; particleIndex <= fluidCell.m_lastIndex; ++particleIndex)
					{
						b3Vector3 particlePos = fluidPosition[particleIndex];
						//particlePos.w = 0.0f;	//	check if necessary(if using vector functions for point-AABB test)
						
						if( expandedRigidAabbMin.x <= particlePos.x && particlePos.x <= expandedRigidAabbMax.x
						 && expandedRigidAabbMin.y <= particlePos.y && particlePos.y <= expandedRigidAabbMax.y
						 && expandedRigidAabbMin.z <= particlePos.z && particlePos.z <= expandedRigidAabbMax.z )
						{
							int pairIndex = atomic_inc(&output[particleIndex].m_numIndicies);
							if(pairIndex < MAX_FLUID_RIGID_PAIRS) output[particleIndex].m_rigidIndicies[pairIndex] = i;
						}
					}
				}
				
			}
	
}
// ////////////////////////////////////////////////////////////////////////////
// Modulo Hash Grid Only
// ////////////////////////////////////////////////////////////////////////////

__kernel void fluidRigidMidphase(__constant b3FluidSphParameters* FP, __global b3Vector3* fluidPosition, 
								__global BodyData* rigidBodies, __global btCollidableGpu* collidables,
								__global b3BvhInfo* bvhInfos, __global btBvhSubtreeInfo* bvhSubtreeInfo, __global btQuantizedBvhNode* bvhNodes,
								__global FluidRigidPairs* midphasePairs, __global FluidRigidPairs* out_pairs, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	b3Scalar particleRadius = FP->m_particleRadius;
	b3Vector3 radiusAabbExtent = (b3Vector3){ particleRadius, particleRadius, particleRadius, 0.0f };
	
	b3Vector3 fluidAabbMin = fluidPosition[i] - radiusAabbExtent;
	b3Vector3 fluidAabbmax = fluidPosition[i] + radiusAabbExtent;
	
	//Convert each entry in midphasePairs into several entries in out_pairs
	int numIndicies = min(midphasePairs[i].m_numIndicies, MAX_FLUID_RIGID_PAIRS);
	midphasePairs[i].m_numIndicies = numIndicies;
	
	for(int n = 0; n < numIndicies; ++n)
	{
		int rigidIndex = midphasePairs[i].m_rigidIndicies[n];
		BodyData rigidBody = rigidBodies[rigidIndex];
		int collidableIndex = rigidBody.m_collidableIdx;
		
		if(collidables[collidableIndex].m_shapeType == SHAPE_CONCAVE_TRIMESH)
		{
			//Traverse the triangle mesh BVH with particle AABB, adding 1 element to out_pairs
			//for each leaf(triangle) that intersects with the particle AABB
			 getIntersectingBvhLeaves(collidables, bvhSubtreeInfo, bvhNodes, bvhInfos, 
										rigidIndex, collidableIndex,
										i, fluidAabbMin, fluidAabbmax, out_pairs);
		}
		else if(collidables[collidableIndex].m_shapeType == SHAPE_COMPOUND_OF_CONVEX_HULLS)
		{
			for(int childShape = 0; childShape < collidables[collidableIndex].m_numChildShapes; ++childShape)
			{
				int childShapeIndex = collidables[collidableIndex].m_shapeIndex + childShape;
				
				int pairIndex = atomic_inc(&out_pairs[i].m_numIndicies);
				if(pairIndex < MAX_FLUID_RIGID_PAIRS) 
				{
					out_pairs[i].m_rigidIndicies[pairIndex] = rigidIndex;
					out_pairs[i].m_rigidSubIndicies[pairIndex] = childShapeIndex;
				}
			}
		}
	}
}

__kernel void fluidRigidNarrowphase(__constant b3FluidSphParameters* FP, 
									__global b3Vector3* fluidPosition, __global FluidRigidPairs* pairs, 
									__global BodyData* rigidBodies, __global btCollidableGpu* collidables,
									__global ConvexPolyhedronCL* convexShapes, __global btGpuFace* faces,
									__global int* convexIndices, __global float4* convexVertices, 
									__global btGpuChildShape* gpuChildShapes, __global FluidRigidContacts* out_contact,
									int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	int numFluidRigidPairs = min(pairs[i].m_numIndicies, MAX_FLUID_RIGID_PAIRS);
	pairs[i].m_numIndicies = numFluidRigidPairs;
	
	FluidRigidPairs currentPairs = pairs[i];
	for(int pair = 0; pair < numFluidRigidPairs; ++pair)
	{
		int rigidIndex = currentPairs.m_rigidIndicies[pair];
	
		BodyData rigidBody = rigidBodies[rigidIndex];
		
		int collidableIndex = rigidBody.m_collidableIdx;
		
		bool isColliding = false;
		float distance;
		float4 normalOnRigid;
		float4 pointOnRigid;
		
		int shapeType = collidables[collidableIndex].m_shapeType;
		int shapeIndex = collidables[collidableIndex].m_shapeIndex;
		float4 rigidPosition = rigidBody.m_pos;
		float4 rigidOrientation = rigidBody.m_quat;
		
		if(shapeType == SHAPE_COMPOUND_OF_CONVEX_HULLS)
		{
			int childShapeIndex = currentPairs.m_rigidSubIndicies[pair];
			int childShapeCollidableIndex = gpuChildShapes[childShapeIndex].m_shapeIndex;
			float4 childWorldPosition = qtRotate(rigidOrientation, gpuChildShapes[childShapeIndex].m_childPosition) + rigidPosition;
			float4 childWorldOrientation = qtMul(rigidOrientation, gpuChildShapes[childShapeIndex].m_childOrientation);
			
			//Replace shape, position, orientation with that of the child shape
			shapeType = collidables[childShapeCollidableIndex].m_shapeType;
			shapeIndex = collidables[childShapeCollidableIndex].m_shapeIndex;
			rigidPosition = childWorldPosition;
			rigidOrientation = childWorldOrientation;
		}
		
		switch(shapeType)
		{
			case SHAPE_CONVEX_HULL:
			{
				isColliding = computeContactSphereConvex( collidableIndex, collidables, convexShapes, convexVertices, convexIndices, faces,
														fluidPosition[i], FP->m_particleRadius, rigidPosition, rigidOrientation,
														&distance, &normalOnRigid, &pointOnRigid );
				normalOnRigid = -normalOnRigid;		//	computeContactSphereConvex() actually returns normal on particle?
			}
				break;
				
			case SHAPE_PLANE:
			{
				float4 rigidPlaneEquation = faces[shapeIndex].m_plane;
		
				isColliding = computeContactSpherePlane(fluidPosition[i], FP->m_particleRadius, 
														rigidPosition, rigidOrientation, rigidPlaneEquation,
														&distance, &normalOnRigid, &pointOnRigid );
			}
				break;
			
			case SHAPE_SPHERE:
			{
				float rigidSphereRadius = collidables[collidableIndex].m_radius;
			
				isColliding = computeContactSphereSphere(fluidPosition[i], FP->m_particleRadius,
														rigidPosition, rigidSphereRadius,
														&distance, &normalOnRigid, &pointOnRigid);
			}
				break;
				
			case SHAPE_CONCAVE_TRIMESH:
			{
				int rigidSubIndex = currentPairs.m_rigidSubIndicies[pair];
			
				float4 triangleVertices[3];
				
				btGpuFace face = faces[ convexShapes[shapeIndex].m_faceOffset + rigidSubIndex ];
		
				for (int j=0;j<3;j++)
				{
					int index = convexIndices[face.m_indexOffset + j];
					float4 vertex = convexVertices[ convexShapes[shapeIndex].m_vertexOffset + index];
					triangleVertices[j] = vertex;
				}
				
				isColliding = computeContactSphereTriangle(triangleVertices, fluidPosition[i], FP->m_particleRadius,
															rigidPosition, rigidOrientation,
															&distance, &normalOnRigid, &pointOnRigid);									
				normalOnRigid = -normalOnRigid;
			}
				break;
			
			default:
				continue;
		}
		
		
		if(isColliding)
		{
			int contactIndex = out_contact[i].m_numContacts;
			
			out_contact[i].m_rigidIndicies[contactIndex] = rigidIndex;
			out_contact[i].m_distances[contactIndex] = distance;
			out_contact[i].m_pointsOnRigid[contactIndex] = pointOnRigid;
			out_contact[i].m_normalsOnRigid[contactIndex] = normalOnRigid;		
			
			++out_contact[i].m_numContacts;
		}
		
		if(out_contact[i].m_numContacts >= MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE) return;
	}
}

__kernel void resolveFluidRigidCollisions(__constant b3FluidSphParameters* FP, 
											__global BodyData* rigidBodies, __global InertiaTensor* rigidInertias,
											__global FluidRigidContacts* contacts, 
											__global b3Vector3* fluidVel, __global b3Vector3* fluidVelEval, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	for(int contactIndex = 0; contactIndex < contacts[i].m_numContacts; ++contactIndex)
	{
		b3Vector3 fluidVelocity = fluidVel[i];
		FluidRigidContacts contact = contacts[i];
	
		int rigidIndex = contact.m_rigidIndicies[contactIndex];
		b3Scalar distance = contact.m_distances[contactIndex];
		b3Vector3 normalOnRigid = contact.m_normalsOnRigid[contactIndex];
		b3Vector3 pointOnRigid = contact.m_pointsOnRigid[contactIndex];
		
		b3Vector3 rigidPosition = rigidBodies[rigidIndex].m_pos;
		b3Vector3 rigidLinearVelocity = rigidBodies[rigidIndex].m_linVel;
		b3Vector3 rigidAngularVelocity = rigidBodies[rigidIndex].m_angVel;
		b3Scalar rigidInvMass = rigidBodies[rigidIndex].m_invMass;
		Matrix3x3 rigidInertiaTensor = rigidInertias[rigidIndex].m_invInertia;
		
		if( distance < 0.0f )
		{
			bool isDynamicRigidBody = (rigidInvMass != 0.0f);
			
			b3Vector3 rigidLocalHitPoint = pointOnRigid - rigidPosition;
			
			b3Vector3 rigidVelocity = (b3Vector3){0.0f, 0.0f, 0.0f, 0.0f};
			if(isDynamicRigidBody) rigidVelocity = rigidLinearVelocity + cross3(rigidAngularVelocity, rigidLocalHitPoint);
			rigidVelocity *= FP->m_simulationScale;
		
			b3Vector3 relativeVelocity = fluidVelocity - rigidVelocity;
			b3Scalar penetratingMagnitude = b3Vector3_dot(relativeVelocity, -normalOnRigid);
			if( penetratingMagnitude < 0.0f ) penetratingMagnitude = 0.0f;
			
			b3Vector3 penetratingVelocity = -normalOnRigid * penetratingMagnitude;
			b3Vector3 tangentialVelocity = relativeVelocity - penetratingVelocity;
			
			penetratingVelocity *= 1.0f + FP->m_boundaryRestitution;
			
			b3Scalar penetration = -distance;
			penetration = (penetration > FP->m_particleMargin) ? penetration : 0.0f;
			b3Scalar positionError = penetration * (FP->m_simulationScale/FP->m_timeStep) * FP->m_boundaryErp;
			
			b3Vector3 particleImpulse = -(penetratingVelocity + (-normalOnRigid*positionError) + tangentialVelocity*FP->m_boundaryFriction);
			
			if(isDynamicRigidBody)
			{
				b3Scalar inertiaParticle = 1.0f / FP->m_particleMass;
				
				b3Vector3 relPosCrossNormal = cross3(rigidLocalHitPoint, normalOnRigid);
				b3Scalar inertiaRigid = rigidInvMass + b3Vector3_dot( mtMul3(relPosCrossNormal, rigidInertiaTensor), relPosCrossNormal );
				
				particleImpulse *= 1.0f / (inertiaParticle + inertiaRigid);
				
				//b3Vector3 worldScaleImpulse = -particleImpulse / FP->m_simulationScale;
				//worldScaleImpulse /= FP->m_timeStep;		//Impulse is accumulated as force
				
				//const b3Vector3& linearFactor = rigidBody->getLinearFactor();
				//accumulatedRigidForce += worldScaleImpulse * linearFactor;
				//accumulatedRigidTorque += rigidLocalHitPoint.cross(worldScaleImpulse * linearFactor) * rigidBody->getAngularFactor();
				
				particleImpulse *= inertiaParticle;
			}
			
			//Leapfrog integration
			b3Vector3 velNext = fluidVelocity + particleImpulse;
			fluidVel[i] = velNext;
		}
	}
}


__kernel void clearRigidFluidContacts(__global RigidFluidContacts* out_rigidFluidContacts, int numRigidBodies)
{
	int rigidIndex = get_global_id(0);
	if(rigidIndex >= numRigidBodies) return;
	
	out_rigidFluidContacts[rigidIndex].m_numContacts = 0;
}

__kernel void mapRigidFluidContacts(__global FluidRigidContacts* fluidRigidContacts,
									__global RigidFluidContacts* out_rigidFluidContacts, int numFluidParticles)
{
	int particleIndex = get_global_id(0);
	if(particleIndex >= numFluidParticles) return;
	
	__global FluidRigidContacts* contacts = &fluidRigidContacts[particleIndex];
	for(int n = 0; n < contacts->m_numContacts; ++n)
	{
		int rigidIndex = contacts->m_rigidIndicies[n];
		
		int pairIndex = atomic_inc(&out_rigidFluidContacts[rigidIndex].m_numContacts);
		if(pairIndex < MAX_FLUID_CONTACTS_PER_DYNAMIC_RIGID) 
		{
			out_rigidFluidContacts[rigidIndex].m_fluidIndicies[pairIndex] = particleIndex;
			out_rigidFluidContacts[rigidIndex].m_contactIndicies[pairIndex] = n;
		}
	}
}

__kernel void resolveRigidFluidCollisions(__constant b3FluidSphParameters* FP, 
											__global BodyData* rigidBodies, __global InertiaTensor* rigidInertias,
											__global FluidRigidContacts* fluidContacts, __global RigidFluidContacts* rigidContacts, 
											__global b3Vector3* fluidVel, int numRigidBodies)
{
	int rigidIndex = get_global_id(0);
	if(rigidIndex >= numRigidBodies) return;
	
	b3Vector3 rigidPosition = rigidBodies[rigidIndex].m_pos;
	b3Vector3 rigidLinearVelocity = rigidBodies[rigidIndex].m_linVel;
	b3Vector3 rigidAngularVelocity = rigidBodies[rigidIndex].m_angVel;
	b3Scalar rigidInvMass = rigidBodies[rigidIndex].m_invMass;
	Matrix3x3 rigidInertiaTensor = rigidInertias[rigidIndex].m_invInertia;
	
	if(rigidInvMass == 0.0f) return;
	
	b3Vector3 accumulatedForce = (b3Vector3){0.0f, 0.0f, 0.0f, 0.0f};
	b3Vector3 accumulatedTorque = (b3Vector3){0.0f, 0.0f, 0.0f, 0.0f};
	
	int numContacts = min(rigidContacts[rigidIndex].m_numContacts, MAX_FLUID_CONTACTS_PER_DYNAMIC_RIGID);
	rigidContacts[rigidIndex].m_numContacts = numContacts;
	
	for(int i = 0; i < numContacts; ++i)
	{
		int particleIndex = rigidContacts[rigidIndex].m_fluidIndicies[i];
		int contactIndex = rigidContacts[rigidIndex].m_contactIndicies[i];
		
		b3Vector3 fluidVelocity = fluidVel[particleIndex];
		
		b3Scalar distance = fluidContacts[particleIndex].m_distances[contactIndex];
		b3Vector3 pointOnRigid = fluidContacts[particleIndex].m_pointsOnRigid[contactIndex];
		b3Vector3 normalOnRigid = fluidContacts[particleIndex].m_normalsOnRigid[contactIndex];
		
		if( distance < 0.0f )
		{
			bool isDynamicRigidBody = (rigidInvMass != 0.0f);
			
			b3Vector3 rigidLocalHitPoint = pointOnRigid - rigidPosition;
			
			b3Vector3 rigidVelocity = (b3Vector3){0.0f, 0.0f, 0.0f, 0.0f};
			if(isDynamicRigidBody) rigidVelocity = rigidLinearVelocity + cross3(rigidAngularVelocity, rigidLocalHitPoint);
			rigidVelocity *= FP->m_simulationScale;
		
			b3Vector3 relativeVelocity = fluidVelocity - rigidVelocity;
			b3Scalar penetratingMagnitude = b3Vector3_dot(relativeVelocity, -normalOnRigid);
			if( penetratingMagnitude < 0.0f ) penetratingMagnitude = 0.0f;
			
			b3Vector3 penetratingVelocity = -normalOnRigid * penetratingMagnitude;
			b3Vector3 tangentialVelocity = relativeVelocity - penetratingVelocity;
			
			penetratingVelocity *= 1.0f + FP->m_boundaryRestitution;
			
			b3Scalar penetration = -distance;
			penetration = (penetration > FP->m_particleMargin) ? penetration : 0.0f;
			b3Scalar positionError = penetration * (FP->m_simulationScale/FP->m_timeStep) * FP->m_boundaryErp;
			
			b3Vector3 particleImpulse = -(penetratingVelocity + (-normalOnRigid*positionError) + tangentialVelocity*FP->m_boundaryFriction);
			
			if(isDynamicRigidBody)
			{
				b3Scalar inertiaParticle = 1.0f / FP->m_particleMass;
				
				b3Vector3 relPosCrossNormal = cross3(rigidLocalHitPoint, normalOnRigid);
				b3Scalar inertiaRigid = rigidInvMass + b3Vector3_dot( mtMul3(relPosCrossNormal, rigidInertiaTensor), relPosCrossNormal );
				
				particleImpulse *= 1.0f / (inertiaParticle + inertiaRigid);
				
				b3Vector3 worldScaleImpulse = -particleImpulse / FP->m_simulationScale;
				worldScaleImpulse /= FP->m_timeStep;		//Impulse is accumulated as force
				
				accumulatedForce += worldScaleImpulse;
				accumulatedTorque += cross3(rigidLocalHitPoint, worldScaleImpulse);
				
				//particleImpulse *= inertiaParticle;
			}
		}	
	}
	
	//Apply accumulated forces
	b3Scalar timeStep = FP->m_timeStep;
		
	rigidLinearVelocity += accumulatedForce * (rigidInvMass * timeStep);
	rigidAngularVelocity += mtMul1(rigidInertiaTensor, accumulatedTorque) * timeStep;
	
	//Limit angular velocity
	float BT_GPU_ANGULAR_MOTION_THRESHOLD = (0.25f * 3.14159254f);
	b3Scalar angVel = sqrt( b3Vector3_dot(rigidAngularVelocity, rigidAngularVelocity) );
	if(angVel*timeStep > BT_GPU_ANGULAR_MOTION_THRESHOLD) rigidAngularVelocity *= (BT_GPU_ANGULAR_MOTION_THRESHOLD/timeStep) / angVel;
	
	rigidBodies[rigidIndex].m_linVel = rigidLinearVelocity;
	rigidBodies[rigidIndex].m_angVel = rigidAngularVelocity;
}
