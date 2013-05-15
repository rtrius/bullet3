
#ifdef cl_amd_printf
	#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

typedef float b3Scalar;
typedef float4 b3Vector3;
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


//Syncronize with 'struct b3FluidSphParametersGlobal' in b3FluidSphParameters.h
typedef struct
{
	b3Scalar m_timeStep;
	b3Scalar m_simulationScale;
	b3Scalar m_speedLimit;
	b3Scalar m_sphSmoothRadius;
	b3Scalar m_sphRadiusSquared;
	b3Scalar m_poly6KernCoeff;
	b3Scalar m_spikyKernGradCoeff;
	b3Scalar m_viscosityKernLapCoeff;
	b3Scalar m_initialSum;
} b3FluidSphParametersGlobal;

//Syncronize with 'struct b3FluidSphParametersLocal' in b3FluidSphParameters.h
typedef struct
{
	b3Vector3 m_aabbBoundaryMin;
	b3Vector3 m_aabbBoundaryMax;
	int m_enableAabbBoundary;
	b3Vector3 m_gravity;
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
} b3FluidSphParametersLocal;


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

///-----------------------------------------------------------------------------
///Bullet3
///-----------------------------------------------------------------------------
#define SHAPE_CONVEX_HULL 3
#define SHAPE_PLANE 4
#define SHAPE_CONCAVE_TRIMESH 5
#define SHAPE_COMPOUND_OF_CONVEX_HULLS 6
#define SHAPE_SPHERE 7

typedef unsigned int u32;

///keep this in sync with btCollidable.h
typedef struct
{
	int m_numChildShapes;
	int blaat2;
	int m_shapeType;
	int m_shapeIndex;
	
} btCollidableGpu;

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
} BodyData;		//b3RigidBodyCL in C++ (Bullet3Collision/NarrowPhaseCollision/b3RigidBodyCL.h)


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
///-----------------------------------------------------------------------------
///Bullet3
///-----------------------------------------------------------------------------




#define MAX_FLUID_RIGID_PAIRS 32
#define MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE 4

///Contains the indicies of rigid bodies whose AABB intersects with that of a single fluid particle
typedef struct
{
	int m_numIndicies;
	int m_rigidIndicies[MAX_FLUID_RIGID_PAIRS];
	
} FluidRigidPairs;

typedef struct
{
	int m_numContacts;
	int m_rigidIndicies[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];
	b3Scalar m_distances[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];
	b3Vector3 m_pointsOnRigid[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];		//World space point
	b3Vector3 m_normalsOnRigid[MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE];		//World space normal
	
} FluidRigidContacts;

__kernel void clearFluidRigidPairsAndContacts(__global FluidRigidPairs* pairs, __global FluidRigidContacts* contacts, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	pairs[i].m_numIndicies = 0;
	contacts[i].m_numContacts = 0;
}

__kernel void fluidRigidBroadphase(__constant b3FluidSphParametersGlobal* FG,  __constant b3FluidSphParametersLocal* FL, 
									__global b3Vector3* fluidPosition, __global b3FluidGridCombinedPos* cellValues, 
									__global b3FluidGridIterator* cellContents, __global btAabbCL* rigidBodyWorldAabbs,
									__global FluidRigidPairs* out_pairs,									
									int numGridCells, int numRigidBodies)
{
	int i = get_global_id(0);
	if(i >= numRigidBodies) return;
	
	b3Scalar gridCellSize = FG->m_sphSmoothRadius / FG->m_simulationScale;
	b3Scalar particleRadius = FL->m_particleRadius;
	b3Vector3 radiusAabbExtent = (b3Vector3){ particleRadius, particleRadius, particleRadius, 0.0f };
	
	btAabbCL rigidAabb = rigidBodyWorldAabbs[i];
	//rigidAabb.m_min.w = 0.0f;	//	check if necessary(if using vector functions for point-AABB test)
	//rigidAabb.m_max.w = 0.0f;
	
	
	b3Scalar maxAabbExtent = gridCellSize * (b3Scalar)B3_FLUID_GRID_COORD_RANGE_HALVED;
	if( fabs(rigidAabb.m_min.x) > maxAabbExtent
	 || fabs(rigidAabb.m_min.y) > maxAabbExtent
	 || fabs(rigidAabb.m_min.z) > maxAabbExtent
	 || fabs(rigidAabb.m_max.x) > maxAabbExtent
	 || fabs(rigidAabb.m_max.y) > maxAabbExtent
	 || fabs(rigidAabb.m_max.z) > maxAabbExtent ) return;
	 
	b3Vector3 expandedRigidAabbMin = rigidAabb.m_min - radiusAabbExtent;
	b3Vector3 expandedRigidAabbMax = rigidAabb.m_max + radiusAabbExtent;
	
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
							int pairIndex = atomic_inc(&out_pairs[particleIndex].m_numIndicies);
							if(pairIndex < MAX_FLUID_RIGID_PAIRS) out_pairs[particleIndex].m_rigidIndicies[pairIndex] = i;
						}
					}
				}
				
			}
	
}

__kernel void fluidRigidNarrowphase(__constant b3FluidSphParametersGlobal* FG, __constant b3FluidSphParametersLocal* FL, 
									__global b3Vector3* fluidPosition, __global FluidRigidPairs* pairs, 
									__global BodyData* rigidBodies, __global btCollidableGpu* collidables,
									__global ConvexPolyhedronCL* convexShapes, __global btGpuFace* faces,
									__global int* convexIndices, __global float4* convexVertices, 
									__global FluidRigidContacts* out_contact,
									int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	FluidRigidPairs currentPairs = pairs[i];
	for(int numRigids = 0; numRigids < currentPairs.m_numIndicies; ++numRigids)
	{
		int rigidIndex = currentPairs.m_rigidIndicies[numRigids];
	
		BodyData rigidBody = rigidBodies[rigidIndex];
		
		int collidableIndex = rigidBody.m_collidableIdx;
		if(collidables[collidableIndex].m_shapeType != SHAPE_CONVEX_HULL) continue;
		
		float distance;
		float4 normalOnRigid;
		float4 pointOnRigid;
		
		bool isColliding = computeContactSphereConvex( collidableIndex,  collidables, convexShapes, convexVertices, convexIndices, faces,
												fluidPosition[i], FL->m_particleRadius, rigidBody.m_pos, rigidBody.m_quat,
												&distance, &normalOnRigid, &pointOnRigid );
										
		if(isColliding)
		{
			int contactIndex = out_contact[i].m_numContacts;
			
			out_contact[i].m_rigidIndicies[contactIndex] = rigidIndex;
			out_contact[i].m_distances[contactIndex] = distance;
			out_contact[i].m_pointsOnRigid[contactIndex] = pointOnRigid;
			out_contact[i].m_normalsOnRigid[contactIndex] = -normalOnRigid;		//computeContactSphereConvex() actually returns normal on particle
			
			++out_contact[i].m_numContacts;
		}
		
		if(out_contact[i].m_numContacts >= MAX_RIGID_CONTACTS_PER_FLUID_PARTICLE) return;
	}
}

__kernel void resolveFluidRigidCollisions(__constant b3FluidSphParametersGlobal* FG, __constant b3FluidSphParametersLocal* FL, 
											__global BodyData* rigidBodies, __global FluidRigidContacts* contacts, 
											__global b3Vector3* fluidVel, __global b3Vector3* fluidVelEval, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;

	//const int NUM_ITERATIONS = 4;
	//for(int iteration = 0; iteration < NUM_ITERATIONS; ++iteration)
	for(int contactIndex = 0; contactIndex < contacts[i].m_numContacts; ++contactIndex)
	{
		b3Vector3 fluidVelocity = fluidVel[i];
		FluidRigidContacts contact = contacts[i];
	
		int rigidIndex = contact.m_rigidIndicies[contactIndex];
		b3Scalar distance = contact.m_distances[contactIndex];
		b3Vector3 normalOnRigid = contact.m_normalsOnRigid[contactIndex];
		b3Vector3 pointOnRigid = contact.m_pointsOnRigid[contactIndex];
		
		b3Vector3 rigidPosition = rigidBodies[rigidIndex].m_pos;
		b3Scalar rigidInvMass = rigidBodies[rigidIndex].m_invMass;
		
		if( distance < 0.0f )
		{
			bool isDynamicRigidBody = (rigidInvMass != 0.0f);
			
			b3Vector3 rigidLocalHitPoint = pointOnRigid - rigidPosition;
			
			b3Vector3 rigidVelocity = (b3Vector3){0.0f, 0.0f, 0.0f, 0.0f};
			//if(isDynamicRigidBody) rigidVelocity = rigidBody->getVelocityInLocalPoint(rigidLocalHitPoint);
			rigidVelocity *= FG->m_simulationScale;
		
			b3Vector3 relativeVelocity = fluidVelocity - rigidVelocity;
			b3Scalar penetratingMagnitude = b3Vector3_dot(relativeVelocity, -normalOnRigid);
			if( penetratingMagnitude < 0.0f ) penetratingMagnitude = 0.0f;
			
			b3Vector3 penetratingVelocity = -normalOnRigid * penetratingMagnitude;
			b3Vector3 tangentialVelocity = relativeVelocity - penetratingVelocity;
			
			penetratingVelocity *= 1.0f + FL->m_boundaryRestitution;
			
			b3Scalar penetration = -distance;
			penetration = (penetration > FL->m_particleMargin) ? penetration : 0.0f;
			b3Scalar positionError = penetration * (FG->m_simulationScale/FG->m_timeStep) * FL->m_boundaryErp;
			//if(iteration != NUM_ITERATIONS - 1 ) positionError = 0;
			
			b3Vector3 particleImpulse = -(penetratingVelocity + (-normalOnRigid*positionError) + tangentialVelocity*FL->m_boundaryFriction);
			
			
			//Leapfrog integration
			b3Vector3 velNext = fluidVelocity + particleImpulse;
			fluidVelEval[i] = (fluidVelocity + velNext) * 0.5f;
			fluidVel[i] = velNext;
		}
	}
}
