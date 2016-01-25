#ifndef B3_GPU_NARROWPHASE_H
#define B3_GPU_NARROWPHASE_H

#include "Bullet3Collision/NarrowPhaseCollision/shared/b3Collidable.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLInclude.h"
#include "Bullet3Common/b3AlignedObjectArray.h"
#include "Bullet3Common/b3Vector3.h"

class b3GpuNarrowPhase
{
protected:

	struct b3GpuNarrowPhaseInternalData*	m_data;
	int m_acceleratedCompanionShapeIndex;
	int	m_static0Index;

	cl_context m_context;
	cl_device_id m_device;
	cl_command_queue m_queue;

	int registerConvexHullShapeInternal(class b3ConvexUtility* convexPtr, b3Collidable& col);
	int registerConcaveMeshShape(b3AlignedObjectArray<b3Vector3>* vertices, b3AlignedObjectArray<int>* indices, b3Collidable& col, const float* scaling);

public:

	


	b3GpuNarrowPhase(cl_context vtx, cl_device_id dev, cl_command_queue q, const struct b3Config& config);

	virtual ~b3GpuNarrowPhase(void);

	int		registerSphereShape(float radius);
	int		registerPlaneShape(const b3Vector3& planeNormal, float planeConstant);

	int registerCompoundShape(b3AlignedObjectArray<b3GpuChildShape>* childShapes);
	int registerFace(const b3Vector3& faceNormal, float faceConstant);
	
	int	registerConcaveMesh(b3AlignedObjectArray<b3Vector3>* vertices, b3AlignedObjectArray<int>* indices,const float* scaling);
	
	//do they need to be merged?
	
	int	registerConvexHullShape(b3ConvexUtility* utilPtr);
	int	registerConvexHullShape(const float* vertices, int strideInBytes, int numVertices, const float* scaling);

	void	writeAllBodiesToGpu();
	void  reset();
	void	readbackAllBodiesToCpu();

	virtual void computeContacts(cl_mem broadphasePairs, int numBroadphasePairs, cl_mem aabbsWorldSpace, int numObjects);
	


	cl_mem	getCollidablesGpu();
	const struct b3Collidable* getCollidablesCpu() const;
	int		getNumCollidablesGpu() const;

	const struct b3SapAabb* getLocalSpaceAabbsCpu() const;

	const struct b3Contact4* getContactsCPU() const;

	cl_mem	getContactsGpu();
	int	getNumContactsGpu() const;

	cl_mem	getAabbLocalSpaceBufferGpu();
	

	int allocateCollidable();

	int getStatic0Index() const
	{
		return m_static0Index;
	}
	b3Collidable& getCollidableCpu(int collidableIndex);
	const b3Collidable& getCollidableCpu(int collidableIndex) const;

	const b3GpuNarrowPhaseInternalData*	getInternalData() const
	{
			return m_data;
	}

	b3GpuNarrowPhaseInternalData*	getInternalData()
	{
			return m_data;
	}

	const struct b3SapAabb& getLocalSpaceAabb(int collidableIndex) const;
};

#endif //B3_GPU_NARROWPHASE_H

