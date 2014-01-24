
#include "b3FluidSphParticleUpdaterCL.h"

#include "Bullet3Common/b3Logging.h"		//B3_PROFILE(name) macro
#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"

#include "Bullet3Fluids/Sph/b3FluidSph.h"
#include "b3FluidSphOpenCL.h"

#include "fluidSphCL.h"

b3FluidSphParticleUpdaterCL::b3FluidSphParticleUpdaterCL(cl_context context, cl_device_id device, cl_command_queue queue) :
	m_createdPosition(context, queue), m_createdVelocity(context, queue),
	m_updatedPositionIndices(context, queue), m_updatedPosition(context, queue),
	m_updatedVelocityIndices(context, queue), m_updatedVelocity(context, queue),
	m_removeSwapSource(context, queue), m_removeSwapTarget(context, queue)
{
	m_context = context;
	m_commandQueue = queue;
	
	//
	const char CL_PROGRAM_PATH[] = "src/Bullet3FluidsOpenCL/fluidSph.cl";
	
	const char* kernelSource = fluidSphCL;	//fluidSphCL.h
	cl_int error;
	char* additionalMacros = 0;
	m_fluidsProgram = b3OpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, 
																	additionalMacros, CL_PROGRAM_PATH);
	b3Assert(m_fluidsProgram);
	
	m_setCreatedParticleAttributesKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "setCreatedParticleAttributes", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_setCreatedParticleAttributesKernel);
	m_applyParticleUpdatesKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "applyParticleUpdates", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_applyParticleUpdatesKernel);
	m_swapRemovedParticlesKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "swapRemovedParticles", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_swapRemovedParticlesKernel);
}

b3FluidSphParticleUpdaterCL::~b3FluidSphParticleUpdaterCL()
{
	clReleaseKernel(m_setCreatedParticleAttributesKernel);
	clReleaseKernel(m_applyParticleUpdatesKernel);
	clReleaseKernel(m_swapRemovedParticlesKernel);
	
	clReleaseProgram(m_fluidsProgram);
}

void b3FluidSphParticleUpdaterCL::createParticlesApplyUpdatesAndRemoveParticles(b3FluidSph* fluid, b3FluidSphOpenCL* fluidDataCL)
{
	b3FluidParticles& particlesCpu = fluid->internalGetParticles();
	b3FluidSphUpdatePacket& updates = fluid->internalGetUpdates();

	//Copy created particle data and updated particle data to GPU
	{
		m_createdPosition.copyFromHost(updates.m_addedParticlePositions, false);
		m_createdVelocity.copyFromHost(updates.m_addedParticleVelocities, false);
		
		m_updatedPositionIndices.copyFromHost(updates.m_updatedPositionsIndex, false);
		m_updatedPosition.copyFromHost(updates.m_updatedPositions, false);
		m_updatedVelocityIndices.copyFromHost(updates.m_updatedVelocitiesIndex, false);
		m_updatedVelocity.copyFromHost(updates.m_updatedVelocities, false);
		
		clFinish(m_commandQueue);
		
		b3Assert( m_createdPosition.size() == m_createdVelocity.size() );
		b3Assert( m_updatedPositionIndices.size() == m_updatedPosition.size() );
		b3Assert( m_updatedVelocityIndices.size() == m_updatedVelocity.size() );
		printf("createdSize: %d, %d\n", m_createdPosition.size(), m_createdVelocity.size() );
		printf("setposSize: %d, %d\n", m_updatedPositionIndices.size(), m_updatedPosition.size() );
		printf("setvelSize: %d, %d\n", m_updatedVelocityIndices.size(), m_updatedVelocity.size() );
	}

	//Create particles and resize arrays, including CPU
	int numCreatedParticles = m_createdPosition.size();
	if(numCreatedParticles)
	{
		int numParticlesBeforeAddingNew = fluid->numParticles();	//This is also the lowest/first index of the newly created particles
	
		int numParticlesAfterAddingNew = numParticlesBeforeAddingNew + numCreatedParticles;
		
		fluidDataCL->resize(numParticlesAfterAddingNew);
		particlesCpu.resize(numParticlesAfterAddingNew);	//Also resize on CPU, to ensure that fluid->numParticles() is correct
	
		b3BufferInfoCL bufferInfo[] = 
		{
			b3BufferInfoCL( m_createdPosition.getBufferCL() ),
			b3BufferInfoCL( m_createdVelocity.getBufferCL() ),
			
			b3BufferInfoCL( fluidDataCL->m_position.getBufferCL() ),
			b3BufferInfoCL( fluidDataCL->m_velocity.getBufferCL() ),
			b3BufferInfoCL( fluidDataCL->m_velocityEval.getBufferCL() )
		};
		
		b3LauncherCL launcher(m_commandQueue, m_setCreatedParticleAttributesKernel);
		launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
		launcher.setConst(numParticlesBeforeAddingNew);
		launcher.setConst(numCreatedParticles);
		
		launcher.launch1D(numCreatedParticles);
	}
	
	//Apply particle position updates
	int numPositionUpdates = m_updatedPositionIndices.size();
	if(numPositionUpdates)
	{
		b3BufferInfoCL bufferInfo[] = 
		{
			b3BufferInfoCL( m_updatedPositionIndices.getBufferCL() ),
			b3BufferInfoCL( m_updatedPosition.getBufferCL() ),
			b3BufferInfoCL( fluidDataCL->m_position.getBufferCL() )
		};
		
		b3LauncherCL launcher(m_commandQueue, m_applyParticleUpdatesKernel);
		launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
		launcher.setConst(numPositionUpdates);
		
		launcher.launch1D(numPositionUpdates);
	}
	
	//Apply particle velocity updates
	int numVelocityUpdates = m_updatedVelocityIndices.size();
	if(numVelocityUpdates)
	{
		{
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( m_updatedVelocityIndices.getBufferCL() ),
				b3BufferInfoCL( m_updatedVelocity.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocity.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_applyParticleUpdatesKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numVelocityUpdates);
			
			launcher.launch1D(numVelocityUpdates);
		}
		
		{
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( m_updatedVelocityIndices.getBufferCL() ),
				b3BufferInfoCL( m_updatedVelocity.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocityEval.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_applyParticleUpdatesKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numVelocityUpdates);
			
			launcher.launch1D(numVelocityUpdates);
		}
	}
	
	clFlush(m_commandQueue);
	
	//Remove marked particles
	b3AlignedObjectArray<int>& removedIndices = updates.m_removedParticleIndices;
		
	int numParticlesBeforeRemove = fluid->numParticles();
	int numRemovedParticles = removedIndices.size();
	int numParticlesAfterRemove = numParticlesBeforeRemove - numRemovedParticles;
	
	if(numRemovedParticles)
	{
		//Process for removing particles:
		//For instance there can be the array of particles:
		// N R N R N N R N N R
		//Where N is an index that is not marked to be removed, and R is marked to be removed.
		//In this case there are 10 particles, with 6 marked N and 4 marked R.
		//
		//Split the array at the last index after particles are removed(6 particles remain, so split occurs at 6th particle).
		// N R N R N N / R N N R
		//The main point is to notice that when the array is divided at the number of particles remaining,
		//the number of particles marked R to the left of the split is always the same as the number of particles marked N to the right.
		//
		//As a result, we can swap the particles marked N on the right with the particles marked R on the left.
		//
		//       |-Swap 3, 8-| 
		//       |           |
		// N R N R N N / R N N R        <-- 2 R on left, 2 N on right
		// 0 1 2 3 4 5   6 7 8 9
		//   |             |
		//   |- Swap 1, 7 -|
		//
		// N N N N N N / R R R R
		//and finally resize/truncate the array to eliminate the removed particles.
		//
		//In order to ensure that the result is deterministic, the array of marked particles
		//is sorted in ascending order and the source and target arrays for the swap are also 
		//sorted in ascending order.
		
		
		//makeUniqueInt() assumes that the array is sorted
		removedIndices.quickSort( b3FluidSphUpdatePacket::AscendingSortPredicate() );
	
		//Remove duplicate indicies
		b3FluidSphUpdatePacket::makeUniqueInt(removedIndices);
		
		//Find indices to the right of the split that are marked N(that is, not marked for remove) and use them as the swap source.
		{
			m_removeSwapSourceCpu.resize(0);
		
			int removedIndex = 0;
			for(int i = numParticlesAfterRemove; i < numParticlesBeforeRemove; ++i) 
			{
				while( i > removedIndices[removedIndex] && removedIndex < removedIndices.size() ) ++removedIndex;
			
				if(i != removedIndices[removedIndex]) m_removeSwapSourceCpu.push_back(i);
			}
		}
		
		//Find indices to the left of the split that are marked R(remove), and use them as the swap target.
		m_removeSwapTargetCpu.resize(0);
		for(int i = 0; i < numRemovedParticles; ++i)
		{
			if(removedIndices[i] < numParticlesAfterRemove) m_removeSwapTargetCpu.push_back(removedIndices[i]);
			else break;
		}
		
		//Copy swap data from CPU to GPU
		{
			m_removeSwapSource.copyFromHost(m_removeSwapSourceCpu, false);
			m_removeSwapTarget.copyFromHost(m_removeSwapTargetCpu, false);
			clFinish(m_commandQueue);
		}
		
		//Perform swap
		{
			int numSwappedParticles = m_removeSwapSource.size();
			b3Assert( m_removeSwapSource.size() == m_removeSwapTarget.size() );
		
			b3BufferInfoCL bufferInfo[] = 
			{
				b3BufferInfoCL( m_removeSwapSource.getBufferCL() ),
				b3BufferInfoCL( m_removeSwapTarget.getBufferCL() ),
				
				b3BufferInfoCL( fluidDataCL->m_position.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocity.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_velocityEval.getBufferCL() )
			};
			
			b3LauncherCL launcher(m_commandQueue, m_swapRemovedParticlesKernel);
			launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(b3BufferInfoCL) );
			launcher.setConst(numSwappedParticles);
			
			launcher.launch1D(numSwappedParticles);
			
			clFinish(m_commandQueue);
		}
		
		//Resize arrays
		fluidDataCL->resize(numParticlesAfterRemove);
		particlesCpu.resize(numParticlesAfterRemove);	//Also resize on CPU, to ensure that fluid->numParticles() is correct
	}
	
	//
	updates.clear();
}



