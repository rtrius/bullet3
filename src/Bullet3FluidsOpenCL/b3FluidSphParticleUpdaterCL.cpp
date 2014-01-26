
#include "b3FluidSphParticleUpdaterCL.h"

#include "Bullet3Common/b3Logging.h"		//B3_PROFILE(name) macro
#include "Bullet3OpenCL/ParallelPrimitives/b3LauncherCL.h"
#include "Bullet3OpenCL/Initialize/b3OpenCLUtils.h"

#include "Bullet3Fluids/Sph/b3FluidSph.h"
#include "b3FluidSphOpenCL.h"

#include "fluidSphCL.h"

b3FluidSphParticleUpdaterCL::b3FluidSphParticleUpdaterCL(cl_context context, cl_device_id device, cl_command_queue queue) :
	m_fill(context, device, queue),
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
	
	m_applyParticleUpdatesKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "applyParticleUpdates", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_applyParticleUpdatesKernel);
	m_swapRemovedParticlesKernel = b3OpenCLUtils::compileCLKernelFromString( context, device, kernelSource, "swapRemovedParticles", &error, m_fluidsProgram, additionalMacros );
	b3Assert(m_swapRemovedParticlesKernel);
}

b3FluidSphParticleUpdaterCL::~b3FluidSphParticleUpdaterCL()
{
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
	}
	//Create particles and resize arrays, including CPU
	int numCreatedParticles = m_createdPosition.size();
	if(numCreatedParticles)
	{
		int numParticlesBeforeAddingNew = fluid->numParticles();	//This is also the lowest/first index of the newly created particles
	
		int numParticlesAfterAddingNew = numParticlesBeforeAddingNew + numCreatedParticles;
		
		fluidDataCL->resize(numParticlesAfterAddingNew);
		particlesCpu.resize(numParticlesAfterAddingNew);	//Also resize on CPU, to ensure that fluid->numParticles() is correct
	
		m_createdPosition.copyToCL(fluidDataCL->m_position.getBufferCL(), numCreatedParticles, 0, numParticlesBeforeAddingNew);
		m_createdVelocity.copyToCL(fluidDataCL->m_velocity.getBufferCL(), numCreatedParticles, 0, numParticlesBeforeAddingNew);
		m_createdVelocity.copyToCL(fluidDataCL->m_velocityEval.getBufferCL(), numCreatedParticles, 0, numParticlesBeforeAddingNew);
		
		//b3FillCL::execute() operates on an array of floats, so there will be issues if b3Scalar == double
		b3Assert( sizeof(b3Scalar) == 4 && sizeof(b3Vector3) == 16 );
		m_fill.execute( reinterpret_cast<b3OpenCLArray<float>&>(fluidDataCL->m_accumulatedForce), 
						b3Scalar(0.0f), numCreatedParticles * 4, numParticlesBeforeAddingNew * 4 );
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
	
	//Remove particles
	int numRemovedParticles = updates.m_removedParticleIndices.size();
	if(numRemovedParticles)
	{
		int numParticlesBeforeRemove = fluid->numParticles();
		int numParticlesAfterRemove = numParticlesBeforeRemove - numRemovedParticles;
		
		//Load indices into updates.m_removeSwapSourceCpu and updates.m_removeSwapTargetCpu
		updates.prepareToRemoveParticles(numParticlesBeforeRemove);
		
		//Copy swap data from CPU to GPU
		{
			m_removeSwapSource.copyFromHost(updates.m_removeSwapSourceCpu, false);
			m_removeSwapTarget.copyFromHost(updates.m_removeSwapTargetCpu, false);
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
				b3BufferInfoCL( fluidDataCL->m_velocityEval.getBufferCL() ),
				b3BufferInfoCL( fluidDataCL->m_accumulatedForce.getBufferCL() )
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



