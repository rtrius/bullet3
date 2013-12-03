/* 
ScreenSpaceFluidRendererGL.h
Copyright (C) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef SCREEN_SPACE_FLUID_RENDERER_GL_H
#define SCREEN_SPACE_FLUID_RENDERER_GL_H

#ifndef __APPLE__
#include <GL/glew.h>
#endif

#include "Bullet3Common/b3Vector3.h"
#include "Bullet3Common/b3AlignedObjectArray.h"

#include "FrameBufferGL.h"


//ScreenSpaceFluidRendererGL constructor may initialize GLEW; will fail if an OpenGL context does not exist
class ScreenSpaceFluidRendererGL
{
	int m_windowWidth;
	int m_windowHeight;

	GLuint m_positionVertexBuffer;
	GLuint m_squareVertexTexcoordBuffer;
	
	GLuint m_generateDepthProgram;
	GLuint m_blurDepthProgram;
	GLuint m_curvatureFlowProgram;
	
	GLuint m_blurThickProgram;
	GLuint m_absorptionAndTransparencyProgram;
	
	GLuint m_generateSurfaceProgram;
	GLuint m_blitProgram;
	
	FrameBufferGL m_frameBuffer;
	
	//Textures are managed(and deleted) by m_frameBuffer
	GLuint m_tempColorTexture;
	GLuint m_tempDepthTexture;
	
	GLuint m_depthTexture;
	GLuint m_blurredDepthTexturePass1;
	GLuint m_blurredDepthTexturePass2;
	
	GLuint m_thickTexture;
	GLuint m_blurredThickTexturePass1;
	GLuint m_blurredThickTexturePass2;
	GLuint m_absorptionAndTransparencyTexture;
	
	GLuint m_surfaceColorTexture;
	GLuint m_surfaceDepthTexture;
	
public:
	ScreenSpaceFluidRendererGL(int screenWidth, int screenHeight);
	~ScreenSpaceFluidRendererGL();

	void render(const float* projectionMatrix, const float* modelviewMatrix, const float* modelviewProjectionMatrix,
				const b3AlignedObjectArray<b3Vector3>& particlePositions, float sphereRadius, 
				float r, float g, float b, float absorptionR, float absorptionG, float absorptionB, bool copyVboFromCpuBuffer);
	
	//Use to load the Vertex Buffer Object(VBO) containing particle positions with OpenCL-OpenGL interop
	//positions are stored as float4, ordered as x,y,z,w
	GLuint getPositionVertexBuffer() { return m_positionVertexBuffer; }
	
	void setWindowResolution(int width, int height)
	{
		if(width <= 0 || height <= 0) return;	//ScreenSpaceFluidRendererGL::render() crashes if dimension is 0
		
		m_windowWidth = width; 
		m_windowHeight = height;
	}
	void setRenderingResolution(int width, int height)
	{
		if(width <= 0 || height <= 0) return;	//ScreenSpaceFluidRendererGL::render() crashes if dimension is 0
		m_frameBuffer.resizeTextures(width, height);
	}
	
	void getWindowResolution(int& out_width, int& out_height)
	{
		out_width = m_windowWidth;
		out_height = m_windowHeight;
	}
	
private:
	void initializeGlew();
	
	void render_stage1_generateDepthTexture(const float* projectionMatrix, const float* modelviewMatrix, 
											const float* modelviewProjectionMatrix,
											int numParticles, float sphereRadius);
	void render_stage2_blurDepthTextureCurvatureFlow(const float* modelviewProjectionMatrix);
	void render_stage2_blurDepthTextureBilateral(const float* modelviewProjectionMatrix);
	void render_stage3_generateThickTexture(const float* modelviewProjectionMatrix, int numParticles, float sphereRadius);
	void render_stage4_blurThickTexture(const float* modelviewProjectionMatrix);
	void render_stage5_generateAbsorptionAndTransparencyTexture(const float* modelviewProjectionMatrix, 
																float absorptionR, float absorptionG, float absorptionB);
	void render_stage6_generateSurfaceTexture(const float* projectionMatrix, const float* modelviewProjectionMatrix, 
											float r, float g, float b, bool useBlurredDepthTexture);
	
	void renderFullScreenTexture(GLuint texture2d_0, GLuint texture2d_1, GLuint texture2d_2);
	
	static GLuint compileProgram(const char* vertexShaderSource, const char* fragmentShaderSource);
};


#endif
