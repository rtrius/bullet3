/* 
ScreenSpaceFluidRendererGL.cpp
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
#include "ScreenSpaceFluidRendererGL.h"

#include <cstdlib> 	//exit()
#include <cstdio>


#define STRINGIFY(A) #A

const char generateDepthVertexShader[] = STRINGIFY(
	#version 330 core \n

	layout (location = 0) in vec4 position;

	uniform mat4 projectionMatrix;
	uniform mat4 modelviewMatrix;
	uniform mat4 modelviewProjectionMatrix;
	
	uniform vec2 screenDimensions;
	uniform float pointRadius;  	//Point size in world space
	out vec3 eyeSpherePosition;
	
	void main()
	{
		eyeSpherePosition = ( modelviewMatrix * vec4(position.xyz, 1.0) ).xyz;
		
		vec4 projectedWidthHeight = projectionMatrix * vec4(pointRadius, pointRadius, eyeSpherePosition.z, 1.0);
		gl_PointSize = screenDimensions.x * projectedWidthHeight.x / projectedWidthHeight.w;
		//gl_PointSize = screenDimensions.y *  projectedWidthHeight.y / projectedWidthHeight.w;
		
		//If USE_RAY_SPHERE_INTERSECTION is true, the square is not large enough to fit the sphere and needs to be expanded
		const float SCALING = 1.6;
		gl_PointSize *= SCALING;
		
		gl_Position = modelviewProjectionMatrix * vec4(position.xyz, 1.0);
	}
);
const char generateDepthFragmentShader[] = STRINGIFY(
	#version 330 core \n
	
	uniform mat4 projectionMatrix;
	uniform float pointRadius;  	//Point size in world space
	in vec3 eyeSpherePosition;      		//Position of sphere center in eye space
	
	out vec4 fragColor;
	void main()
	{
		const float INV_SCALING_SQUARED = (1.0/1.6) * (1.0/1.6);	//If gl_PointSize is scaled by N, then this should be (1/N)^2
		const int USE_RAY_SPHERE_INTERSECTION = 1;
		
		vec3 normalAtPointOnSphere;
		normalAtPointOnSphere.xy = gl_PointCoord.xy*2.0 - 1.0;
		float r2 = dot(normalAtPointOnSphere.xy, normalAtPointOnSphere.xy);
		if( (USE_RAY_SPHERE_INTERSECTION == 0) && r2 > 1.0*INV_SCALING_SQUARED) discard;
		
		float nearness = (USE_RAY_SPHERE_INTERSECTION) ? sqrt(1.0 - r2) : sqrt(1.0*INV_SCALING_SQUARED - r2);
		normalAtPointOnSphere.z = nearness;
		normalAtPointOnSphere = normalize(normalAtPointOnSphere);
		
		vec3 pointOnSphere = eyeSpherePosition + normalAtPointOnSphere*pointRadius;
		
		if(USE_RAY_SPHERE_INTERSECTION)
		{
			vec3 rayDirection = normalize( eyeSpherePosition + vec3(gl_PointCoord.xy*2.0 - 1.0, 0.0)*pointRadius );
			//vec3 rayDirection = normalize(pointOnSphere.xyz);
			vec3 sphereCenter = eyeSpherePosition;
			
			float b = -2.0 * dot(rayDirection, sphereCenter);
			float c = dot(sphereCenter, sphereCenter) - pointRadius * pointRadius * INV_SCALING_SQUARED;	//	determine cause of scaling
			
			float discriminant = b*b - 4.0 * c;
			if(discriminant < 0.0) discard;
			
			float discriminant_sqrt = sqrt(discriminant);
			float t1 = (-b + discriminant_sqrt) * 0.5;
			float t2 = (-b - discriminant_sqrt) * 0.5;
			
			pointOnSphere = rayDirection * min(t1, t2);
			normalAtPointOnSphere = normalize(pointOnSphere - eyeSpherePosition);
		}
		
		vec4 clipSpacePosition = projectionMatrix * vec4(pointOnSphere, 1.0);
		
		float thickness = 1.0;	//nearness * 0.03;	
		float depth = clipSpacePosition.z / clipSpacePosition.w;
		
		fragColor = vec4( vec3(1.0), thickness );
		gl_FragDepth = depth;
	}
);

const char fullScreenTextureVertexShader[] = STRINGIFY(
	#version 330 core \n
	
	layout (location = 0) in vec2 positionAndTexcoord;

	uniform mat4 modelviewProjectionMatrix;
	
	out vec2 texcoord;
	void main()
	{
		gl_Position = modelviewProjectionMatrix * vec4(positionAndTexcoord, 0, 1);
		texcoord = positionAndTexcoord;
	}
);

const char curvatureFlowShader[] = STRINGIFY(
	#version 330 core \n
	
	uniform vec2 focalLength;
	uniform vec2 texelSize;
	
	uniform float timeStep;
	uniform sampler2D depthTexture;
	
	in vec2 texcoord;
	out vec4 fragColor;
	
	//See "Screen Space Fluid Rendering with Curvature Flow"
	//by W. J. Van Der Laan, S. Green, M. Sainz. 
	void main()
	{
		float depth = texture(depthTexture, texcoord).x;
		
		const bool WRITE_CURVATURE_TO_COLOR_TEXTURE = false;
		if(WRITE_CURVATURE_TO_COLOR_TEXTURE && depth == gl_DepthRange.far)
		{
			gl_FragDepth = depth;
			fragColor = vec4(0.0, 0.0, 0.0, 1.0);
			return;
		}
		
		if(depth == gl_DepthRange.far) discard;
		
		vec2 upTexel = texcoord + vec2(0, texelSize.y);
		vec2 downTexel = texcoord + vec2(0, -texelSize.y);
		vec2 leftTexel = texcoord + vec2(texelSize.x, 0);
		vec2 rightTexel = texcoord + vec2(-texelSize.x, 0);

		float upDepth = texture(depthTexture, upTexel).x;
		float downDepth = texture(depthTexture, downTexel).x;
		float leftDepth = texture(depthTexture, leftTexel).x;
		float rightDepth = texture(depthTexture, rightTexel).x;
		
		//First order partial derivatives
		float dz_dx = (leftDepth + depth)*0.5 - (depth + rightDepth)*0.5;
		float dz_dy = (upDepth + depth)*0.5 - (depth + downDepth)*0.5;
		
		//Second order partial derivatives
		float d2z_dx2 = (leftDepth - 2.0*depth + rightDepth);	//d2z_dx2 == d^2z / dx^2
		float d2z_dy2 = (upDepth - 2.0*depth + downDepth);	
		
		float Cx = (2.0 / focalLength.x) * texelSize.x;		//texelSize.x == (1.0 / screenWidthInPixels)
		float Cy = (2.0 / focalLength.y) * texelSize.y;		//texelSize.y == (1.0 / screenHeightInPixels)
		float CxSquared = Cx*Cx;
		float CySquared = Cy*Cy;
		
		float D = CySquared*dz_dx*dz_dx + CxSquared*dz_dy*dz_dy + CxSquared*CySquared*depth*depth;
		
		float dD_dx;
		float dD_dy;
		{
			vec2 upLeftTexel = texcoord + vec2(texelSize.x, texelSize.y);
			vec2 upRightTexel = texcoord + vec2(-texelSize.x, texelSize.y);
			vec2 downLeftTexel = texcoord + vec2(texelSize.x, -texelSize.y);
			vec2 downRightTexel = texcoord + vec2(-texelSize.x, -texelSize.y);
		
			float upLeftDepth = texture(depthTexture, upLeftTexel).x;
			float upRightDepth = texture(depthTexture, upRightTexel).x;
			float downLeftDepth = texture(depthTexture, downLeftTexel).x;
			float downRightDepth = texture(depthTexture, downRightTexel).x;
			
			//Mixed partial derivatives
			float d2z_dxdy = (upLeftDepth - downLeftDepth - upRightDepth + downRightDepth)*0.25;
			float d2z_dydx = d2z_dxdy;
		
			dD_dx = 2.0*CySquared*d2z_dx2*dz_dx + 2.0*CxSquared*d2z_dxdy*dz_dy + 2.0*CxSquared*CySquared*dz_dx*depth;	
			dD_dy = 2.0*CySquared*d2z_dydx*dz_dx + 2.0*CxSquared*d2z_dy2*dz_dy + 2.0*CxSquared*CySquared*dz_dy*depth;	
		}
		
		float Ex = 0.5*dz_dx*dD_dx - d2z_dx2*D;
		float Ey = 0.5*dz_dy*dD_dy - d2z_dy2*D;
		
		float doubleH = (Cy*Ex + Cx*Ey) / pow(D, 1.5);
		float H = doubleH*0.5;	//H should be in [-1, 1]

		//If abs(H) is high, neighboring depths are likely part of another surface(discontinuous depth)
		const float H_THRESHOLD = 1.0;
		gl_FragDepth = ( abs(H) < H_THRESHOLD ) ? depth - H*timeStep : depth;
		//gl_FragDepth = ( abs(H) < H_THRESHOLD ) ? (upDepth + downDepth + leftDepth + rightDepth + depth)*0.20 : depth;
		
		if(WRITE_CURVATURE_TO_COLOR_TEXTURE)
		{
			if( abs(H) < H_THRESHOLD )
			{
				float curvature = H;
				float r = (curvature < 0) ? abs(curvature) : 0.0;
				float g = (curvature > 0) ? abs(curvature) : 0.0;
					
				fragColor = vec4( vec3(r, g, 0.0), 1.0 );
			}
			else fragColor = vec4(0.0, 0.0, 0.0, 1.0);
		}
	}
);

const char bilateralFilter1dFragmentShader_depth[] = STRINGIFY(
	#version 330 core \n
	
	uniform float texelSize;
	uniform float filterRadiusPixels;
	uniform float blurScale;			//blurScale: Lower values increase blur
	uniform float blurDepthFalloff;		//blurDepthFalloff: Higher values decrease blurring between pixels of differing intensity
	uniform vec2 blurDirection;
	uniform sampler2D depthTexture;
	
	in vec2 texcoord;
	void main()
	{
		float depth = texture(depthTexture, texcoord).x;
		if(depth == gl_DepthRange.far) discard;
		
		float sum = 0.0;
		float wsum = 0.0;
		for(float x = -filterRadiusPixels; x <= filterRadiusPixels; x += 1.0)
		{	
			float neighborDepth = texture(depthTexture, texcoord + blurDirection*texelSize*x).x;

			//Spatial domain
			float r = x * blurScale;
			float w = exp(-r*r);

			//Range domain
			float r2 = (neighborDepth - depth) * blurDepthFalloff;
			float g = exp(-r2*r2);
			
			sum += neighborDepth * w * g;
			wsum += w * g;
		}

		if(wsum > 0.0) sum /= wsum;
		gl_FragDepth = sum;
	}
);
const char bilateralFilter1dFragmentShader_alpha[] = STRINGIFY(
	#version 330 core \n
	
	uniform float texelSize;
	uniform float filterRadiusPixels;
	uniform float blurScale;			//blurScale: Lower values increase blur
	uniform float blurDepthFalloff;		//blurDepthFalloff: Higher values decrease blurring between pixels of differing intensity
	uniform vec2 blurDirection;
	uniform sampler2D alphaTexture;
	
	in vec2 texcoord;
	out vec4 fragColor;
	void main()
	{
		float depth = texture(alphaTexture, texcoord).a;
		
		float sum = 0.0;
		float wsum = 0.0;
		for(float x = -filterRadiusPixels; x <= filterRadiusPixels; x += 1.0)
		{	
			float neighborDepth = texture(alphaTexture, texcoord + blurDirection*texelSize*x).a;

			//Spatial domain
			float r = x * blurScale;
			float w = exp(-r*r);

			//Range domain
			float r2 = (neighborDepth - depth) * blurDepthFalloff;
			float g = exp(-r2*r2);
			
			sum += neighborDepth * w * g;
			wsum += w * g;
		}

		if(wsum > 0.0) sum /= wsum;
		fragColor = vec4( vec3(1.0), sum );
	}
);

const char absorptionAndTransparencyFragmentShader[] = STRINGIFY(
	#version 330 core \n

	//Beer's law / absorption constants(xyz == rgb)
	//Controls the darkening of the fluid's color based on its thickness
	//For a constant k, (k > 1) == darkens faster; (k < 1) == darkens slower; (k == 0) == disable
	uniform vec3 absorption;
	uniform sampler2D thicknessTexture;
	
	in vec2 texcoord;
	out vec4 fragColor;
	void main()
	{
		float thickness = texture(thicknessTexture, texcoord).a;
		fragColor = vec4( exp(-absorption.x * thickness), exp(-absorption.y * thickness), exp(-absorption.z * thickness), thickness );
	}
);

const char generateSurfaceFragmentShader[] = STRINGIFY(
	#version 330 core \n
	
	vec3 getEyePos(sampler2D depthTexture, vec2 texCoord, mat4 projectionMatrix)
	{
		float depth = texture(depthTexture, texCoord).x;
		
		vec4 unprojectedPosition = inverse(projectionMatrix) * vec4(texCoord.xy*2.0 - 1.0, depth, 1.0);
		
		vec3 eyePosition = unprojectedPosition.xyz / unprojectedPosition.w;
		return eyePosition;
	}
	
	uniform mat4 depthProjectionMatrix;		//Projection matrix used to generate depth values
	uniform vec2 texelSize;
	uniform vec4 baseColor;
	uniform sampler2D depthTexture;
	uniform sampler2D absorptionAndTransparencyTexture;
	
	in vec2 texcoord;
	out vec4 fragColor;
	void main()
	{
		float depth = texture(depthTexture, texcoord).x;
		if(depth == gl_DepthRange.far) discard;
		
		//Calculate normal using texCoords, depth, and projection matrix 
		vec3 eyePosition = getEyePos(depthTexture, texcoord, depthProjectionMatrix);

		vec3 ddx = getEyePos(depthTexture, texcoord + vec2(texelSize.x, 0), depthProjectionMatrix) - eyePosition;
		vec3 ddx2 = eyePosition - getEyePos(depthTexture, texcoord + vec2(-texelSize.x, 0), depthProjectionMatrix);
		if( abs(ddx.z) > abs(ddx2.z) ) ddx = ddx2;

		vec3 ddy = getEyePos(depthTexture, texcoord + vec2(0, texelSize.y), depthProjectionMatrix) - eyePosition;
		vec3 ddy2 = eyePosition - getEyePos(depthTexture, texcoord + vec2(0, -texelSize.y), depthProjectionMatrix);
		if( abs(ddy.z) > abs(ddy2.z) ) ddy = ddy2;
		
		vec3 normal = normalize( cross(ddx, ddy) );
		
		//
		const vec3 LIGHT_DIRECTION = vec3(0.577, 0.577, 0.577);
		const float SHININESS = 40.0;
		
		float diffuse = max( 0.0, dot(normal, LIGHT_DIRECTION) );
		
		vec3 v = normalize(-eyePosition);			//Normalized vector pointing at camera/eye
		vec3 h = normalize(LIGHT_DIRECTION + v);	//Normalized vector halfway between LIGHT_DIRECTION and v
		float specular = pow( max(0.0, dot(normal, h)), SHININESS );
		
		//
		vec4 absorptionAndTransparency = texture(absorptionAndTransparencyTexture, texcoord);
		
		const float MINIMUM_ALPHA = 0.70;
		vec3 color = baseColor.xyz * absorptionAndTransparency.xyz * diffuse + specular;
		float alpha = MINIMUM_ALPHA + absorptionAndTransparency.w * (1.0 - MINIMUM_ALPHA);
		
		//fragColor = vec4(color, alpha);
		fragColor = vec4(baseColor.xyz * diffuse + specular, 1.0);
		
		//Convert depth from Normalized Device Coordinates(NDC) to Window/Screen coordinates
		gl_FragDepth = (gl_DepthRange.diff*depth + gl_DepthRange.near + gl_DepthRange.far) * 0.5;
		
		//
		const int DISPLAY_FLUID = 0;
		const int DISPLAY_DEPTH = 1;
		const int DISPLAY_LINEAR_DEPTH = 2;
		const int DISPLAY_NORMAL = 3;
		const int DISPLAY_THICKNESS = 4;
		const int DISPLAY_ABSORPTION = 5;
		const int DISPLAY_CURVATURE = 6;
		
		const int DISPLAY_MODE = DISPLAY_FLUID;
		switch(DISPLAY_MODE)
		{
			case DISPLAY_DEPTH:
				fragColor = vec4( vec3(depth), 1.0 );
				break;
			case DISPLAY_LINEAR_DEPTH:	
				float linearDepth = min(1.0, -eyePosition.z / 10000.0);		//10000.0 == far z value
				fragColor = vec4( vec3(linearDepth), 1.0 );
				break;
			case DISPLAY_NORMAL:
				fragColor = vec4( (normal + 1.0) * 0.5, 1.0 );
				break;
			case DISPLAY_THICKNESS:
				float thickness = absorptionAndTransparency.a;
				fragColor = vec4( vec3(thickness), 1.0 );
				break;
			case DISPLAY_ABSORPTION:
				fragColor = vec4(absorptionAndTransparency.xyz, 1.0);
				break;
			case DISPLAY_CURVATURE:
				//Use m_tempColorTexture for absorptionAndTransparencyTexture when enabling this
				fragColor = texture(absorptionAndTransparencyTexture, texcoord);
				gl_FragDepth = 0.0;
				break;
				
			case DISPLAY_FLUID:
			default:
				//Fluid color is set above
				break;
		}
	}
);

const char blitFragmentShader[] = STRINGIFY(
	#version 330 core \n
	
	uniform sampler2D rgbaTexture;
	uniform sampler2D depthTexture;
	
	in vec2 texcoord;
	out vec4 fragColor;
	void main()
	{
		fragColor = texture(rgbaTexture, texcoord);
		gl_FragDepth = texture(depthTexture, texcoord).x;
	}
);


ScreenSpaceFluidRendererGL::ScreenSpaceFluidRendererGL(int screenWidth, int screenHeight)
{
	initializeGlew();

	//
	m_windowWidth = screenWidth;
	m_windowHeight = screenHeight;
	
	//
	m_generateDepthProgram = compileProgram(generateDepthVertexShader, generateDepthFragmentShader);
	m_blurDepthProgram = compileProgram(fullScreenTextureVertexShader, bilateralFilter1dFragmentShader_depth);
	m_curvatureFlowProgram = compileProgram(fullScreenTextureVertexShader, curvatureFlowShader);
	
	m_blurThickProgram = compileProgram(fullScreenTextureVertexShader, bilateralFilter1dFragmentShader_alpha);
	m_absorptionAndTransparencyProgram = compileProgram(fullScreenTextureVertexShader, absorptionAndTransparencyFragmentShader);
	
	m_generateSurfaceProgram = compileProgram(fullScreenTextureVertexShader, generateSurfaceFragmentShader);
	m_blitProgram = compileProgram(fullScreenTextureVertexShader, blitFragmentShader);
	
	//Generate vertex buffer for particle positions
	glGenBuffers(1, &m_positionVertexBuffer);
	
	//Generate vertex buffer for rendering full screen rectangle
	{
		//Arranged for GL_TRIANGLE_STRIP; order: Upper Left, Upper Right, Lower Left, Lower Right
		const GLfloat SQUARE_VERTICES[4 * 2] = { 0.0f,1.0f, 1.0f,1.0f, 0.0f,0.0f, 1.0f,0.0f };
		
		glGenBuffers(1, &m_squareVertexTexcoordBuffer);
		
		glBindBuffer(GL_ARRAY_BUFFER, m_squareVertexTexcoordBuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(SQUARE_VERTICES), SQUARE_VERTICES, GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
	
	//
	m_frameBuffer.initialize(screenWidth, screenHeight);
	m_tempColorTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_RGBA_TEXTURE);
	m_tempDepthTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	
	m_depthTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	m_blurredDepthTexturePass1 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	m_blurredDepthTexturePass2 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	
	m_thickTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_ALPHA_TEXTURE);
	m_blurredThickTexturePass1 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_ALPHA_TEXTURE);
	m_blurredThickTexturePass2 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_ALPHA_TEXTURE);
	m_absorptionAndTransparencyTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_RGBA_TEXTURE);
	
	m_surfaceColorTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_RGBA_TEXTURE);
	m_surfaceDepthTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	
}
ScreenSpaceFluidRendererGL::~ScreenSpaceFluidRendererGL()
{
	glDeleteProgram(m_generateDepthProgram);
	
	glDeleteProgram(m_blurDepthProgram);
	glDeleteProgram(m_curvatureFlowProgram);
	
	glDeleteProgram(m_blurThickProgram);
	glDeleteProgram(m_absorptionAndTransparencyProgram);
	
	glDeleteProgram(m_generateSurfaceProgram);
	glDeleteProgram(m_blitProgram);
	
	//
	glDeleteBuffers(1, &m_positionVertexBuffer);
	glDeleteBuffers(1, &m_squareVertexTexcoordBuffer);
	
	//
	m_frameBuffer.deactivate();
}
	
void ScreenSpaceFluidRendererGL::render(const float* projectionMatrix, const float* modelviewMatrix, const float* modelviewProjectionMatrix,
										const b3AlignedObjectArray<b3Vector3>& particlePositions, float sphereRadius, 
										float r, float g, float b, float absorptionR, float absorptionG, float absorptionB, bool copyVboFromCpuBuffer)
{		
	//Column major order; this is the modelview projection matrix
	//that results from calling(OpenGL 2.0):
	//	glMatrixMode(GL_MODELVIEW);	
	//	glLoadIdentity();
	//	glMatrixMode(GL_PROJECTION);
	//	glLoadIdentity();
	//	glOrtho(0, 1, 0, 1, -1, 1);
	const GLfloat renderFullScreenTextureMatrix[16] = {  2,  0,  0,  0,
														 0,  2,  0,  0,
														 0,  0, -1,  0,
														-1, -1,  0,  1	};

	b3Assert( sizeof(b3Vector3) == 16 );
	
	int numParticles = particlePositions.size();
	if(!numParticles) return;
	
	//Load particle positions
	if(copyVboFromCpuBuffer)
	{
		glBindBuffer(GL_ARRAY_BUFFER, m_positionVertexBuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(b3Vector3) * numParticles, &particlePositions[0], GL_DYNAMIC_DRAW);
	}
	
	//
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);
	
	bool renderingResolutionDiffers = ( m_windowWidth != m_frameBuffer.getWidth() || m_windowHeight != m_frameBuffer.getHeight() );
	if(renderingResolutionDiffers) glViewport( 0, 0, m_frameBuffer.getWidth(), m_frameBuffer.getHeight() );
	
	render_stage1_generateDepthTexture( projectionMatrix, modelviewMatrix, modelviewProjectionMatrix, numParticles, sphereRadius);
	
	const bool BLUR_DEPTH_TEXTURE = 1;
	if(BLUR_DEPTH_TEXTURE)
	{
		const bool USE_CURVATURE_FLOW = 0;
		if(USE_CURVATURE_FLOW) render_stage2_blurDepthTextureCurvatureFlow(renderFullScreenTextureMatrix);
		else render_stage2_blurDepthTextureBilateral(renderFullScreenTextureMatrix);
	}
	
	render_stage3_generateThickTexture(renderFullScreenTextureMatrix, numParticles, sphereRadius);
	
	render_stage4_blurThickTexture(renderFullScreenTextureMatrix);
	
	render_stage5_generateAbsorptionAndTransparencyTexture(renderFullScreenTextureMatrix, absorptionR, absorptionG, absorptionB);
	
	render_stage6_generateSurfaceTexture(projectionMatrix, renderFullScreenTextureMatrix, r, g, b, BLUR_DEPTH_TEXTURE);
	
	//Blit results to the main/window frame buffer
	if(renderingResolutionDiffers) glViewport(0, 0, m_windowWidth, m_windowHeight);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glUseProgram(m_blitProgram);
	glUniformMatrix4fv( glGetUniformLocation(m_blitProgram, "modelviewProjectionMatrix"), 1, false, renderFullScreenTextureMatrix );
	glUniform1i( glGetUniformLocation(m_blitProgram, "rgbaTexture"), 0 );
	glUniform1i( glGetUniformLocation(m_blitProgram, "depthTexture"), 1 );
		//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		renderFullScreenTexture(m_surfaceColorTexture, m_surfaceDepthTexture, 0);
	glDisable(GL_BLEND);
	glUseProgram(0);
	
	//Default clear color for Bullet demos
	glClearColor(1.0f, 1.0, 1.0f, 1.0f);
}
	
void ScreenSpaceFluidRendererGL::initializeGlew()
{
#ifndef __APPLE__
	GLenum errorCode = glewInit();
	if(errorCode != GLEW_OK)
	{
		printf( "GLEW error: %s(%d) \n", glewGetErrorString(errorCode), errorCode );
	}
	
	const int NUM_REQUIRED_EXTENSIONS = 6;
	const char* requiredExtensions[NUM_REQUIRED_EXTENSIONS] =
	{
		//"GL_VERSION_4_2",
		"GL_VERSION_3_3",
		"GL_VERSION_2_0",
		"GL_ARB_multitexture",
		"GL_ARB_vertex_buffer_object",
		"GL_ARB_framebuffer_object",
		"GL_ARB_point_sprite"
	};
	
	bool areRequiredExtensionsMissing = false;
	for(int i = 0; i < NUM_REQUIRED_EXTENSIONS; ++i)
	{
		if( !glewIsSupported(requiredExtensions[i]) )
		{
			printf("Required extension: %s is missing.\n", requiredExtensions[i]);
			areRequiredExtensionsMissing = true;
		}
	}
	
	if(areRequiredExtensionsMissing) 
	{
		printf("ScreenSpaceFluidRendererGL: Required OpenGL extensions missing.\n");
		printf("Press enter to exit.\n");
		getchar();
		exit(-1);
	}
#endif
}

void ScreenSpaceFluidRendererGL::render_stage1_generateDepthTexture(const float* projectionMatrix, const float* modelviewMatrix, 
																	const float* modelviewProjectionMatrix,
																	int numParticles, float sphereRadius)
{
	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_POINT_SPRITE);
	glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	glPointParameteri(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);

	glUseProgram(m_generateDepthProgram);
	
	{
		float screenWidth = static_cast<float>( m_frameBuffer.getWidth() );
		float screenHeight = static_cast<float>( m_frameBuffer.getHeight() );
		glUniformMatrix4fv( glGetUniformLocation(m_generateDepthProgram, "projectionMatrix"), 1, false, projectionMatrix );
		glUniformMatrix4fv( glGetUniformLocation(m_generateDepthProgram, "modelviewMatrix"), 1, false, modelviewMatrix );
		glUniformMatrix4fv( glGetUniformLocation(m_generateDepthProgram, "modelviewProjectionMatrix"), 1, false, modelviewProjectionMatrix );
		glUniform2f( glGetUniformLocation(m_generateDepthProgram, "screenDimensions"), screenWidth, screenHeight );
		glUniform1f( glGetUniformLocation(m_generateDepthProgram, "pointRadius"), sphereRadius );
		
		glBindBuffer(GL_ARRAY_BUFFER, m_positionVertexBuffer);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
		
		m_frameBuffer.attachAndSetRenderTargets(m_tempColorTexture, m_depthTexture);
			glDrawArrays(GL_POINTS, 0, numParticles);
		m_frameBuffer.detachAndUseDefaultFrameBuffer();
		
		glDisableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
	
	glUseProgram(0);
	glPointParameteri(GL_POINT_SPRITE_COORD_ORIGIN, GL_UPPER_LEFT);	//Reset to default, GL_UPPER_LEFT
	glDisable(GL_POINT_SPRITE);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void ScreenSpaceFluidRendererGL::render_stage2_blurDepthTextureCurvatureFlow(const float* modelviewProjectionMatrix)
{
	glDepthFunc(GL_ALWAYS);
	glUseProgram(m_curvatureFlowProgram);
	
	//const float HORIZONTAL_FOV_RADIANS = B3_RADS_PER_DEG * 90.0f;
	//const float VERTICAL_FOV_RADIANS = B3_RADS_PER_DEG * 75.0f;

	//const float SCREEN_WIDTH = 1.0f;
	//const float SCREEN_HEIGHT = 0.75f;
	
	//	focal length probably incorrect
	//const float FOCAL_LENGTH_X = (SCREEN_WIDTH * 0.5f) / b3Tan(HORIZONTAL_FOV_RADIANS * 0.5f);
	//const float FOCAL_LENGTH_Y = (SCREEN_HEIGHT * 0.5f) / b3Tan(VERTICAL_FOV_RADIANS * 0.5f);
	
	const float FOCAL_LENGTH_X = 0.5;
	const float FOCAL_LENGTH_Y = 0.5;
	
	float texelSizeX = 1.0f / static_cast<float>( m_frameBuffer.getWidth() );
	float texelSizeY = 1.0f / static_cast<float>( m_frameBuffer.getHeight() );
	
	glUniformMatrix4fv( glGetUniformLocation(m_curvatureFlowProgram, "modelviewProjectionMatrix"), 1, false, modelviewProjectionMatrix );
	glUniform2f( glGetUniformLocation(m_curvatureFlowProgram, "focalLength"), FOCAL_LENGTH_X, FOCAL_LENGTH_Y );
	glUniform2f( glGetUniformLocation(m_curvatureFlowProgram, "texelSize"), texelSizeX, texelSizeY );
	glUniform1f( glGetUniformLocation(m_curvatureFlowProgram, "timeStep"), 0.001f );
	glUniform1i( glGetUniformLocation(m_curvatureFlowProgram, "depthTexture"), 0 );

	GLuint depthTextures[2] = { m_blurredDepthTexturePass2, m_depthTexture };
	
	//Since textures are swapped, apply an odd number of iterations to ensure that the final result is in depthTextures[0]
	const int ITERATIONS = 120;
	const int ACTUAL_ITERATIONS = (ITERATIONS % 2) ? ITERATIONS : ITERATIONS + 1;
	for(int i = 0; i < ACTUAL_ITERATIONS; ++i)
	{
		m_frameBuffer.attachAndSetRenderTargets(m_tempColorTexture, depthTextures[0]);
			renderFullScreenTexture(depthTextures[1], 0, 0);		
		m_frameBuffer.detachAndUseDefaultFrameBuffer();
		
		//Final output will be in depthTextures[0] == m_blurredDepthTexturePass2
		GLuint swap = depthTextures[0];
		depthTextures[0] = depthTextures[1];
		depthTextures[1] = swap;
	}

	glDepthFunc(GL_LESS);
	glUseProgram(0);
}
void ScreenSpaceFluidRendererGL::render_stage2_blurDepthTextureBilateral(const float* modelviewProjectionMatrix)
{
	glDepthFunc(GL_ALWAYS);
	glUseProgram(m_blurDepthProgram);
	
	glUniformMatrix4fv( glGetUniformLocation(m_blurDepthProgram, "modelviewProjectionMatrix"), 1, false, modelviewProjectionMatrix );
		
	//First pass blurs along the x-axis
	glUniform1f( glGetUniformLocation(m_blurDepthProgram, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getWidth() ) );
	glUniform1f( glGetUniformLocation(m_blurDepthProgram, "filterRadiusPixels"), 64.0f );
	glUniform1f( glGetUniformLocation(m_blurDepthProgram, "blurScale"), 0.3f );
	glUniform1f( glGetUniformLocation(m_blurDepthProgram, "blurDepthFalloff"), 600.0f );
	glUniform2f( glGetUniformLocation(m_blurDepthProgram, "blurDirection"), 1.0f, 0.0f );
	glUniform1i( glGetUniformLocation(m_blurDepthProgram, "depthTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_tempColorTexture, m_blurredDepthTexturePass1);
		renderFullScreenTexture(m_depthTexture, 0, 0);		
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	//Second pass blurs along the y-axis
	glUniform1f( glGetUniformLocation(m_blurDepthProgram, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getHeight() ) );
	glUniform1f( glGetUniformLocation(m_blurDepthProgram, "filterRadiusPixels"), 64.0f );
	glUniform1f( glGetUniformLocation(m_blurDepthProgram, "blurScale"), 0.3f );
	glUniform1f( glGetUniformLocation(m_blurDepthProgram, "blurDepthFalloff"), 600.0f );
	glUniform2f( glGetUniformLocation(m_blurDepthProgram, "blurDirection"), 0.0f, 1.0f );
	glUniform1i( glGetUniformLocation(m_blurDepthProgram, "depthTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_tempColorTexture, m_blurredDepthTexturePass2);
		renderFullScreenTexture(m_blurredDepthTexturePass1, 0, 0);
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glDepthFunc(GL_LESS);
	glUseProgram(0);
}
void ScreenSpaceFluidRendererGL::render_stage3_generateThickTexture(const float* modelviewProjectionMatrix, int numParticles, float sphereRadius)
{	
	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_POINT_SPRITE);
	glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	glUseProgram(m_generateDepthProgram);

	glPointParameteri(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);
	
	glDepthFunc(GL_ALWAYS);
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE);
	{
		float screenWidth = static_cast<float>( m_frameBuffer.getWidth() );
		float screenHeight = static_cast<float>( m_frameBuffer.getHeight() );
		float lesserDistance = (screenWidth > screenHeight) ?  screenHeight : screenWidth;
		
		glUniformMatrix4fv( glGetUniformLocation(m_generateDepthProgram, "modelviewProjectionMatrix"), 1, false, modelviewProjectionMatrix );
		glUniform1f( glGetUniformLocation(m_generateDepthProgram, "pointScale"), lesserDistance );
		glUniform1f( glGetUniformLocation(m_generateDepthProgram, "pointRadius"), sphereRadius );
		
		glBindBuffer(GL_ARRAY_BUFFER, m_positionVertexBuffer);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
		
		m_frameBuffer.attachAndSetRenderTargets(m_thickTexture, m_tempDepthTexture);
			glDrawArrays(GL_POINTS, 0, numParticles);
		m_frameBuffer.detachAndUseDefaultFrameBuffer();
		
		glDisableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
	glDepthFunc(GL_LESS);
	glDisable(GL_BLEND);
	
	glUseProgram(0);
	glPointParameteri(GL_POINT_SPRITE_COORD_ORIGIN, GL_UPPER_LEFT);	//Reset to default, GL_UPPER_LEFT
	glDisable(GL_POINT_SPRITE);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}
void ScreenSpaceFluidRendererGL::render_stage4_blurThickTexture(const float* modelviewProjectionMatrix)
{	
	glDepthFunc(GL_ALWAYS);
	glUseProgram(m_blurThickProgram);
	
	glUniformMatrix4fv( glGetUniformLocation(m_blurThickProgram, "modelviewProjectionMatrix"), 1, false, modelviewProjectionMatrix );
	
	//First pass blurs along the x-axis
	glUniform1f( glGetUniformLocation(m_blurThickProgram, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getWidth() ) );
	glUniform1f( glGetUniformLocation(m_blurThickProgram, "filterRadiusPixels"), 32.0f );
	glUniform1f( glGetUniformLocation(m_blurThickProgram, "blurScale"), 0.1f );
	glUniform1f( glGetUniformLocation(m_blurThickProgram, "blurDepthFalloff"), 6.0f );
	glUniform2f( glGetUniformLocation(m_blurThickProgram, "blurDirection"), 1.0f, 0.0f );
	glUniform1i( glGetUniformLocation(m_blurThickProgram, "alphaTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_blurredThickTexturePass1, m_tempDepthTexture);
		renderFullScreenTexture(m_thickTexture, 0, 0);		
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	//Second pass blurs along the y-axis
	glUniform1f( glGetUniformLocation(m_blurThickProgram, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getHeight() ) );
	glUniform1f( glGetUniformLocation(m_blurThickProgram, "filterRadiusPixels"), 32.0f );
	glUniform1f( glGetUniformLocation(m_blurThickProgram, "blurScale"), 0.1f );
	glUniform1f( glGetUniformLocation(m_blurThickProgram, "blurDepthFalloff"), 6.0f );
	glUniform2f( glGetUniformLocation(m_blurThickProgram, "blurDirection"), 0.0f, 1.0f );
	glUniform1i( glGetUniformLocation(m_blurThickProgram, "alphaTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_blurredThickTexturePass2, m_tempDepthTexture);
		renderFullScreenTexture(m_blurredThickTexturePass1, 0, 0);
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glDepthFunc(GL_LESS);
	glUseProgram(0);
}
void ScreenSpaceFluidRendererGL::render_stage5_generateAbsorptionAndTransparencyTexture(const float* modelviewProjectionMatrix, 
																						float absorptionR, float absorptionG, float absorptionB)
{
	glDepthFunc(GL_ALWAYS);
	glUseProgram(m_absorptionAndTransparencyProgram);
	
	glUniformMatrix4fv( glGetUniformLocation(m_absorptionAndTransparencyProgram, "modelviewProjectionMatrix"), 1, false, modelviewProjectionMatrix );
	glUniform3f( glGetUniformLocation(m_absorptionAndTransparencyProgram, "absorption"), absorptionR, absorptionG, absorptionB);
	glUniform1i( glGetUniformLocation(m_absorptionAndTransparencyProgram, "thicknessTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_absorptionAndTransparencyTexture, m_tempDepthTexture);
		renderFullScreenTexture(m_blurredThickTexturePass2, 0, 0);
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glDepthFunc(GL_LESS);
	glUseProgram(0);
}

void ScreenSpaceFluidRendererGL::render_stage6_generateSurfaceTexture(const float* projectionMatrix, const float* modelviewProjectionMatrix, 
																		float r, float g, float b, bool useBlurredDepthTexture)
{
	float texelSize_x = 1.0f / static_cast<float>( m_frameBuffer.getWidth() );
	float texelSize_y = 1.0f / static_cast<float>( m_frameBuffer.getHeight() );

	glUseProgram(m_generateSurfaceProgram);
	glUniformMatrix4fv( glGetUniformLocation(m_generateSurfaceProgram, "depthProjectionMatrix"), 1, false, projectionMatrix );
	glUniformMatrix4fv( glGetUniformLocation(m_generateSurfaceProgram, "modelviewProjectionMatrix"), 1, false, modelviewProjectionMatrix );
	glUniform2f( glGetUniformLocation(m_generateSurfaceProgram, "texelSize"), texelSize_x, texelSize_y );
	glUniform4f( glGetUniformLocation(m_generateSurfaceProgram, "baseColor"), r, g, b, 1.0f );
	glUniform1i( glGetUniformLocation(m_generateSurfaceProgram, "depthTextureBlurred"), 0 );
	glUniform1i( glGetUniformLocation(m_generateSurfaceProgram, "absorptionAndTransparencyTexture"), 1 );
	
	GLuint depthTexture = (useBlurredDepthTexture) ? m_blurredDepthTexturePass2 : m_depthTexture;
	
	m_frameBuffer.attachAndSetRenderTargets(m_surfaceColorTexture, m_surfaceDepthTexture);
		renderFullScreenTexture(depthTexture, m_absorptionAndTransparencyTexture, 0);
		//renderFullScreenTexture(depthTexture, m_tempColorTexture, 0);		//For display of curvature
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glUseProgram(0);
}

void ScreenSpaceFluidRendererGL::renderFullScreenTexture(GLuint texture2d_0, GLuint texture2d_1, GLuint texture2d_2)
{
	//Enable states
	glEnable(GL_TEXTURE_2D);
	
	//
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, texture2d_2);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, texture2d_1);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texture2d_0);
	
	//Render
	{
		glBindBuffer(GL_ARRAY_BUFFER, m_squareVertexTexcoordBuffer);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
		
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		
		glDisableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
	
	//Disable states
	glDisable(GL_TEXTURE_2D);
	
	//
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, 0);
}

 GLuint ScreenSpaceFluidRendererGL::compileProgram(const char* vertexShaderSource, const char* fragmentShaderSource)
{
	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

	glShaderSource(vertexShader, 1, &vertexShaderSource, 0);
	glShaderSource(fragmentShader, 1, &fragmentShaderSource, 0);
	
	glCompileShader(vertexShader);
	glCompileShader(fragmentShader);

	//
	GLuint program = glCreateProgram();

	glAttachShader(program, vertexShader);
	glAttachShader(program, fragmentShader);

	glLinkProgram(program);

	GLint success = 0;
	glGetProgramiv(program, GL_LINK_STATUS, &success);

	if(!success) 
	{
		const int MAX_STRING_LENGTH = 65536;
		b3AlignedObjectArray<char> string;
		string.resize(MAX_STRING_LENGTH);
		char* stringStart = &string[0];
		
		glGetProgramInfoLog(program, MAX_STRING_LENGTH, 0, stringStart);
		printf("GL Program Build Log:\n");
		printf("%s\n", stringStart);
		
		glGetShaderInfoLog(vertexShader, MAX_STRING_LENGTH, 0, stringStart);
		printf("Vertex Shader Build Log:\n");
		printf("%s\n", stringStart);
		
		glGetShaderInfoLog(fragmentShader, MAX_STRING_LENGTH, 0, stringStart);
		printf("Fragment Shader Build Log:\n");
		printf("%s\n", stringStart);
		
		glDeleteProgram(program);
		program = 0;
	}
	
	//If program was compiled successfully, marks shaders for deletion(not deleted until program is deleted)
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
	
	return program;
}


