
#version 330 core

layout (location = 0) in vec4 position;

uniform mat4 projectionMatrix;
uniform mat4 modelviewMatrix;
uniform vec2 screenDimensions;
uniform float sphereRadius;  	//Sphere size in world space

out vec3 eyeSpherePosition;

void main()
{
	eyeSpherePosition = ( modelviewMatrix * vec4(position.xyz, 1.0) ).xyz;
	
	vec4 projectedWidthHeight = projectionMatrix * vec4(sphereRadius, sphereRadius, eyeSpherePosition.z, 1.0);
	gl_PointSize = screenDimensions.x * projectedWidthHeight.x / projectedWidthHeight.w;
	//gl_PointSize = screenDimensions.y *  projectedWidthHeight.y / projectedWidthHeight.w;
	
	//Perspective projection distorts the sphere, so the square is not large enough and needs to be expanded
	const float SCALING = 1.6;
	gl_PointSize *= SCALING;
	
	mat4 modelviewProjectionMatrix = projectionMatrix * modelviewMatrix;
	gl_Position = modelviewProjectionMatrix * vec4(position.xyz, 1.0);
}

