
#version 330 core

uniform mat4 projectionMatrix;
uniform float sphereRadius;  	//Sphere size in world space
uniform vec4 sphereColor;
in vec3 eyeSpherePosition;     	//Position of sphere center in eye space

out vec4 fragColor;
void main()
{
	const float SCALING = 1.6;		//Make sure this matches SCALING in raycastSphereVS.glsl
	const float INV_SCALING_SQUARED = (1.0/SCALING) * (1.0/SCALING);	//If gl_PointSize is scaled by N, then this should be (1/N)^2
	
	vec3 pointOnSphere;
	vec3 normalAtPointOnSphere;
	
	//Perform ray-sphere intersection test
	{
		vec3 rayDirection = normalize( eyeSpherePosition + vec3(gl_PointCoord.xy*2.0 - 1.0, 0.0)*sphereRadius );
		vec3 sphereCenter = eyeSpherePosition;
		
		float b = -2.0 * dot(rayDirection, sphereCenter);
		float c = dot(sphereCenter, sphereCenter) - sphereRadius * sphereRadius * INV_SCALING_SQUARED;	//	determine cause of scaling
		
		float discriminant = b*b - 4.0 * c;
		if(discriminant < 0.0) discard;		//Negative discriminant == no intersection
		
		float discriminant_sqrt = sqrt(discriminant);
		float t1 = (-b + discriminant_sqrt) * 0.5;
		float t2 = (-b - discriminant_sqrt) * 0.5;
		
		pointOnSphere = rayDirection * min(t1, t2);
		normalAtPointOnSphere = normalize(pointOnSphere - eyeSpherePosition);
	}
	
	vec4 clipSpacePosition = projectionMatrix * vec4(pointOnSphere, 1.0);
	float depth = clipSpacePosition.z / clipSpacePosition.w;
	
	const vec3 LIGHT_DIRECTION = vec3(0.577, 0.577, 0.577);
	const float SHININESS = 40.0;
	
	float diffuse = max( 0.0, dot(normalAtPointOnSphere, LIGHT_DIRECTION) );
	
	vec3 v = normalize(-eyeSpherePosition);		//Normalized vector pointing at camera/eye
	vec3 h = normalize(LIGHT_DIRECTION + v);	//Normalized vector halfway between LIGHT_DIRECTION and v
	float specular = pow( max(0.0, dot(normalAtPointOnSphere, h)), SHININESS );
	
	vec3 color = sphereColor.xyz * diffuse + specular;
	fragColor = vec4(color, 1.0);
	
	//Convert depth from Normalized Device Coordinates(NDC) to Window/Screen coordinates
	gl_FragDepth = (gl_DepthRange.diff*depth + gl_DepthRange.near + gl_DepthRange.far) * 0.5;
}
