//this file is autogenerated using stringify.bat (premake --stringify) in the build folder of this project
static const char* raycastSphereVertexShader= \
"#version 330 core\n"
"layout (location = 0) in vec4 position;\n"
"uniform mat4 projectionMatrix;\n"
"uniform mat4 modelviewMatrix;\n"
"uniform vec2 screenDimensions;\n"
"uniform float sphereRadius;  	//Sphere size in world space\n"
"out vec3 eyeSpherePosition;\n"
"void main()\n"
"{\n"
"	eyeSpherePosition = ( modelviewMatrix * vec4(position.xyz, 1.0) ).xyz;\n"
"	\n"
"	vec4 projectedWidthHeight = projectionMatrix * vec4(sphereRadius, sphereRadius, eyeSpherePosition.z, 1.0);\n"
"	gl_PointSize = screenDimensions.x * projectedWidthHeight.x / projectedWidthHeight.w;\n"
"	//gl_PointSize = screenDimensions.y *  projectedWidthHeight.y / projectedWidthHeight.w;\n"
"	\n"
"	//Perspective projection distorts the sphere, so the square is not large enough and needs to be expanded\n"
"	const float SCALING = 1.6;\n"
"	gl_PointSize *= SCALING;\n"
"	\n"
"	mat4 modelviewProjectionMatrix = projectionMatrix * modelviewMatrix;\n"
"	gl_Position = modelviewProjectionMatrix * vec4(position.xyz, 1.0);\n"
"}\n"
;
