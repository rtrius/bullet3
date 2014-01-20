	project "Bullet3Fluids"

	language "C++"
				
	kind "StaticLib"

	includedirs {
		".."
	}		
	targetdir "../../bin"

	files {
		"**.cpp",
		"**.h"
	}