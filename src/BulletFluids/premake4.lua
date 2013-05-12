	project "BulletFluids"

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