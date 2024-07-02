project "ptp"
	kind "StaticLib"
	language "C++"
	cppdialect "C++latest"
	staticruntime "on"

	targetdir ("build/bin/" .. outputdir .. "%{prj.name}")
	objdir ("build/bin-int/" .. outputdir .. "%{prj.name}")

	files
	{
		"src/**.h",
		"src/**.hpp",
		"src/**.cpp",
	}

	includedirs
	{
		"src/",
		"../third_party",
	}

	filter "configurations:Debug"
		runtime "Debug"
		symbols "on"
		defines {
			"PTP_DEBUG"
		}

	filter "configurations:Release"
		runtime "Release"
		optimize "on"
		defines {
			"PTP_RELEASE"
		}

	filter "configurations:Distribute"
		runtime "Release"
		optimize "on"
		defines {
			"PTP_DISTRIBUTE"
		}
	