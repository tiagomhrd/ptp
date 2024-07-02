workspace "ptp"
	configurations {"Debug", "Release", "Distribute"}
	platforms{"x64"}
	startproject "tests"
	filter { "platforms:x64" }
	architecture "x86_64"

outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"


include "third_party/Catch2"
include "tests"
include "ptp"