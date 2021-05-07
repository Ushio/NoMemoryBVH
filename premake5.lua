include "libs/PrLib"

workspace "SampleApp"
    location "build"
    configurations { "Debug", "Release" }
    startproject "main"

architecture "x86_64"

externalproject "prlib"
	location "libs/PrLib/build" 
    kind "StaticLib"
    language "C++"

project "main"
    kind "ConsoleApp"
    language "C++"
    targetdir "bin/"
    systemversion "latest"
    flags { "MultiProcessorCompile", "NoPCH" }

    -- Src
    files { "main.cpp" }

    -- UTF8
    postbuildcommands { 
        "mt.exe -manifest ../utf8.manifest -outputresource:$(TargetDir)$(TargetName).exe -nologo"
    }

    -- USD
    defines { "NOMINMAX", "__TBB_NO_IMPLICIT_LINKAGE" }
    includedirs { "libs/USDvs2019/include" }
    libdirs { "libs/USDvs2019/lib" }

    filter {"Debug"}
        libdirs { "libs/USDvs2019/libd" }
    filter {"Release"}
        libdirs { "libs/USDvs2019/lib" }
    filter{}

    links { 
        "usdGeom",
        "usdHydra",
        "usdLux",
        "usdMedia",
        "usdRender",
        "usdRi",
        "usdShade",
        "usdSkel",
        "usdUI",
        "usdUtils",
        "usdVol",
        "vt",
        "work",
        "ar",
        "arch",
        "gf",
        "js",
        "kind",
        "ndr",
        "pcp",
        "plug",
        "sdf",
        "sdr",
        "tf",
        "trace",
        "usd",
        "tbb",
    }

    
    filter {"Debug"}
        debugenvs { "PATH=%PATH%;$(ProjectDir)../libs/USDvs2019/libd" }
    filter {"Release"}
        debugenvs { "PATH=%PATH%;$(ProjectDir)../libs/USDvs2019/lib" }
    filter{}
    -- prlib
    -- setup command
    -- git submodule add https://github.com/Ushio/prlib libs/prlib
    -- premake5 vs2017
    dependson { "prlib" }
    includedirs { "libs/prlib/src" }
    libdirs { "libs/prlib/bin" }
    filter {"Debug"}
        links { "prlib_d" }
    filter {"Release"}
        links { "prlib" }
    filter{}

    symbols "On"

    filter {"Debug"}
        runtime "Debug"
        targetname ("Main_Debug")
        optimize "Off"
    filter {"Release"}
        runtime "Release"
        targetname ("Main")
        optimize "Full"
    filter{}
