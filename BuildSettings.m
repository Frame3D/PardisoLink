(* ::Package:: *)

(*
One can also link against the MKL shipped with Mathematica. 
For this, just use the alternative value of "LibraryDirectories". 
See also https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor 
for more details on compiler and linker options.
*)
Switch[$OperatingSystem,
  "MacOSX", (* Compilation settings for OS X *)
  $buildSettings = With[{MKLROOT="/opt/intel/mkl"},
  {
	"CompileOptions" -> {" -DMKL_ILP64"}
	, "LinkerOptions"->{" -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl -liomp5"}
    , "IncludeDirectories" -> {}
    (*, "LibraryDirectories" -> {FileNameJoin[{MKLROOT,"lib"}]}*)
    , "LibraryDirectories" -> {FileNameJoin[{$InstallationDirectory,"SystemFiles","Libraries",$SystemID}]}
 }],

  "Unix", (* Compilation settings for Linux *)
  $buildSettings = {
	"CompileOptions" -> {" -DMKL_ILP64"}
	, "LinkerOptions"->{" -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"}
    , "IncludeDirectories" -> {}
    , "LibraryDirectories" -> {FileNameJoin[{$InstallationDirectory,"SystemFiles","Libraries",$SystemID}]}
  },

  "Windows", (* Compilation settings for Windows *)
  $buildSettings = {
    "CompileOptions" -> {"/EHsc", "/wd4244", "/DNOMINMAX", " /DMKL_ILP64"}
	, "LinkerOptions"->{" mkl_intel_ilp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib"}
    , "IncludeDirectories" -> {}
    , "LibraryDirectories" -> {FileNameJoin[{$InstallationDirectory,"SystemFiles","Libraries",$SystemID}]}
  }
]
