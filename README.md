# PardisoLink
A Mathematica interface for Intel MKL Pardiso

Download the zip archive, unzip it, and copy the resulting folder `PardisoLink` into

    FileNameJoin[{$UserBaseDirectory, "Applications"}]
    
Then the package can be loaded with

    Needs["PardisoLink`"]
    
After loading the package,you may visit 

    FileNameJoin[{PardisoLink`Private`$packageDirectory, "Documentation", "Examples", "Example.nb"}]

for a short tutorial.

See also the discussion on
https://mathematica.stackexchange.com/a/219940


Many thanks to Szabolcs Horvát whose package ``"LTemplate`"`` is used to wrestle down LibraryLink.
