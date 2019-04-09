# Installation 

So far, the package has not been registered with Julia.

Steps:

   1. Create a directory for Julia packages which have not been registered, let us denote this
      as PKG_DIR
    
   2. cd to PKG_DIR and clone this repository:

      Either from WIAS rhodecode server
````
     git clone  https://repos.wias-berlin.de/users/fuhrmann/projects/julia-packages/VoronoiFVM
````
     Or from github:
````
     git clone  https://github.com/j-fu/VoronoiFVM.jl VoronoiFVM
````
   
   3. Add the following line to  the file `.julia/config/startup.jl` in your  home directory (create this file it does not exist). PKG_DIR must be the full path name of the directory.
````
     push!(LOAD_PATH, "PKG_DIR")
````

Now, `import VoronoiFVM` should work in Julia scripts.
