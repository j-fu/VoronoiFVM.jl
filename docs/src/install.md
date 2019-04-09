# Installation 

So far, the package has not been registered with the Julia central, however Julia provides
the possibility to install packages from any accessible git repo.

## Canonical Julia method

Go into package manager  (`]` key) and add the package via the repo name:

````
     pkg> add https://github.com/j-fu/VoronoiFVM.jl
````

or 

````
     pkg> add https://repos.wias-berlin.de/users/fuhrmann/projects/julia-packages/VoronoiFVM.jl
````

This also will consistently pull in any necessary package dependencies.

You can run the test suite

````
     pkg> test VoronoiFVM
````

If you don't like this anymore (or want to switch between remote repos), you can remove it:

````
     pkg> remove VoronoiFVM
````

## Old method (which still works)
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

Now, `using VoronoiFVM` should work in Julia scripts.
