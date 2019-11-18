#=
This file has been created after the corresponding example in TriangleMesh.jl
=#
using BinDeps

isfile("deps.jl") && rm("deps.jl")

@BinDeps.setup
libtriangle = library_dependency("libtriangle", aliases = ["libtriangle.dylib"], runtime = true)

rootdir = BinDeps.depsdir(libtriangle)
srcdir = joinpath(rootdir, "src")
prefix = joinpath(rootdir, "usr")
libdir = joinpath(prefix, "lib")
mkpath(libdir)

if Sys.iswindows()
    # not checked
    libfile = joinpath(libdir, "libtriangle.dll")
    arch = "x86"
    if Sys.WORD_SIZE == 64
        arch = "x64"
    end
    @build_steps begin
        FileRule(libfile, @build_steps begin
                 BinDeps.run(@build_steps begin
                             ChangeDirectory(srcdir)
                             `cmd /c compile.bat all $arch`
                             `cmd /c copy libtriangle.dll $libfile`
                             `cmd /c copy triangle.h $headerdir`
                             `cmd /c copy tesselate.h $headerdir`
                             `cmd /c copy commondefine.h $headerdir`
                             `cmd /c compile.bat clean $arch`
                             end) end) end
    
    provides(Binaries, URI(libfile), libtriangle)
else
    libname = "libtriangle.so"
    ldflags = "--shared"
    if Sys.isapple()
        libname = "libtriangle.dylib"
#        ldflags="-dynamiclib -undefined suppress -flat_namespace"
    end
    libfile = joinpath(libdir, libname)
    provides(BinDeps.BuildProcess, (@build_steps begin
                                    FileRule(libfile, @build_steps begin
                                             BinDeps.ChangeDirectory(srcdir)
                                             `cc  -DTRILIBRARY -fPIC -DNDEBUG -DNO_TIMER $(ldflags) -o $(libfile) triangle.c `
                                             end)
                                    end), libtriangle)
end

@BinDeps.install Dict(:libtriangle => :libtriangle)
