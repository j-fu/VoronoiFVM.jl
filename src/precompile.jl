## Attempt on https://discourse.julialang.org/t/22-seconds-to-3-and-now-more-lets-fix-all-of-the-differentialequations-jl-universe-compile-times/66313
let
    while false
        # Problem similar to:  https://github.com/JuliaLang/julia/issues/35972
        # If we enable this, code in the tests starts allocating
        # https://github.com/jkrumbiegel/GridLayoutBase.jl/pull/21#discussion_r684676916
        function lin1()
            n=5
            X=0:1.0/n:1
            grid=simplexgrid(X,X)
            ispec=1
            
            function diffusion_flux(D::T) where T
                (y,u,edge)-> y[1]=D(u[1,1]+u[1,2])*(u[1,1]-u[1,2])
            end
            
            flux!(y,u,edge)= y[1]= u[1,1]-u[1,2]
            reaction!(y,u,node)=  y[1]= u[1,1]^3-10
            function bc!(args...)
                boundary_dirichlet!(args...,region=1,value=0)
                boundary_dirichlet!(args...,region=3,value=1)
            end
            
            sys=VoronoiFVM.System(grid; species=1,flux=flux!,reaction=reaction!,breaction=bc!,check_allocs=false)
            solution=solve(sys; verbose=false)
            return sum(solution)
        end
        lin1()
        break
    end
end
