#=

# 115: 1D heterogeneous catalysis
 ([source code](SOURCE_URL))

Let $\Omega=(0,1)$, $\Gamma_1=\{0\}$, $\Gamma_2=\{1\}$
Regard a system of three species: $A,B,C$ and let 
$u_A=[A]$, $u_B=[B]$ and $u_C=[C]$ be their corresponding concentrations.

Species $A$ and $B$ exist in the interior of the domain, species $C$
lives a the boundary $\Gamma_1$.  We assume a heterogeneous reaction scheme
where $A$ reacts to $C$ and $C$ reacts to $B$:

```math
\begin{aligned}
      A &\leftrightarrow C\\
      C &\leftrightarrow B 
\end{aligned}
```
with reaction constants $k_{AC}^\pm$ and k_{BC}^\pm$.

In $\Omega$, both $A$ and $B$ are transported through diffusion:

```math
\begin{aligned}
\partial_t u_B - \nabla\cdot D_A \nabla u_A & = f_A\\
\partial_t u_B - \nabla\cdot D_B \nabla u_B & = 0\\
\end{aligned}
```
Here, $f(x)$ is a source term creating $A$.
On $\Gamma_2$, we set boundary conditions
```math
\begin{aligned}
D_A \nabla u_A & = 0\\
u_B&=0
\end{aligned}
```
describing no normal flux for $A$ and zero concentration of $B$.
On $\Gamma_1$, we use the mass action law to describe the boundary reaction and
the evolution of the boundary concentration $C$. We assume that there is a limited
amount of surface sites $S$ for species C, so in fact A has to react with a free
surface site in order to become $C$ which reflected by the factor $1-u_C$. The same
is true for $B$.
```math
\begin{aligned}
R_{AC}(u_A, u_C)&=k_{AC}^+ u_A(1-u_C) - k_{AC}^-u_C\\
R_{BC}(u_C, u_B)&=k_{BC}^+ u_B(1-u_C) - k_{BC}^-u_C\\
- D_A \nabla u_A  + S R_{AC}(u_A, u_C)& =0 \\
- D_B \nabla u_B  + S R_{BC}(u_B, u_C)& =0 \\
\partial_t C  - R_{AC}(u_A, u_C) - R_{BC}(u_B, u_C) &=0
\end{aligned}
```

=#

# Module header
module Example115_HeterogeneousCatalysis1D
using Printf
using VoronoiFVM

# Main function
function main(;n=10,Plotter=nothing,verbose=false,tend=1, dense=false)
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    N=length(X)
    
    grid=VoronoiFVM.Grid(X)
    ## By default, \Gamma_1 at X[1] and \Gamma_2 is at X[end]
    
    ## Species numbers
    iA=1
    iB=2
    iC=3


    ## Diffusion flux for species A and B
    D_A=1.0
    D_B=1.0e-2
    function flux!(f,u,edge)
        uk=viewK(edge,u)
        ul=viewL(edge,u)
        f[iA]=D_A*(uk[iA]-ul[iA])
        f[iB]=D_B*(uk[iB]-ul[iB])
    end

    ## Storage term of species A and B
    function storage!(f,u,node)
        f[iA]=u[iA]
        f[iB]=u[iB]
    end

    ## Source term for species a around 0.5
    function source!(f,node)
        x1=node[1]-0.5
        f[iA]=exp(-100*x1^2)
    end

    ## Reaction constants (p = + , m = -)
    ## Choosen to prefer path A-> C -> B
    ## More over, A reacts faster than to C than C to B
    ## leading to "catalyst poisoning", i.e. C taking up most of the
    ## available catalyst sites 
    kp_AC=100.0
    km_AC=1.0

    kp_BC=0.1
    km_BC=1.0

    S=0.01
    
    R_AC(u_A, u_C)=kp_AC*u_A*(1-u_C) - km_AC*u_C
    R_BC(u_B, u_C)=kp_BC*u_B*(1-u_C) - km_BC*u_C
    
    function breaction!(f,u,node)
        if  node.region==1
            f[iA]=S*R_AC(u[iA], u[iC])
            f[iB]=S*R_BC(u[iB], u[iC])
            f[iC]=-R_BC(u[iB], u[iC])-R_AC(u[iA], u[iC])
        end
    end

    ## This is for the term \partial_t u_C at the boundary
    function bstorage!(f,u,node)
        if  node.region==1
            f[iC]=u[iC]
        end
    end
    
    physics=VoronoiFVM.Physics(
        num_species=3,
        breaction=breaction!,
        bstorage=bstorage!,
        flux=flux!,
        storage=storage!,
        source=source!
    )
    
    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end

    ## Enable species in bulk resp
    enable_species!(sys,iA,[1])
    enable_species!(sys,iB,[1])

    ## Enable surface species
    enable_boundary_species!(sys,iC,[1])
    
    ## Set Dirichlet bc for species B on \Gamma_2
    boundary_dirichlet!(sys,iB,2,0.0)

    ## Initial values
    inival=unknowns(sys)
    inival.=0.0
    U=unknowns(sys)

    tstep=0.01
    time=0.0

    ## Data to store surface concentration vs time
    T=zeros(0)
    u_C=zeros(0)

    while time<tend
        time=time+tstep
        solve!(U,inival,sys,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        ## Record  boundary species
        push!(T,time)
        push!(u_C,U[iC,1])

        if isplots(Plotter)
            Plots=Plotter
            p1=Plots.plot(grid.coord[1,:],U[iA,:], grid=true, label="[A]")
            Plots.plot!(p1,grid.coord[1,:],U[iB,:], label="[B]",
                        title=@sprintf("max_A=%.5f max_B=%.5f u_C=%.5f\n",maximum(U[iA,:]),maximum(U[iB,:]),u_C[end]),
                        ylabel="[A], [B]", xlabel="x",legend=:topright,framestyle=:full)
            p2=Plots.plot(T,u_C,ylabel="[C]",xlabel="t",framestyle=:full, label="[C]")
            p=Plots.plot(p1,p2,layout=(2,1))
            Plots.gui(p)
        end
    end
    return U[iC,1]
end

# Test evaluation
function test()
    testval=0.87544440641274
    main(dense=false) ≈ testval && main(dense=true) ≈ testval
end

# Module end
end
