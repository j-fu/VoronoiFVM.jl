module test_transientsol
using VoronoiFVM

function make_transientsol(;n=10,M=5,N=100,in_memory=true)
    makevec(k)=[ k+i*j for i=1:M,j=1:N]
    sol=TransientSolution(0,makevec(0),in_memory=in_memory)
    for k=1:n
        append!(sol,k,makevec(k))
    end
    sol
end

function test()
    msol=make_transientsol(in_memory=true)
    dsol=make_transientsol(in_memory=false)
    msol==dsol
end

end
