### A Pluto.jl notebook ###
# v0.17.0

using Markdown
using InteractiveUtils

# ╔═╡ 80263ff2-a7c5-42c1-bd3e-ee7799622f76
begin
    using Pkg
    Pkg.activate(".testenv")
    Pkg.add("Revise")
    using Revise
    Pkg.develop("VoronoiFVM")
    Pkg.add(["PlutoVista", "PlutoUI"])
end

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
    using PlutoVista
    using VoronoiFVM
	using PlutoUI
end

# ╔═╡ 5e13b3db-570c-4159-939a-7e2268f0a102
md"""
# Bernoulli function test

We test the implementation of the Bernoulli function in VoronoiFVM against the evaluation
with BigFloat. This allows to optimize thresholds for switching between evaluation expressions.
"""

# ╔═╡ 20ce5908-7fe7-4e8f-baf8-c858d2d12ce5
TableOfContents(title="")

# ╔═╡ 1774defe-1939-45f9-95e5-6a9fff70f3f2
md"""
## The candidates
"""

# ╔═╡ b47b781b-ec11-485a-9de6-061ad0957f46
"""
Benchmark implementation  using BigFloat
"""
function B_Big(x)
	bx=BigFloat(x)
	Float64(bx/(exp(bx)-one(bx)))
end

# ╔═╡ 01dd7870-45a7-42e5-ad26-7933bdec0c60
"""
Naive implementation using FLoat64
"""
B_Float64(x)=x/(exp(x)-1)

# ╔═╡ 36158a4c-4a4b-44c8-9226-448713304335
"""
Approximation for large positive and negative arguments
"""
B_large(x)= x>0 ? 0.0 : -x

# ╔═╡ 54de85ff-8d61-4d15-b7c6-939ce0ea10fc
"""
Horner scheme approximation for small  positive and negative arguments
"""
B_horner(x)=VoronoiFVM.bernoulli_horner(x)

# ╔═╡ 8e77f24b-b3fd-4069-82ba-799dc6aa602a
"""
Implementation  in VoronoiFVM
"""
B_vfvm(x)=VoronoiFVM.fbernoulli(x)

# ╔═╡ 984bd562-271b-4e7e-b289-95d95e92ec2b
md"""
This implementation uses `B_horner` if  `|exp(x)-1|<` $(VoronoiFVM.bernoulli_small_threshold) and 
`B_large` if `|x|>` $(VoronoiFVM.bernoulli_large_threshold).
"""

# ╔═╡ c3fd0ff2-7111-4165-ad93-d6d7257301fa
md"""
## Approximation for small x

For small values of x, the horner scheme approximation is used, as the exponential expression runs into errors in the vicinity of zero. As as long as its error is large than that of the Taylor approximation calculated with the Horner scheme, we should use the later one. 
"""

# ╔═╡ 56ff3f5c-6fe9-4d44-a5ae-449c42efca62
smallX=collect(-0.5:1.0e-4+1.0e-8:0.5);

# ╔═╡ 6e7c197b-8ad2-4b9d-a0bc-40a48db32387
let
	p=PlutoVistaPlot(resolution=(500,200),legend=:ct,xlabel="x",ylabel="error")
    plot!(p,smallX,abs.(B_Big.(smallX).-B_Float64.(smallX)),yscale=:log,label="|B_Big(x)-B_Float64(x)|")
    plot!(p,smallX,abs.(B_Big.(smallX).- B_horner.(smallX)),label="|B_Big(x)-B_horner(x)|")
end

# ╔═╡ 65aecf78-f88b-4399-abac-717c3c62a285
md"""
## Approximation for large x


Here, an important aspect is to prevent floating point overflow for large x.
"""

# ╔═╡ 26cdb920-291a-4b54-963f-fd9bd610662f
largeX=-100:1.00001e-3:100;

# ╔═╡ 8633969b-33cb-486a-8731-2ec7dcb881d7
let
	p=PlutoVistaPlot(resolution=(500,200),legend=:lt)
    plot!(p,largeX,abs.(B_Big.(largeX).-B_Float64.(largeX)) ,yscale=:log,label="|B_Big(x)-B_Float64(x)|")
    plot!(p,largeX,abs.(B_Big.(largeX).-B_large.(largeX)),label="|B_Big(x)-B_large(x)|")
end

# ╔═╡ 5a293797-beb9-493e-af12-d978c50d6148
md"""
## Test of the Implementation in VoronoiFVM
"""

# ╔═╡ ed8f172e-1d74-4514-b1e6-815ec9c87ae5
maxerror(X)=maximum(abs.(B_Big.(X).-B_vfvm.(X)));

# ╔═╡ feb21ce6-0ddc-45fb-90f4-1e46261a9110
plot(smallX,abs.(B_Big.(smallX).-B_vfvm.(smallX)),label="|B_Big(x)-B_vfvm(x)|",legend=:lt,yscale=:log,resolution=(500,200),color=:red,xlabel="x",ylabel="error", title="Max error for small x: $(maxerror(smallX))" )

# ╔═╡ 91d3b907-9053-4467-a8ab-be9c5597741a
plot(largeX,abs.(B_Big.(largeX).-B_vfvm.(largeX)),label="|B_Big(x)-B_vfvm(x)|",legend=:lb,yscale=:log,resolution=(500,200),color=:red,xlabel="x",ylabel="error", title="Max error for large x: $(maxerror(largeX))" )

# ╔═╡ 12f268d2-baa6-4c1b-bed5-e9df53b469fc
html"<hr>"

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─5e13b3db-570c-4159-939a-7e2268f0a102
# ╠═20ce5908-7fe7-4e8f-baf8-c858d2d12ce5
# ╠═1774defe-1939-45f9-95e5-6a9fff70f3f2
# ╠═b47b781b-ec11-485a-9de6-061ad0957f46
# ╠═01dd7870-45a7-42e5-ad26-7933bdec0c60
# ╠═36158a4c-4a4b-44c8-9226-448713304335
# ╠═54de85ff-8d61-4d15-b7c6-939ce0ea10fc
# ╠═8e77f24b-b3fd-4069-82ba-799dc6aa602a
# ╠═984bd562-271b-4e7e-b289-95d95e92ec2b
# ╠═c3fd0ff2-7111-4165-ad93-d6d7257301fa
# ╠═56ff3f5c-6fe9-4d44-a5ae-449c42efca62
# ╠═6e7c197b-8ad2-4b9d-a0bc-40a48db32387
# ╠═65aecf78-f88b-4399-abac-717c3c62a285
# ╠═26cdb920-291a-4b54-963f-fd9bd610662f
# ╠═8633969b-33cb-486a-8731-2ec7dcb881d7
# ╠═5a293797-beb9-493e-af12-d978c50d6148
# ╠═ed8f172e-1d74-4514-b1e6-815ec9c87ae5
# ╠═feb21ce6-0ddc-45fb-90f4-1e46261a9110
# ╠═91d3b907-9053-4467-a8ab-be9c5597741a
# ╟─12f268d2-baa6-4c1b-bed5-e9df53b469fc
# ╠═80263ff2-a7c5-42c1-bd3e-ee7799622f76
