using Test

@time begin
    @time begin
        print("             iliq:")
        include("../examples/iliq.jl")
        @test run_iliq(n=20,pyplot=false,dlcap=false,verbose=false) ≈ 0.9999546021312723
        @test run_iliq(n=20,pyplot=false,dlcap=true,verbose=false) ≈ .010759276468375045
    end

    @time begin
        print("  nlpoisson-1spec:")
        include("../examples/nlpoisson-1spec.jl")
        @test run_nlpoisson_1spec(n=10,verbose=false) ≈ 0.3371249631439964
    end
    @time begin
        print("  nlpoisson-2spec:")
        include("../examples/nlpoisson-2spec.jl")
        @test run_nlpoisson_2spec(n=10,verbose=false) ≈ 0.7117546972922056
    end
    @time begin
        print("           test2d:")
        include("../examples/test2d.jl")
        @test run_test2d(n=10,verbose=false) ≈ 0.3554284760906605
    end
    @time begin
        print("       test2d-rea:")
        include("../examples/test2d-rea.jl")
        @test run_test2d_rea(n=10,verbose=false) ≈ 0.014566189535134827
    end
    @time begin
        print("      test2d-brea:")
        include("../examples/test2d-brea.jl")
        @test run_test2d_brea(n=10,verbose=false) ≈ 0.008761335823958986
    end
    @time begin
        print("test2d-brea-bspec:")
        include("../examples/test2d-brea-bspec.jl")
        @test run_test2d_brea_bspec(n=10,verbose=false) ≈ 0.0020781361856598
    end
print("            all:")
end







