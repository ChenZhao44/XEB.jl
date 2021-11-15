using Test
using XEB

@testset "Rule-1" begin
    g = google_layout_53(2, 20)
    add_noise!(g, 1, 7)
    add_noise!(g, 1, 9)
    update_noise!(g)
    @test all(isequal(:noise), gates(g, 1)) && all(isequal(:noise), gates(g, 2))
end

@testset "Rule-2" begin
    g = google_layout_53(2, 20)
    add_noise!(g, 1, 7)
    add_noise!(g, 2, 7)
    update_noise!(g)
    @test all(isequal(:noise), gates(g, 1)) && all(isequal(:noise), gates(g, 2))
end

g = google_layout_53(2, 20)
add_noise!(g, 1, 7)
add_noise!(g, 1, 9)
update_noise!(g)
@test all(isequal(:noise), gates(g, 1)) && all(isequal(:noise), gates(g, 2))
XEB.update_id!(g)
gates(g, 2)