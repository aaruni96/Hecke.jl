
@testset "interval intersections" begin
  # fully disjoint
  i1 = interval(0, 1)
  i2 = interval(2, 3)
  @test intersect(i1, i2) == (false,)

  # mirror image disjoint
  i1 = interval(2, 3)
  i2 = interval(0, 1)
  @test intersect(i1, i2) == (false,)

  # simple intersections
  i1 = interval(0, 2)
  i2 = interval(1, 3)
  @test intersect(i1, i2) == (true, interval(1, 2))

  # other way round
  i1 = interval(1, 3)
  i2 = interval(0, 2)
  @test intersect(i1, i2) == (true, interval(1, 2))

  # one inside the other
  i1 = interval(0, 3)
  i2 = interval(1, 2)
  @test intersect(i1, i2) == (true, interval(1, 2))

  # other way round
  i1 = interval(1, 2)
  i2 = interval(0, 3)
  @test intersect(i1, i2) == (true, interval(1, 2))

end
