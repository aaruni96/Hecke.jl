# using Hecke

"""
Aaruni's special intervals for packing rectangles in squares
"""
struct interval
  a::Rational
  b::Rational
  interval(a, b) = a > b ? error("out of order") : new(a, b)
end

struct rectangle
  i1::interval    # interval on X axis
  i2::interval    # interval on Y axis

  function rectangle(
    a::Rational,
    b::Rational,
    c::Rational,
    d::Rational
  )
    @assert a<b "out of order"
    @assert c<d "out of order"
    return new(interval(a,b), interval(c,d))
  end

  function rectangle(
    k::Integer,
    x::Rational=0//1,
    y::Rational=0//1,
    r::Bool=false
  )
    i1 = interval(0, 1//(k+1))
    i2 = interval(0, 1//k)
    if r
      i1, i2 = i2, i1
    end
    i1 = i1 + x
    i2 = i2 + y
    return new(i1, i2)
  end

#  function rectangle(k::Integer, x::Rational=0, y::Rational=0, r::Bool=false)
#    println(k)
#    println(x)
#    println(y)
#    println(r)
#    return new(k, x//1, y//1, r)
#  end

  function rectangle(k::Integer, r::Bool=false)
    return rectangle(k, 0//1, 0//1, r)
  end

  rectangle(i1::interval, i2::interval) = new(i1, i2)
end

function Base.:+(i::interval, x::Rational)
  return interval(i.a+x, i.b+x)
end

#function Base.:+(i::interval, x::Rational)
#  return i+x
#end

function Base.:+(x::Rational, i::interval)
  return i + x
end

#function Base.:+(x::Rational, i::interval)
#  return i + x
#end

function midpoint(i::interval)
  return (i.a + i.b ) // 2
end

"""
test docstring
"""
function length(i::interval)
  return i.b - i.a
end

function coordinates(r::rectangle)
  return [(r.i1.a,r.i2.a), (r.i1.b, r.i2.a), (r.i1.a, r.i2.b), (r.i1.b, r.i2.b)]
end

"""
test test test
"""
function intersect(a::interval, b::interval)
  try
    i = interval(maximum([a.a, b.a]), minimum([a.b, b.b]))
    return (true, i)
  catch
    return (false, nothing)
  end
end

function intersect(r1::rectangle, r2::rectangle)
  i1 = intersect(r1.i1, r2.i1)
  i2 = intersect(r1.i2, r2.i2)
  print(i1)
  print(i2)
  if !i1[1] || !i2[1]
    return (false, nothing)
  end
  return (true, rectangle(i1[2], i2[2]))
end

function rectarray()
  return [
    rectangle(1),
    rectangle(2, 1//2, 0//1, true),
    rectangle(3),
    rectangle(4)
  ]
end

#=
if :Plots in names(Main, imported=true)

  # plot(bar([1//2,1//4, 3//4], [1,1,1], bar_width=[1,0.5, 0.5], color=collect(colors), fillto = [0,0,0.5], alpha = [0.1, 0.8, 0.8]))

  colors = collect(keys(Plots.Colors.color_names))

  plot(
    bar(
      Rational.([midpoint(r.i1) for r in rectarray]), #midpoint of bar
      Rational.([r.i2.b for r in rectarray]), # height of bar
      bar_width = Rational.([length(r.i1) for r in rectarray]),
      color = colors,
      fillto = Rational.([r.i2.a for r in rectarray]),
      alpha = [0.8 for r in rectarray]
    )
  )

end
=#
