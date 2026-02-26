using Hecke

"""
Aaruni's special intervals for packing rectangles in squares
"""
struct interval
  a::QQFieldElem
  b::QQFieldElem
  interval(a, b) = a > b ? error("out of order") : new(a, b)
end

struct rectangle
  i1::interval    # interval on X axis
  i2::interval    # interval on Y axis

  function rectangle(
    a::QQFieldElem,
    b::QQFieldElem,
    c::QQFieldElem,
    d::QQFieldElem
  )
    @assert a<b "out of order"
    @assert c<d "out of order"
    return new(interval(a,b), interval(c,d))
  end

  function rectangle(
    a::Rational,
    b::Rational,
    c::Rational,
    d::Rational
  )
    return rectangle(QQ(a), QQ(b), QQ(c), QQ(d))
  end

  function rectangle(     ##
    k::Integer,
    x::QQFieldElem=QQ(0),
    y::QQFieldElem=QQ(0),
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

  #=
  function rectangle(     ##
    k::Integer,
    x::Rational=0//1,
    y::Rational=0//1,
    r::Bool=false
  )
    println(k)
    println(x)
    println(y)
    println(r)
    return rectangle(k, QQ(x), QQ(y), r)
  end
  =#

#  function rectangle(k::Integer, r::Bool=false)
#    return rectangle(k, QQ(0), QQ(0), r)
#  end

  rectangle(i1::interval, i2::interval) = new(i1, i2)
end

function Base.:+(i::interval, x::QQFieldElem)
  return interval(i.a+x, i.b+x)
end

function Base.:+(i::interval, x::Rational)
  return i+QQ(x)
end

function Base.:+(x::Rational, i::interval)
  return i + QQ(x)
end

function Base.:+(x::QQFieldElem, i::interval)
  return i + x
end

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
function intersect(a::interval, b::interval) :: Union{interval, Bool}
  try
    i = interval(maximum([a.a, b.a]), minimum([a.b, b.b]))
    return i
  catch
    return false
  end
end

function intersect(r1::rectangle, r2::rectangle) :: Union{rectangle, Bool}
  i1 = intersect(r1.i1, r2.i1)
  i2 = intersect(r1.i2, r2.i2)
  if typeof(i1) == Bool || typeof(i2) == Bool
    return false
  end
  return rectangle(i1, i2)
end
#export interval

function rectarray()
  return [
    rectangle(1),
    rectangle(2, QQ(1//2), QQ(0), true),
    rectangle(3),
    rectangle(4)
  ]
end

export interval
export rectangle
export intersect
export rectarray

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
