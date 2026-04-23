module MyModule
using Hecke
using Plots

"""
One dimensional intervals in QQ.
"""
struct interval
  a::QQFieldElem
  b::QQFieldElem
  interval(a, b) = a > b ? error("out of order") : new(a, b)

  function interval(a, b)
    if a > b
      println("Problem while contsructing interval")
      @show a
      @show b
      error("Out of order error! a must always be smaller than b!")
    else
      return new(a,b)
    end
  end
end # end struct

function Base.show(io::IO, mime::MIME"text/plain", i::interval)
  print(io, "[$(i.a) $(i.b)]")
end

function Base.:(==)(i1::interval, i2::interval)
  if i1.a == i2.a && i1.b == i2.b
    return true
  end
  return false
end

"""
A rectangle, expressed as a product of two intervals. Interval i1 is the interval on X-axis.
Interval i2 is the interval on y-axis.
"""
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

  """
  Returns a rectangle of sides 1/k, and 1/(k+1), shifted by x,y.
  r = false is a verticle rectangle
  r = true is a horizontal rectangle
  """
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

    function rectangle(
      k::Integer,
      x::Union{Rational, Integer},
      y::Union{Rational, Integer},
      r::Bool = false
    )
      return rectangle(k, QQ(x), QQ(y), r)
    end

  rectangle(i1::interval, i2::interval) = new(i1, i2)
end #end of rectangle type

function Base.show(io::IO, r::rectangle)
  print(io, "[$(r.i1.a) $(r.i1.b)][$(r.i2.a) $(r.i2.b)]")
end

function Base.show(io::IO, mime::MIME"text/plain", r::rectangle)
    print(io, "[$(r.i1.a) $(r.i1.b)][$(r.i2.a) $(r.i2.b)]")
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

function Base.:+(r::rectangle, xy::Tuple{Number, Number})
  return rectangle(r.i1 + QQ(xy[1]), r.i2 + QQ(xy[2]))
end

function Base.:+(r::rectangle, xy::Tuple{QQFieldElem, QQFieldElem})
  return rectangle(r.i1 + xy[1], r.i2 + xy[2])
end

function Base.:+(r1::rectangle, r2::rectangle )
  return r1 + (r2.i1.a, r2.i2.a)
end

function midpoint(i::interval)
  return (i.a + i.b ) // 2
end

"""
Checks if i2 sits perfectly inside i1, edge overlapping allowed.
"""
function contains(i1::interval, i2::interval)
  if i2.a >= i1.a
    if i2.b <= i1.b
      return true
    end
  end
  return false
end

"""
Checks if r2 sits perfectly inside r1, edge overlapping allowed.
"""
function contains(r1::rectangle, r2::rectangle)
  if contains(r1.i1, r2.i1)
    if contains(r1.i2, r2.i2)
      return true
    end
  end
  return false
end

"""
test docstring
"""
function Base.length(i::interval)
  return i.b - i.a
end

function Base.length(r::rectangle)
  return maximum(length.([r.i1, r.i2]))
end

function width(r::rectangle)
  return minimum(length.([r.i1, r.i2]))
end

function coordinates(r::rectangle)
  return [(r.i1.a,r.i2.a), (r.i1.b, r.i2.a), (r.i1.a, r.i2.b), (r.i1.b, r.i2.b)]
end

"""
test test test
"""
function intersect(a::interval, b::interval)
  il = maximum([a.a, b.a])
  iu = minimum([a.b, b.b])
  if il <= iu
    i = interval(maximum([a.a, b.a]), minimum([a.b, b.b]))
    if (length(i) == 0)
      return (false, )
    end
    return (true, i)
  else
    return (false, )
  end
end

function area(r::rectangle)
  return length(r.i1) * length(r.i2)
end

"""
changed docstring
"""
function intersect(r1::rectangle, r2::rectangle)
  i1 = intersect(r1.i1, r2.i1)
  i2 = intersect(r1.i2, r2.i2)
  if ! i1[1] || ! i2[1]
    return (false, )
  end
  return (true, rectangle(i1[2], i2[2]))
end

function rectarray()
  return [
    rectangle(1),
    rectangle(2, QQ(1//2), QQ(0), true),
    rectangle(3),
    rectangle(4)
  ]
end


"""
r1 is a big rectangle
r2 is a small rectangle

cut out an r2 shaped hole into r1
"""
function Base.:-(r1::rectangle, r2::rectangle)
  i = intersect(r1, r2)
  if !i[1]
    error("Cannot subtract non intersecting rectanlges!")
  end
  if i[2].i1 == r1.i1
    ni1 = r1.i1
    ni2 = interval(0, 0)
    try
      ni2 = interval(r2.i2.b, r1.i2.b)
    catch
      error("Can only subtract a small rectangle from a larger rectangle")
    end
    return rectangle(ni1, ni2)
  elseif i[2].i2 == r1.i2
    ni2 = r1.i2
    ni1 = interval(0,0)
    try
      ni1 = interval(r2.i1.b, r1.i1.b)
    catch
      error("Can only subtract a small rectangle from a larger rectangle")
    end
    return rectangle(ni1, ni2)
  else
    error(
      "for our purposes, one of the above two must always be true (even though this is not true",
      "in general)"
    )
  end
end

function add_to_list!(newrect::rectangle, rectlist::Vector{rectangle})
  for r in rectlist
    if intersect(newrect, r)[1]
      error("could not tile without intersection!")
    end
  end
  push!(rectlist, newrect)
end

function place_next(rects::Vector{rectangle}, boxes::Vector{rectangle})
    k = length(rects) + 1

    rn = rectangle(k, QQ(0), QQ(0), false)

    t = rectangle(0)
    j = -1
    # bad rule 1: place in first box\
    for i in 1:length(boxes)
      #check if rectangle fits in boxes
      rnk = rn+(boxes[i].i1.a, boxes[i].i2.a)
      if contains(boxes[i], rnk)
        t = boxes[i]
        j = i
        break
      end
    end
    if j == -1
      println("=====================")
      println("Could not fit rectnew in a box")
      println("=====================")
      error("fail here")
      return rects, boxes
      #error("rectangle didn't fit in any box")
    end
    newrects = deepcopy(rects)
    rnt = rn+t
    add_to_list!(rnt, newrects)

    # questionable rule2: update the boxlist
    # split the box into two at some point in some direction
    newboxes = deepcopy(boxes)
    c = coordinates(rn+t)[4]
    d1 = 1-c[1]
    d2 = 1-c[2]
    if d1 < d2 # split across x axis
      t1, t2 = rectangle(t.i1, interval(t.i2.a, c[2])), rectangle(t.i1, interval(c[2], t.i2.b))
    else  #split across y axis
      t1, t2 = rectangle(interval(t.i1.a, c[1]), t.i2), rectangle(interval(c[1], t.i1.b), t.i2)
    end

    # subtract the latest rectangle from one of the two new boxes
    if intersect(rnt, t1)[1]
      t1 = t1 - rnt
    elseif intersect(rnt, t2)[1]
      t2 = t2 - rnt
    else
      error("impossible case is an error")
    end
    newboxes[j] = t1
    push!(newboxes, t2)

    # delete any empty boxes
    deleteat!(newboxes, area.(newboxes) .== 0)

    return newrects, newboxes
end

function place_first()
    rects = [
        rectangle(1, 0, 0, false),
        rectangle(2, 1//2, 0, true)
        #rectangle(2, 1//2, 0, false)
        ]

    boxes = [
        rectangle(1//2, 1//1, 1//3, 1//1)
        #rectangle(1//2+1//3, 1//1, 0//1, 1//2),
        #rectangle(1//2, 1//1, 1//2, 1//1)
    ]
    return (rects, boxes)
end

"""
test
"""
function plot_stuff(rects)

  if :Plots in names(Main, imported=true)
    println("plotting")

    # plot(bar([1//2,1//4, 3//4], [1,1,1], bar_width=[1,0.5, 0.5], color=collect(colors), fillto = [0,0,0.5], alpha = [0.1, 0.8, 0.8]))

    colors = collect(keys(Plots.Colors.color_names))
    prects = deepcopy(rects)
    push!(prects, rectangle(0//1,1//1,0//1,1//1))
    plot(
      bar(
        Rational.([midpoint(r.i1) for r in rects]), #midpoint of bar
        Rational.([r.i2.b for r in rects]), # height of bar
        bar_width = Rational.([length(r.i1) for r in rects]),
        color = colors,
        fillto = Rational.([r.i2.a for r in rects]),
        alpha = [0.8 for r in rects]
      ),
      show=true,
      reuse=true
    )
  else
    println("Not plotting")
  end

end


export rectangle
export interval
export place_first
export place_next
export plot_stuff
export area



end

# using .MyModule

#=

rects, boxes = place_first()
newrects, newboxes = deepcopy(rects), deepcopy(boxes)

i = 3
while true
  try
    newrects, newboxes = place_next(newrects, newboxes)
    plot_stuff(newrects)
    i = i+1
    println(i)
    sleep(0.1)
  catch
    println("Reached a point where we can no longer stupidly fit things!")
    @show i
    break
  end
end




=====================
Could not fit rectnew in a box
=====================
Reached a point where we can no longer stupidly fit things!
i = 757

=#

