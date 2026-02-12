using Hecke

export interval

"""
Aaruni's special intervals for packing rectangles in squares
"""
struct interval
  a::QQFieldElem
  b::QQFieldElem
end

struct rectangle
  i1::interval    # interval on X axis
  i2::interval    # interval on Y axis
end

r1 = rectangle(interval(0, 3//10), interval(8//10, 1))


"""
test docstring
"""
function length(i::interval)
  return i.b - i.a
end

function product(a::interval, b::interval)

end


function intersect(a::interval, b::interval)

end
