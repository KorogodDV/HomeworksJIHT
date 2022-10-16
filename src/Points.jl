module Points

import LinearAlgebra
import Base: *, +, /

export Point, neighbors, Circle, Square, center

struct Point
    x::Real
    y::Real
end

Base.:+(a::Point, b::Point) = Point(a.x + b.x, a.y + b.y)
Base.:-(a::Point) = Point(-a.x, -a.y)
Base.:-(a::Point, b::Point) = a + (-b)
Base.:*(a::Point, b::Real) = Point(a.x * b, a.y * b)
Base.:/(a::Point, b::Real) = a * (1 / b)
Base.:*(b::Real, a::Point) = a * b

LinearAlgebra.dot(a::Point, b::Point) = a.x * b.x + a.y * b.y
LinearAlgebra.norm(a::Point) = sqrt(LinearAlgebra.dot(a, a))
inf_norm(p::Point) = max(abs(p.x), abs(p.y))


center(points) = sum(points) / length(points)

function neighbors(points, origin::Point, k::Integer)
    if k <= 0
        return Point[]
    elseif k >= length(points)
        return deleteat!(vec(points), findall(x -> x == origin, vec(points)))
    else
        distances = []
        res = Point[]
        for p in points
            if p != origin
                append!(distances, [LinearAlgebra.norm(p-origin)])
                append!(res, [p])
            end
        end
    end
    return res[sortperm(distances)[1:k]]
end

struct Circle
    o::Point
    radius::Real
end

Base.:in(p::Point, circ::Circle) = (LinearAlgebra.norm(p - circ.o) <= circ.radius)

struct Square
    o::Point
    side::Real
end

Base.:in(p::Point, sq::Square) = (inf_norm(p - sq.o) <= sq.side / 2)

function center(points, area)
    sum = Point(0, 0)
    points_in_area = 0
    for p in points
        if p in area
            sum += p
            points_in_area += 1
        end
    end
    return sum / points_in_area
end

end # module