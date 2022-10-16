using ScoredTests
using Points
using Test
import LinearAlgebra


const AWARD1 = 1
const AWARD2 = 2

ScoredTests.DefaultScoring.award = AWARD1
ScoredTests.DefaultScoring.penalty = 0


maints = ScoredTestSet("Реализация модуля Points")

# Конструкторы struct Point
maints *= let ts = ScoredTestSet("Конструкторы struct Point")
    ts *= @scoredtest Point(1, 2) isa Point name="Point(::Int, ::Int)"
    ts *= @scoredtest Point(1.0, 2.0) isa Point name="Point(::Float64, ::Float64)"
    ts *= @scoredtest Point(1, 2.5) isa Point name="Point(::Int, ::Float64)"
    ts *= @scoredtest Point(1.0, 2) isa Point name="Point(::Float64, ::Int)"
    ts *= @scoredtest Point(1//1, 2) isa Point name="Point(::Rational, Int)" award=AWARD2
end

# Линейные операции
maints *= let ts = ScoredTestSet("Линейные операции")
    ts *= @scoredtest Point(1, 2) + Point(4, 5) == Point(5, 7) name="p + q"
    ts *= @scoredtest Point(1, 2) - Point(4, 5) == Point(-3, -3) name="p - q"
    ts *= @scoredtest -Point(1, 2) == Point(-1, -2) name="-p"
    ts *= @scoredtest 2.5 * Point(1, 2) == Point(2.5, 5.0) name="α * p"
    ts *= @scoredtest -2.5 * Point(1, 2) == Point(-2.5, -5.0) name="α * p"
    ts *= @scoredtest Point(1, 2) * 2.5 == Point(2.5, 5.0) name="p * α"
    ts *= @scoredtest Point(1, 2) / 2 == Point(0.5, 1.0) name="p / α"
end

maints *= let ts = ScoredTestSet("Расширения для LinearAlgebra")
    ts *= @scoredtest LinearAlgebra.norm(Point(1, 1)) ≈ √2 name="Норма"
    ts *= @scoredtest LinearAlgebra.dot(Point(1, 2), Point(2, 3)) == 8 name="Скалярное произведение"
end

maints *= let ts = ScoredTestSet("Конструкторы Circle")
    ts *= @scoredtest Circle(Point(0, 0), 2) isa Circle name="Circle(::Point{Int}, ::Int)"
    ts *= @scoredtest Circle(Point(0, 0), 2.0) isa Circle name="Circle(::Point{Int}, ::Float64)" award=AWARD2
end

maints *= let ts = ScoredTestSet("Точка в круге (расширение Base)")
    ts *= @scoredtest Point(1, 1) in Circle(Point(0, 0), 2)
    ts *= @scoredtest Point(1, 0) in Circle(Point(2.0, 0.0), nextfloat(1.0))
end

maints *= let ts = ScoredTestSet("Конструкторы Square")
    ts *= @scoredtest Square(Point(0, 0), 2) isa Square name="Square(::Point{Int}, ::Int)"
    ts *= @scoredtest Square(Point(0, 0), 2.0) isa Square name="Square(::Point{Int}, ::Float64)" award=AWARD2
end

maints *= let ts = ScoredTestSet("Точка в квадрате (расширение Base)")
    ts *= @scoredtest Point(1, 1) in Square(Point(0, 0), 2)
    ts *= @scoredtest !(Point(1, 1) in Square(Point(2.0, 2.0), prevfloat(2.0)))
end

maints *= let ts = ScoredTestSet("Центр масс")
    ts *= @scoredtest center(i * Point(1, 1) for i in 1:10) == Point(5.5, 5.5) name="Без указания области" award=AWARD2

    try
        points = [Point(x, y) for x in -2:2, y in -2:2]
        s = Square(Point(0, 0) + Point(0.5, 0.5), 2.0)
        ts *= @scoredtest center(points, Circle(Point(1, 2), 1)) == Point(1.0, 1.75) name="Лежащие в круге" award=AWARD2
        ts *= @scoredtest center(points, s) == Point(0.5, 0.5) name="Лежащие в квадрате" award=AWARD2
    catch
        ts *= @scoredtest error("Не реализованы операции") name="Лежащие в круге"
        ts *= @scoredtest error("Не реализованы операции") name="Лежащие в квадрате"
    end
    ts
end

maints *= let ts = ScoredTestSet("k-Ближайших соседей")
    try
        points = [Point(x, y) for x in -2:2, y in -2:2]
        origin = Point(0, 0)
        ts *= @scoredtest isempty(neighbors(points, origin, 0)) award=AWARD2
        ts *= @scoredtest isempty(neighbors(points, origin, -1)) award=AWARD2

        n4 = Set(Point(x, y) for x in -1:1, y in -1:1 if  abs(x + y) == 1)
        ts *= @scoredtest first(neighbors(points, origin, 1)) in n4 award=AWARD2
        ts *= @scoredtest all(x -> x in n4, neighbors(points, origin, 4)) award=AWARD2

        ts *= @scoredtest all(x -> x in setdiff(Set(points[2:4, 2:4]), [origin]), neighbors(points, origin, 8)) award=AWARD2
        ts *= @scoredtest Set(neighbors(points, origin, 2 * length(points))) == setdiff(Set(points), [origin]) award=AWARD2
    catch
        for _ in 1:6
            ts *= @scoredtest error("Не реализованы операции") award=AWARD2
        end
    end
    ts
end

maints *= let ts = ScoredTestSet("Дополнительно")
    try
        points = [Point(x, y) for x in -2:2, y in -2:2]
        origin = points[3, 3]
        try
            for k in 1:length(points)+1
                @inferred neighbors(points, origin, k)
            end
            ts *= @scoredtest true award=4 name="Стабильность функции neighbors"
        catch
            ts *= @scoredtest false name="Стабильность функции neighbors"
        end

        try
            @inferred center(points, Circle(Point(1, 2), 1))
            ts *= @scoredtest true award=AWARD2 name="Стабильность функции center для круга"
        catch
            ts *= @scoredtest false name="Стабильность функции center для круга"
        end

        try
            @inferred center(points, Square(Point(0, 0) + Point(0.5, 0.5), 2.0))
            ts *= @scoredtest true award=AWARD2 name="Стабильность функции center для квадрата"
        catch
            ts *= @scoredtest false name="Стабильность функции center для квадрата"
        end
    catch
        ts *= @scoredtest error("Не реализованы операции") award=4
        ts *= @scoredtest error("Не реализованы операции") award=AWARD2
        ts *= @scoredtest error("Не реализованы операции") award=AWARD2
    end
    ts
end

printsummary(maints)
