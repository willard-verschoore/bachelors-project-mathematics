from solver import find_solutions

K.<x> = FunctionField(QQbar)
O = K.maximal_order()
O_inf = K.maximal_order_infinite()
places =  [O.ideal(x).place(), O.ideal(x - 1).place(), O_inf.ideal(1/x).place()]
height = len(places) + K.genus() - 2

print("Algebraic field example:")
for unit_1, unit_2 in find_solutions(places, height):
    print(f"f = {unit_1}, g = {unit_2}")
print()

K.<x> = FunctionField(GF(5))
O = K.maximal_order()
O_inf = K.maximal_order_infinite()
places = [O.ideal(x).place(), O.ideal(x + 2).place(), O.ideal(x^2 + x + 1).place(), O_inf.ideal(1/x).place()]
height = sum([place.degree() for place in places]) + K.genus() - 2

print("Finite field example:")
for unit_1, unit_2 in find_solutions(places, height):
    print(f"f = {unit_1}, g = {unit_2}")
print()

K.<x> = FunctionField(QQbar)
_.<y> = K[]
L.<y> = K.extension(y^3 + x^3 * y + x)

places_0 = [place for (place, value) in L(y).divisor().list()]
places_1 = [place for (place, value) in L(y + 1).divisor().list()]
places = list(set(places_0 + places_1))
height = 3

print("Extended algebraic field example:")
for unit_1, unit_2 in find_solutions(places, height):
    print(f"f = {unit_1}, g = {unit_2}")
print()

K.<x> = FunctionField(GF(5))
_.<y> = K[]
L.<y> = K.extension(y^3 + x^3 * y + x)

places_0 = [place for (place, value) in L(y).divisor().list()]
places_1 = [place for (place, value) in L(y + 1).divisor().list()]
places = list(set(places_0 + places_1))
height = 3

print("Extended finite field example:")
for unit_1, unit_2 in find_solutions(places, height):
    print(f"f = {unit_1}, g = {unit_2}")
print()
