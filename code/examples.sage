from solver import find_solutions

def get_height_bound(field, places):
    return sum([place.degree() for place in places]) + field.genus() - 2

def print_example(field, places, height):
    print(f"field: {field}")
    print(f"places: {places}")
    print(f"height: {height}")
    print("------------------------------------------------")

    for unit_1, unit_2 in find_solutions(places, height):
        print(f"{unit_1} + {unit_2} = 1")
    print("------------------------------------------------")

K.<x> = FunctionField(QQbar)
O = K.maximal_order()
O_inf = K.maximal_order_infinite()
places =  [O.ideal(x).place(), O.ideal(x - 1).place(), O_inf.ideal(1/x).place()]
height = get_height_bound(K, places)

print_example(K, places, height)
print()

K.<x> = FunctionField(GF(5))
O = K.maximal_order()
O_inf = K.maximal_order_infinite()
places = [O.ideal(x).place(), O.ideal(x + 2).place(), O.ideal(x^2 + x + 1).place(), O_inf.ideal(1/x).place()]
height = get_height_bound(K, places)

print_example(K, places, height)
print()

K.<x> = FunctionField(GF(5))
_.<y> = K[]
L.<y> = K.extension(y^3 + x^3 * y + x)
places_0 = [place for (place, value) in L(y).divisor().list()]
places_1 = [place for (place, value) in L(y + 1).divisor().list()]
places = list(set(places_0 + places_1))
height = get_height_bound(L, places)

print_example(L, places, height)
