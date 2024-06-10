class BoundedValues:
    def __init__(self, bound, value_count):
        self._bound = bound
        self._value_count = value_count

        self._current_bound = 0
        self._positive_partitions = Partitions(0)
        self._positive_partition = 0
        self._negative_partitions = Partitions(0)
        self._negative_partition = 0
        self._permutations = Permutations([])
        self._permutation = 0

    def __iter__(self):
        return self

    def __next__(self):
        self.__next_permutation()
        return self._permutations[self._permutation]

    def __next_bound(self):
        if self._current_bound + 1 > self._bound:
            raise StopIteration

        self._current_bound += 1
        self.__reset_positive_partitions()

    def __next_positive_partition(self):
        if self._positive_partition + 1 >= self._positive_partitions.cardinality():
            self.__next_bound()
            return

        self._positive_partition += 1
        self.__reset_negative_partitions()

    def __next_negative_partition(self):
        if self._negative_partition + 1 >= self._negative_partitions.cardinality():
            self.__next_positive_partition()
            return

        self._negative_partition += 1
        self.__reset_permutations()

    def __next_permutation(self):
        if self._permutation + 1 >= self._permutations.cardinality():
            self.__next_negative_partition()
            return

        self._permutation += 1

    def __reset_positive_partitions(self):
        self._positive_partitions = Partitions(self._current_bound, max_length=self._value_count-1)
        self._positive_partition = 0
        self.__reset_negative_partitions()

    def __reset_negative_partitions(self):
        remaining = self._value_count - len(self._positive_partitions[self._positive_partition])
        self._negative_partitions = Partitions(self._current_bound, max_length=remaining)
        self._negative_partition = 0
        self.__reset_permutations()

    def __reset_permutations(self):
        positive_values = self._positive_partitions[self._positive_partition]
        negative_values = [-v for v in self._negative_partitions[self._negative_partition]]
        remaining = self._value_count - len(positive_values) - len(negative_values)
        values = positive_values + negative_values + [0] * remaining
        self._permutations = Permutations(values)
        self._permutation = 0

def compute_field_element(values, places):
    divisor = sum([value // place.degree() * place for value, place in zip(values, places)])
    basis = divisor.basis_function_space()
    if len(basis) == 1: return 1 / basis[0]

def generate_candidates(places, height):
    for values in BoundedValues(height, len(places)):
        for value, place in zip(values, places):
            if value % place.degree() != 0: break
        else:
            candidate = compute_field_element(values, places)
            if candidate is not None:
                yield candidate

def find_solutions(places, height):
    canditates = generate_candidates(places, height)
    for [x, y] in Combinations(canditates, 2):
        ratio = x.derivative() / y.derivative()
        if not ratio.degree() == 0: continue

        determinant = x * y.derivative() - y * x.derivative()
        if determinant == 0: continue

        c_1 = y.derivative() / determinant
        c_2 = -x.derivative() / determinant
        if not c_1.degree() == 0 or not c_2.degree() == 0: continue

        unit_1 = c_1 * x
        unit_2 = c_2 * y

        yield unit_1, unit_2

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
