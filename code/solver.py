from sage.all import *


class BoundedValues:
    """
    A generator for all sequences of values for the valuations of an S-unit.
    """

    def __init__(self, height_bound, value_count):
        """
        Initialize the generator.

        :param bound: The maximum height of the S-unit.
        :param value_count: The number of places in the set S.
        """

        self._height_bound = height_bound
        self._value_count = value_count

        self._positive_partitions = Partitions(height_bound, max_length=value_count - 1)
        self._positive_partition = 0
        self._negative_partitions = Partitions(height_bound, max_length=value_count - 1)
        self._negative_partition = 0
        self._permutations = Permutations(
            [height_bound, -height_bound] + [Integer(0)] * (value_count - 2)
        )
        self._permutation = -1

    def __iter__(self):
        return self

    def __next__(self):
        self.__next_permutation()
        return [x for x in self._permutations[self._permutation]]

    def __next_positive_partition(self):
        if self._positive_partition + 1 >= self._positive_partitions.cardinality():
            raise StopIteration

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

    def __reset_negative_partitions(self):
        remaining = self._value_count - len(
            self._positive_partitions[self._positive_partition]
        )
        self._negative_partitions = Partitions(self._height_bound, max_length=remaining)
        self._negative_partition = 0
        self.__reset_permutations()

    def __reset_permutations(self):
        positive_values = self._positive_partitions[self._positive_partition]
        negative_values = [
            -v for v in self._negative_partitions[self._negative_partition]
        ]
        remaining = self._value_count - len(positive_values) - len(negative_values)
        values = positive_values + negative_values + [Integer(0)] * remaining
        self._permutations = Permutations(values)
        self._permutation = 0


def compute_field_element(places, values):
    """
    Compute a field element corresponding to the divisor given by a set of places and
    coefficients. If the divisor is not principal, return None.

    :param places: The places in the support of the divisor.
    :param values: The coefficients of the divisor, multiplied by the degree of the
                   corresponding place.
    """

    divisor = sum(
        [value // place.degree() * place for value, place in zip(values, places)]
    )
    basis = divisor.basis_function_space()
    if len(basis) == 1:
        return 1 / basis[0]


def generate_candidates(places, height):
    """
    Generate all S-units of a given height.

    :param places: The places in the set S.
    :param height: The height of the S-units.
    """

    for values in BoundedValues(height, len(places)):
        for value, place in zip(values, places):
            if value % place.degree() != 0:
                break
        else:
            candidate = compute_field_element(places, values)
            if candidate is not None:
                yield candidate, values


def compare(vx, vy):
    """
    Compare two sequences of integers in lexicographic order.

    :param vx: The first sequence.
    :param vy: The second sequence.
    """

    for i in range(len(vx)):
        if vx[i] < vy[i]:
            return -1
        if vx[i] > vy[i]:
            return 1
    return 0


def is_first_pair(vx, vy):
    """
    Check if a pair of S-units is the first of its group of conjugates. Ordering is
    based on the lexicographic order of the sequence of values v_P * deg(P) for each
    place P in the ordered set S.

    :param vx: The sequence of values of the first S-unit.
    :param vy: The sequence of values of the second S-unit.
    """

    if compare(vx, vy) > 0:  # Ensure consistent ordering of pairs.
        vx, vy = vy, vx

    # Compare to first alternative.
    vx_alt = [-vx[i] for i in range(len(vx))]
    vy_alt = [vy[i] - vx[i] for i in range(len(vx))]
    if compare(vx_alt, vy_alt) > 0:  # Ensure consistent ordering of pairs.
        vx_alt, vy_alt = vy_alt, vx_alt
    if compare(vx + vy, vx_alt + vy_alt) > 0:
        return False

    # Compare to second alternative.
    vx_alt = [-vy[i] for i in range(len(vx))]
    vy_alt = [vx[i] - vy[i] for i in range(len(vx))]
    if compare(vx_alt, vy_alt) > 0:  # Ensure consistent ordering of pairs.
        vx_alt, vy_alt = vy_alt, vx_alt
    if compare(vx + vy, vx_alt + vy_alt) > 0:
        return False

    return True


def find_solutions(places, height_bound):
    """
    Find all pairs of S-units up to a given height that solve the S-unit equation.

    :param places: The places in the set S.
    :param height_bound: The maximum height of the S-units.
    """

    for height in range(1, height_bound + 1):
        canditates = generate_candidates(places, height)
        for [(x, vx), (y, vy)] in Combinations(canditates, 2):
            # Only consider one of every group of six conjugate pairs.
            if not is_first_pair(vx, vy):
                continue

            # In order to solve the linear equation the ratio of the derivatives must be
            # a constant. Improve performance by filtering out failing pairs early.
            ratio = x.derivative() / y.derivative()
            if not ratio.degree() == 0:
                continue

            determinant = x * y.derivative() - y * x.derivative()
            if determinant == 0:
                continue

            c_1 = y.derivative() / determinant
            c_2 = -x.derivative() / determinant
            if not c_1.degree() == 0 or not c_2.degree() == 0:
                continue

            unit_1 = c_1 * x
            unit_2 = c_2 * y

            yield unit_1, unit_2
