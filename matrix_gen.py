from numpy.linalg import eig, eigh
import numpy as np
import math
from utils import matrix_inverse, dot


def matrix_gen(eig, size):  # it's decorator
    def mg(f):
        f.eig = eig
        f.size = size
        return f
    return mg


@matrix_gen(eigh, 1)
def chain_matrix(B, H, J=1):
    K = J * B
    h = H * B
    m = np.zeros((2, 2))
    m[0, 0]           = np.exp(K + h)
    m[1, 0] = m[0, 1] = np.exp(-K)
    m[1, 1]           = np.exp(K - h)
    return m


@matrix_gen(eigh, 2)
def ladder_matrix(B, H, J=1):
    K = J * B
    h = H * B
    m = np.zeros((4, 4))
    m[0, 0]           = np.exp(3 * K + 2 * h)  # np.exp(3 * B + 2 * H * B)
    m[1, 0] = m[0, 1] = np.exp(h)  # np.exp(H * B)
    m[2, 0] = m[0, 2] = np.exp(h)  # np.exp(H * B)
    m[3, 0] = m[0, 3] = np.exp(K)  # np.exp(B)
    m[1, 1]           = np.exp(K)  # np.exp(B)
    m[1, 2] = m[2, 1] = np.exp(-3 * K)  # np.exp(-3 * B)
    m[1, 3] = m[3, 1] = np.exp(-h)  # np.exp(-H * B)
    m[2, 2]           = np.exp(-K)  # np.exp(B)
    m[2, 3] = m[3, 2] = np.exp(-h)  # np.exp(-H * B)
    m[3, 3]           = np.exp(3 * K - 2 * h)  # np.exp(3 * B - 2 * H * B)
    return m


# Generacja macierzy dla dowolnego rozmiaru
## Funkcje pomocnicze
def int_to_sigmas(sigmas_asint, size):
    """
    wstep: kazdy uklad spinow moze byc zapisany jako liczba binarna,
    ktora moze byc przekonwertowana do postaci dziesietnej

    ta funkcja robi odwrotnie - zamienia liczbe w zapisie dziesietnym
    na konfiguracje spinow w postaci listy 1 i -1
    """
    binstr = bin(sigmas_asint)[2:].zfill(size)
    binlist = [1 if x == '1' else -1 for x in binstr]
    return binlist


def sum_products_of_cyclic_pairs(sigmas):
    """oblicza
    \sum_i \sigma_i * \sigma_{i+1}
    """
    prev_sigmas = sigmas[:-1]
    next_sigmas = sigmas[1:]  # + [sigmas[0]]
    result = sum([ps * ns for ps, ns in zip(prev_sigmas, next_sigmas)])
    return result


def beta_coef(sigmas, sigmasprim, *args):
    """
    oblicza wspoczynnik przy beta
        a = \sum_i \sigma_i * \sigma_i'
        b = \sum_i \sigma_i * \sigma_{i+1}
        c = \sum_i \sigma_i' * \sigma_{i+1}'
        wsp = a + b + c
    """
    a = sum([s * spr for s, spr in zip(sigmas, sigmasprim)])
    b = sum_products_of_cyclic_pairs(sigmas) / 2.
    c = sum_products_of_cyclic_pairs(sigmasprim) / 2.
    return a + b + c


def sum_products_of_s1_s2p(sigmas, sigmasprim):
    """
        \sum_i \sigma_i * \sigma_{i+1}'
    """
    s1 = sigmas
    s2 = sigmasprim[1:] + [sigmasprim[0]]
    return sum([s * spr for s, spr in zip(s1, s2)])


def beta_coef_triangle(sigmas, sigmasprim, *args):
    """
    oblicza wspoczynnik przy beta
        a = \sum_i \sigma_i * \sigma_i'
        b = \sum_i \sigma_i * \sigma_{i+1}
        c = \sum_i \sigma_i' * \sigma_{i+1}'
        d = \sum_i \sigma_i * \sigma_{i+1}''
        wsp = a + b + c + d
    """
    abc = beta_coef(sigmas, sigmasprim)
    d = sum_products_of_s1_s2p(sigmas, sigmasprim)
    return abc + d


def beta_coef_cube(edges):
    """
    oblicza wspoczynnik przy beta
        a = \sum_i \sigma_i * \sigma_i'
        b = \sum_i \sigma_i * \sigma_{i+1}
        c = \sum_i \sigma_i' * \sigma_{i+1}'

        wsp = a + b + c + d

        edges to `dodatkowe` krawdzie pomiedzy spinami -
        to znaczy takie, ktorych nie ma w siatce 2d

        conns to suma iloczynow takich spinow, ktore sa polaczone
        `dodatkowa` krawedzia

        s2---|
        |    |
        s1   |
        |    |
        s0---|

        w reprezentacji siatki kubicznej 2x2 jako siatki 2d
        s0 i s1, s1 i s2 sa poloczone `standardowo`
        a s0 i s2 `dodatkowo`
    """
    def f(sigmas, sigmasprim):
        abc = beta_coef(sigmas, sigmasprim)
        conns = sum([sigmas[v1] * sigmas[v2] * p for v1, v2, p in edges])
        return abc + conns
    return f


def H_coef(sigmas, sigmasprim):
    """
    oblicza wspolczynnik przy H
        sum_i (\sigma_i + \sigma_i') / 2
    """
    return (sum(sigmas) + sum(sigmasprim)) / 2.


def sum_pm(sigmas, start=1):
    i = start
    sum_ = 0
    for sigma in sigmas:
        sum_ += i * sigma
        i *= -1
    return sum_


start = 1
def H_coef_field1(sigmas, sigmasprim):
    """
    oblicza wspolczynnik przy H
        ale tak, ze na siatce kwadratowej na zmiane sa plusy
        i minusy, z takim wspolczynnikiem jest brana sigma
        -+-+-+-
        +-+-+-+
        -+-+-+-
        +-+-+-+
    """
    start = 1
    r = sum_pm(sigmas, start=start) + sum_pm(sigmasprim, start=-start)
    # print r
    return r


def H_coef_field2(sigmas, sigmasprim):
    """
    oblicza wspolczynnik przy H
        ale tak, ze na siatce kwadratowej na zmiane sa plusy
        i minusy, z takim wspolczynnikiem jest brana sigma
        -+-+-+-
        +-+-+-+
        -+-+-+-
        +-+-+-+
    """
    start = -1
    r = sum_pm(sigmas, start=start) + sum_pm(sigmasprim, start=-start)
    # print r
    return r


def coefs(sigmas_asint, sigmasprim_asint, size, H_coef=H_coef, beta_coef=beta_coef):
    """
    liczy wartosc wspolczynnikow przy beta i H
    dla komorki danej jako para liczb w zapisie dziesietnym
    """
    sigmas = int_to_sigmas(sigmas_asint, size)
    sigmasprim = int_to_sigmas(sigmasprim_asint, size)
    bc = beta_coef(sigmas, sigmasprim)
    Hc = H_coef(sigmas, sigmasprim)
    return (bc, Hc)


def coefs_field1(sigmas_asint, sigmasprim_asint, size):
    return coefs(sigmas_asint, sigmasprim_asint, size, H_coef=H_coef_field1)


def coefs_field2(sigmas_asint, sigmasprim_asint, size):
    return coefs(sigmas_asint, sigmasprim_asint, size, H_coef=H_coef_field2)


def coefs_triangle(sigmas_asint, sigmasprim_asint, size):
    return coefs(sigmas_asint, sigmasprim_asint, size, beta_coef=beta_coef_triangle)


def coefs_cube(edges):
    def f(sigmas_asint, sigmasprim_asint, size):
        return coefs(sigmas_asint, sigmasprim_asint, size, beta_coef=beta_coef_cube(edges))
    return f


def gen_matrix_coefs(size, coefs=coefs):
    """
    tworzy macierz par,
    pierwsza wartosc to wspolczynnik przy beta w danej komorce,
    druga wartosc to wspolczynnik przy H w danej komorce.
    """
    matrix = []
    msize = (2 ** size)
    size_rev_range = range(msize - 1, -1, -1)
    for row_i in size_rev_range:
        row = []
        for col_i in size_rev_range:
            row.append((coefs(row_i, col_i, size)))
        matrix.append(row)
    return matrix


def gen_matrix_coefs_triangle(size):
    return gen_matrix_coefs(size, coefs=coefs_triangle)


def gen_matrix_coefs_cube(size, edges):
    return gen_matrix_coefs(size, coefs=coefs_cube(edges))


def gen_matrix_coefs_field1(size):
    return gen_matrix_coefs(size, coefs=coefs_field1)


def gen_matrix_coefs_field2(size):
    return gen_matrix_coefs(size, coefs=coefs_field2)


def gen_matrix_coefs_field(size):
    m1 = gen_matrix_coefs(size, coefs=coefs_field1)
    m2 = gen_matrix_coefs(size, coefs=coefs_field2)
    m2i = matrix_inverse(m2)
    m = []
    for row in m1:
        m.append([])
        for col in m2i:
            r = []
            for x1, x2 in zip(row, col):
                r.append((x1[0] + x2[0], x1[1] + x2[1]))
            m[-1].append(r)
    return m


def gen_rectangle_eval(size):
    """
    tworzy funkcje zwracajaca macierz przejscia o podanym
    rozmiarze, ktorej argumentami sa beta i H
    """
    coefs_matrix = gen_matrix_coefs(size)

    @matrix_gen(eigh, size)
    def rectangle_grid_matrix(bet, H):
        m = []
        for row in coefs_matrix:
            m.append([])
            for bi, hi in row:
                m[-1].append(math.exp(bi * bet + hi * bet * H))
        return m
    return rectangle_grid_matrix

gr = gen_rectangle_eval


def gen_rectangle_eval_af(size):
    """
    tworzy funkcje zwracajaca macierz przejscia o podanym
    rozmiarze, ktorej argumentami sa beta i H
    """
    coefs_matrix = gen_matrix_coefs(size)

    @matrix_gen(eigh, size)
    def rectangle_grid_matrix_af(bet, H):
        m = []
        for row in coefs_matrix:
            m.append([])
            for bi, hi in row:
                m[-1].append(math.exp(-bi * bet + hi * bet * H))
        return m
    return rectangle_grid_matrix_af

ga = gen_rectangle_eval_af


def gen_triangle_eval(size):
    """
    tworzy funkcje zwracajaca macierz przejscia o podanym
    rozmiarze, ktorej argumentami sa beta i H
    """
    coefs_matrix = gen_matrix_coefs_triangle(size)

    @matrix_gen(eig, size)
    def triangle_grid_matrix(bet, H):
        m = []
        for row in coefs_matrix:
            m.append([])
            for bi, hi in row:
                m[-1].append(math.exp(bi * bet + hi * bet * H))
        return m
    return triangle_grid_matrix

gt = gen_triangle_eval


def gen_triangle_eval_af(size):
    """
    tworzy funkcje zwracajaca macierz przejscia o podanym
    rozmiarze, ktorej argumentami sa beta i H
    """
    coefs_matrix = gen_matrix_coefs_triangle(size)

    @matrix_gen(eig, size)
    def triangle_grid_matrix(bet, H):
        m = []
        for row in coefs_matrix:
            m.append([])
            for bi, hi in row:
                m[-1].append(math.exp(-bi * bet + hi * bet * H))
        return m
    return triangle_grid_matrix

gtaf = gen_triangle_eval_af


def gen_field_eval(size):
    """
    tworzy funkcje zwracajaca macierz przejscia o podanym
    rozmiarze, ktorej argumentami sa beta i H
    """
    coefs_matrix = gen_matrix_coefs_field(size)

    @matrix_gen(eig, size)
    def field_grid_matrix(bet, H):
        m = []
        for row in coefs_matrix:
            m.append([])
            for cs in row:
                r = 0
                for bi, hi in cs:
                    r += math.exp(bi * bet + hi * bet * H)
                m[-1].append(r)
        return m
    return field_grid_matrix

gf = gen_field_eval


def gen_field_eval_af(size):
    """
    tworzy funkcje zwracajaca macierz przejscia o podanym
    rozmiarze, ktorej argumentami sa beta i H
    """
    coefs_matrix = gen_matrix_coefs_field(size)

    @matrix_gen(eig, size)
    def field_grid_matrix(bet, H):
        m = []
        for row in coefs_matrix:
            m.append([])
            for cs in row:
                r = 0
                for bi, hi in cs:
                    r += math.exp(-bi * bet + hi * bet * H)
                m[-1].append(r)
        return m
    return field_grid_matrix

gfa = gen_field_eval_af


print '2x2'
_gc2x2 = gen_matrix_coefs_cube(4, [(0, 3, 1)])


@matrix_gen(eigh, 4)
def cube2x2_grid_matrix(bet, H):
    m = []
    for row in _gc2x2:
        m.append([])
        for bi, hi in row:
            m[-1].append(math.exp(bi * bet + hi * bet * H))
    return m

c22 = cube2x2_grid_matrix

print '2x3'
_gc2x3 = gen_matrix_coefs_cube(6, [(0, 5, 1), (1, 4, 1)])


@matrix_gen(eigh, 6)
def cube2x3_grid_matrix(bet, H):
    m = []
    for row in _gc2x3:
        m.append([])
        for bi, hi in row:
            m[-1].append(math.exp(bi * bet + hi * bet * H))
    return m

c23 = cube2x3_grid_matrix

# print '3x3'
# _gc3x3 = gen_matrix_coefs_cube(9, [(0, 5, 1), (1, 4, 1), (3, 8, 1), (4, 7, 1)])


# @matrix_gen(eigh, 9)
# def cube3x3_grid_matrix(bet, H):
#     m = []
#     for row in _gc3x3:
#         m.append([])
#         for bi, hi in row:
#             m[-1].append(math.exp(bi * bet + hi * bet * H))
#     return m

# c33 = cube3x3_grid_matrix

# print '3x4'
# _gc3x4 = gen_matrix_coefs_cube(12, [(0, 7, 1), (1, 6, 1), (2, 5, 1), (4, 11, 1), (5, 10, 1), (6, 9, 1)])
# print 'generated 3x4'


# @matrix_gen(eigh, 12)
# def cube3x4_grid_matrix(bet, H):
#     m = []
#     for row in _gc3x4:
#         m.append([])
#         for bi, hi in row:
#             m[-1].append(math.exp(bi * bet + hi * bet * H))
#     return m
