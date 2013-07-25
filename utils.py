def max_snd_eig_with_vec(eigs):
    """Zwraca najwieksza i druga z koleji wartosc wlasna i ich wektory"""
    max_eig = 0
    snd_eig = 0
    max_vec_i = None
    snd_vec_i = None
    for i, eigv in enumerate(eigs[0]):
        if eigv > max_eig:
            snd_eig = max_eig
            snd_vec_i = max_vec_i
            max_eig = eigv
            max_vec_i = i
        elif eigv > snd_eig:
            snd_eig = eigv
            snd_vec_i = i
    max_vec = eigs[1][:,max_vec_i]
    snd_vec = eigs[1][:,max_vec_i]
    return max_eig, max_vec, snd_eig, snd_vec


def matrix_inverse(m):
    r = []
    n = len(m)
    for i in range(0, n):
        r.append([0] * n)
    for i in range(0, n):
        for j in range(0, n):
            r[i][j] = m[j][i]
    return r


def dot(v1, v2):
    bi = 0
    hi = 0
    for x, y in zip(v1, v2):
        bi1, hi1 = x
        bi2, hi2 = y
        bi += bi1 + bi2
        hi += hi1 + hi2
    return bi, hi
