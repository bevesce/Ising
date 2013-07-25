def frange(start, stop, step=None, howmany=None):
    if not step and not howmany:
        howmany = 100

    range_list = []
    if not step and howmany:
        step = (stop - start) / float(howmany)
    av = start

    while av < stop:
        range_list.append(av)
        av += step
    return range_list

def_Bs = frange(0., 2, howmany=50)
def_Hs_for_B = [0., 0.5, 1.]

def_Bs_for_H = [0.1, 0.7, 1.5]
def_Hs = frange(-4, 4, howmany=50)

Bs_af = [1., 1.5, 2., 3., 4.]
Hs_af = frange(-6, 6, step=0.1)

dr = frange(-3.0, 3.0, howmany=500)
dl = [0.3, 0.44068, 3.0]
drb = frange(0.00001, 1., howmany=100)

# beta = 5, 3, 2, 1, 0.5, 0.2
