import pickle
from function import *
from utils import max_snd_eig_with_vec
from datetime import datetime
from notify import notify
import frange
import matrix_gen
import matplotlib.pyplot as plt


results_path = 'results/'

ls_eval_counter = 0
MX_B = True


class Ising(object):
    def __init__(self, matrix_gen=None, Hs=frange.dr, HBs=frange.dl, Bs=frange.drb, BHs=frange.dl):
        if not matrix_gen:
            return

        start = datetime.now()
        self.matrix_gen = matrix_gen
        self.Hs = Hs
        self.HBs = HBs
        self.Bs = Bs
        self.BHs = BHs

        self.ls_of_H, self.vec_of_H, self.snd_ls_of_H, self.snd_vec_of_H = ls_of('H', matrix_gen, Hs, HBs)

        ls_of_H_m, _, _, _ = ls_of('H', matrix_gen, Hs, [b - dshift for b in HBs])
        ls_of_H_p, _, _, _ = ls_of('H', matrix_gen, Hs, [b + dshift for b in HBs])
        self.ls_of_H.set_label(xlabel='H', ylabel='l')
        self.snd_ls_of_H.set_label(xlabel='H', ylabel='l')

        self.fs_of_H = f_H(self.ls_of_H)
        self.us_of_H = u_H(self.ls_of_H, ls_of_H_m, ls_of_H_p)
        self.Cs_of_H = C_H(self.ls_of_H, ls_of_H_m, ls_of_H_p)
        self.Ms_of_H = M_H(self.ls_of_H)
        self.Xs_of_H = X_H(self.ls_of_H)

        self.ls_of_B, self.vec_of_B, self.snd_ls_of_B, self.snd_vec_of_B = ls_of('B', matrix_gen, Bs, BHs)
        if MX_B:
            ls_of_B_m, _, _, _ = ls_of('B', matrix_gen, Bs, [h - dshift for h in BHs])
        if MX_B:
            ls_of_B_p, _, _, _ = ls_of('B', matrix_gen, Bs, [h + dshift for h in BHs])
        self.ls_of_B.set_label(xlabel='B', ylabel='l')
        self.snd_ls_of_B.set_label(xlabel='B', ylabel='l')
        self.fs_of_B = f_B(self.ls_of_B)
        self.us_of_B = u_B(self.ls_of_B)
        self.Cs_of_B = C_B(self.ls_of_B)
        if MX_B:
            self.Ms_of_B = M_B(self.ls_of_B, ls_of_B_m, ls_of_B_p)
        if MX_B:
            # print 'jej'
            self.Xs_of_B = X_B(self.ls_of_B, ls_of_B_m, ls_of_B_p)
        end = datetime.now()
        diff = str(end - start)
        notify(matrix_gen.__name__, diff)
        global ls_eval_counter
        ls_eval_counter = 0

    def to_dict(self):
        d = {}
        d["Hs"]           = self.Hs
        d["HBs"]          = self.HBs
        d["Bs"]           = self.Bs
        d["BHs"]          = self.BHs
        d["ls_of_H"]      = self.ls_of_H.to_dict()
        d["snd_ls_of_H"]  = self.snd_ls_of_H.to_dict()
        d["vec_of_H"]     = self.vec_of_H.to_dict()
        d["snd_vec_of_H"] = self.snd_vec_of_H.to_dict()
        d["fs_of_H"]      = self.fs_of_H.to_dict()
        d["us_of_H"]      = self.us_of_H.to_dict()
        d["Cs_of_H"]      = self.Cs_of_H.to_dict()
        d["Ms_of_H"]      = self.Ms_of_H.to_dict()
        d["Xs_of_H"]      = self.Xs_of_H.to_dict()
        d["ls_of_B"]      = self.ls_of_B.to_dict()
        d["snd_ls_of_B"]  = self.snd_ls_of_B.to_dict()
        d["vec_of_B"]     = self.vec_of_B.to_dict()
        d["snd_vec_of_B"] = self.snd_vec_of_B.to_dict()
        d["fs_of_B"]      = self.fs_of_B.to_dict()
        d["us_of_B"]      = self.us_of_B.to_dict()
        d["Cs_of_B"]      = self.Cs_of_B.to_dict()
        if MX_B: d["Ms_of_B"]      = self.Ms_of_B.to_dict()
        if MX_B: d["Xs_of_B"]      = self.Xs_of_B.to_dict()
        return d

    @staticmethod
    def from_dict(d):
        ni = Ising()
        ni.Hs           = d["Hs"]
        ni.HBs          = d["HBs"]
        ni.Bs           = d["Bs"]
        ni.BHs          = d["BHs"]
        ni.ls_of_H      = FunctionBundle.from_dict(d["ls_of_H"])
        ni.snd_ls_of_H  = FunctionBundle.from_dict(d["snd_ls_of_H"])
        ni.vec_of_H     = FunctionBundle.from_dict(d["vec_of_H"])
        ni.snd_vec_of_H = FunctionBundle.from_dict(d["snd_vec_of_H"])
        ni.fs_of_H      = FunctionBundle.from_dict(d["fs_of_H"])
        ni.us_of_H      = FunctionBundle.from_dict(d["us_of_H"])
        ni.Cs_of_H      = FunctionBundle.from_dict(d["Cs_of_H"])
        ni.Ms_of_H      = FunctionBundle.from_dict(d["Ms_of_H"])
        ni.Xs_of_H      = FunctionBundle.from_dict(d["Xs_of_H"])
        ni.ls_of_B      = FunctionBundle.from_dict(d["ls_of_B"])
        ni.snd_ls_of_B  = FunctionBundle.from_dict(d["snd_ls_of_B"])
        ni.vec_of_B     = FunctionBundle.from_dict(d["vec_of_B"])
        ni.snd_vec_of_B = FunctionBundle.from_dict(d["snd_vec_of_B"])
        ni.fs_of_B      = FunctionBundle.from_dict(d["fs_of_B"])
        ni.us_of_B      = FunctionBundle.from_dict(d["us_of_B"])
        ni.Cs_of_B      = FunctionBundle.from_dict(d["Cs_of_B"])
        if MX_B: ni.Ms_of_B = FunctionBundle.from_dict(d["Ms_of_B"])
        if MX_B: ni.Xs_of_B = FunctionBundle.from_dict(d["Xs_of_B"])
        return ni

    @staticmethod
    def load(filepath):
        with open(filepath, 'r') as f:
            d = pickle.load(f)
            return Ising.from_dict(d)

    def save(self, filename_prefix=''):
        if filename_prefix:
            filename_prefix += '_'
        filepath = results_path + filename_prefix + self.filename_sufix()
        with open(filepath, 'w') as f:
            pickle.dump(self.to_dict(), f)
        return filepath

    def filename_sufix(self):
        return '{mgen}_{size}_h_{hmin}-{hmax}-{hb}_b_{bmin}-{bmax}-{bh}-{date}'.format(
            mgen=self.matrix_gen.__name__,
            size=self.matrix_gen.size,
            hmin=self.Hs[0],
            hmax=self.Hs[-1],
            bmin=self.Bs[0],
            bmax=self.Hs[-1],
            hb=','.join([str(x) for x in self.HBs]),
            bh=','.join([str(x) for x in self.BHs]),
            date=datetime.now().isoformat()
            )

    def plot(self, draw_legend=True, show=False):
        self.plot_of_B(draw_legend)
        self.plot_of_H(draw_legend)
        if show:
            plt.show()

    def plot_of_H(self, draw_legend=True):
        plt.figure()
        set_label_start(r'$\beta$ = ')
        self.ls_of_H.plot(draw_legend)
        plt.figure()
        self.fs_of_H.plot(draw_legend)
        plt.figure()
        self.us_of_H.plot(draw_legend)
        plt.figure()
        self.Cs_of_H.plot(draw_legend)
        plt.figure()
        self.Ms_of_H.plot(draw_legend)
        plt.figure()
        self.Xs_of_H.plot(draw_legend)

    def plot_of_B(self, draw_legend=True):
        plt.figure()
        set_label_start('H = ')
        self.ls_of_B.plot(draw_legend)
        plt.figure()
        self.fs_of_B.plot(draw_legend)
        plt.figure()
        self.us_of_B.plot(draw_legend)
        plt.figure()
        self.Cs_of_B.plot(draw_legend)
        if MX_B: plt.figure()
        if MX_B: self.Ms_of_B.plot(draw_legend)
        if MX_B: plt.figure()
        if MX_B: self.Xs_of_B.plot(draw_legend)

    def plot_vec(self, how_many=11, ks=None):
        self.plot_vec_of_B(how_many, ks)
        self.plot_vec_of_H(how_many, ks)

    def plot_vec_of_B(self, how_many=11, ks=None):
        if ks is None:
            ks = self.vec_of_H.functions.keys()
        if not (isinstance(ks, list) or isinstance(ks, tuple)):
            ks = [ks]
        for k in ks:
            # print k
            self.vec_of_B.functions[k].plot_vec(how_many=how_many)
            plt.figure()

    def plot_vec_of_H(self, how_many=11, ks=None):
        if ks is None:
            ks = self.vec_of_H.functions.keys()
        if not (isinstance(ks, list) or isinstance(ks, tuple)):
            ks = [ks]
        for k in ks:
            # print k
            self.vec_of_H.functions[k].plot_vec(how_many=how_many)
            plt.figure()


def ls_of(of_what, matrix_gen, xs, params):
    global ls_eval_counter
    max_eig_functions = []
    snd_eig_functions = []
    max_vec_functions = []
    snd_vec_functions = []
    for param in params:
        max_eigs, max_vecs, snd_eigs, snd_vecs = [], [], [], []
        for x in xs:
            ls_eval_counter += 1
            # print ls_eval_counter
            if of_what == 'H':
                eigs = matrix_gen.eig(matrix_gen(param, x))
            elif of_what == 'B':
                eigs = matrix_gen.eig(matrix_gen(x, param))
            max_eig, max_vec, snd_eig, snd_vec = max_snd_eig_with_vec(eigs)
            max_eigs.append(max_eig)
            max_vecs.append(max_vec)
            snd_eigs.append(snd_eig)
            snd_vecs.append(snd_vec)
        f_max_eig = Function(xs=xs, ys=max_eigs, label=str(param))
        f_max_eig.param = param
        max_eig_functions.append(f_max_eig)

        f_max_vec = Function(xs=xs, ys=max_vecs, label=str(param))
        f_max_vec.param = param
        max_vec_functions.append(f_max_vec)
        f_snd_eig = Function(xs=xs, ys=snd_eigs, label=str(param))
        f_snd_eig.param = param
        snd_eig_functions.append(f_snd_eig)

        f_snd_vec = Function(xs=xs, ys=snd_vecs, label=str(param))
        f_snd_vec.param = param
        snd_vec_functions.append(f_snd_vec)
    fb_max_eig = FunctionBundle(max_eig_functions, params)
    fb_max_vec = FunctionBundle(max_vec_functions, params)
    fb_snd_eig = FunctionBundle(snd_eig_functions, params)
    fb_snd_vec = FunctionBundle(snd_vec_functions, params)
    return fb_max_eig, fb_max_vec, fb_snd_eig, fb_snd_vec


class Isings(object):
    def __init__(self, sizes=[], matrix_gen=None, Hs=frange.dr, HBs=frange.dl, Bs=frange.drb, BHs=frange.dl):
        self.isings = {}
        for size in sizes:
            # print "SIZE ", size
            self.isings[size] = Ising(matrix_gen(size), Hs, HBs, Bs, BHs)

    def save(self, filename_prefix=''):
        names = {}
        for size, i in self.isings.items():
            name = i.save()
            names[size] = name
        filepath = results_path + filename_prefix + self.filename_sufix()
        with open(filepath, 'w') as f:
            pickle.dump(names, f)
        return filepath

    @staticmethod
    def load(filepath):
        with open(filepath, 'r') as f:
            i = Isings()
            d = pickle.load(f)
            for size, filename in d.items():
                i.isings[size] = Ising.load(filename)
            return i

    def filename_sufix(self):
        k = self.isings.keys()[0]
        return '{mgen}_{size}_{date}'.format(
            mgen=self.isings[k].matrix_gen.__name__,
            size=len(self.isings),
            date=datetime.now().isoformat()
        )

    def keys(self):
        k = self.isings.keys()[0]
        # print 'B:', self.isings[k].ls_of_H.functions.keys()
        k2 = self.isings[k].ls_of_H.functions.keys()[0]
        # print 'H:', self.isings[k].ls_of_B.functions.keys()
        k3 = self.isings[k].ls_of_B.functions.keys()[0]
        # print self.isings[k].ls_of_H.functions[k2].xs
        # print self.isings[k].ls_of_B.functions[k3].xs

    def plot(self, what='', H=None, B=None):
        """what = l, f, u, C, M, X, v"""
        set_label_start('N = ')
        if not what:
            # print 'what!?'
            return
        if H and B:
            # print 'H albo B1'
            return
        if H is None and B is None:
            # print 'H albo B2'
            return

        what += 's_of_'
        if H:
            what += 'H'
        else:
            what += 'B'

        param = H
        if H is None:
            param = B
        for size, i in self.isings.items():
            i.__getattribute__(what).functions[param].plot(label=str(size))
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    def plot_vec_of_H(self, H=None, B=None, how_many=11):
        if H is None:
            H = self.isings[2].vec_of_H.functions.keys()[0]

        for size, i in self.isings.items():
            i.vec_of_H.functions[B].plot_vec(label=str(size), param=H, how_many=how_many)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    def plot_vec_of_B(self, H=None, B=None, how_many=11):
        if B is None:
            B = self.isings[2].vec_of_B.functions.keys()[0]

        for size, i in self.isings.items():
            i.vec_of_B.functions[H].plot_vec(label=str(size), param=B, how_many=how_many)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


def plot_vec(mg, B, H):
        plt.xlabel('nr.')
        plt.ylabel('$V^2$')
        e = mg.eig(mg(B, H))
        # print e[0]
        # print e[1]
        max_eig, max_vec, snd_eig, snd_vec = max_snd_eig_with_vec(
            e
        )
        # print max_eig, max_vec
        ys = [y ** 2 for y in max_vec]
        plt.axis([-10, len(max_vec) + 9, 0, max(ys)])
        plt.vlines(range(0, len(ys)), 0, ys, color=color())
        # return ys

        # plt.plot(range(0, len(max_vec)), [y ** 2 for y in max_vec], label='N = ' + str(N), color=color())
