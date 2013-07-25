import matplotlib.pyplot as plt
import numpy as np

dshift = 0.01

colors = [
    'green', 'blue', 'red', 'magenta', 'cyan'
]

color_idx = 0


def color():
    global color_idx
    color_idx = (color_idx + 1) % len(colors)
    return colors[color_idx]


def restart_color():
    global color_idx
    color_idx = 0

label_start = r'$\beta$ = '


def set_label_start(s):
    global label_start
    label_start = s


class Function(object):
    def __init__(self, xs=[], ys=[], xlabel="x", ylabel="y", label=""):
        self.xs = xs
        self.ys = ys
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.label = label
        self.param = None

    def to_dict(self):
        return {
            "xs": self.xs,
            "ys": self.ys,
            "xlabel": self.xlabel,
            "ylabel": self.ylabel,
            "label": self.label,
            "param": self.param,
        }

    @staticmethod
    def from_dict(d):
        new        = Function()
        new.xs     = d["xs"]
        new.ys     = d["ys"]
        new.xlabel = d["xlabel"]
        new.ylabel = d["ylabel"]
        new.label  = d["label"]
        new.param  = d["param"]
        return new

    def set_label(self, xlabel=None, ylabel=None):
        if not xlabel is None:
            self.xlabel = xlabel
        if not ylabel is None:
            self.ylabel = ylabel

    def plot(self, label=None):
        if not label:
            label = self.label

        xlabel, ylabel = get_true_label(self.xlabel, self.ylabel)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.plot(self.xs, self.ys, label=label_start + str(label), color=color())

    def plot_vec(self, draw_legend=True, how_many=None, label=None, param=None):
        if how_many is None:
            how_many = len(self.xs)
        plt.xlabel(' ')
        plt.ylabel(' ')

        if not param:
            if how_many:
                fxs = self.xs[::len(self.xs)/how_many]
                fys = self.ys[::len(self.ys)/how_many]
            else:
                fxs = self.xs
                fys = self.ys
        else:
            fxs = [param]
            fys = [self.ys[self.xs.index(param)]]
        for x, ys in zip(fxs, fys):
            plt.plot(range(0, len(ys)), [y ** 2 for y in ys], label=(str(x) if (label is None) else label), color=color())
        if draw_legend:
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    def diff(self, ylabel=None, label=None):
        res_xs = []
        res_ys = []
        prev_x, prev_y = self.xs[0], self.ys[0]
        for x, y in zip(self.xs[1:], self.ys[1:]):
            try:
                dx = (x - prev_x) / 2 + prev_x
                dy = (y - prev_y) / (x - prev_x)
                res_xs.append(dx)
                res_ys.append(dy)
            except Exception as e:
                print e
                pass
            prev_x, prev_y = x, y

        nf = Function(
            xs=res_xs,
            ys=res_ys,
            xlabel=self.xlabel,
            ylabel=ylabel or (self.ylabel + "'"),
            label=label or self.label
        )
        nf.param = self.param
        return nf

    def snd_diff(self, ylabel=None, label=None):
        res_xs = []
        res_ys = []
        prevprev_x, prevprev_y = self.xs[0], self.ys[0]
        prev_x, prev_y = self.xs[1], self.ys[1]
        for x, y in zip(self.xs[2:], self.ys[2:]):
            try:
                dx = prev_x
                h = (x - prevprev_x) / 2
                dy = (y - 2 * prev_y + prevprev_y) / (h ** 2)
                res_xs.append(dx)
                res_ys.append(dy)
            except Exception as e:
                print e
            prevprev_x, prevprev_y = prev_x, prev_y
            prev_x, prev_y = x, y

        nf = Function(
            xs=res_xs,
            ys=res_ys,
            xlabel=self.xlabel,
            ylabel=ylabel or (self.ylabel + "''"),
            label=label or self.label
        )
        nf.param = self.param
        return nf

    def param_apply(self, f, ylabel=None, label=None):
        if not self.param is None:
            nf = Function(
                xs=self.xs,
                ys=[f(y, self.param) for y in self.ys],
                xlabel=self.xlabel,
                ylabel=ylabel or ("{1}({0})".format(self.xlabel, f.__name__)),
                label=label or self.label
            )
            nf.param = self.param
            return nf
        else:
            raise Exception('no param!')

    def apply(self, f, ylabel=None, label=None):
        nf = Function(
            xs=self.xs,
            ys=[f(y) for y in self.ys],
            xlabel=self.xlabel,
            ylabel=ylabel or ("{1}({0})".format(self.xlabel, f.__name__)),
            label=label or self.label
        )
        nf.param = self.param
        return nf

    def apply2(self, f, ylabel=None, label=None):
        nf = Function(
            xs=self.xs,
            ys=[f(x, y) for x, y in zip(self.xs, self.ys)],
            xlabel=self.xlabel,
            ylabel=ylabel or ("{1}({0})".format(self.xlabel, f.__name__)),
            label=label or self.label
            )
        nf.param = self.param
        return nf

    def __add__(self, other):
        f1 = [(x, y) for x, y in zip(self.xs, self.ys)]
        f2 = {x:y for x, y in zip(other.xs, other.ys)}
        ys = []
        for x, y in f1:
            ys.append(y + f2[x])
        nf = Function(
            xs=self.xs,
            ys=ys,
            xlabel=self.xlabel,
            ylabel=self.ylabel + ' + ' + other.ylabel,
            label=self.label
            )
        nf.param = self.param
        return nf

    def __sub__(self, other):
        f1 = [(x, y) for x, y in zip(self.xs, self.ys)]
        f2 = {x:y for x, y in zip(other.xs, other.ys)}
        ys = []
        for x, y in f1:
            ys.append(y - f2[x])
        nf = Function(
            xs=self.xs,
            ys=ys,
            xlabel=self.xlabel,
            ylabel=self.ylabel + ' - ' + other.ylabel,
            label=self.label
            )
        nf.param = self.param
        return nf


from collections import OrderedDict


class FunctionBundle(object):
    def __init__(self, functions_list=[], labels_list=None):
        if not labels_list:
            labels_list = [str(x) for x in range(0, len(functions_list))]

        self.functions = OrderedDict()
        for labe, fun in  zip(labels_list, functions_list):
            self.functions[labe] = fun

    def to_dict(self):
        return {
            label: self.functions[label].to_dict() for label in self.functions
        }

    @staticmethod
    def from_dict(d):
        nfb = FunctionBundle()
        nfb.functions = {
            k: Function.from_dict(d[k]) for k in d
        }
        return nfb

    def plot(self, draw_legend=True):
        for label in self.functions:
            self.functions[label].label = str(label)
            self.functions[label].plot()
        if draw_legend:
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    def diff(self):
        new_functions = [
            self.functions[label].diff() for label in self.functions
        ]
        return FunctionBundle(new_functions, self.labels())

    def snd_diff(self):
        new_functions = [
            self.functions[label].snd_diff() for label in self.functions
        ]
        return FunctionBundle(new_functions, self.labels())

    def param_apply(self, f):
        new_functions = [
            self.functions[label].param_apply(f) for label in self.functions
        ]
        return FunctionBundle(new_functions, self.labels())

    def apply(self, f):
        new_functions = [
            self.functions[label].apply(f) for label in self.functions
        ]
        return FunctionBundle(new_functions, self.labels())

    def apply2(self, f):
        new_functions = [
            self.functions[label].apply2(f) for label in self.functions
        ]
        return FunctionBundle(new_functions, self.labels())

    def set_label(self, xlabel=None, ylabel=None):
        for function in self.functions.values():
            function.set_label(xlabel, ylabel)

    def __getitem__(self, label):
        return self.functions[label]

    def labels(self):
        return self.functions.keys()

    def __add__(self, other):
        nfs = []
        for label in self.functions:
            nfs.append(self.functions[label] + other.functions[label])
        return FunctionBundle(nfs, self.labels())

    def __sub__(self, other):
        nfs = []
        for label in self.functions:
            nfs.append(self.functions[label] - other.functions[label])
        return FunctionBundle(nfs, self.labels())


def f_H(ls):
    result = ls.apply(lambda x: np.log(x)).param_apply(lambda x, b: x / b)
    result.set_label(ylabel='f')
    return result


delta = 0.05


def u_H(ls, ls_m, ls_p):
    ls_m = ls.apply(lambda y: np.log(y))
    ls = ls.apply(lambda y: np.log(y))
    ls_p = ls.apply(lambda y: np.log(y))
    result = (ls_p - ls_m).apply(lambda x: x / (2 * delta)).param_apply(lambda x, p: x / p)
    # czy na pewno tu jest dzielenie przez B?
    result.set_label(ylabel='u')
    return result


def C_H(ls, ls_m, ls_p):
    ls_m = ls.apply(lambda y: np.log(y))
    ls = ls.apply(lambda y: np.log(y))
    ls_p = ls.apply(lambda y: np.log(y))
    result = (ls_m + ls_p - ls - ls).param_apply(lambda x, p: x / delta * p ** 2)
    result.set_label(ylabel='C')
    return result


def M_H(ls):
    #           to dzielenie z zamiany h na H <- h = HB
    result = ls.apply(lambda y: (np.log(y))).diff().param_apply(lambda x, p: x / p)
    result.set_label(ylabel='M')
    return result


def X_H(ls):
    result = ls.param_apply(lambda x, p: x / p).apply(lambda y: (np.log(y))).snd_diff()
    result.set_label(ylabel='X')
    return result


def f_B(ls):
    result = ls.apply2(lambda x, y: np.log(y) / (x + 0.))
    result.set_label(ylabel='f')
    return result


def u_B(ls):
    result = ls.apply(lambda y: np.log(y)).diff()
    result.set_label(ylabel='u')
    return result


def C_B(ls):
    result = ls.apply(
        lambda y: (np.log(y))
    ).snd_diff().apply2(
        lambda x, y: x * x * y
    )
    result.set_label(ylabel='C')
    return result


def M_B(ls, ls_m, ls_p):
    ls_m = ls.apply(lambda y: np.log(y)).param_apply(lambda x, p: x / p)
    ls = ls.apply(lambda y: np.log(y)).param_apply(lambda x, p: x / p)
    ls_p = ls.apply(lambda y: np.log(y)).param_apply(lambda x, p: x / p)
    result = (ls_p - ls_m).apply(lambda x: x / (2 * delta)).param_apply(lambda x, p: x / p)
    result.set_label(ylabel='M')
    return result


def X_B(ls, ls_m, ls_p):
    ls_m = ls_m.apply(lambda y: np.log(y)).apply2(lambda x, y: y/x)
    ls = ls.apply(lambda y: np.log(y)).apply2(lambda x, y: y/x)
    ls_p = ls_p.apply(lambda y: np.log(y)).apply2(lambda x, y: y/x)

    old_k = ls_m.functions.keys()
    for k in old_k:
        ls_m.functions[k + dshift] = ls_m.functions[k]
        del ls_m.functions[k]

    old_k = ls_p.functions.keys()
    for k in old_k:
        ls_p.functions[k - dshift] = ls_p.functions[k]
        del ls_p.functions[k]

    result = (ls_p + ls_m - ls - ls).apply(lambda x: x / (delta**2))
    result = result.apply2(lambda x, y: y/x)
    # result = result.diff()
    result.set_label(ylabel='$\chi$')
    result.set_label(xlabel=r'$\beta$')

    return result


_xs = np.linspace(-3, 3)
sin = Function(xs=_xs, ys=[np.sin(x) for x in _xs], xlabel='x', ylabel='sin(x)')


def get_true_label(xlabel, ylabel):
    ys = {
        'x': r'$\chi$',
        'X': r'$\chi$',
        'l': r'$\lambda$'
    }
    xs = {
        'B': r'$\beta$'
    }
    return xs.get(xlabel, xlabel), ys.get(ylabel, ylabel)
