from __future__ import division

import numpy as np
import scipy
from scipy.misc import (factorial,
                        comb as nchoosek,
                        )

from .helpers import nest, autocat

import logging
LOG = logging.getLogger(__name__)

#import decorator

#@decorator.decorator
#def logcall(fn, *args, **dargs):
#    LOG.debug(fn.__name__)
#    return fn(*args, **dargs)

IMAG_CONST = scipy.sqrt(-1)
PI_CONST = np.pi
NAN_CONST = np.NaN


class Pipeline(object):

    def geometric_moments_approx(self, points_array, faces_array, N):
        raise NotImplementedError()

    def geometric_moments_exact(self, points_array, faces_array, N):
        raise NotImplementedError()


class SerialPipeline(Pipeline):

    def geometric_moments_exact(self, points_array, faces_array, N):
        n_facets, n_vertices = faces_array.shape[:2]
        assert n_vertices == 3
        moments_array = np.zeros([N + 1, N + 1, N + 1])
        monomial_array = self.monomial_precalc(points_array, N)
        for face in faces_array:
            vertex_list = [points_array[_i, ...] for _i in face]
            Cf_list = [monomial_array[_i, ...] for _i in face]
            Vf = self.facet_volume(vertex_list)
            moments_array += Vf * self.term_Sijk(Cf_list, N)
        return self.factorial_scalar(N) * moments_array

    def factorial_scalar(self, N):
        i, j, k = np.mgrid[0:N + 1, 0:N + 1, 0:N + 1]
        return factorial(i) * factorial(j) * factorial(k) / (factorial(i + j + k + 2) * (i + j + k + 3))

    def monomial_precalc(self, points_array, N):
        n_points = points_array.shape[0]
        monomial_array = np.zeros([n_points, N + 1, N + 1, N + 1])
        tri_array = self.trinomial_precalc(N)
        for point_indx, point in enumerate(points_array):
            monomial_array[point_indx, ...] = self.mon_comb(
                point, tri_array, N)
        return monomial_array

    def mon_comb(self, vertex, tri_array, N, out=None):
        x, y, z = vertex
        c = np.zeros([N + 1, N + 1, N + 1])
        for i, j, k in nest(lambda: range(N + 1),
                            lambda _i: range(N - _i + 1),
                            lambda _i, _j: range(N - _i - _j + 1),
                            ):
            c[i, j, k] = tri_array[i, j, k] * \
                np.power(x, i) * np.power(y, j) * np.power(z, k)
        return c

    def term_Sijk(self, Cf_list, N):
        S = np.zeros([N + 1, N + 1, N + 1])
        C0, C1, C2 = Cf_list
        Dabc = self.term_Dabc(C1, C2, N)
        for i, j, k, ii, jj, kk in nest(lambda: range(N + 1),
                                        lambda _i: range(N - _i + 1),
                                        lambda _i, _j: range(N - _i - _j + 1),
                                        lambda _i, _j, _k: range(_i + 1),
                                        lambda _i, _j, _k, _ii: range(_j + 1),
                                        lambda _i, _j, _k, _ii, _jj: range(
                                            _k + 1),
                                        ):
            S[i, j, k] += C0[ii, jj, kk] * Dabc[i - ii, j - jj, k - kk]
        return S

    def trinomial_precalc(self, N):
        tri_array = np.zeros([N + 1, N + 1, N + 1])
        for i, j, k in nest(lambda: range(N + 1),
                            lambda _i: range(N - _i + 1),
                            lambda _i, _j: range(N - _i - _j + 1)
                            ):
            tri_array[i, j, k] = self.trinomial(i, j, k)
        return tri_array

    def trinomial(self, i, j, k):
        return factorial(i + j + k) / (factorial(i) * factorial(j) * factorial(k))

    def facet_volume(self, vertex_list):
        return np.linalg.det(autocat(vertex_list, axis=1))

    def term_Dabc(self, C1, C2, N):
        D = np.zeros([N + 1, N + 1, N + 1])
        for i, j, k, ii, jj, kk in nest(lambda: range(N + 1),
                                        lambda _i: range(N + 1),
                                        lambda _i, _j: range(N + 1),
                                        lambda _i, _j, _k: range(_i + 1),
                                        lambda _i, _j, _k, _ii: range(_j + 1),
                                        lambda _i, _j, _k, _ii, _jj: range(
                                            _k + 1)
                                        ):
            D[i, j, k] += C1[ii, jj, kk] * C2[i - ii, j - jj, k - kk]
        return D

    def zernike(self, G, N):
        V = np.zeros([N + 1, N + 1, N + 1], dtype=complex)
        for a, b, c, alpha in nest(lambda: range(int(N / 2) + 1),
                                   lambda _a: range(N - 2 * _a + 1),
                                   lambda _a, _b: range(N - 2 * _a - _b + 1),
                                   lambda _a, _b, _c: range(_a + _c + 1),
                                   ):
            V[a, b, c] += np.power(IMAG_CONST, alpha) * \
                nchoosek(a + c, alpha) * G[2 * a + c - alpha, alpha, b]

        W = np.zeros([N + 1, N + 1, N + 1], dtype=complex)
        for a, b, c, alpha in nest(lambda: range(int(N / 2) + 1),
                                   lambda _a: range(N - 2 * _a + 1),
                                   lambda _a, _b: range(N - 2 * _a - _b + 1),
                                   lambda _a, _b, _c: range(_a + 1),
                                   ):
            W[a, b, c] += np.power(-1, alpha) * np.power(2, a - alpha) * \
                nchoosek(a, alpha) * V[a - alpha, b, c + 2 * alpha]

        X = np.zeros([N + 1, N + 1, N + 1], dtype=complex)
        for a, b, c, alpha in nest(lambda: range(int(N / 2) + 1),
                                   lambda _a: range(N - 2 * _a + 1),
                                   lambda _a, _b: range(N - 2 * _a - _b + 1),
                                   lambda _a, _b, _c: range(_a + 1),
                                   ):
            X[a, b, c] += nchoosek(a, alpha) * W[a - alpha, b + 2 * alpha, c]

        Y = np.zeros([N + 1, N + 1, N + 1], dtype=complex)
        for l, nu, m, j in nest(lambda: range(N + 1),
                                lambda _l: range(int((N - _l) / 2) + 1),
                                lambda _l, _nu: range(_l + 1),
                                lambda _l, _nu, _m: range(int((_l - _m) / 2) + 1),
                                ):
            Y[l, nu, m] += self.Yljm(l, j, m) * X[nu + j, l - m - 2 * j, m]

        Z = np.zeros([N + 1, N + 1, N + 1], dtype=complex)
        for n, l, m, nu, in nest(lambda: range(N + 1),
                                 lambda _n: range(_n + 1),
                                 # there's an if...mod missing in this but it
                                 # still works?
                                 lambda _n, _l: range(_l + 1),
                                 lambda _n, _l, _m: range(int((_n - _l) / 2) + 1),
                                 ):
            k = (n - l) / 2
            Z[n, l, m] += (3 / (4 * PI_CONST)) * \
                self.Qklnu(k, l, nu) * np.conj(Y[l, nu, m])

        for n, l, m in nest(lambda: range(N + 1),
                            lambda _n: range(n + 1),
                            lambda _n, _l: range(l + 1),
                            ):
            if np.mod(np.sum([n, l, m]), 2) == 0:
                Z[n, l, m] = np.real(
                    Z[n, l, m]) - np.imag(Z[n, l, m]) * IMAG_CONST
            else:
                Z[n, l, m] = -np.real(Z[n, l, m]) + \
                    np.imag(Z[n, l, m]) * IMAG_CONST

        return Z

    def Yljm(self, l, j, m):
        aux_1 = np.power(-1, j) * (np.sqrt(2 * l + 1) / np.power(2, l))
        aux_2 = self.trinomial(
            m, j, l - m - 2 * j) * nchoosek(2 * (l - j), l - j)
        aux_3 = np.sqrt(self.trinomial(m, m, l - m))
        y = (aux_1 * aux_2) / aux_3
        return y

    def Qklnu(self, k, l, nu):
        aux_1 = np.power(-1, k + nu) / np.power(4.0, k)
        aux_2 = np.sqrt((2 * l + 4 * k + 3) / 3.0)
        aux_3 = self.trinomial(
            nu, k - nu, l + nu + 1) * nchoosek(2 * (l + nu + 1 + k), l + nu + 1 + k)
        aux_4 = nchoosek(2.0 * (l + nu + 1), l + nu + 1)
        return (aux_1 * aux_2 * aux_3) / aux_4

    def feature_extraction(self, Z, N):
        F = np.zeros([N + 1, N + 1]) - 1  # +NAN_CONST
        for n in range(N + 1):
            for l in range(n + 1):
                if np.mod(n - l, 2) != 0:
                    continue
                aux_1 = Z[n, l, 0:(l + 1)]
                if l > 0:
                    aux_2 = np.conj(aux_1[1:(l + 1)])
                    for m in range(0, l):
                        aux_2[m] = aux_2[m] * np.power(-1, m + 1)
                    aux_2 = np.flipud(aux_2)
                    aux_1 = np.concatenate([aux_2, aux_1])
                F[n, l] = np.linalg.norm(aux_1, ord=2)
        F = F.transpose()
        return F[F >= 0]

import multiprocessing as mp


def _mp_geometric_moments_exact_worker(pipeline, vertex_list, Cf_list, N):
    Vf = pipeline.facet_volume(vertex_list)  # volume of the whole face
    return Vf * pipeline.term_Sijk(Cf_list, N)


def _mp_mon_comb_worker(pipeline, *args, **dargs):
    return pipeline.mon_comb(*args, **dargs)


class MultiprocPipeline(SerialPipeline):

    def geometric_moments_exact(self, points_array, faces_array, N):
        n_facets, n_vertices = faces_array.shape[:2]
        assert n_vertices == 3
        moments_array = np.zeros([N + 1, N + 1, N + 1])
        monomial_array = self.monomial_precalc(points_array, N)
        process_pool = mp.Pool()
        for face in faces_array:
            vertex_list = [points_array[_i, ...] for _i in face]
            monomial_list = [monomial_array[_i, ...] for _i in face]
            process_pool.apply_async(_mp_geometric_moments_exact_worker,
                                     args=(self, vertex_list, monomial_list, N),
                                     callback=moments_array.__iadd__,
                                     )
        process_pool.close()
        process_pool.join()
        return self.factorial_scalar(N) * moments_array

    def monomial_precalc(self, points_array, N):
        n_points = points_array.shape[0]
        monomial_array = np.zeros([n_points, N + 1, N + 1, N + 1])
        tri_array = self.trinomial_precalc(N)
        process_pool = mp.Pool()
        for point_indx, point in enumerate(points_array):
            def get_callback(_i):
                def __callback(result):
                    monomial_array[_i, ...] = result
                return __callback
            process_pool.apply_async(_mp_mon_comb_worker,
                                     args=(self, point, tri_array, N),
                                     callback=get_callback(point_indx),
                                     )
        process_pool.close()
        process_pool.join()
        return monomial_array

import itertools as it


def threeD_reversed(C):
    return C[::-1, ::-1, ::-1]


class NumpyOptimizations(Pipeline):

    def term_Dabc(self, C1, C2, N):
        D = np.zeros_like(C1)
        for a, b, c in it.product(range(N + 1), repeat=3):
            c1 = C1[:a + 1, :b + 1, :c + 1]
            c2 = threeD_reversed(C2[:a + 1, :b + 1, :c + 1])
            D[a, b, c] = np.sum(c1 * c2)
        return D

    def term_Sijk(self, Cf_list, N):
        S = np.zeros([N + 1, N + 1, N + 1])
        C0, C1, C2 = Cf_list
        Dabc = self.term_Dabc(C1, C2, N)
        for i, j, k in nest(lambda: range(N + 1),
                            lambda _i: range(N - _i + 1),
                            lambda _i, _j: range(N - _i - _j + 1),
                            ):
            C_ijk = C0[:i + 1, :j + 1, :k + 1]
            D_ijk = threeD_reversed(Dabc[:i + 1, :j + 1, :k + 1])
            S[i, j, k] += np.sum(C_ijk * D_ijk)
        return S

    def trinomial_precalc(self, N):
        i, k, j = np.mgrid[0:N + 1, 0:N + 1, 0:N + 1]
        return factorial(i + j + k) / (factorial(i) * factorial(j) * factorial(k))

    def mon_comb(self, vertex, tri_array, N):
        i, j, k = np.mgrid[0:N + 1, 0:N + 1, 0:N + 1]
        x, y, z = vertex
        return tri_array * (x ** i) * (y ** j) * (z ** k)

class KoehlOptimizations(Pipeline):

    def geometric_moments_exact(self, points_array, faces_array, N):
        n_facets, n_vertices = faces_array.shape[:2]
        assert n_vertices == 3
        moments_array = np.zeros([N+1, N+1, N+1])
        for face in faces_array:
            vertex_list = [points_array[_i, ...] for _i in face]
            moments_array += self.facet_contribution(vertex_list, N)
        return self.factorial_scalar(N)*moments_array

    def facet_contribution(self, vertex_list, N):
        Vf = self.facet_volume(vertex_list)
        Cf = self.term_Cijk(vertex_list[2], N)
        Df = self.term_Dijk(vertex_list[1], N, Cf)
        return Vf*self.term_Sijk(vertex_list[0], N, Df)

    def term_Cijk(self, vertex, N):
        return self.work_loop(vertex, N)

    def term_Dijk(self, vertex, N, Cijk):
        return self.work_loop(vertex, N, Cijk)

    def term_Sijk(self, vertex, N, Dijk):
        return self.work_loop(vertex, N, Dijk)

    def work_loop(self, vertex, N, prev=None):
        R = prev
        if R is None:
            R = np.zeros([N+1, N+1, N+1])
        Q = np.zeros([N+1, N+1, N+1])
        Q[0, 0, 0] = 1.0

        #recursion_term = lambda _X, (x, y, z), mask: \
        recursion_term = lambda _X, x, y, z, mask: \
            np.roll(_X, 1, axis=0)[mask]*x + \
            np.roll(_X, 1, axis=1)[mask]*y + \
            np.roll(_X, 1, axis=2)[mask]*z
        i, j, k = np.mgrid[:N+1, :N+1, :N+1]
        order = (i+j+k)
        for n in range(N):
            mask = (order==n+1)
            _Q = recursion_term(Q, vertex, mask)
            Q[mask] = _Q + R[mask]
        return Q

def _kmp_geometric_moments_exact_worker(self, vertex_list, N):
    return self.facet_contribution(vertex_list, N)

class KoehlMultiproc(KoehlOptimizations):
    def geometric_moments_exact(self, points_array, faces_array, N):
        n_facets, n_vertices = faces_array.shape[:2]
        assert n_vertices == 3
        moments_array = np.zeros([N+1, N+1, N+1])
        process_pool = mp.Pool()
        for face in faces_array:
            vertex_list = [points_array[_i, ...] for _i in face]
            process_pool.apply_async(_kmp_geometric_moments_exact_worker,
                                     args=(self, vertex_list, N),
                                     callback=moments_array.__iadd__,
                                     )
        process_pool.close()
        process_pool.join()
        return self.factorial_scalar(N)*moments_array


#DefaultPipeline = type('DefaultPipeline', (SerialPipeline,), {})
#DefaultPipeline = type(
#     'DefaultPipeline', (NumpyOptimizations, MultiprocPipeline,), {})
#DefaultPipeline = type(
#    'DefaultPipeline', (KoehlOptimizations, SerialPipeline), {})
DefaultPipeline = type(
    'DefaultPipeline', (KoehlMultiproc, SerialPipeline), {})
