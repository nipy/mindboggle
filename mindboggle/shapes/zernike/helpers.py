import numpy as np

def nested_loop(stack, args):
    if len(stack) != 0:
        fn = stack.pop()
        for i in fn(*args):
            for j in nested_loop(stack, args+[i]):
                yield (i,)+j
        stack.append(fn)
    else:
        yield tuple()

def nest(*_stack):
    return nested_loop(list(reversed(_stack)), [])

def autocat(arrs, **dargs):
    axis = dargs.pop('axis', None)
    if axis is None:
        return np.concatenate(arrs, **dargs)
    ndim = arrs[0].ndim
    assert all([ a.ndim == ndim for a in arrs])
    if axis >= ndim:
        arrs = tuple([ np.expand_dims(a,axis) for a in arrs ])
    return np.concatenate(arrs, axis=axis)

def main():

    def zero(): return list(range(0, 3))
    def one(_x): return list(range(0, _x))
    def two(_x, _y): return list(range(0, _y))

    for i in nest(zero, one, two):
        print('{0} {1}'.format('#'*5, i))

if __name__ == '__main__':
    main()
