import numpy as np
import numpy.linalg as lin
from numpy.linalg import norm
from itertools import takewhile


def iterate(func, start):
    while True:
        yield start
        start = func(*start)

def steep(A: np, b, alfa, dir, x0=None, itmax: int=100, toll: float=1e-10):
    err = lambda num, den: norm(num)/norm(den)

    if x0 is None:
        x0 = np.zeros_like(b)
        
    def stop_criteria(i, args):
        res, _, _ = args
        return i < itmax and err(res, b) >= toll
    
    def iteration(res, x, p):
        Ap = A @ p
        a = alfa(A, res, p)
        r = res + a*Ap
        return r, x + a*p, dir(r, res, p)
    
    res = A @ x0 - b
    
    iter = iterate(iteration, (res, x0, -res))
    arr = list(dict(takewhile(lambda param: stop_criteria(*param), enumerate(iter))).values()) + \
        [next(iter)] 
    
    x_arr = [x for res, x, p in arr]
    r_arr = [err(res, b) for res, x, p in arr]
    return x_arr[-1], len(x_arr), r_arr 