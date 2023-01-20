
"""
Criado en 05/08/21 as 22:00
@autor: Valdomiro
"""
import numpy as np    
def secante(f, x0, x1, TOL, MAX_IT, e, M):
    x2 = 0
    i = 2
    assert(x1 != x0), "intervalos dados são iguais"
    
    while(np.abs(x1-x0) > TOL):
        x2 = (x0*f(M,x1,e)-x1*f(M,x0,e))/(f(M,x1,e)-f(M,x0,e))
        x0 = x1
        x1 = x2
        i += 1
        assert(i < MAX_IT), "maximo de iterações atingidas"
    return x2
