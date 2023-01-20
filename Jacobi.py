# -*- coding: utf-8 -*-
"""
Criado en 05/08/21 as 22:00
@autor: Valdomiro
"""
# import numpy as np


def jacobi(H,w,solu_vector,num_iteracoes, TOL = 1e-5):
    n = 0 
    soma_colunas = 0
    while n < num_iteracoes:
        for i in range(len(H)):
            x = w[i]
            for j in range(len(H[0])):
                if i != j:         
                    soma_colunas += H[i][j]*solu_vector[j]
                    
            solu_vector[i] = (x-soma_colunas)/H[i][i]
            soma_colunas = 0
        n += 1
    return solu_vector
    