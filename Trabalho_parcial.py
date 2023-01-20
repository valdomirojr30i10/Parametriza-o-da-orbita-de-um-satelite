"""
Criado en 05/08/21 as 22:00
@autor: Valdomiro
"""
import numpy as np
import matplotlib.pyplot as plt
import Jacobi as jacobi
import secante as metodo


# Pontos recebidos em kilometros informando 3 posicoes em que o satelite esteve
x = [20621.3, 34642.3, 21168.5]
y = [5214.2009061, 3201.7095080, 5193.6201775]

# Cria um sistea linear na forma: Ax[i]**2 + Bx[i] + C = y[i]**2 
def SL(x,y):
    # Conversao da unidade de kilometros para metros:
    x = [i*1e3 for i in x]
    y = [i*1e3 for i in y]
    
    H = np.array([[x[0]**2,x[0],1],
                  [x[1]**2,x[1],1],
                  [x[2]**2,x[2],1]], dtype = np.float64)    
    w = np.array([y[0]**2,y[1]**2,y[2]**2], dtype = np.float64)
    
    return H,w

# apos obter os valores de A,B,C atraves do metodo de jacobi para 
# solucao de sistemas lineares, os parametros sao encontrados na
# funcao abaixo 
def gera_parametros(A,B,C):
    a = ((((-1)*B)/(2*A))**2 - (C/A))**(1/2)
    b = a*(-A)**(1/2)
    xc = ((-1)*B)/(2*A)
    return a,b,xc
    
# Desenho feito a partir dos parametros 'a' e 'b' da elipse:
def desenho_mov_orbital(a,b,xc,f):
    t = np.linspace(0,2*np.pi,100000)
    X = [a*np.cos(i) + xc for i in t]
    Y = [b*np.sin(i)      for i in t]
    plt.plot(X,Y)
    plt.show()
    
# Calculo de anomalias medias para 30 pontos diferentes:
def anomalias_medias(a,b,xc,G,massa_terra):
    n = 30
    P = 2*np.pi*(a**3/(G*massa_terra))
    p_sobre_n = P/n
    t = p_sobre_n
    m = np.zeros(n)
    for i in range(n):
        m[i] = ((2*np.pi)/P)*t
        t += p_sobre_n
    return m

# Funcao a ser usada no metodo da secante: f(E) = M + esinE - E = 0
def Exc(M,E,e):
    return M + (e*np.sin(E)) - E


def ponto_fixo(Exc,i0,T,N):
    i = i0 + 1
    n = 0 
    while(i-i0 > T):
        i = (Exc(i0))
        i0 = i
        n += 1
        assert (n < N), "excedeu o numero maximo de iteracoes"
    return i
    
def anomalias_excentricas(m,a,b,e):
    Tam = len(m)
    E = np.zeros(Tam)
    for i in range(Tam):
        E[i] = metodo.secante(Exc,0,3,1e-5,500, e, m[i])
    return E,e

def anomalias_verdadeiras(E,e):
    Tam = len(E)
    V = np.zeros(Tam)
    for i in range(Tam):
        V[i] = ((1+e)/(1-e))**(1/2) * np.tan(E[i]/2)
    return V

def pos_satelite(V,a,e,xc,f):
    Tam = len(V)
    r = np.zeros(Tam)
    ang = np.zeros(Tam)
    x_final = np.zeros(Tam)
    y_final = np.zeros(Tam)
    
    for i in range(Tam):
        ang[i] = 2*np.arctan(V[i])
        r[i] = a*((1-e**2)/(1+e*np.cos(ang[i])))
        x_final[i] = r[i]*np.cos(ang[i]) + xc + f
        y_final[i] = r[i]*np.sin(ang[i])
        
    return x_final, y_final

def desenho_dos_pontos(x_final, y_final):
    plt.plot(x_final,y_final, 'o')
    plt.show()
def imp_dados(a,b,xc,A,B,C):
    print("valores de A,B,C encontrados: \n{}\n{}\n{}".format( A,B,C))
    print("\nParametros encontrados [a,b,xc]: \n{}\n{}\n{}".format(a,b,xc))
    print("excentricidade encontrada: ", (1- (b**2/a**2))**(1/2))
def t_parcial(x,y):
    H,w = SL(x,y)
    A,B,C = jacobi.jacobi(H,w,np.zeros(len(w)),10000)
    a,b,xc = gera_parametros(A, B, C)
    e = (1- (b**2/a**2))**(1/2)
    desenho_mov_orbital(a, b, xc, a*e)
    m = anomalias_medias(a,b,xc,6.67408e-11,5.972e24)
    E, e = anomalias_excentricas(m, a, b, e)
    V = anomalias_verdadeiras(E,e)
    x_final, y_final = pos_satelite(V, a, e, xc, a*e)
    desenho_dos_pontos(x_final,y_final)
    imp_dados(a, b, xc, A, B, C)
t_parcial(x, y) 
