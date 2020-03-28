import numpy as np
import pandas as pd
import os


class F():
    def __init__(self, value):
        self.value = value

    def __add__(self, other):
        if isinstance(other, F):
            return self.value + other.value
        return self.value + other

    def __sub__(self, other):
        return self.value - other

    def __matmul__(self, other):
        return self.value * other


#Dados experimentais
alfa = 8.414499908
beta = 0.1969820003
gama = -15.31551146
I0 = 121
mort = 2/100 #livre para variar
R0 = 2.79 #livre para variar
S0 = 210000000
epsilon = (alfa*(R0 -1))/(I0 * mort) #Fator de correção do número de casos documentados no dia 0

#Constantes
kc = (beta * R0)/(S0*(R0 - 1))
fat = (alfa * beta)/(epsilon*I0)
kr = kc*S0 - beta - fat

#Estado Inicial
I = F(I0*epsilon)
R = F(-842.9092681)
M = F(-6.901011552E+00)
S = F(S0 - I.value - R.value - M.value)

#Operador diferencial com relação ao tempo

def d_dt(f):
    if f == S:
        df = -kc * S.value * I.value
    elif f == I:
        df = kc * S.value * I.value - kr*I.value - fat*I.value
    elif f == R:
        df = kr*I.value
    else:
        df = fat*I.value
    return df

sufixo = "_R2-79_mort0-02"
step_division = 100000
dt = 1/step_division
table = []
dias = 15

d = -1                                                                                    #modo até acabar infecção
while I.value > 1:                                                                        #modo até acabar infecção
#for d in range(0, dias):                                                                    #modo de dias específicos
    d += 1                                                                                #modo até acabar infecção
    # Joga os valores na tabela
    print("dia {}\n"
          "S:{:.3f}\n"
          "I:{:.3f}\n"
          "R:{:.3f}\n"
          "M:{:.3f}\n"
          "Soma: {:.3f}"
          "\n".format(d, S.value, I.value, R.value, M.value,
                      S.value + I.value + R.value + M.value))
    table.append([f.value for f in (S, I, R, M)])

    # Impede o último loop de executar desnecessariamente
#    if d == dias - 1:                                                                    #modo de dias específicos
#        break                                                                            #modo de dias específicos

    for i in range(0, step_division):
        #Utilizando as diferenciais para obter novos valores das funções
        new = {}
        for f in S, I, R, M:
            new[f] = f.value + d_dt(f) * dt

        #Atualizando valores das funções
        for f in S, I, R, M:
            f.value = new[f]

data = pd.DataFrame(table, index=list(range(d+1)), columns=list("SIRM"))
data.to_csv(os.getcwd() + "\\output" + sufixo + ".csv")