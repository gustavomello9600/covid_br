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


#Estado Inicial
I = F(84.43875125)
S = F(-(I - 210000000))
R = F(0)
M = F(0)

#Operador diferencial com relação ao tempo

def d_dt(f):
    if f == S:
        df = -1.45092774571429E-09 * S.value * I.value
    elif f == I:
        df = 1.4509277E-09 * S.value * I.value - 1.04027002481127E-01*I.value - 5.18261278769005E-03*I.value
    elif f == R:
        df = 1.04027002481127E-01*I.value
    else:
        df = 5.18261278769005E-03*I.value
    return df

interval = "_seconds"
step_division = 24*60*60
dt = 1/step_division
table = []
dias = 15

for d in range(0, dias):
    for i in range(0, step_division):
        #Utilizando as diferenciais para obter novos valores das funções
        new = {}
        for f in S, I, R, M:
            new[f] = f.value + d_dt(f) * dt

        #Atualizando valores das funções
        for f in S, I, R, M:
            f.value = new[f]

    print("dia {}\n"
          "S:{:.3f}\n"
          "I:{:.3f}\n"
          "R:{:.3f}\n"
          "M:{:.3f}\n"
          "Soma: {:.3f}"
          "\n".format(d, S.value, I.value, R.value, M.value,
                      S.value+I.value+R.value+M.value))
    table.append([f.value for f in (S,I,R,M)])

data = pd.DataFrame(table, index=list(range(dias)), columns=list("SIRM"))
data.to_csv(os.getcwd() + "\\output" + interval + ".csv")