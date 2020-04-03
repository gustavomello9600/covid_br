import pandas as pd
import os


class StartFunction(dict):
    def __init__(self, value, family="S"):
        self[0] = value
        self.family = family

    def __hash__(self):
        return hash(self.family)


#Dados experimentais
alfa = 8.414499908
beta = 0.1969820003
gama = -15.31551146
mort = 0.823/100 #livre para variar
R0 = 2.79 #livre para variar
S0 = 210000000

#Constantes
M0 = alfa + gama
kc = (beta * R0)/(S0*(R0 - 1))
fat = mort * (beta/(R0 - 1))
kr = (1 - mort) * (beta/(R0 - 1))
I0 = alfa*beta/fat

#Estado Inicial
S, I, R, M = {},{},{},{}
I = StartFunction(I0, "I")
R = StartFunction(((1-mort)/mort) * M0, "R")
M = StartFunction(M0, "M")
S = StartFunction(S0 - I[0] - R[0] - M[0], "S")

#Realidade dos leitos brasileiros
taxa_de_hospitalização = 2.67/100
total_de_leitos = 42630

#Diferenciais
def d_dt(f):
    if f == S:
        return lambda s, i, r, m, t: -kc*s*i
    elif f == I:
        return lambda s, i, r, m, t: (kc*s - kr - fat)*i
    elif f == R:
        return lambda s, i, r, m, t: kr*i
    else:
        return lambda s, i, r, m, t: fat * i

#Condições da simulação
leitos_limitados = True
lockdown_vertical = False
steps = 10
dt = 1/steps
dias = 15

#Variáveis da iteração
total_de_mortos_sem_leito = 0
table = []
d = -1

while I[d + 1] > 1:
    d += 1
    # Joga os valores na tabela
    print("dia {}\n"
          "S:{:.3f}\n"
          "I:{:.3f}\n"
          "R:{:.3f}\n"
          "M:{:.3f}\n"
          "Soma: {:.3f}"
          "\n".format(d, S[d], I[d], R[d], M[d],
                      S[d] + I[d] + R[d] + M[d]))
    table.append([f[d] for f in (S, I, R, M)])

    R_ant = R[d]

    #Implementação do algoritmo Runge-Kutta
    #Tem o efeito de atualizar o valor de cada f para o dia d
    for i in range(0, steps):
        #try:
        t = 10000*i + d
        k = dict()
        #Determinação dos k1s
        for f in (S, I, R, M):
            k[(f, 1)] = dt*(d_dt(f)(S[t], I[t], R[t], M[t], t))
        #Determinação dos k2s
        for f in (S, I, R, M):
            k[(f, 2)] = dt*(d_dt(f)(S[t] + k[(S, 1)]/2,
                                    I[t] + k[(I, 1)]/2,
                                    R[t] + k[(R, 1)]/2,
                                    M[t] + k[(M, 1)]/2,
                                    t + dt/2))
        #Determinação dos k3s
        for f in (S, I, R, M):
            k[(f, 3)] = dt*(d_dt(f)(S[t] + k[(S, 2)]/2,
                                    I[t] + k[(I, 2)]/2,
                                    R[t] + k[(R, 2)]/2,
                                    M[t] + k[(M, 2)]/2,
                                    t + dt/2))
        #Determinação dos k4s
        for f in (S, I, R, M):
            k[(f, 4)] = dt * (d_dt(f)(S[t] + k[(S, 3)],
                                      I[t] + k[(I, 3)],
                                      R[t] + k[(R, 3)],
                                      M[t] + k[(M, 3)],
                                      t + dt))
        #Atualização dos valores das funções
        for f in (S, I, R, M):
            f[t + 10000] = f.pop(t) + (k.pop((f, 1)) + 2*k.pop((f, 2)) + 2*k.pop((f, 3)) + k.pop((f, 4)))/6
        #except:
         #   breakpoint()

    for f in (S, I, R, M):
        f[d + 1] = f.pop(t + 10000)

    if I[d + 1] * taxa_de_hospitalização > total_de_leitos and leitos_limitados:  # Roda se leitos não forem suficientes
        R[d] = R_ant
        Recuperariam = R[d + 1] - R[d]  # Número de pessoas que se recuperariam se houvessem leitos suficientes
        if Recuperariam > total_de_leitos / taxa_de_hospitalização:
            # Limita o número de recuperados dentre os que precisariam ser internados ao número de internações possíveis
            Recuperam = (1 - taxa_de_hospitalização) * Recuperariam + total_de_leitos
            M_semleito = Recuperariam - Recuperam
            total_de_mortos_sem_leito += M_semleito
            print("No dia {}, morreram {:.0f} pessoas sem leito".format(d + 1, M_semleito))
            R[d + 1] = Recuperam + R[d]
            M[d + 1] += M_semleito

    if d == 18 and lockdown_vertical: kc = kc * 1.5

print("{} pessoas morreram sem leito".format(total_de_mortos_sem_leito))

#Exporta os dados como um arquivo .csv
sufixo = "RK_R{}_mort{}_{}leitos".format(*["{:.2f}".format(n).replace(".", "-") for n in (R0, 100*mort)],total_de_leitos)
data = pd.DataFrame(table, index=list(range(d+1)), columns=list("SIRM"))
data.to_csv(os.getcwd() + "\\output" + sufixo + ".csv")