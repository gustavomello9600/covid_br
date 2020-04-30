import pandas as pd
import os


class StartFunction(dict):
    def __init__(self, value, family="S"):
        self[0] = value
        self.family = family

    def __hash__(self):
        return hash(self.family)

def contagem_de_mortes_se(r0=2.79, mortalidade=1.305, TL=42630, hosp=2.67, dL=0, CL=1.5):
    """
    Simula o número de mortes ao final de uma pandemia situada em território brasileiro para um dado cenário.

    Modelo da Simulação: SIRM
    Método de Resolução das Equações Diferenciais Ordinárias: Runge-Kutta
    """

    #Dados experimentais
    alfa = 8.414499908
    beta = 0.1969820003
    gama = -15.31551146
    mort = mortalidade/100
    R0 = r0
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
    taxa_de_hospitalização = hosp/100
    total_de_leitos = TL

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
    leitos_limitados = bool(TL)
    lockdown_vertical = bool(dL)
    steps = 10
    dt = 1/steps
    dias = 15

    #Variáveis da iteração
    total_de_mortos_sem_leito = 0
    d = -1

    while I[d + 1] > 1:
        d += 1

        #Salvando número de recuperados no dia anterior
        R_ant = R[d]

        #Implementação do algoritmo Runge-Kutta
        for i in range(0, steps):
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

        for f in (S, I, R, M):
            f[d + 1] = f.pop(t + 10000)

        #condição de leitos limitados
        if I[d + 1] * taxa_de_hospitalização > total_de_leitos and leitos_limitados:  # Roda se leitos não forem suficientes
            R[d] = R_ant
            Recuperariam = R[d + 1] - R[d]  # Número de pessoas que se recuperariam se houvessem leitos suficientes
            if Recuperariam > total_de_leitos / taxa_de_hospitalização:
                # Limita o número de recuperados dentre os que precisariam ser internados ao número de internações possíveis
                Recuperam = (1 - taxa_de_hospitalização) * Recuperariam + total_de_leitos
                M_semleito = Recuperariam - Recuperam
                total_de_mortos_sem_leito += M_semleito
                R[d + 1] = Recuperam + R[d]
                M[d + 1] += M_semleito

        #condição de lockdown vertical no dia 1 de Abril
        if d == dL and lockdown_vertical: kc = kc * CL

    return "{:.2f} milhões/ {:.2f} milhões".format(M[d + 1]/1e6,
                                                               total_de_mortos_sem_leito/1e6)

R0s = [1.5, 1.7, 1.8, 2, 2.5, 2.79, 3, 3.5]
morts = [0.2, 0.5, 0.7, 0.832, 1, 1.2, 1.5]
cenário = {}
cenário["A"] = {"TL": 0,
            "dL": 0}
cenário["B"] = {"TL": 42630,
            "hosp": 2.67,
            "dL": 0}
cenário["C"] = {"TL": 42630,
            "hosp": 1,
            "dL": 0}
cenário["D"] = {"TL": 42630/10,
            "hosp": 2.67,
            "dL": 0}
cenário["E"] = {"TL": 42630/10,
            "hosp": 1,
            "dL": 0}
cenário["F"] = {"TL": 42630,
            "hosp": 1,
            "dL": 18}
cenário["G"] = {"TL": 42630,
            "hosp": 1,
            "dL": 48}
cenário["H"] = {"TL": 42630,
            "hosp": 1,
            "dL": 78}
cenário["I"] = {"TL": 42630,
            "hosp": 1,
            "dL": 18,
            "CL": 2}
cenário["J"] = {"TL": 42630,
            "hosp": 1,
            "dL": 48,
            "CL": 2}
cenário["K"] = {"TL": 42630,
            "hosp": 1,
            "dL": 78,
            "CL": 2}
cenário["L"] = {"TL": 1e5,
            "hosp": 1,
            "dL": 0}

for índice in cenário:
    table = []
    for m in morts:
        line = []
        for r in R0s:
            line.append(contagem_de_mortes_se(r, m, **cenário[índice]))
        table.append(line)

    # Exporta os dados como um arquivo .csv
    data = pd.DataFrame(table, index=morts, columns=R0s)
    data.to_csv(os.getcwd() + "\\cenário" + índice + ".csv")
    print("Cenário {} pronto".format(índice))