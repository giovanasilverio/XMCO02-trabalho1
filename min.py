import numpy as np

TOL = 1e-9
MAX_ITER = 10000

qtd_vertices = None
qtd_arestas = None

s = None
t = None

arestas = []
capacidades = []

def ler_instancia(i):
    global qtd_vertices, qtd_arestas, s, t, arestas, capacidades

    with open(f"instance{i}.min", "r", encoding="utf-8") as f:
        for linha in f:
            linha = linha.strip()

            if not linha or linha.startswith("c"):
                continue

            partes = linha.split()

            if partes[0] == "p":
                qtd_vertices = int(partes[2])

            elif partes[0] == "n":
                if partes[2] == "s":
                    s = int(partes[1])
                else:
                    t = int(partes[1])

            else:
                arestas.append((int(partes[1]), int(partes[2])))
                capacidades.append(float(partes[3]))

    qtd_arestas = len(arestas)


def construir_modelo():
    qtd_variaveis = qtd_vertices + qtd_arestas

    A = []
    B = []
    C = np.zeros(qtd_variaveis)

    for e, (u, v) in enumerate(arestas):

        # cortek - du + dv >= 0
        linha = np.zeros(qtd_variaveis)
        linha[u - 1] -= 1.0
        linha[v - 1] += 1.0
        linha[qtd_vertices + e] = 1.0
        A.append(linha)
        B.append(0.0)

        # cortek + du - dv >= 0
        linha = np.zeros(qtd_variaveis)
        linha[u - 1] += 1.0
        linha[v - 1] -= 1.0
        linha[qtd_vertices + e] = 1.0
        A.append(linha)
        B.append(0.0)

        C[qtd_vertices + e] = capacidades[e]

    # ds = 0
    linha = [0.0] * qtd_variaveis
    linha[s - 1] = 1.0
    A.append(linha)
    B.append(0.0)

    # dt = 1
    linha = [0.0] * qtd_variaveis
    linha[t - 1] = 1.0
    A.append(linha)
    B.append(1.0)

    return A, B, C

def custos_reduzidos(A, C, y):
    return [
        C[j] - sum(y[i] * A[i][j] for i in range(len(y)))
        for j in range(len(C))
    ]

def simplex(A, B, C, base):
    m = len(B)
    n = len(C)
    tableau = [A[i][:] + [B[i]] for i in range(m)]

    for _ in range(MAX_ITER):
        reduzidos = [
            C[j] - sum(C[base[i]] * tableau[i][j] for i in range(m))
            for j in range(n)
        ]

        entra = next(
            (j for j in range(n) if reduzidos[j] > TOL),
            None
        )

        if entra is None:
            break

        possiveis = [
            i for i in range(m)
            if tableau[i][entra] > TOL
        ]

        sai = min(
            possiveis,
            key=lambda i: (
                tableau[i][-1] / tableau[i][entra],
                base[i]
            )
        )

        pivo = tableau[sai][entra]
        tableau[sai] = [valor / pivo for valor in tableau[sai]]

        for i in range(m):
            if i != sai:
                fator = tableau[i][entra]

                if abs(fator) > TOL:
                    tableau[i] = [
                        tableau[i][j] - fator * tableau[sai][j]
                        for j in range(n + 1)
                    ]

        base[sai] = entra

    valores = np.zeros(n)

    for i in range(m):
        if base[i] < n:
            valores[base[i]] = tableau[i][-1]

    valor = sum(C[j] * valores[j] for j in range(n))

    return valor, valores, base, tableau

def resolver_rsp(A, B, y, J):
    linhas = []
    rhs = []
    tipos = []
    qtd_desigualdades = 2 * qtd_arestas

    for i in range(len(B)):
        linha = [A[i][j] for j in J]

        if i >= qtd_desigualdades or y[i] > TOL:
            linhas.append(linha)
            rhs.append(B[i])
            tipos.append("=")
        else:
            linhas.append([-valor for valor in linha])
            rhs.append(-B[i])
            tipos.append("<=")

    m = len(rhs)
    n = len(J)
    matriz = [linha[:] for linha in linhas]
    base = []
    artificiais = []

    for i in range(m):
        indice = len(matriz[0])
        for k in range(m):
            matriz[k].append(1.0 if k == i else 0.0)

        base.append(indice)

        if tipos[i] == "=":
            artificiais.append(indice)

    custos = np.zeros(len(matriz[0]))

    for j in artificiais:
        custos[j] = -1.0

    valor, valores, base, tableau = simplex(matriz, rhs, custos, base)

    if -valor <= TOL:
        x_restrito = np.zeros(n)

        for i in range(m):
            if base[i] < n:
                x_restrito[base[i]] = tableau[i][-1]

        return {
            "factivel": True,
            "w": 0.0,
            "x_restrito": x_restrito,
            "pi": None
        }

    mapa = []

    for i in range(len(B)):
        if i >= qtd_desigualdades or y[i] > TOL:
            mapa.append((i, 1.0))
            mapa.append((i, -1.0))
        else:
            mapa.append((i, 1.0))

    linhas = []

    for j in J:
        linhas.append([
            sinal * A[i][j]
            for i, sinal in mapa
        ])

    linhas.append([1.0] * len(mapa))
    rhs = [0.0] * len(J) + [1.0]

    n_pi = len(mapa)
    n_cert = len(rhs)
    matriz = []

    for i in range(n_cert):
        linha = linhas[i][:] + [0.0] * n_cert
        linha[n_pi + i] = 1.0
        matriz.append(linha)

    custos = [
        B[i] * sinal
        for i, sinal in mapa
    ] + [0.0] * n_cert

    base = [n_pi + i for i in range(n_cert)]

    w, valores, _, _ = simplex(matriz, rhs, custos, base)

    pi = np.zeros(len(B))

    for valor_variavel, (i, sinal) in zip(valores, mapa):
        pi[i] += sinal * valor_variavel

    return {
        "factivel": False,
        "w": w,
        "x_restrito": None,
        "pi": pi
    }

def resolver_primal_dual(A, B, C):
    qtd_variaveis = len(C)
    qtd_restricoes = len(B)

    y = np.zeros(qtd_restricoes)
    historico = []

    for iteracao in range(1, MAX_ITER + 1):
        d = custos_reduzidos(A, C, y)

        J = [
            j for j in range(qtd_variaveis)
            if abs(d[j]) <= TOL
        ]

        rsp = resolver_rsp(A, B, y, J)

        valor_dual = sum(B[i] * y[i] for i in range(qtd_restricoes))

        historico.append({
            "iteracao": iteracao,
            "tamanho_J": len(J),
            "w": rsp["w"],
            "valor_dual": valor_dual
        })

        if rsp["factivel"]:
            x = np.zeros(qtd_variaveis)

            for posicao, j in enumerate(J):
                x[j] = rsp["x_restrito"][posicao]

            return {
                "x": x,
                "y": y,
                "objetivo_primal": sum(C[j] * x[j] for j in range(qtd_variaveis)),
                "objetivo_dual": valor_dual,
                "iteracoes": iteracao,
                "historico": historico,
                "J": J
            }

        pi = rsp["pi"]
        passos = []

        for j in range(qtd_variaveis):
            alfa = sum(A[i][j] * pi[i] for i in range(qtd_restricoes))

            if d[j] > TOL and alfa > TOL:
                passos.append(d[j] / alfa)

        for i in range(2 * qtd_arestas):
            if pi[i] < -TOL:
                passos.append(y[i] / (-pi[i]))

        theta = min(passos)
        y = [y[i] + theta * pi[i] for i in range(qtd_restricoes)]

if __name__ == "__main__":
    instancia = input("Escolha uma instancia para execucao (1 - 5): ")

    ler_instancia(instancia)
    A, B, C = construir_modelo()

    resultado = resolver_primal_dual(A, B, C)

    print(resultado["x"][:qtd_vertices])
    print(resultado["x"][qtd_vertices:])
    print(resultado["y"])
    print(resultado["objetivo_primal"])
    print(resultado["objetivo_dual"])
    print(resultado["iteracoes"])
    print(resultado["historico"])