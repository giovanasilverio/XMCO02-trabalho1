import numpy as np

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
        linha = [0.0] * qtd_variaveis
        linha[u - 1] -= 1.0
        linha[v - 1] += 1.0
        linha[qtd_vertices + e] = 1.0
        A.append(linha)
        B.append(0.0)

        # cortek + du - dv >= 0
        linha = [0.0] * qtd_variaveis
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

    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)

    print("A:\n", A)
    print("B:\n", B)
    print("C:\n", C)

if __name__ == "__main__":
    instancia = input("Escolha uma instancia para execucao (1 - 5): ")

    ler_instancia(instancia)

    np.set_printoptions(linewidth=300)
    construir_modelo()