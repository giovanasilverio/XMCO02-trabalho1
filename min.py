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
                qtd_arestas = int(partes[3])

            elif partes[0] == "n":
                if partes[2] == "s":
                    s = int(partes[1])
                else:
                    t = int(partes[1])

            else:
                arestas.append((int(partes[1]), int(partes[2])))
                capacidades.append(float(partes[3]))

    print(qtd_vertices, qtd_arestas, s, t)

    for i in range(len(arestas)):
        print(arestas[i], capacidades[i])

if __name__ == "__main__":
    instancia = input("Escolha uma instancia para execucao (1 - 5): ")

    ler_instancia(instancia)