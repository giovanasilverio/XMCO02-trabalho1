import pulp

def ler_instancia(caminho):
    qtd_vertices = None
    s = None
    t = None
    arestas = []
    capacidades = []

    with open(caminho, "r", encoding="utf-8") as f:
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

    return {
        "qtd_vertices": qtd_vertices,
        "s": s,
        "t": t,
        "arestas": arestas,
        "capacidades": capacidades,
    }

def resolver_com_pulp(instancia):
    qtd_vertices = instancia["qtd_vertices"]
    s = instancia["s"]
    t = instancia["t"]
    arestas = instancia["arestas"]
    capacidades = instancia["capacidades"]

    problema = pulp.LpProblem("corte_minimo", pulp.LpMinimize)

    d = {
        v: pulp.LpVariable(f"d_{v}", lowBound=0)
        for v in range(1, qtd_vertices + 1)
    }

    x = {
        e: pulp.LpVariable(f"x_{e+1}", lowBound=0)
        for e in range(len(arestas))
    }

    problema += pulp.lpSum(capacidades[e] * x[e] for e in range(len(arestas)))

    for e, (u, v) in enumerate(arestas):
        problema += x[e] - d[u] + d[v] >= 0, f"corte_pos_{e+1}"
        problema += x[e] + d[u] - d[v] >= 0, f"corte_neg_{e+1}"

    problema += d[s] == 0, "fixa_s"
    problema += d[t] == 1, "fixa_t"

    problema.solve(pulp.PULP_CBC_CMD(msg=False))

    if pulp.LpStatus[problema.status] != "Optimal":
        raise RuntimeError(f"PuLP: Status {pulp.LpStatus[problema.status]}")

    valores_d = {v: d[v].value() or 0.0 for v in range(1, qtd_vertices + 1)}
    valores_x = {e: x[e].value() or 0.0 for e in range(len(arestas))}

    # Construção da partição S e T
    S = [v for v in range(1, qtd_vertices + 1) if valores_d[v] <= 0.5]
    T = [v for v in range(1, qtd_vertices + 1) if valores_d[v] > 0.5]

    return {
        "objetivo": pulp.value(problema.objective),
        "d": valores_d,
        "x": valores_x,
        "S": S,
        "T": T,
    }


if __name__ == "__main__":
    instancia_num = input("Escolha uma instancia para execucao (1 - 5): ")
    caminho = f"instance{instancia_num}.min"

    instancia = ler_instancia(caminho)
    solucao = resolver_com_pulp(instancia)

    print("\n--- Resultado PuLP/CBC ---")
    print(f"Valor Otimizacao (Corte Minimo): {solucao['objetivo']:.10g}")
    print(f"Conjunto S: {solucao['S']}")
    print(f"Conjunto T: {solucao['T']}")