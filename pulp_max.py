import pulp
from max import ler_instancia


def resolver_com_pulp(instancia, variaveis=None):
    problema = pulp.LpProblem("fluxo_maximo_multiplas_mercadorias", pulp.LpMaximize)
    arestas = instancia["arestas"]
    chaves = [
        (mercadoria, indice)
        for mercadoria in range(1, instancia["num_mercadorias"] + 1)
        for indice in range(len(arestas))
    ]
    x = {
        (mercadoria, indice): pulp.LpVariable(f"x_{mercadoria}_{indice}", lowBound=0)
        for mercadoria, indice in chaves
    }

    problema += pulp.lpSum(
        instancia["premios"][mercadoria] * x[(mercadoria, indice)]
        for mercadoria, indice in chaves
        if arestas[indice][1] == instancia["destinos"][mercadoria]
    )

    for indice, (_, _, capacidade) in enumerate(arestas):
        problema += pulp.lpSum(
            x[(mercadoria, indice)]
            for mercadoria in range(1, instancia["num_mercadorias"] + 1)
        ) <= capacidade

    for mercadoria in range(1, instancia["num_mercadorias"] + 1):
        fonte = instancia["fontes"][mercadoria]
        destino = instancia["destinos"][mercadoria]
        for vertice in range(1, instancia["num_vertices"] + 1):
            if vertice in (fonte, destino):
                continue
            sai = pulp.lpSum(
                x[(mercadoria, indice)]
                for indice, (origem, _, _) in enumerate(arestas) if origem == vertice
            )
            entra = pulp.lpSum(
                x[(mercadoria, indice)]
                for indice, (_, destino_aresta, _) in enumerate(arestas)
                if destino_aresta == vertice
            )
            problema += sai - entra == 0

    problema.solve(pulp.PULP_CBC_CMD(msg=False))
    if pulp.LpStatus[problema.status] != "Optimal":
        raise RuntimeError(f"PuLP: {pulp.LpStatus[problema.status]}")

    valores = {chave: (variavel.value() or 0.0) for chave, variavel in x.items()}
    return pulp.value(problema.objective), valores


def comparar_com_pulp(instancia, variaveis, solucao):
    objetivo_pulp, valores_pulp = resolver_com_pulp(instancia, variaveis)
    print(f"Objetivo primal-dual: {solucao['objetivo_primal']:.4f}")
    print(f"Objetivo PuLP/CBC:    {objetivo_pulp:.4f}")
    return objetivo_pulp, valores_pulp


def main():
    caminho = "mc_instance5.max"
    instancia = ler_instancia(caminho)
    objetivo, valores = resolver_com_pulp(instancia)
    print(f"Valor objetivo (PuLP/CBC): {objetivo:.4f}")


if __name__ == "__main__":
    main()
