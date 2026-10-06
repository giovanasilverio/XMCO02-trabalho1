from fractions import Fraction

def ler_instancia(caminho):
    num_vertices = 0
    num_arestas = 0
    num_mercadorias = 0

    arestas = []
    fontes = {}
    destinos = {}
    premios = {}

    with open(caminho, "r", encoding="utf-8") as arquivo:
        for linha in arquivo:
            linha = linha.strip()

            if linha == "":
                continue

            if linha[0] == "c":
                continue

            partes = linha.split()

            if linha[0] == "p":
                num_vertices = int(partes[2])
                num_arestas = int(partes[3])
                num_mercadorias = int(partes[4]) 
                
            elif linha[0] == "n":
                vertice = int(partes[1])
                tipo_no = partes[2]
                id_mercadoria = int(partes[3])

                if tipo_no == "s":
                    fontes[id_mercadoria] = vertice

                elif tipo_no == "t":
                    destinos[id_mercadoria] = vertice

            elif linha[0] == "w":
                id_mercadoria = int(partes[1])
                premio = int(partes[2])

                premios[id_mercadoria] = premio

            elif linha[0] == "a":
                origem = int(partes[1])
                destino = int(partes[2])
                capacidade = int(partes[3])

                arestas.append(
                    (origem, destino, capacidade)
                )


    instancia = {
        "num_vertices": num_vertices,
        "num_arestas": num_arestas,
        "num_mercadorias": num_mercadorias,
        "arestas": arestas,
        "fontes": fontes,
        "destinos": destinos,
        "premios": premios
    }

    return instancia


def construir_variavel(instancia):
    arestas = instancia["arestas"]
    num_mercadorias = instancia["num_mercadorias"]

    variaveis = {}
    coluna = 0

    for mercadoria in range(1, num_mercadorias + 1):
        for origem, destino, _ in arestas:
            variaveis[
                (mercadoria, origem, destino)
            ] = coluna
            coluna += 1

    return variaveis


def construir_objetivo(instancia, variaveis):
    destinos = instancia["destinos"]
    premios = instancia["premios"]

    vetor = [0] * len(variaveis)

    for chave, coluna in variaveis.items():
        mercadoria, origem, destino_aresta = chave

        destino_mercadoria = destinos[mercadoria]

        if destino_aresta == destino_mercadoria:
            vetor[coluna] = premios[mercadoria]

    return vetor

def construir_capacidades(instancia, variaveis):
    arestas = instancia["arestas"]
    linhas = []
    lados_direitos = []

    for origem, destino, capacidade in arestas:
        linha = [0] * len(variaveis)

        for mercadoria in range(1, instancia["num_mercadorias"] + 1):
            coluna = variaveis[(mercadoria, origem, destino)]
            linha[coluna] = 1

        linhas.append(linha)
        lados_direitos.append(capacidade)

    return linhas, lados_direitos


def construir_conservacao(instancia, variaveis):
    arestas = instancia["arestas"]
    linhas = []
    lados_direitos = []
    
    for mercadoria in range(1, instancia["num_mercadorias"] + 1):
        for vertice in range(1, instancia["num_vertices"] + 1):
            if vertice == instancia["fontes"][mercadoria] or vertice == instancia["destinos"][mercadoria]:
                continue
            linha = [0] * len(variaveis)
            for origem, destino, capacidade in arestas:
                coluna = variaveis[(mercadoria, origem, destino)]
                if origem == vertice:
                    linha[coluna] = +1
                if destino == vertice:
                    linha[coluna] = -1   
                                
            linhas.append(linha)
            lados_direitos.append(0)    
                
    return linhas, lados_direitos

def construir_modelo_fluxo(instancia):
    variaveis = construir_variavel(instancia)
    objetivo = construir_objetivo(instancia, variaveis)
    linhas_capacidade, limites_capacidade = construir_capacidades(instancia, variaveis)
    linhas_conservacao, zeros_conservacao = construir_conservacao(instancia, variaveis)
    return variaveis, objetivo, linhas_capacidade, limites_capacidade, linhas_conservacao, zeros_conservacao

def _resolver_rsp(colunas, b, max_pivos=10000, tol=1e-9):
    """Minimiza a soma das artificiais com pivoteamento Simplex manual."""
    m = len(b)
    n = len(colunas)
    sinais = [1 if valor >= 0 else -1 for valor in b]
    tableau = []
    for i in range(m):
        linha = [sinais[i] * coluna[i] for coluna in colunas]
        linha += [1.0 if i == j else 0.0 for j in range(m)]
        linha.append(abs(b[i]))
        tableau.append(linha)

    base = [n + i for i in range(m)]
    custos = [0.0] * n + [1.0] * m
    total_colunas = n + m

    for _ in range(max_pivos):
        reduzidos = [
            custos[j] - sum(custos[base[i]] * tableau[i][j] for i in range(m))
            for j in range(total_colunas)
        ]
        # Regra de Bland: primeira variável com custo reduzido negativo.
        entra = next((j for j in range(total_colunas) if reduzidos[j] < -tol), None)
        if entra is None:
            break

        possiveis = [i for i in range(m) if tableau[i][entra] > tol]
        if not possiveis:
            raise RuntimeError("Subproblema restrito ilimitado")
        sai = min(possiveis, key=lambda i: (tableau[i][-1] / tableau[i][entra], base[i]))

        pivo = tableau[sai][entra]
        tableau[sai] = [valor / pivo for valor in tableau[sai]]
        for i in range(m):
            if i != sai and abs(tableau[i][entra]) > tol:
                fator = tableau[i][entra]
                tableau[i] = [tableau[i][j] - fator * tableau[sai][j]
                              for j in range(total_colunas + 1)]
        base[sai] = entra
    else:
        raise RuntimeError("Limite de pivôs do subproblema atingido")

    w_otimo = sum(custos[base[i]] * tableau[i][-1] for i in range(m))
    valores = [0.0] * n
    for i in range(m):
        if base[i] < n:
            valores[base[i]] = tableau[i][-1]

    # Multiplicadores do RSP: custo reduzido da artificial = 1 - pi'_i.
    pi = [sinais[i] * (1 - reduzidos[n + i]) for i in range(m)]
    return w_otimo, valores, pi


def resolver_primal_dual(modelo, max_iteracoes=1000, tol=1e-8):
    (variaveis, objetivo, linhas_capacidade, limites_capacidade,
     linhas_conservacao, zeros_conservacao) = modelo
    n = len(objetivo)
    m1 = len(linhas_capacidade)
    m2 = len(linhas_conservacao)
    m = m1 + m2

    # min -c'x, com capacidades escritas como -Ax >= -b.
    A = [[-float(v) for v in linha] for linha in linhas_capacidade]
    A += [[float(v) for v in linha] for linha in linhas_conservacao]
    b = [-float(v) for v in limites_capacidade] + [float(v) for v in zeros_conservacao]
    custos = [-float(v) for v in objetivo] + [0.0] * m1

    # Uma coluna para cada fluxo, seguida pelas colunas das variáveis de excesso.
    colunas = [[A[i][j] for i in range(m)] for j in range(n)]
    for i in range(m1):
        coluna = [0.0] * m
        coluna[i] = -1.0
        colunas.append(coluna)

    def custos_reduzidos(y):
        return [custos[j] - sum(y[i] * coluna[i] for i in range(m))
                for j, coluna in enumerate(colunas)]

    # y = 0 não serve se há prêmios positivos. Este M garante dual-factibilidade.
    M = max([0.0] + [float(v) for v in objetivo])
    y = [M] * m1 + [0.0] * m2
    historico = []
    status = "max_iteracoes_atingido"
    J = []
    valores = []

    for iteracao in range(1, max_iteracoes + 1):
        d = custos_reduzidos(y)
        J = [j for j in range(len(colunas)) if abs(d[j]) <= tol]
        w_otimo, valores, pi = _resolver_rsp([colunas[j] for j in J], b)
        historico.append((iteracao, len(J), w_otimo,
                          sum(limites_capacidade[i] * y[i] for i in range(m1))))
        if w_otimo <= tol:
            status = "otimo"
            break

        # Maior passo que preserva os custos reduzidos não negativos.
        passos = []
        for j, coluna in enumerate(colunas):
            alfa = sum(pi[i] * coluna[i] for i in range(m))
            if d[j] > tol and alfa > tol:
                passos.append(d[j] / alfa)
        if not passos:
            status = "infactivel"
            break
        theta = min(passos)
        y = [y[i] + theta * pi[i] for i in range(m)]

    x = [0.0] * n
    if status == "otimo":
        for posicao, j in enumerate(J):
            if j < n:
                x[j] = valores[posicao]

    y_capacidade = y[:m1]
    y_conservacao = [-valor for valor in y[m1:]]
    return {
        "x": x,
        "y_capacidade": y_capacidade,
        "y_conservacao": y_conservacao,
        "objetivo_primal": sum(objetivo[j] * x[j] for j in range(n)),
        "objetivo_dual": sum(limites_capacidade[i] * y_capacidade[i]
                              for i in range(m1)),
        "status": status,
        "iteracoes": len(historico),
        "historico": historico,
    }

def validar_solucao(instancia, variaveis, solucao, tol=1e-5):
    x = solucao["x"]
    arestas = instancia["arestas"]
    num_mercadorias = instancia["num_mercadorias"]

    for coluna in variaveis.values():
        if x[coluna] < -tol:
            return False

    for origem, destino, capacidade in arestas:
        total = sum(
            x[variaveis[(mercadoria, origem, destino)]]
            for mercadoria in range(1, num_mercadorias + 1)
        )
        if total > capacidade + tol:
            return False

    for mercadoria in range(1, num_mercadorias + 1):
        fonte = instancia["fontes"][mercadoria]
        destino_merc = instancia["destinos"][mercadoria]
        for vertice in range(1, instancia["num_vertices"] + 1):
            if vertice in (fonte, destino_merc):
                continue
            entra = sum(
                x[variaveis[(mercadoria, o, d)]]
                for o, d, _ in arestas if d == vertice
            )
            sai = sum(
                x[variaveis[(mercadoria, o, d)]]
                for o, d, _ in arestas if o == vertice
            )
            if abs(sai - entra) > tol:
                return False

    objetivo_recalculado = 0.0
    for mercadoria in range(1, num_mercadorias + 1):
        destino_merc = instancia["destinos"][mercadoria]
        premio = instancia["premios"][mercadoria]
        for origem, destino, _ in arestas:
            if destino == destino_merc:
                objetivo_recalculado += premio * x[variaveis[(mercadoria, origem, destino)]]

    if abs(objetivo_recalculado - solucao["objetivo_primal"]) > 1e-3:
        return False

    return True


def exibir_resultado(instancia, variaveis, solucao, valido):
    x = solucao["x"]
 
    print("=" * 70)
    print("SOLUÇÃO — MÉTODO PRIMAL-DUAL DO SIMPLEX")
    print("=" * 70)
    print(f"Status: {solucao['status']}   (iterações: {solucao['iteracoes']})")
 
    print("\n-- Variáveis primais não-nulas (fluxo por mercadoria e aresta) --")
    algum = False
    for (mercadoria, origem, destino), coluna in sorted(variaveis.items()):
        valor = x[coluna]
        if abs(valor) > 1e-6:
            algum = True
            print(f"  f[mercadoria={mercadoria}, {origem} -> {destino}] = {valor:.4f}")
    if not algum:
        print("  (nenhuma variável não-nula)")
 
    print("\n-- Variáveis duais das restrições de capacidade (y_uv >= 0) --")
    algum = False
    for i, (origem, destino, capacidade) in enumerate(instancia["arestas"]):
        y = solucao["y_capacidade"][i]
        if abs(y) > 1e-6:
            algum = True
            print(f"  y[{origem} -> {destino}] (capacidade {capacidade}) = {y:.4f}")
    if not algum:
        print("  (nenhuma variável dual de capacidade ativa)")
 
    print("\n-- Variáveis duais das restrições de conservação (livres) --")
    indice = 0
    algum = False
    for mercadoria in range(1, instancia["num_mercadorias"] + 1):
        fonte = instancia["fontes"][mercadoria]
        destino_merc = instancia["destinos"][mercadoria]
        for vertice in range(1, instancia["num_vertices"] + 1):
            if vertice in (fonte, destino_merc):
                continue
            valor = solucao["y_conservacao"][indice]
            if abs(valor) > 1e-6:
                algum = True
                print(f"  pi[mercadoria={mercadoria}, vértice={vertice}] = {valor:.4f}")
            indice += 1
    if not algum:
        print("  (nenhuma variável dual de conservação não-nula)")
 
    print(f"\nValor ótimo (objetivo primal): {solucao['objetivo_primal']:.4f}")
    print(f"Valor ótimo (objetivo dual):   {solucao['objetivo_dual']:.4f}")
 
    print("\nSolução válida" if valido else "\nSolução inválida")

def main():
    instancia = ler_instancia("mc_instance3.max")
    modelo = construir_modelo_fluxo(instancia)
    variaveis = modelo[0]

    print(f"  vértices={instancia['num_vertices']}  arestas={instancia['num_arestas']}  "
          f"mercadorias={instancia['num_mercadorias']}")
    print(f"  fontes={instancia['fontes']}  destinos={instancia['destinos']}  "
          f"premios={instancia['premios']}")
    print(f"  quantidade de variáveis do modelo: {len(variaveis)}")
 
    solucao = resolver_primal_dual(modelo)
    valido = validar_solucao(instancia, variaveis, solucao)
    exibir_resultado(instancia, variaveis, solucao, valido)

if __name__ == "__main__":
    main()