import numpy as np

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
            if linha == "" or linha[0] == "c":
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
                premios[int(partes[1])] = float(partes[2])

            elif linha[0] == "a":
                arestas.append((int(partes[1]), int(partes[2]), float(partes[3])))

    return {
        "num_vertices": num_vertices,
        "num_arestas": len(arestas),
        "num_mercadorias": num_mercadorias,
        "arestas": arestas,
        "fontes": fontes,
        "destinos": destinos,
        "premios": premios,
    }


def construir_variavel(instancia):
    variaveis = {}
    coluna = 0
    for mercadoria in range(1, instancia["num_mercadorias"] + 1):
        for indice_aresta in range(len(instancia["arestas"])):
            variaveis[(mercadoria, indice_aresta)] = coluna
            coluna += 1
    return variaveis


def construir_objetivo(instancia, variaveis):
    vetor = np.zeros(len(variaveis))
    for (mercadoria, indice_aresta), coluna in variaveis.items():
        if instancia["arestas"][indice_aresta][1] == instancia["destinos"][mercadoria]:
            vetor[coluna] = instancia["premios"][mercadoria]
    return vetor


def construir_modelo_fluxo(instancia):
    variaveis = construir_variavel(instancia)
    lucro = construir_objetivo(instancia, variaveis) 
    arestas = instancia["arestas"]
    K = instancia["num_mercadorias"]
    n = len(variaveis)
    m1 = len(arestas)

    # linhas de conservação: um par (mercadoria, vértice) para cada v != s_i, t_i
    linhas_cons = {}
    for mercadoria in range(1, K + 1):
        for v in range(1, instancia["num_vertices"] + 1):
            if v in (instancia["fontes"][mercadoria], instancia["destinos"][mercadoria]):
                continue
            linhas_cons[(mercadoria, v)] = m1 + len(linhas_cons)
    m = m1 + len(linhas_cons)

    A = np.zeros((m, n + m1))
    b = np.zeros(m)
    for e, (origem, destino, cap) in enumerate(arestas):
        b[e] = -cap                              #  -sum f >= -cap
        A[e, n + e] = -1.0                       #  excesso:  -sum f - s = -cap
    for (mercadoria, e), col in variaveis.items():
        origem, destino, _ = arestas[e]
        A[e, col] = -1.0
        if (mercadoria, origem) in linhas_cons:
            A[linhas_cons[(mercadoria, origem)], col] += 1.0     # sai de 'origem'
        if (mercadoria, destino) in linhas_cons:
            A[linhas_cons[(mercadoria, destino)], col] -= 1.0    # entra em 'destino'

    custo = np.r_[-lucro, np.zeros(m1)]
    return {
        "variaveis": variaveis, "lucro": lucro, "custo": custo, "A": A, "b": b,
        "n": n, "m1": m1, "m": m, "linhas_cons": linhas_cons,
    }

class SubproblemaRestrito:
    def __init__(self, A, b, tol=1e-9, refatorar_a_cada=100):
        self.A, self.b = A, b
        self.m, self.N = A.shape
        self.tol = tol
        self.refatorar_a_cada = refatorar_a_cada
        self.sinais = np.where(b >= 0, 1.0, -1.0)
        self.base = np.arange(self.N, self.N + self.m)        # artificiais
        self.Binv = np.diag(self.sinais)
        self.xB = np.abs(b).astype(float)
        self.pivos = 0
        self.pivos_desde_refat = 0
        # colunas esparsas (no máximo 3 não nulos por coluna): índices e valores,
        # preenchidos com a linha fictícia m (valor 0) para vetorizar o produto A'pi
        k = int(max(1, (A != 0).sum(axis=0).max()))
        self.idx = np.full((self.N, k), self.m, dtype=np.int64)
        self.val = np.zeros((self.N, k))
        for j in range(self.N):
            nz = np.flatnonzero(A[:, j])
            self.idx[j, :len(nz)] = nz
            self.val[j, :len(nz)] = A[nz, j]

    def produto_transposto(self, pi):
        """A' pi usando a estrutura esparsa."""
        pi_ext = np.append(pi, 0.0)
        return (self.val * pi_ext[self.idx]).sum(axis=1)

    def direcao(self, j):
        """B^-1 A_j usando só os não nulos da coluna."""
        nz = self.idx[j] < self.m
        return self.Binv[:, self.idx[j][nz]] @ self.val[j][nz]

    def _matriz_base(self):
        B = np.empty((self.m, self.m))
        for pos, j in enumerate(self.base):
            if j < self.N:
                B[:, pos] = self.A[:, j]
            else:
                B[:, pos] = 0.0
                B[j - self.N, pos] = self.sinais[j - self.N]
        return B

    def refatorar(self):
        self.Binv = np.linalg.inv(self._matriz_base())
        self.xB = self.Binv @ self.b
        self.xB[np.abs(self.xB) < 1e-11] = 0.0
        self.pivos_desde_refat = 0

    def multiplicadores(self):
        """pi = c_B' B^-1 (custo 1 nas artificiais da base, 0 nas demais)."""
        art = self.base >= self.N
        return self.Binv[art].sum(axis=0)

    def valor(self):
        return float(self.xB[self.base >= self.N].sum())

    def _escolher_saida(self, dB):
        """Teste da razão; empates (tolerância) -> menor índice de variável (Bland)."""
        linhas = np.flatnonzero(dB > 1e-9)
        if len(linhas) == 0:
            return None, None
        razoes = np.maximum(self.xB[linhas], 0.0) / dB[linhas]
        rmin = razoes.min()
        empatadas = linhas[razoes <= rmin + 1e-9 * (1.0 + abs(rmin))]
        sai = int(empatadas[np.argmin(self.base[empatadas])])
        return sai, float(max(self.xB[sai], 0.0) / dB[sai])

    def resolver(self, permitidas):
        visitadas = set()
        while True:
            pi = self.multiplicadores()
            alfa = self.produto_transposto(pi)        # pi'A_j ; custo reduzido no RSP = -alfa
            elegiveis = permitidas & (alfa > self.tol)
            elegiveis[self.base[self.base < self.N]] = False
            candidatos = np.flatnonzero(elegiveis)
            if len(candidatos) == 0:
                break                                # RSP ótimo
            entra = int(candidatos[0])               # Bland

            dB = self.direcao(entra)
            sai, theta = self._escolher_saida(dB)
            if sai is None:
                raise RuntimeError("RSP ilimitado (impossível: o RSP é limitado por zero)")

            if theta > 1e-12:
                visitadas.clear()                    # w diminuiu: novo trecho
            else:
                chave = (tuple(sorted(self.base.tolist())), entra, sai)
                if chave in visitadas:
                    raise RuntimeError("Ciclo detectado no RSP (Bland)")
                visitadas.add(chave)

            self.xB -= theta * dB
            self.xB[sai] = theta
            self.xB[np.abs(self.xB) < 1e-11] = 0.0
            linha_pivo = self.Binv[sai] / dB[sai]
            afetadas = np.flatnonzero(np.abs(dB) > 1e-14)      # só as linhas com dB_i != 0
            afetadas = afetadas[afetadas != sai]
            self.Binv[afetadas] -= np.outer(dB[afetadas], linha_pivo)
            self.Binv[sai] = linha_pivo
            self.base[sai] = entra

            self.pivos += 1
            self.pivos_desde_refat += 1
            if self.pivos_desde_refat >= self.refatorar_a_cada:
                self.refatorar()
        self.refatorar()                          
        # após refatorar, a base pode ainda admitir melhora numérica; repete se necessário
        pi = self.multiplicadores()
        alfa = self.produto_transposto(pi)
        el = permitidas & (alfa > 10 * self.tol)
        el[self.base[self.base < self.N]] = False
        if el.any():
            return self.resolver(permitidas)
        return self.valor(), pi

    def solucao(self):
        x = np.zeros(self.N)
        mask = self.base < self.N
        x[self.base[mask]] = self.xB[mask]
        return x


def resolver_primal_dual(modelo, max_iteracoes=10000, tol=1e-9, tol_w=1e-7):
    A, b, c = modelo["A"], modelo["b"], modelo["custo"]
    n, m1, m = modelo["n"], modelo["m1"], modelo["m"]
    arestas_cap = np.array([-bi for bi in b[:m1]])
    N = n + m1
    y = np.zeros(m)
    lucro = modelo["lucro"]
    for (mercadoria, e), col in modelo["variaveis"].items():
        y[e] = max(y[e], lucro[col])
    d = c - A.T @ y
    assert d.min() >= -tol, "y inicial não é dual-factível"

    rsp = SubproblemaRestrito(A, b, tol=tol)
    historico = []
    status = "max_iteracoes_atingido"
    fluxo_ok = False

    for iteracao in range(1, max_iteracoes + 1):
        d = c - A.T @ y
        J = d <= tol                                  
        J[rsp.base[rsp.base < N]] = True            
        w, pi = rsp.resolver(J)                     
        obj_dual = float(arestas_cap @ y[:m1])
        historico.append((iteracao, int(J.sum()), w, obj_dual, rsp.pivos))
        if w <= tol_w:
            status = "otimo"
            fluxo_ok = True
            break

        alfa = A.T @ pi
        possiveis = (~J) & (alfa > tol)
        if not possiveis.any():
            status = "infactivel"
            break
        theta = float(np.min(d[possiveis] / alfa[possiveis]))
        y = y + theta * pi

    x = rsp.solucao()[:n] if fluxo_ok else np.zeros(n)
    d = c - A.T @ y
    return {
        "x": x,
        "y": y,
        "y_capacidade": y[:m1].copy(),
        "y_conservacao": y[m1:].copy(),
        "objetivo_primal": float(modelo["lucro"] @ x),
        "objetivo_dual": float(arestas_cap @ y[:m1]),
        "custo_reduzido_min": float(d.min()),
        "status": status,
        "iteracoes": len(historico),
        "pivos_rsp": rsp.pivos,
        "historico": historico,
    }

def validar_solucao(instancia, modelo, solucao, tol=1e-5):
    if solucao["status"] != "otimo":
        return False
    x = solucao["x"]
    variaveis = modelo["variaveis"]
    arestas = instancia["arestas"]
    K, V = instancia["num_mercadorias"], instancia["num_vertices"]
    origem = np.array([a[0] for a in arestas])
    destino = np.array([a[1] for a in arestas])
    cap = np.array([a[2] for a in arestas])

    F = np.zeros((K, len(arestas)))
    for (mercadoria, e), col in variaveis.items():
        F[mercadoria - 1, e] = x[col]

    if F.min() < -tol:
        return False
    if np.any(F.sum(axis=0) > cap + tol):
        return False

    obj = 0.0
    for k in range(1, K + 1):
        s, t = instancia["fontes"][k], instancia["destinos"][k]
        saldo = np.zeros(V + 1)
        np.add.at(saldo, origem, F[k - 1])
        np.add.at(saldo, destino, -F[k - 1])
        for v in range(1, V + 1):
            if v not in (s, t) and abs(saldo[v]) > tol:
                return False
        obj += instancia["premios"][k] * F[k - 1][destino == t].sum()
    if abs(obj - solucao["objetivo_primal"]) > 1e-3:
        return False

    # certificado de otimalidade: factibilidade dual e gap nulo
    if solucao["custo_reduzido_min"] < -tol:
        return False
    if abs(solucao["objetivo_primal"] - solucao["objetivo_dual"]) > 1e-4 * (1 + abs(obj)):
        return False
    return True

def exibir_resultado(instancia, modelo, solucao, valido):
    x = solucao["x"]
    
    print("\n-- Variáveis primais não-nulas (fluxo por mercadoria e aresta) --")
    algum = False
    for (mercadoria, e), col in sorted(modelo["variaveis"].items()):
        o, d, _ = instancia["arestas"][e]
        if abs(x[col]) > 1e-6:
            algum = True
            print(f"  f[mercadoria={mercadoria}, aresta={e + 1}, {o} -> {d}] = {x[col]:.4f}")
    if not algum:
        print("  (nenhuma variável não-nula)")
        
    print("\n-- Variáveis duais das restrições de capacidade (y_uv >= 0) --")
    algum = False
    for e, (o, d, cap) in enumerate(instancia["arestas"]):
        yv = solucao["y_capacidade"][e]
        if abs(yv) > 1e-6:
            algum = True
            print(f"  y[aresta={e + 1}, {o} -> {d}] (capacidade {cap:g}) = {yv:.4f}")
    if not algum:
        print("  (nenhuma variável dual de capacidade ativa)")
 
    print("\n-- Variáveis duais das restrições de conservação (livres) --")
    algum = False
    for (mercadoria, v), linha in sorted(modelo["linhas_cons"].items()):
        yv = solucao["y"][linha]
        if abs(yv) > 1e-6:
            algum = True
            print(f"  pi[mercadoria={mercadoria}, vértice={v}] = {yv:.4f}")
    if not algum:
        print("  (nenhuma variável dual de conservação não-nula)")

    print(f"\nValor ótimo (objetivo primal): {solucao['objetivo_primal']:.4f}")
    print(f"Valor ótimo (objetivo dual):   {solucao['objetivo_dual']:.4f}")
    print("\nSolução válida" if valido else "\nSolução inválida")

def main():
    instancia = ler_instancia("mc_instance1.max")
    modelo = construir_modelo_fluxo(instancia)

    solucao = resolver_primal_dual(modelo)
    valido = validar_solucao(instancia, modelo, solucao)
    exibir_resultado(instancia, modelo, solucao, valido)


if __name__ == "__main__":
    main()