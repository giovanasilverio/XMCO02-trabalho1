"""
Fluxo máximo com múltiplas mercadorias — método PRIMAL-DUAL (Dantzig–Ford–Fulkerson)
conforme a Aula 14 (Métodos Exatos / XMCO02).  Somente Python + NumPy.

Uso:  python teste.py caminho/da/instancia.max

Modelo (forma "min c'x, Ax >= b, x >= 0", como nos slides)
---------------------------------------------------------
  variáveis  : f[i,e] >= 0  (fluxo da mercadoria i na aresta e; arestas paralelas
               são colunas distintas, identificadas pelo ÍNDICE da linha 'a')
  min        : - sum_i p_i * sum_{e=(v,t_i)} f[i,e]
  capacidade : - sum_i f[i,e] >= - cap_e            (dual u_e >= 0)
  conservação:   saída - entrada = 0   em v != s_i, t_i   (dual pi[i,v] livre)

Para ter um subproblema em igualdade, cada restrição de capacidade recebe uma
variável de excesso s_e >= 0 (coluna -e_e, custo 0).  O vetor dual completo é
y = (u, pi) e o custo reduzido de uma coluna j é d_j = c_j - y'A_j.

Algoritmo primal-dual
---------------------
 Passo 0  y dual-factível (NÃO é y = 0: o custo dos arcos que entram em t_i é
          negativo).  Usa-se u_e = maior prêmio entre as mercadorias que terminam
          no vértice final de e, e pi = 0.
 Passo 1  J = { j : d_j = 0 }.
 Passo 2  RSP:  min sum w   s.a.  A_J x_J + (artificiais) = b,  x_J >= 0.
 Passo 3  w* = 0  =>  x do RSP é primal-factível e, por folgas complementares
          (x_j > 0 só para j em J), ótimo.  Parar.
 Passo 4  w* > 0  =>  pi* = multiplicadores duais do RSP;
          theta = min { d_j / (pi*'A_j) : j fora de J, pi*'A_j > 0 };  y += theta*pi*.
 Passo 5  Novo J; voltar ao Passo 1.

Sobre o RSP (causa dos ciclos do código anterior)
-------------------------------------------------
O RSP é resolvido por simplex sobre as colunas de J, com a MESMA base de uma
iteração externa para a seguinte (partida a quente).  Isso é essencial: as
colunas básicas têm pi*'A_j = 0, logo seu custo reduzido não muda em
y + theta*pi* e elas continuam em J; a base anterior continua válida e as colunas
que acabam de entrar em J têm custo reduzido negativo no RSP.  O objetivo dual
cresce exatamente theta*w* > 0 por iteração externa e o RSP prossegue do ponto onde
parou, o que garante o progresso.  Reiniciar o RSP do zero a cada passo (como antes)
permite escolher, entre os muitos multiplicadores ótimos degenerados, um pi*
diferente a cada vez; nesse caso o passo theta pode encolher indefinidamente
(0,5; 0,25; 0,125...) com w* > 0 constante.
"""

import sys
import time

import numpy as np


# --------------------------------------------------------------------------- #
# Leitura da instância
# --------------------------------------------------------------------------- #
def ler_instancia(caminho):
    num_vertices = 0
    num_arestas = 0
    num_mercadorias = 0

    arestas = []     # (origem, destino, capacidade); o ÍNDICE na lista identifica a aresta
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


# --------------------------------------------------------------------------- #
# Construção do modelo
# --------------------------------------------------------------------------- #
def construir_variavel(instancia):
    """(mercadoria, índice da aresta) -> coluna.  Arestas paralelas ficam distintas."""
    variaveis = {}
    coluna = 0
    for mercadoria in range(1, instancia["num_mercadorias"] + 1):
        for indice_aresta in range(len(instancia["arestas"])):
            variaveis[(mercadoria, indice_aresta)] = coluna
            coluna += 1
    return variaveis


def construir_objetivo(instancia, variaveis):
    """Prêmio (lucro) de cada variável de fluxo: p_i se a aresta entra em t_i."""
    vetor = np.zeros(len(variaveis))
    for (mercadoria, indice_aresta), coluna in variaveis.items():
        if instancia["arestas"][indice_aresta][1] == instancia["destinos"][mercadoria]:
            vetor[coluna] = instancia["premios"][mercadoria]
    return vetor


def construir_modelo_fluxo(instancia):
    """
    Devolve um dicionário com a matriz do PL na forma de igualdade
        A_ext x_ext = b,   x_ext = (fluxos, excessos de capacidade) >= 0,
    linhas = [capacidades (m1) ; conservações (m2)].
    """
    variaveis = construir_variavel(instancia)
    lucro = construir_objetivo(instancia, variaveis)       # a PL minimiza -lucro
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


# --------------------------------------------------------------------------- #
# Subproblema primal restrito (simplex revisado com partida a quente)
# --------------------------------------------------------------------------- #
class SubproblemaRestrito:
    """
    RSP:  min soma(artificiais)  s.a.  A x + S w = b,  x_j >= 0 (j em J), x_j = 0 (j fora de J),
    S = diag(sinais) deixa as artificiais iniciais com valor |b_i| >= 0.
    A base é mantida entre chamadas; só o conjunto J (colunas permitidas) muda.
    Regra de Bland com empates tratados por tolerância (ver resolver).
    """

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
        """
        Simplex sobre as colunas 'permitidas' (máscara booleana de tamanho N).

        Regra de Bland: entra a coluna elegível de menor índice (custo reduzido
        negativo no RSP) e, no teste da razão, empates (com tolerância) saem pelo
        menor índice de variável.  O RSP é massivamente degenerado (muitos b_i = 0);
        o tratamento do empate com tolerância é o que mantém a garantia de Bland
        em ponto flutuante.  Base repetida só pode ocorrer dentro de um trecho de
        pivôs degenerados (fora dele w decresce), então a detecção de ciclo guarda
        apenas as bases do trecho atual e sinaliza erro se uma se repetir.
        """
        visitadas = set()
        while True:
            pi = self.multiplicadores()
            alfa = self.A.T @ pi                      # pi'A_j ; custo reduzido no RSP = -alfa
            elegiveis = permitidas & (alfa > self.tol)
            elegiveis[self.base[self.base < self.N]] = False
            candidatos = np.flatnonzero(elegiveis)
            if len(candidatos) == 0:
                break                                # RSP ótimo
            entra = int(candidatos[0])               # Bland

            dB = self.Binv @ self.A[:, entra]
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
            dB_outras = dB.copy()
            dB_outras[sai] = 0.0
            self.Binv -= np.outer(dB_outras, linha_pivo)
            self.Binv[sai] = linha_pivo
            self.base[sai] = entra

            self.pivos += 1
            self.pivos_desde_refat += 1
            if self.pivos_desde_refat >= self.refatorar_a_cada:
                self.refatorar()
        self.refatorar()                              # confirma valor e multiplicadores sem deriva
        # após refatorar, a base pode ainda admitir melhora numérica; repete se necessário
        pi = self.multiplicadores()
        alfa = self.A.T @ pi
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


# --------------------------------------------------------------------------- #
# Laço primal-dual
# --------------------------------------------------------------------------- #
def resolver_primal_dual(modelo, max_iteracoes=100000, tol=1e-9, tol_w=1e-7):
    A, b, c = modelo["A"], modelo["b"], modelo["custo"]
    n, m1, m = modelo["n"], modelo["m1"], modelo["m"]
    arestas_cap = np.array([-bi for bi in b[:m1]])
    N = n + m1

    # Passo 0: y dual-factível.  u_e = maior prêmio entre as mercadorias cujo destino
    # é a cabeça de e; pi = 0.  Verifica-se d >= 0 (em particular o excesso exige u >= 0).
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
        J = d <= tol                                  # Passo 1 (custo reduzido zero)
        J[rsp.base[rsp.base < N]] = True              # básicas permanecem em J (pi'A_j = 0)
        w, pi = rsp.resolver(J)                       # Passos 2-3
        obj_dual = float(arestas_cap @ y[:m1])
        historico.append((iteracao, int(J.sum()), w, obj_dual, rsp.pivos))
        if w <= tol_w:
            status = "otimo"
            fluxo_ok = True
            break

        # Passo 4: maior passo que preserva d >= 0
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


# --------------------------------------------------------------------------- #
# Validação (independente da matriz do modelo: usa só os dados da instância)
# --------------------------------------------------------------------------- #
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


# --------------------------------------------------------------------------- #
# Saída
# --------------------------------------------------------------------------- #
def exibir_resultado(instancia, modelo, solucao, valido, tempo=None):
    x = solucao["x"]
    print("=" * 70)
    print("SOLUÇÃO — MÉTODO PRIMAL-DUAL DO SIMPLEX")
    print("=" * 70)
    extra = f", tempo: {tempo:.2f}s" if tempo is not None else ""
    print(f"Status: {solucao['status']}   (iterações primal-duais: {solucao['iteracoes']}, "
          f"pivôs do RSP: {solucao['pivos_rsp']}{extra})")

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
    instancia = ler_instancia("mc_instance5.max")
    modelo = construir_modelo_fluxo(instancia)

    print(f"  vértices={instancia['num_vertices']}  arestas={instancia['num_arestas']}  "
          f"mercadorias={instancia['num_mercadorias']}")
    print(f"  fontes={instancia['fontes']}  destinos={instancia['destinos']}  "
          f"premios={instancia['premios']}")
    print(f"  quantidade de variáveis do modelo: {len(modelo['variaveis'])}")

    inicio = time.perf_counter()
    solucao = resolver_primal_dual(modelo)
    tempo = time.perf_counter() - inicio
    valido = validar_solucao(instancia, modelo, solucao)
    exibir_resultado(instancia, modelo, solucao, valido, tempo)


if __name__ == "__main__":
    main()