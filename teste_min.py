"""
Corte mínimo s-t — método PRIMAL-DUAL (Aula 14, Métodos Exatos / XMCO02).
Somente Python + NumPy.   Uso:  python corte_minimo_primal_dual.py instance3.min

Modelo (P1 dos slides; arestas tratadas como não direcionadas, como no exemplo):
    min  sum_e c_e x_e
    s.a. x_e - d_u + d_v >= 0        (dual y+_e >= 0)
         x_e + d_u - d_v >= 0        (dual y-_e >= 0)
         d_s = 0                     (dual lambda_s livre)
         d_t = 1                     (dual lambda_t livre)
         x, d >= 0
Dual:
    max  lambda_t
    s.a. y+_e + y-_e <= c_e                                  (coluna x_e)
         sum_{e=(v,.)} (-y+_e + y-_e) + sum_{e=(.,v)} (y+_e - y-_e) + [lambda_s | lambda_t] <= 0   (coluna d_v)
         y+, y- >= 0

Cada desigualdade ganha uma variável de excesso (coluna -e_i, custo 0) para que o
subproblema restrito trabalhe em igualdade.

Algoritmo (passos dos slides)
  Passo 0  y = 0 é dual-factível aqui, pois c >= 0 (o custo de d_v, dos excessos é 0).
  Passo 1  J = { j : c_j - y'A_j = 0 }.
  Passo 2  RSP:  min soma(w)  s.a.  A_J x_J + S w = b,  x_J >= 0,  x_j = 0 fora de J.
  Passo 3  w* = 0  =>  x do RSP é primal-factível e, por folgas complementares, ótimo.
  Passo 4  w* > 0  =>  pi* = multiplicadores duais do RSP;
           theta = min { (c_j - y'A_j)/(pi*'A_j) : j fora de J, pi*'A_j > 0 };  y += theta*pi*.
  Passo 5  novo J, volta ao Passo 1.
O RSP é resolvido por simplex revisado (Bland, empates por tolerância) com partida a
quente: a base é mantida entre iterações externas (as colunas básicas têm pi'A_j = 0 e
permanecem em J), o que garante o progresso do laço externo.
"""

import sys
import time

import numpy as np


# --------------------------------------------------------------------------- #
# Leitura
# --------------------------------------------------------------------------- #
def ler_instancia(caminho):
    n_declarado = m_declarado = None
    s = t = None
    arestas, capacidades = [], []
    with open(caminho, "r", encoding="utf-8") as f:
        for linha in f:
            linha = linha.strip()
            if not linha or linha[0] == "c":
                continue
            p = linha.split()
            if p[0] == "p":
                n_declarado, m_declarado = int(p[2]), int(p[3])
            elif p[0] == "n":
                if p[2] == "s":
                    s = int(p[1])
                elif p[2] == "t":
                    t = int(p[1])
            elif p[0] == "a":
                arestas.append((int(p[1]), int(p[2])))
                capacidades.append(float(p[3]))
    if n_declarado is None:
        raise ValueError("Instância inválida: linha 'p' não encontrada.")
    if s is None or t is None:
        raise ValueError("Instância inválida: fonte ou sumidouro ausente.")
    aviso = None
    if m_declarado != len(arestas):
        # usa as arestas realmente lidas (são elas que definem o problema)
        aviso = (f"AVISO: o cabeçalho declara {m_declarado} arestas, "
                 f"mas o arquivo contém {len(arestas)}; usando {len(arestas)}.")
    return {"V": n_declarado, "E": len(arestas), "s": s, "t": t,
            "arestas": arestas, "cap": np.array(capacidades), "aviso": aviso}


# --------------------------------------------------------------------------- #
# Modelo:  A x = b,  x = (d, x_arestas, excessos) >= 0
# --------------------------------------------------------------------------- #
def construir_modelo(inst):
    V, E, s, t = inst["V"], inst["E"], inst["s"], inst["t"]
    m = 2 * E + 2
    N = V + E + 2 * E
    A = np.zeros((m, N))
    b = np.zeros(m)
    for e, (u, v) in enumerate(inst["arestas"]):
        A[2 * e, V + e] += 1.0            #  x_e - d_u + d_v - s = 0
        A[2 * e, u - 1] -= 1.0
        A[2 * e, v - 1] += 1.0
        A[2 * e, V + E + 2 * e] = -1.0
        A[2 * e + 1, V + e] += 1.0        #  x_e + d_u - d_v - s = 0
        A[2 * e + 1, u - 1] += 1.0
        A[2 * e + 1, v - 1] -= 1.0
        A[2 * e + 1, V + E + 2 * e + 1] = -1.0
    A[2 * E, s - 1] = 1.0                 #  d_s = 0
    A[2 * E + 1, t - 1] = 1.0             #  d_t = 1
    b[2 * E + 1] = 1.0
    c = np.zeros(N)
    c[V:V + E] = inst["cap"]
    return {"A": A, "b": b, "c": c, "V": V, "E": E, "N": N, "m": m}


# --------------------------------------------------------------------------- #
# Subproblema primal restrito (simplex revisado, partida a quente)
# --------------------------------------------------------------------------- #
class SubproblemaRestrito:
    def __init__(self, A, b, tol=1e-9, refatorar_a_cada=100):
        self.A, self.b = A, b
        self.m, self.N = A.shape
        self.tol = tol
        self.refatorar_a_cada = refatorar_a_cada
        self.sinais = np.where(b >= 0, 1.0, -1.0)
        self.base = np.arange(self.N, self.N + self.m)       # artificiais
        self.Binv = np.diag(self.sinais)
        self.xB = np.abs(b).astype(float)
        self.pivos = 0
        self.pivos_desde_refat = 0
        k = int(max(1, (A != 0).sum(axis=0).max()))          # colunas esparsas
        self.idx = np.full((self.N, k), self.m, dtype=np.int64)
        self.val = np.zeros((self.N, k))
        for j in range(self.N):
            nz = np.flatnonzero(A[:, j])
            self.idx[j, :len(nz)] = nz
            self.val[j, :len(nz)] = A[nz, j]

    def produto_transposto(self, pi):
        return (self.val * np.append(pi, 0.0)[self.idx]).sum(axis=1)

    def direcao(self, j):
        nz = self.idx[j] < self.m
        return self.Binv[:, self.idx[j][nz]] @ self.val[j][nz]

    def refatorar(self):
        B = np.zeros((self.m, self.m))
        for pos, j in enumerate(self.base):
            if j < self.N:
                B[:, pos] = self.A[:, j]
            else:
                B[j - self.N, pos] = self.sinais[j - self.N]
        self.Binv = np.linalg.inv(B)
        self.xB = self.Binv @ self.b
        self.xB[np.abs(self.xB) < 1e-11] = 0.0
        self.pivos_desde_refat = 0

    def multiplicadores(self):
        """pi = c_B' B^-1  (custo 1 nas artificiais básicas, 0 nas demais)."""
        return self.Binv[self.base >= self.N].sum(axis=0)

    def valor(self):
        return float(self.xB[self.base >= self.N].sum())

    def _escolher_saida(self, dB):
        linhas = np.flatnonzero(dB > 1e-9)
        if len(linhas) == 0:
            return None, None
        razoes = np.maximum(self.xB[linhas], 0.0) / dB[linhas]
        rmin = razoes.min()
        emp = linhas[razoes <= rmin + 1e-9 * (1.0 + abs(rmin))]
        sai = int(emp[np.argmin(self.base[emp])])             # Bland
        return sai, float(max(self.xB[sai], 0.0) / dB[sai])

    def resolver(self, permitidas):
        visitadas = set()                                     # bases do trecho degenerado atual
        while True:
            pi = self.multiplicadores()
            alfa = self.produto_transposto(pi)                # custo reduzido no RSP = -alfa
            el = permitidas & (alfa > self.tol)
            el[self.base[self.base < self.N]] = False
            cand = np.flatnonzero(el)
            if len(cand) == 0:
                break
            entra = int(cand[0])                              # Bland
            dB = self.direcao(entra)
            sai, theta = self._escolher_saida(dB)
            if sai is None:
                raise RuntimeError("RSP ilimitado (impossível: limitado por zero)")
            if theta > 1e-12:
                visitadas.clear()
            else:
                chave = (tuple(sorted(self.base.tolist())), entra, sai)
                if chave in visitadas:
                    raise RuntimeError("Ciclo detectado no RSP (Bland)")
                visitadas.add(chave)
            self.xB -= theta * dB
            self.xB[sai] = theta
            self.xB[np.abs(self.xB) < 1e-11] = 0.0
            linha_pivo = self.Binv[sai] / dB[sai]
            afet = np.flatnonzero(np.abs(dB) > 1e-14)
            afet = afet[afet != sai]
            self.Binv[afet] -= np.outer(dB[afet], linha_pivo)
            self.Binv[sai] = linha_pivo
            self.base[sai] = entra
            self.pivos += 1
            self.pivos_desde_refat += 1
            if self.pivos_desde_refat >= self.refatorar_a_cada:
                self.refatorar()
        # refatora (elimina deriva) só se houve pivôs suficientes ou se w está perto de 0,
        # quando o valor/multiplicadores serão usados para decidir a otimalidade
        if self.pivos_desde_refat > 0 and (self.pivos_desde_refat >= 25 or self.valor() <= 1e-5):
            self.refatorar()
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


# --------------------------------------------------------------------------- #
# Laço primal-dual
# --------------------------------------------------------------------------- #
def resolver_primal_dual(modelo, max_iteracoes=100000, tol=1e-9, tol_w=1e-7):
    A, b, c = modelo["A"], modelo["b"], modelo["c"]
    N = modelo["N"]

    y = np.zeros(modelo["m"])                 # Passo 0: y = 0 é dual-factível pois c >= 0
    assert (c - A.T @ y).min() >= -tol, "y inicial não é dual-factível"

    rsp = SubproblemaRestrito(A, b, tol=tol)
    historico = []
    status = "max_iteracoes_atingido"
    for it in range(1, max_iteracoes + 1):
        d = c - A.T @ y
        J = d <= tol                                          # Passo 1
        J[rsp.base[rsp.base < N]] = True
        w, pi = rsp.resolver(J)                               # Passos 2-3
        historico.append((it, int(J.sum()), w, float(b @ y), rsp.pivos))
        if w <= tol_w:
            status = "otimo"
            break
        alfa = A.T @ pi                                       # Passo 4
        possiveis = (~J) & (alfa > tol)
        if not possiveis.any():
            status = "infactivel"
            break
        theta = float(np.min(d[possiveis] / alfa[possiveis]))
        y = y + theta * pi

    xc = rsp.solucao() if status == "otimo" else np.zeros(N)
    return {
        "xcompleto": xc, "y": y, "status": status,
        "objetivo_primal": float(c @ xc),
        "objetivo_dual": float(b @ y),
        "custo_reduzido_min": float((c - A.T @ y).min()),
        "iteracoes": len(historico), "pivos_rsp": rsp.pivos, "historico": historico,
    }


# --------------------------------------------------------------------------- #
# Recuperação, validação e saída
# --------------------------------------------------------------------------- #
def recuperar_solucao(inst, modelo, sol):
    V, E = inst["V"], inst["E"]
    xc, y = sol["xcompleto"], sol["y"]
    d = xc[:V]
    x = xc[V:V + E]
    S = [v for v in range(1, V + 1) if d[v - 1] <= 0.5]
    T = [v for v in range(1, V + 1) if d[v - 1] > 0.5]
    Sset = set(S)
    valor_particao = sum(inst["cap"][e] for e, (u, v) in enumerate(inst["arestas"])
                         if (u in Sset) != (v in Sset))
    return {"d": d, "x": x, "yMais": y[0:2 * E:2], "yMenos": y[1:2 * E:2],
            "lambdaS": y[2 * E], "lambdaT": y[2 * E + 1],
            "S": S, "T": T, "valorParticao": float(valor_particao)}


def validar_solucao(inst, modelo, sol, rec, tol=1e-5):
    """Só devolve True/False.  Primal (independente da matriz), partição, dual e gap."""
    if sol["status"] != "otimo":
        return False
    d, x = rec["d"], rec["x"]
    if d.min() < -tol or x.min() < -tol:
        return False
    if abs(d[inst["s"] - 1]) > tol or abs(d[inst["t"] - 1] - 1) > tol:
        return False
    for e, (u, v) in enumerate(inst["arestas"]):
        if x[e] < abs(d[u - 1] - d[v - 1]) - tol:
            return False
    if abs(float(inst["cap"] @ x) - sol["objetivo_primal"]) > 1e-6:
        return False
    if inst["s"] not in rec["S"] or inst["t"] not in rec["T"]:
        return False
    if abs(rec["valorParticao"] - sol["objetivo_primal"]) > 1e-6:       # corte realmente induzido
        return False
    if sol["custo_reduzido_min"] < -tol:                                 # dual factível
        return False
    if abs(sol["objetivo_primal"] - sol["objetivo_dual"]) > 1e-6 * (1 + abs(sol["objetivo_primal"])):
        return False
    return True


def imprimir_resultados(inst, sol, rec, valido, tempo):
    print("=" * 70)
    print("CORTE MÍNIMO s-t — MÉTODO PRIMAL-DUAL DO SIMPLEX")
    print("=" * 70)
    print(f"Vértices: {inst['V']}  Arestas: {inst['E']}  Fonte: {inst['s']}  Sumidouro: {inst['t']}")
    print(f"Status: {sol['status']}   (iterações primal-duais: {sol['iteracoes']}, "
          f"pivôs do RSP: {sol['pivos_rsp']}, tempo: {tempo:.2f}s)")

    print("\n-- Variáveis primais não-nulas --")
    for v, val in enumerate(rec["d"], start=1):
        if abs(val) > 1e-6:
            print(f"  d[{v}] = {val:.6g}")
    for e, val in enumerate(rec["x"]):
        if abs(val) > 1e-6:
            u, v = inst["arestas"][e]
            print(f"  x[aresta={e + 1}, {u}-{v}] = {val:.6g}  (capacidade {inst['cap'][e]:g})")

    print("\n-- Variáveis duais não-nulas --")
    for e, (u, v) in enumerate(inst["arestas"]):
        if abs(rec["yMais"][e]) > 1e-6:
            print(f"  y+[aresta={e + 1}, {u}-{v}] = {rec['yMais'][e]:.6g}")
        if abs(rec["yMenos"][e]) > 1e-6:
            print(f"  y-[aresta={e + 1}, {u}-{v}] = {rec['yMenos'][e]:.6g}")
    print(f"  lambda_s = {rec['lambdaS']:.6g}")
    print(f"  lambda_t = {rec['lambdaT']:.6g}")

    print(f"\nS = {rec['S']}")
    print(f"T = {rec['T']}")
    print(f"\nValor ótimo (objetivo primal): {sol['objetivo_primal']:.6f}")
    print(f"Valor ótimo (objetivo dual):   {sol['objetivo_dual']:.6f}")
    print(f"Valor do corte pela partição:  {rec['valorParticao']:.6f}")
    print("\nSolução válida" if valido else "\nSolução inválida")


def main():
    inst = ler_instancia("instance5.min")
    if inst["aviso"]:
        print(inst["aviso"])
    modelo = construir_modelo(inst)
    t0 = time.perf_counter()
    sol = resolver_primal_dual(modelo)
    tempo = time.perf_counter() - t0
    rec = recuperar_solucao(inst, modelo, sol)
    valido = validar_solucao(inst, modelo, sol, rec)
    imprimir_resultados(inst, sol, rec, valido, tempo)


if __name__ == "__main__":
    main()