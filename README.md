# Trabalho 1 — Métodos Exatos

Trabalho da disciplina **Métodos Exatos (XMAC02)**, do curso de Ciência da Computação da **Universidade Federal de Itajubá (UNIFEI)**.

- **Professor:** Sávio S. Dias
- **Período:** 2º semestre de 2026
- **Linguagem utilizada:** Python

## Objetivo

Implementar o **algoritmo primal-dual do Simplex** para resolver dois problemas de otimização em grafos:

1. Corte mínimo s-t.
2. Fluxo máximo com múltiplas mercadorias.

Os resultados devem ser comparados com os obtidos utilizando **PuLP (COIN-OR)**.

## Organização dos arquivos

| Arquivo | Descrição |
| --- | --- |
| `min.py` | Implementação referente ao corte mínimo s-t. |
| `pulp_min.py` | Resolução do corte mínimo utilizando PuLP. |
| `instance1.min` a `instance5.min` | Instâncias de teste do corte mínimo. |
| `max.py` | Implementação referente ao fluxo máximo com múltiplas mercadorias. |
| `pulp_max.py` | Resolução do fluxo com múltiplas mercadorias utilizando PuLP. |
| `mc_instance1.max` a `mc_instance5.max` | Instâncias de teste do fluxo com múltiplas mercadorias. |
| `README.md` | Documentação do trabalho. |

## 1. Corte mínimo s-t

Dado um grafo $G=(V,E)$, com capacidades $c_{uv}\geq 0$ para cada aresta e dois vértices distintos $s$ e $t$, encontrar um conjunto de arestas cuja remoção desconecte $s$ de $t$, minimizando a soma de suas capacidades.

### Modelo matemático

As variáveis $x_{uv}$ estão associadas às arestas, e as variáveis $d_u$ aos vértices.

$$
\min \sum_{(u,v)\in E} c_{uv}x_{uv}
$$

Sujeito a:

$$
x_{uv}\geq d_u-d_v
\qquad \forall (u,v)\in E
$$

$$
x_{uv}\geq d_v-d_u
\qquad \forall (u,v)\in E
$$

$$
d_s=0,\qquad d_t=1
$$

$$
x_{uv}\geq 0
\qquad \forall (u,v)\in E
$$

$$
d_u\geq 0
\qquad \forall u\in V
$$

## 2. Fluxo máximo com múltiplas mercadorias

Dado um grafo $G=(V,E)$, com capacidades $c_{uv}\geq 0$, e $k$ pares fonte-sumidouro $(s_i,t_i)$, cada um com prêmio $p_i$, determinar os fluxos enviados simultaneamente, respeitando as capacidades compartilhadas das arestas.

### Modelo matemático

A variável $f_{i,uv}$ representa o fluxo da mercadoria $i$ na aresta $(u,v)$.

$$
\max \sum_{i=1}^{k} p_i
\sum_{(v,t_i)\in E} f_{i,vt_i}
$$

Sujeito a:

**Capacidade compartilhada:**

$$
\sum_{i=1}^{k} f_{i,uv}\leq c_{uv}
\qquad \forall (u,v)\in E
$$

**Conservação do fluxo nos vértices intermediários:**

$$
\sum_{(u,v)\in E} f_{i,uv}
-
\sum_{(v,u)\in E} f_{i,vu}
=0
\qquad
\forall u\in V\setminus\{s_i,t_i\},
\quad \forall i=1,\ldots,k
$$

**Não negatividade:**

$$
f_{i,uv}\geq 0
\qquad \forall (u,v)\in E,
\quad \forall i=1,\ldots,k
$$

O objetivo maximiza a soma dos fluxos que chegam ao destino de cada mercadoria, ponderados pelos respectivos prêmios.

## Comparação com PuLP

Para cada instância, devem ser comparados os valores da função objetivo e a viabilidade das soluções obtidas.

Quando houver múltiplas soluções ótimas, os valores individuais das variáveis podem diferir entre as implementações, mesmo que o valor ótimo da função objetivo seja igual.

Diferenças numéricas devem ser analisadas considerando a precisão utilizada.

## Dependências

Para a comparação com o solver, é necessário instalar PuLP:

```bash
python -m pip install pulp
```

As demais dependências e a forma de selecionar as instâncias devem ser verificadas nos respectivos scripts.


A apresentação será realizada por grupo, com execução de instâncias e perguntas sobre a implementação.

Conforme o enunciado, é expressamente proibido o uso de LLMs para geração de código neste trabalho. A identificação desse uso implica não originalidade e anulação da nota.
