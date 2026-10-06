Trabalho 1 — Métodos Exatos

Trabalho da disciplina Métodos Exatos (XMAC02), do Bacharelado em Ciência da Computação da Universidade Federal de Itajubá (UNIFEI), referente ao 2º semestre de 2026, ministrada pelo professor Sávio S. Dias.

Objetivo

Implementar o algoritmo primal-dual do Simplex para resolver dois problemas de otimização em grafos: corte mínimo s-t e fluxo máximo com múltiplas mercadorias. Os resultados devem ser comparados com os obtidos utilizando PuLP (COIN-OR).

Organização dos arquivos

Arquivo

Finalidade

README.md

Documentação e especificações do trabalho.

min.py

Código referente ao problema de corte mínimo s-t.

pulp_min.py

Resolução do corte mínimo com PuLP para comparação.

instance1.min a instance5.min

Instâncias de teste do corte mínimo.

mc.py

Código referente ao fluxo máximo com múltiplas mercadorias.

resolver_pulp.py

Resolução do fluxo com múltiplas mercadorias com PuLP para comparação.

mc_instance1.max a mc_instance5.max

Instâncias de teste do fluxo com múltiplas mercadorias.

__pycache__/

Arquivos de cache gerados automaticamente pelo Python.

A associação dos scripts aos problemas segue os nomes apresentados no repositório. Os argumentos de execução e as dependências adicionais devem ser consultados nos respectivos códigos.

1. Corte mínimo s-t

Dado um grafo $G=(V,E)$, com capacidades $c_{uv}\geq 0$ e dois vértices distintos $s$ e $t$, encontrar um conjunto de arestas cuja remoção desconecte $s$ de $t$, minimizando a soma de suas capacidades.

Modelo matemático

As variáveis são $x_{uv}$, associadas às arestas, e $d_u$, associadas aos vértices.

$$
\min \sum_{(u,v)\in E} c_{uv}x_{uv}
$$

Sujeito a:

$$
x_{uv}\geq d_u-d_v \qquad \forall (u,v)\in E
$$

$$
x_{uv}\geq d_v-d_u \qquad \forall (u,v)\in E
$$

$$
d_s=0,\qquad d_t=1
$$

$$
x_{uv}\geq 0 \qquad \forall (u,v)\in E
$$

$$
d_u\geq 0 \qquad \forall u\in V
$$

As duas desigualdades por aresta impõem $x_{uv}\geq |d_u-d_v|$. O objetivo minimiza a soma das capacidades ponderadas por essas variáveis.

2. Fluxo máximo com múltiplas mercadorias

Dado um grafo $G=(V,E)$, com capacidades $c_{uv}\geq 0$, e $k$ pares fonte-sumidouro $(s_i,t_i)$, cada um com prêmio $p_i$, determinar os fluxos enviados simultaneamente, respeitando as capacidades compartilhadas das arestas.

Modelo matemático

A variável $f_{i,uv}\geq 0$ representa o fluxo da mercadoria $i$ na aresta $(u,v)$.

$$
\max \sum_{i=1}^{k} p_i\sum_{(v,t_i)\in E} f_{i,vt_i}
$$

Sujeito a:

Capacidade compartilhada:

$$
\sum_{i=1}^{k}f_{i,uv}\leq c_{uv}
\qquad \forall (u,v)\in E
$$

Conservação do fluxo nos vértices intermediários:

$$
\sum_{(u,v)\in E} f_{i,uv}
-\sum_{(v,u)\in E} f_{i,vu}=0
\qquad \forall u\in V\setminus{s_i,t_i},\quad \forall i=1,\ldots,k
$$

Não negatividade:

$$
f_{i,uv}\geq 0
\qquad \forall (u,v)\in E,\quad \forall i=1,\ldots,k
$$

O objetivo soma os fluxos que chegam ao destino de cada mercadoria, ponderados pelos respectivos prêmios. A conservação exige que, nos vértices intermediários, o fluxo que entra seja igual ao fluxo que sai para cada mercadoria.

Requisitos do trabalho

Implementar o algoritmo primal-dual do Simplex para os dois problemas.

Ler e resolver as instâncias fornecidas em arquivos de texto, no formato dicionário/modelo indicado no enunciado.

Apresentar os valores das variáveis primais e duais não nulas na solução ótima, além do valor da função objetivo.

Resolver os mesmos modelos utilizando PuLP e comparar os resultados.

Apresentar a implementação ao professor, que executará instâncias e fará perguntas sobre o código.

Atividade bônus

Implementar o algoritmo de Stoer-Wagner para o corte mínimo global e comparar seu resultado com a solução do problema de corte mínimo s-t.

O corte global não fixa os vértices $s$ e $t$. Por isso, seu valor pode diferir do valor do corte mínimo s-t, devendo essa diferença ser considerada na comparação.

Comparação dos resultados

Para cada instância, a comparação deve registrar:

Valor da função objetivo obtido pela implementação primal-dual.

Valor da função objetivo obtido com PuLP.

Variáveis primais e duais não nulas da solução ótima.

Eventuais diferenças numéricas, considerando a precisão utilizada.

Soluções ótimas podem ter valores diferentes para as variáveis quando houver múltiplos ótimos. Assim, a comparação deve considerar também a viabilidade das soluções e o valor da função objetivo.

Para instalar PuLP:

python -m pip install pulp

Outras dependências e a forma de seleção das instâncias devem ser verificadas nos scripts antes da execução.
