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

