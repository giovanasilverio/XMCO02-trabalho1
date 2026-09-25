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

            # Fonte ou destino de uma mercadoria:
            # n <vertice> <s/t> <mercadoria>
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
            # Sentido origem -> destino
            variaveis[
                (mercadoria, origem, destino)
            ] = coluna
            coluna += 1

            # Sentido destino -> origem
            variaveis[
                (mercadoria, destino, origem)
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


def main():
    instancia = ler_instancia("mc_instance1.max")

    variaveis = construir_variavel(instancia)

    objetivo = construir_objetivo(
        instancia,
        variaveis
    )

    print("Instância:")
    print(instancia)

    print("\nQuantidade de variáveis:")
    print(len(variaveis))

    print("\nVariáveis e respectivas colunas:")
    for variavel, coluna in variaveis.items():
        print(variavel, "-> coluna", coluna)

    print("\nCoeficientes não nulos da função objetivo:")
    for variavel, coluna in variaveis.items():
        if objetivo[coluna] != 0:
            print(
                variavel,
                "->",
                objetivo[coluna]
            )


if __name__ == "__main__":
    main()