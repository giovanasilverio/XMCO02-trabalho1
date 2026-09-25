def ler_instancia(caminho):
    num_vertices = 0
    num_arestas = 0
    num_mercadorias = 0
    id_mercadoria = 0 

    arestas = []
    fontes = {}
    destinos = {}

    with open(caminho, "r", encoding="utf-8") as arquivo:
        for linha in arquivo:
            linha = linha.strip()

            if linha == "":
                continue

            if linha[0] == "c":
                continue

            if linha[0] == "p":
                partes = linha.split()
                num_vertices = int(partes[2])
                num_arestas = int(partes[3])
                num_mercadorias = int(partes[4])
                
            if linha[0] == "n":
                partes = linha.split()
                vertice = int(partes[1])
                tipo = str(partes[2])
                id_mercadoria = int(partes[3])
                
                if tipo == "s":
                    fontes[id_mercadoria] = vertice
                    
                if tipo == "t":
                    destinos[id_mercadoria] = vertice
            
            if linha[0] == "a":
                partes = linha.split()
                origem = int(partes[1])
                destino = int(partes[2])
                capacidade = int(partes[3])
                arestas.append((origem, destino, capacidade))
                

    instancia = {
        "num_vertices": num_vertices,
        "num_arestas": num_arestas,
        "num_mercadorias": num_mercadorias,
        "arestas": arestas,
        "fontes": fontes,
        "destinos": destinos
    }

    return instancia


def main():
    instancia = ler_instancia("mc_instance1.max")
    print(instancia)


if __name__ == "__main__":
    main()