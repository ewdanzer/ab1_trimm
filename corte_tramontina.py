import os
import csv
from Bio import SeqIO

def corte_rapido_tramontina(pasta_entrada, pasta_saida, nota_de_corte):
    # Verifica se a pasta existe, se não, ela cria
    caminho_saida = os.path.join(pasta_saida, "trimm")
    os.makedirs(caminho_saida, exist_ok=True)

    # Obter a lista de arquivos com a extensão .ab1 na pasta de entrada
    arquivo_ab1 = [file for file in os.listdir(pasta_entrada) if file.endswith(".ab1")]

    # Criar uma lista para armazenar os nomes dos arquivos .ab1
    nomes = []

    for file in arquivo_ab1:
        # Obter o caminho completo para o arquivo de entrada
        arquivo_entrada = os.path.join(pasta_entrada, file)

        # Criar o nome do arquivo de saída
        arquivo_saida = os.path.join(caminho_saida, f"{os.path.splitext(file)[0]}_trim.fasta")

        # Abrir o arquivo de entrada no formato .ab1
        sequencia_ab1 = SeqIO.parse(arquivo_entrada, "abi")

        # Lista para armazenar as sequências cortadas
        sequencias_cortadas = []

        for sequence in sequencia_ab1:
            # Obter a qualidade de cada base na sequência
            base_qualities = sequence.letter_annotations["phred_quality"]

            # Encontrar o índice da base inicial com qualidade acima do limiar
            inicio_sequencia = next(i for i, qual in enumerate(base_qualities) if qual >= nota_de_corte)

            # Encontrar o índice da base final com qualidade acima do limiar
            end_index = next(i for i, qual in enumerate(base_qualities[::-1]) if qual >= nota_de_corte)
            end_index = len(sequence) - end_index

            # Cortar a sequência com base nos índices encontrados
            trimmed_sequence = sequence[inicio_sequencia:end_index]

            # Adicionar a sequência cortada à lista
            sequencias_cortadas.append(trimmed_sequence)

        # Salvar as sequências cortadas no arquivo de saída
        SeqIO.write(sequencias_cortadas, arquivo_saida, "fasta")

        # Adicionar o nome do arquivo .ab1 à lista
        nomes.append(file)

    print("Pontas com qualidade ruim foram cortadas com sucesso!")

    # Criar a planilha CSV na pasta "trimm"
    csv_file = os.path.join(caminho_saida, "nomes.csv")
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Nome do Arquivo .ab1"])
        writer.writerows([[name] for name in nomes])

    print("Planilha CSV foi gerada com sucesso!")


pasta_entrada = "."  # Pasta atual, mas você pode alterar para uma pasta qualquer
pasta_saida = "."  # Pasta atual, ou, pode alterar para a pasta que voce quiser
nota_de_corte = 20  # Nota mínima para que seja considerado um par de base bom

corte_rapido_tramontina(pasta_entrada, pasta_saida, nota_de_corte)

