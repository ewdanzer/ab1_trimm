import os
import csv
from Bio import SeqIO

def corte_rapido_tramontina(pasta_entrada, pasta_saida, nota_de_corte):
    # Verifica se a pasta existe, se não, ela cria
    caminho_saida = os.path.join(pasta_saida, "trimm")
    os.makedirs(caminho_saida, exist_ok=True)

    # Verifica e cria uma lista com todos os arquivos com extensão .ab1 presentes na pasta
    arquivo_ab1 = [file for file in os.listdir(pasta_entrada) if file.endswith(".ab1")]

    # Cria uma lista para armazenar os nomes dos arquivos .ab1
    nomes = []

    for file in arquivo_ab1:
        # Fixa o caminho da pasta
        arquivo_entrada = os.path.join(pasta_entrada, file)

        # Fixa o nome do arquivo de saída
        arquivo_saida = os.path.join(caminho_saida, f"{os.path.splitext(file)[0]}_trim.fasta")

        # Faz o seqIO abrir o arquivo
        sequencia_ab1 = SeqIO.parse(arquivo_entrada, "abi")

        # armazena as sequências cortadas
        sequencias_cortadas = []

        for sequence in sequencia_ab1:
            # Fixa a qualidade de cada base na sequência baseado no numero descrito mais abaixo
            base_qualities = sequence.letter_annotations["phred_quality"]

            # Encontra o início da seuencia de base com qualidade acima da nota de corte
            inicio_sequencia = next(i for i, qual in enumerate(base_qualities) if qual >= nota_de_corte)

            # Encontra a ultima base final com qualidade acima da nota de corte
            end_index = next(i for i, qual in enumerate(base_qualities[::-1]) if qual >= nota_de_corte)
            end_index = len(sequence) - end_index

            # copia a sequência com base no inicio e final verificado anteriormente
            trimmed_sequence = sequence[inicio_sequencia:end_index]

            # Adiciona a sequência cortada à lista
            sequencias_cortadas.append(trimmed_sequence)

        # Salva as sequências cortadas no arquivo de saída
        SeqIO.write(sequencias_cortadas, arquivo_saida, "fasta")

        # Adiciona o nome do arquivo .ab1 à lista
        nomes.append(file)
        
    # Aviso no terminal
    print("Pontas com qualidade ruim foram cortadas com sucesso!")

    # Cria a planilha CSV na pasta "trimm"
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

