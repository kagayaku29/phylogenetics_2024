# Where Do We Come From? What Are We? Where Are We Going?

Task1: https://docs.google.com/document/d/1uV-HjQk-K1WCGMnkjM_lkI4p4QMIK8H3qZJAHX2oQVI/edit

Data: https://figshare.com/ndownloader/files/30768763

Neanderthal samples: https://figshare.com/ndownloader/files/30768766

Denisovan samples: https://figshare.com/ndownloader/files/30768775


## Step1

Создание окружения, установка туллов.
```
  conda create -n phylo_hw
  conda activate phylo_hw
  conda install -c bioconda mafft
  conda install -c bioconda iqtree2
  conda install -c bioconda fasttree
```

Запуск из текущей папки, объединенние всех скачанных файлов с расширением .fasta
`cat *.fasta > all_sequences.fasta`

Используем MAFFT для множественного выравнивания последовательносте и записываем результат в файл aligned_sequences.fasta.
`mafft --auto all_sequences.fasta > aligned_sequences.fasta `

Запуск fasttree для получения файла .newick
`FastTree -gtr -nt -boot 1000 aligned_sequences.fasta > tree_with_bootstrap_1.newick`

Для визуализации можно имкользовать код на Pytone
```
from Bio import Phylo
import matplotlib.pyplot as plt

tree = Phylo.read("tree_with_bootstrap_1.newick", "newick")

Phylo.draw(tree)
plt.show()
```
![image](https://github.com/user-attachments/assets/d1e32baf-e767-4d80-8d0f-5dbef6541495)

Так же можно воспользоваться ITOl (https://itol.embl.de/ ) ,у сервиса есть графическая веб версия, бесплатный функционал широкий и он делает красивую картинку. 
РЕЗУЛЬТАТЫ ПОЛУЧЕННЫЕ С ПОМОШЬЮ ITOl - tree.pdf, tree2.pdf, tree3.pdf

## Step2. Расчет возраста митохондриальной Евы

Для мтДНК человека скорость мутаций может составлять примерно 1 мутация на 5000 лет. Возраст = Количество мутаций между исследуемыми последовательностями/2
```
from Bio import AlignIO

alignment = AlignIO.read("aligned_sequences.fasta", "fasta")

def get_sequence_by_id(alignment, seq_id):
    for record in alignment:
        if record.id == seq_id:
            return record.seq
    return None

african_id = "FJ713601.1"
non_african_id = "KY934476.1"

african_seq = get_sequence_by_id(alignment, african_id)
non_african_seq = get_sequence_by_id(alignment, non_african_id)

if african_seq is None or non_african_seq is None:
    print("Одна или обе последовательности не найдены.")
    exit()

def count_mutations(seq1, seq2):
    mutations = 0
    for a, b in zip(seq1, seq2):
        if a != b and a != '-' and b != '-':  
            mutations += 1
    return mutations

mutations = count_mutations(african_seq, non_african_seq)

years_per_mutation = 5000
age = (mutations * years_per_mutation) / 2 

print(f"Количество мутаций между {african_id} и {non_african_id}: {mutations}")
print(f"Возраст митохондриальной Евы: {age} лет")
```
Таким образом: Количество мутаций между FJ713601.1 и KY934476.1: 80
Возраст митохондриальной Евы: 200000.0 лет

## Step3. Неандертальцы и Денисовцы
Этапы повторяются. Результат см. tree2.pdf

1. Возраст самого позднего предка неандертальца — современного человека:
```
from Bio import AlignIO

alignment = AlignIO.read("aligned_new_sequences.fasta", "fasta")

def get_sequence_by_id(alignment, seq_id):
    for record in alignment:
        if record.id == seq_id:
            return record.seq
    return None

human_id = "KY934476.1"
neanderthal_id = "KX198084.1"

human_seq = get_sequence_by_id(alignment, human_id)
neanderthal_seq = get_sequence_by_id(alignment, neanderthal_id)

if human_seq is None or neanderthal_seq is None:
    print("Одна или обе последовательности не найдены.")
    exit()

def count_mutations(seq1, seq2):
    mutations = 0
    for a, b in zip(seq1, seq2):
        if a != b and a != '-' and b != '-': 
            mutations += 1
    return mutations

mutations = count_mutations(human_seq, neanderthal_seq)

years_per_mutation = 5000
age = (mutations * years_per_mutation) / 2  

print(f"Количество мутаций между {human_id} и {neanderthal_id}: {mutations}")
print(f"Возраст общего предка: {age} лет")
```
Количество мутаций между KY934476.1 и KX198084.1: 212
Возраст общего предка: 530000.0 лет

2. Возраст расхождения Неандертальцев и Денисовцев
```
    from Bio import AlignIO

alignment = AlignIO.read("aligned_new_sequences.fasta", "fasta")
def get_sequence_by_id(alignment, seq_id):
    for record in alignment:
        if record.id == seq_id:
            return record.seq
    return None
neanderthal_id = "KX198084.1"
denisovan_ids = ["FN673705.1", "FR695060.1", "KT780370.1"]

neanderthal_seq = get_sequence_by_id(alignment, neanderthal_id)
if neanderthal_seq is None:
    print(f"Не найдена последовательность неандертальца с ID {neanderthal_id}")
    exit()
def count_mutations(seq1, seq2):
    mutations = 0
    for a, b in zip(seq1, seq2):
        if a != b and a != '-' and b != '-':  
            mutations += 1
    return mutations
total_mutations = 0
mutation_counts = []

for denisovan_id in denisovan_ids:
    denisovan_seq = get_sequence_by_id(alignment, denisovan_id)

    if denisovan_seq is None:
        print(f"Не найдена последовательность денисовца с ID {denisovan_id}")
        continue

    mutations = count_mutations(neanderthal_seq, denisovan_seq)
    mutation_counts.append(mutations)
    total_mutations += mutations

    print(f"Количество мутаций между {neanderthal_id} и {denisovan_id}: {mutations}")

if len(mutation_counts) > 0:
    average_mutations = total_mutations / len(mutation_counts)
else:
    print("Нет данных для вычисления среднего количества мутаций.")
    exit()
years_per_mutation = 5000
age = (average_mutations * years_per_mutation) / 2 

print(f"\nСреднее количество мутаций между неандертальцем и денисовцами: {average_mutations}")
print(f"Возраст общего предка денисовцев и неандертальцев (с учётом деления на 2): {age} лет")
```
Количество мутаций между KX198084.1 и FN673705.1: 391
Количество мутаций между KX198084.1 и FR695060.1: 391
Количество мутаций между KX198084.1 и KT780370.1: 369

Среднее количество мутаций между неандертальцем и денисовцами: 383.6666666666667
Возраст общего предка денисовцев и неандертальцев (с учётом деления на 2): 959166.6666666667 лет

Вывод: Около 959166 лет назад произошло расхождение неандертальцев и денисовцев
см tree3.pdf
