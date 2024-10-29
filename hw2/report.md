# Отчет по анализу филогенетического дерева гена RAG2

## Введение
Эволюционные отношения челюстноротых (челюстных позвоночных), которые включают хрящевых рыб (хрящевых рыб), лопастеперых рыб (целакантов и двоякодышащих рыб), четвероногих и лучепёрых рыб (лучеперых рыб), обсуждаются уже почти столетие. Целью данного анализа является построение и сравнение филогенетических деревьев для гена RAG2 у различных видов. Мы используем методы построения деревьев, включая выравнивание последовательностей, построение филогенетического дерева, и сравнение деревьев с использованием Robinson-Foulds distance.

## Данные
Исходный набор данных состоит из последовательностей гена RAG2 у различных видов. Данные включают последовательности с идентификаторами:

```bash
AF369089.1Trionyx sinensis RAG2 (RAG2) gene, partial cds
AF369088.1Typhlonectes natans RAG2 (RAG2) gene, partial cds
AF369087.1Latimeria menadoensis RAG2 (RAG2) gene, partial cds
AF369086.1Protopterus sp. IMCB-2001 RAG2 (RAG2) gene, partial cds
AF369085.1Pachytriton sp. IMCB-2001 RAG2 (RAG2) gene, partial cds
AF369084.1Torpedo californica RAG2 (RAG2) gene, partial cds
AF369083.1Triakis sp. IMCB-2001 RAG2 (RAG2) gene, partial cds
AF369082.1Chiloscyllium punctatum RAG2 (RAG2) gene, partial cds
AF369081.1Osteoglossum sp. IMCB-2001 RAG2 (RAG2) gene, partial cds
AF369080.1Amia calva RAG2 (RAG2) gene, partial cds
AF369079.1Amia calva RAG2 (RAG2) gene, partial cds
AF369078.1Lepisosteus osseus RAG2 (RAG2) gene, partial cds
AF369077.1Lepisosteus osseus RAG2 (RAG2) gene, partial cds
AF369076.1Polyodon spathula RAG2 (RAG2) gene, partial cds
AF369075.1Polyodon spathula RAG2 (RAG2) gene, partial cds
AF369074.1Acipenser sp. IMCB-2001 RAG2 (RAG2) gene, partial cds
AF369073.1Acipenser sp. IMCB-2001 RAG2 (RAG2) gene, partial cds
AF369072.1Polypterus sp. IMCB-2001 RAG2 (RAG2) gene, partial cds
```

## Подготовка среды
Для выполнения анализа мы создали виртуальную среду и установили необходимые инструменты:

```bash
conda create --name phylogen
conda activate phylogen
conda install bioconda::fasttree -y
```
## Сбор и объединение данных
Мы объединили все файлы с последовательностями RAG2 в один файл rag2.fasta для дальнейшего выравнивания.

```bash
cat data/rag2/*.fasta >> rag2.fasta
```

## Множественное выравнивание последовательностей
Мы использовали MAFFT для множественного выравнивания последовательностей:

```bash
mafft --auto rag2.fasta > rag2_align.fasta
```

- `mafft --auto`: Запускает MAFFT с автоматическим выбором стратегии выравнивания в зависимости от размера и сложности входных данных.
- `rag2.fasta`: Это входной файл, содержащий наши последовательности.
- `> rag2_align.fasta`: Результат выравнивания будет сохранен в файл `rag2_align.fasta`.

## Построение филогенетического дерева
 После выравнивания мы использовали FastTree для построения филогенетического дерева в формате Newick:
  
```bash
FastTree -gtr -boot 1000 -quote -nt rag2_align.fasta > rag2_tree.newick
```

## Визуализация iTol:
![Филогенетическое дерево](https://raw.githubusercontent.com/kagayaku29/phylogenetics_2024/04124eaa0792633d3fc1ff472cb9e7b065b99c36/hw2/iTol1.svg)

## Форматирование данных для MrBayes и выбор аутгруппы
Для байесовского анализа дерева использовали MrBayes, преобразовав выровненные данные в формат Nexus:

```bash
seqmagick convert --output-format nexus --alphabet dna rag2_align.fasta rag2_align.nex
```

Выбрала Torpedo californica и Triakis sp. в качестве аутгруппы, исходя из аналогии со статьей.

## Байесовский анализ с MrBayes
Для запуска MrBayes использовали следующие параметры:

```bash
execute formatted_rag2_align.nex
outgroup AF369084.1_Torpedo_californica AF369083.1_Triakis_sp
mcmcp ngen=1000000 nruns=2 nchains=2 samplefreq=200 burninfrac=0.2
mcmc
sump
sumt
```

## Визуализация iTol:
![Филогенетическое дерево](https://raw.githubusercontent.com/kagayaku29/phylogenetics_2024/1f312d906b7c76e018b301e977cdbf681ed45a18/hw2/iTol2.svg)

## Сравнение деревьев с использованием Robinson-Foulds Distance
Для сравнения деревьев использовали DendroPy, который позволяет рассчитать расстояние Robinson-Foulds:

```python
# Установка DendroPy
!pip install dendropy

# Импорт библиотек
import dendropy
from dendropy.calculate import treecompare

# Загрузка деревьев из файлов
tree1_path = "tree1.tre"
tree2_path = "tree2.newick"
tree1 = dendropy.Tree.get(path=tree1_path, schema="newick")
tree2 = dendropy.Tree.get(path=tree2_path, schema="newick")

# Вычисление RF distance
rf_distance = treecompare.symmetric_difference(tree1, tree2)
print("Robinson-Foulds Distance:", rf_distance)
```
Значение RF distance получилось 62, что свидетельствует о значительных различиях в структуре двух деревьев.

## Заключение
В результате анализа мы построили филогенетические деревья для гена RAG2 и рассчитали степень различий между деревьями с помощью Robinson-Foulds distance. Значительное значение RF distance указывает на то, что деревья имеют различную топологию, что может отражать различия в эволюционной истории у исследуемых видов.
