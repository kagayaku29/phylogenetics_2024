# Отчет по анализу филогенетического дерева гена RAG2

## Введение
Целью данного анализа является построение и сравнение филогенетических деревьев для гена RAG2 у различных видов. Мы используем методы построения деревьев, включая выравнивание последовательностей, построение филогенетического дерева, и сравнение деревьев с использованием Robinson-Foulds distance.

## Данные
Исходный набор данных состоит из последовательностей гена RAG2 у различных видов. Данные включают последовательности с идентификаторами:

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

