# Where Do We Come From? What Are We? Where Are We Going?

Task1: https://docs.google.com/document/d/1uV-HjQk-K1WCGMnkjM_lkI4p4QMIK8H3qZJAHX2oQVI/edit

Data: https://figshare.com/ndownloader/files/30768763

Neanderthal samples: https://figshare.com/ndownloader/files/30768766

Denisovan samples: https://figshare.com/ndownloader/files/30768775


## Step1

Создание окружения, установка туллов.
```conda create -n phylo_hw
  conda activate phylo_hw
  conda install -c bioconda mafft
  conda install -c bioconda iqtree2
  conda install -c bioconda fasttree```

Запуск из текущей папки, объединенние всех скачанных файлов с расширением .fasta
`cat *.fasta > all_sequences.fasta`

Используем MAFFT для множественного выравнивания последовательносте и записываем результат в файл aligned_sequences.fasta.
`mafft --auto all_sequences.fasta > aligned_sequences.fasta `

Запуск fasttree для получения файла .newick
`FastTree -gtr -nt -boot 1000 aligned_sequences.fasta > tree_with_bootstrap_1.newick`
