Author:     Hong Wenjie

Email:      xiaotinghua@zju.edu.cn

institute:  Key laboratory of Microbiol technology and Bioinformatics of Zhejiang Province

This program is designed for annotation of antimicrobal resistance(AMR), insertion sequence(IS), transposon(Tn) and integron(In) in bacteria.

### Install:
Bacant is a python3.X script, running on linux. 
You should install BLAST and add it in environment variable, you can download from `https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/`. BLAST version is 2.7.1 in bacant.

* One:
  You can download from github by `git clone https://github.com/xthua/bacant.git`. Then execute `python setup.py install`.
* Two:
  You can install BacAnt from [PyPI](https://pypi.org/project/BacAnt) by `pip install BacAnt`.
* Three:
  When installation, maybe packages confict exists. Use conda to fix it. `conda create -n bacant python=3.6` create a virtual environment, `pip install BacAnt`.

### Run:
BacAnt can accept FASTA and GENBANK format file. Attention on GENBANK format file, it should follow standard format.
There are three input parameter, "-n" means FASTA, "-g" means GENBANK, "-D" means input dir contains FASTA or GENBANK.
* Simply, you can just run:
```
bacant -n FASTA -o outdir
bacant -g GENBANK -o outdir
bacant -D input_dir -o outdir
```
* For more parameter, you can run:
```
bacant -h
```
* Here are some import parameter:

parameter  | description
---- | -----
--nucleotide(-n) | FASTA file
--genbank(-g) | GENBANK file
--indir(-D) | input dirname
--resultdir(-o) | output dirname
--coverages(-c) | filtering coverage, default is "60\|60\|60\|60", four numbers represents AMR,IS,In,Tn in turn
--identities(-i) | filtering identity, default is "90\|90\|90\|90", four numbers represents AMR,IS,In,Tn in turn

      
### Output:

filename  | description
---- | -----
*.gb | GENBANK format annotation
AMR.tsv | filtered resistance annotation
AMR.possible.tsv | all possible resistance annotation
replicon.tsv | replicon annotation
integron.filter.tsv | most like integron
integron.detail.tsv | integron_finder result,detail descripton of integron structure
ISfinder.filter.tsv | filtered insertion sequences
ISfinder.tsv | all possible insertion sequences
transposon.filter.tsv | transposon element after overlap screen
transposon.possible.tsv | all possible transposon element
