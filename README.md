Author:     Hong Wenjie

Email:      xiaotinghua@zju.edu.cn

institute:  Key laboratory of Microbiol technology and Bioinformatics of Zhejiang Province

This program is designed for annotation of antimicrobal resistance(AMR), transposon(Tn) and integron(In) in bacteria.

### Install:
Bacant is a python3.X script, running on linux. 
You should install BLAST and add it in environment variable, you can download from `https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/`. BLAST version is 2.7.1 in bacant.

* One:
  You can download from github by `git clone https://github.com/xthua/bacant.git`. Then execute `python setup.py install`.
* Two:
  You can install BacAnt from [PyPI](https://pypi.org/project/BacAnt) by `pip install BacAnt`.
* Three:
  You can install BacAnt from [anaconda](https://anaconda.org/bacant/bacant).
  First,`conda create -n bacant` create a virtual environment.Then `conda install -c bacant -c conda-forge -c bioconda bacant ` install bacant.

### Run:
BacAnt can accept FASTA and GENBANK format file(single or multi sequences in one file). Attention on GENBANK format file, it should follow standard format.
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
--coverages(-c) | filtering coverage, default is "60\|\|60\|60", four numbers represents AMR,In,Tn in turn
--identities(-i) | filtering identity, default is "90\|\|90\|90", four numbers represents AMR,In,Tn in turn

### Databases:
We have updated database to v2.0(2021.05.11). You can download from [here](http://bacant.net/static/database/v2.0/bacant-db-v2.0.tar.gz).

User can define their custom databases, and when run bacant ,just add parameter -p(--path) for databases dirname.

Here are databases structure:

.
├── IntegronDB
│   ├── Integron.fasta    Integron reference sequences in FASTA format,sequence id must be description|accession,eg: In0|PAU49101
│   ├── Integron.nhr
│   ├── Integron.nin
│   └── Integron.nsq
├── ResDB
│   ├── Res.fasta         Resistance gene reference sequences in FASTA format,sequence id must be database name~~~gene~~~accession~~~description,eg: ncbi~~~1567214_ble~~~NG_047553.1~~~BLEOMYCIN BLMA family bleomycin binding protein
│   ├── Res.nhr
│   ├── Res.nin
│   └── Res.nsq
└── TransposonDB
    ├── Transposon.fasta  Transposon reference sequences in FASTA format,sequence id must be description|accession,eg: Tn2009|CP001937
    ├── Transposon.nhr
    ├── Transposon.nin
    └── Transposon.nsq
      
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
annotation.html | output visualization
