Author:     Hong Wenjie at 20191126

Email:      xiaotinghua@zju.edu.cn

institute:  Key laboratory of Microbiol technology and Bioinformatics of Zhejiang Province

This program is designed for annotation of antimicrobal resistance(AMR), insertion sequence(IS), transposon(Tn) and integron(In) in bacteria.

### Install:
Bacant is a python3.X script, running on linux. You can download from github by `git clone https://github.com/xthua/bacant.git`.
* First:  
  You should install BLAST and add it in environment variable, you can download from `https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/`.
  BLAST version is 2.7.1 in bacant.
* Second:  
  Read requirements.txt , check python module dependecies or just run: pip install -r requirements.txt.
  Be care of module version, incompatible version will report bug.
* Last:  
  In bacant/Integron_finder/software, there are three exceute file require Execute Permission.
  That means `cd bacant/Integron_finder/software && chmod +x cmsearch && chmod +x hmmsearch && chmod +x prodigal`


### Run:
Bacant can accept FASTA and GENBANK format file. Attention on GENBANK format file, it should follow standard format.
There are three input parameter, "-n" means FASTA, "-g" means GENBANK, "-D" means input dir contains FASTA or GENBANK.
* Simply, you can just run:
```
python main.py -n FASTA -o outdir
python main.py -g GENBANK -o outdir
python main.py -D input_dir -o outdir
```
* For more parameter, you can run:
```
python main.py -h
```
* Here are some import parameter:

parameter  | description
---- | -----
--nucleotide(-n) | FASTA file
--genbank(-g) | GENBANK file
--indir(-D) | input dirname
--resultdir(-o) | output dirname
--coverages(-c) | filtering coverage, default is "60|60|60|60", four numbers represents AMR,IS,In,Tn in turn
--identities(-i) | filtering identity, default is "90|90|90|90", four numbers represents AMR,IS,In,Tn in turn

      
### Output:

filename  | description
---- | -----
annotation.gb | GENBANK format annotation
AMR.xls | resistance annotation
integron.filter.xls | most like integron
integron.detail.xls | integron_finder result,detail descripton of integron structure
ISfinder.filter.xls | insertion sequence
transposon.filter.xls | transposon element
