Author:     Hong Wenjie at 20191126

Email:      xiaotinghua@zju.edu.cn

institute:  Key laboratory of Microbiol technology and Bioinformatics of Zhejiang Province

This program is designed for annotation of resistance gene, insertion sequence, transposon and integron in bacteria.


Run:

    bacant is a python3.X script, running on linux.
    1.BLAST is required, it shoud be added in environment variable.
      Read requirements.txt , check your dependecies or just run: pip install -r requirements.txt.
      Otherwise,in bacant/Integron_finder/software, there are three exceute file require Execute Permission.
      That means  "cd bacant/Integron_finder/software && chmod +x cmsearch && chmod +x hmmsearch && chmod +x prodigal"
    2.You can input FASTA or GENBANK format file:
      python main.py -n FASTA -o outdir
      python main.py -g GENBANK -o outdir
      python main.py -D input_dir -o outdir
      
Output:

    annotation.gb           GENBANK format annotation
    AMR.xls                 resistance annotation
    integron.filter.xls     most like integron
    integron.detail.xls     integron_finder result,detail descripton of integron structure
    ISfinder.filter.xls     insertion sequence
    transposon.filter.xls   transposon element
