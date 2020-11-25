from Bio import SeqIO

f = open("integron_ann")
ann = {}
for lines in f:
    lines = lines.strip()
    if lines:
        line = lines.split("\t")
        if len(line)>1:
            ann[line[0]] = line[1]
f.close()

for record in SeqIO.parse("Integron.fasta","fasta"):
    id = record.id
    if id in ann:
        id = id+'|'+ann[id]
    seq = record.seq
    with open("out.fasta","a") as w:
        w.write('>'+id+'\n'+str(seq)+'\n')
    #break
