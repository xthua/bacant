from Bio import SeqIO

f = open("transposon_ann")
ann = {}
for lines in f:
    lines = lines.strip()
    if lines:
        line = lines.split("\t")
        accs = line[-1].split(",")
        type = line[0]
        if accs != "NA":
            for acc in accs:
                acc = acc.split(".")[0]
                ann[acc] = type
f.close()
#print(ann)
for record in SeqIO.parse("Transposon.fasta","fasta"):
    id = record.id
    if id in ann:
        id = id+'|'+ann[id]
    seq = record.seq
    with open("out.fasta","a") as w:
        w.write('>'+id+'\n'+str(seq)+'\n')
    #break
