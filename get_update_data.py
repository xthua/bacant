import os
import sys
from Bio import SeqIO
import re


info_file = os.path.abspath(sys.argv[1])
new_file = os.path.abspath(sys.argv[2])
old_file = os.path.abspath(sys.argv[3])
outfile = os.path.abspath(sys.argv[4])
tempfile = outfile+".temp"

d = {}
f=open(info_file)
lines=f.readlines()
f.close()
for line in lines[1:]:
	data = line.strip('\n').split('\t')
	gene_name = data[0]
	if not gene_name:
		gene_name = data[1]
	product_name = data[3]
	class_ = data[7].strip()
	type_ = data[5]

	ref_nucl_acc = data[10]
	nucl_acc = data[13]
	d[ref_nucl_acc] = "ncbi~~~%s~~~%s~~~%s %s"%(gene_name,nucl_acc,class_,product_name)
	d[nucl_acc] = "ncbi~~~%s~~~%s~~~%s %s"%(gene_name,nucl_acc,class_,product_name)

seqs = {}
for rec in SeqIO.parse(old_file,'fasta'):
	acc = re.split("~+",rec.id)[2]
	seqs[acc] = rec

add_items = 0
w=open(tempfile,'w+')
for rec in SeqIO.parse(new_file,'fasta'):
	acc = rec.id.split('|')[2]
	if acc not in seqs:
		add_items +=1
		w.write(">"+d[acc]+'\n'+str(rec.seq)+'\n')
w.close()

os.system("cat %s %s >%s"%(old_file,tempfile,outfile))
if os.path.exists(tempfile):
	os.remove(tempfile)
outdir = os.path.dirname(outfile)
os.system('makeblastdb -in %s -out %s/Res -dbtype nucl'%(outfile,outdir))
print("add %s items"%add_items)