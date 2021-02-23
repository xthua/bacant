#coding:utf-8

import sys,os,re,argparse
import matplotlib
from matplotlib import pyplot as plt 
import numpy as np
from Bio import SeqIO
import shutil
from Bio.Alphabet import generic_dna
from Bio import Seq
from BCBio import GFF
import re
from reportlab.lib import colors
from reportlab.lib.units import cm
import subprocess
import time
import warnings
import codecs
from threading import Thread
from BacAnt.Integron_Finder.scripts.finder import run_integron_finder
warnings.filterwarnings('ignore')

#hongwj 20191126
#usage: python3.7 main.py -o result -n test.fasta -c "90|90|90|90" -i "60|60|60|60" -d "resDB|ISDatabase|IntegronDB|TransposonDB" or python main.py -o result -g genbank.gb or python main.py -o result -D indir

#arg parse
def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nucleotide","-n",help="nucleotide file")
    parser.add_argument("--resultdir","-o",help="resultdir")
    parser.add_argument("--databases","-d",help="multi databases : default=resDB|ISDatabase|IntegronDB|TransposonDB",default="resDB|ISDatabase|IntegronDB|TransposonDB")
    parser.add_argument("--genbank","-g",help="genbank file")
    parser.add_argument("--coverages","-c",help="filtering coverage",default="60|60|60|60")
    parser.add_argument("--identities","-i",help="filtering identity",default="90|90|90|90")
    parser.add_argument("--indir","-D",help="input dir")
    args = parser.parse_args()
    nucleotide = args.nucleotide
    genbank = args.genbank
    resultdir = args.resultdir.rstrip("/")
    databases = args.databases
    coverages = args.coverages
    identities = args.identities
    indir = args.indir
    return nucleotide,genbank,resultdir,databases,coverages,identities,indir

#conver gff to genbank
def convert(gff_file, fasta_file,resultdir,gb_file):
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(check_gff(fix_ncbi_id(gff_iter)), gb_file, "genbank")

def fix_ncbi_id(fasta_iter):
    """GenBank identifiers can only be 16 characters; try to shorten NCBI.
    """
    for rec in fasta_iter:
        if len(rec.name) > 16 and rec.name.find("|") > 0:
            new_id = [x for x in rec.name.split("|") if x][-1]
            print("Warning: shortening NCBI name %s to %s" % (rec.id, new_id))
            rec.id = new_id
            rec.name = new_id
        yield rec

def check_gff(gff_iterator):
    """Check GFF files before feeding to SeqIO to be sure they have sequences.
    """
    for rec in gff_iterator:
        if isinstance(rec.seq, Seq.UnknownSeq):
            print("Warning: FASTA sequence not found for '%s' in GFF file" % (
                    rec.id))
            rec.seq.alphabet = generic_dna
        yield flatten_features(rec)

def flatten_features(rec):
    """Make sub_features in an input rec flat for output.
    GenBank does not handle nested features, so we want to make
    everything top level.
    """
    out = []
    for f in rec.features:
        cur = [f]
        while len(cur) > 0:
            nextf = []
            for curf in cur:
                out.append(curf)
                if len(curf.sub_features) > 0:
                    nextf.extend(curf.sub_features)
            cur = nextf
    rec.features = out
    return rec

#check software
def check_dependencies():
    blastn_path = shutil.which("blastn")
    if not blastn_path:
        print("blastn not found")
        
    blastp_path = shutil.which("blastp")
    if not blastp_path:
        print("blastp not found")
        
    return blastn_path,blastp_path

#parse input file
def format_fasta(annot,seq,num):
    format_seq=""
    for i, char in enumerate(seq):
        format_seq += char
        if (i + 1) % num == 0:
            format_seq+="\n"
    return annot + format_seq + "\n"

def parse_genbank(genbank,resultdir):
    gb_seqs=SeqIO.parse(genbank,"gb")
    for gb_seq in gb_seqs:
        complete_seq=str(gb_seq.seq)
        if "accessions" in gb_seq.annotations and gb_seq.description != "":
            complete_annot=">"+gb_seq.annotations["accessions"][0] + " " + gb_seq.description + "\n"
        elif "accessions" not in gb_seq.annotations and gb_seq.description != "":
            complete_annot=">"+ gb_seq.description + "\n"
        else:
            complete_annot=">"+'sequence'+"\n"
        complete_fasta = format_fasta(complete_annot, complete_seq, 60)
        index = 0
        feature_list = []
        for key in gb_seq.features:
            feature_list.append(key.type)
        if 'CDS' in feature_list:
            for key in gb_seq.features:
                if key.type == "CDS":
                    if 'locus_tag' in key.qualifiers:
                        with open(resultdir+'/gb_location.txt','a') as w:
                            w.write(str(key.location).split(":")[0].split("[")[1].strip('<')+':'+str(key.location).split(":")[1].split("]")[0].strip('>')+'\t'+str(key.qualifiers['locus_tag'][0])+'\t0\n')
                    else:
                        index = index + 1
                        with open(resultdir+'/gb_location.txt','a') as w:
                            w.write(str(key.location).split(":")[0].split("[")[1].strip('<')+':'+str(key.location).split(":")[1].split("]")[0].strip('>')+'\t'+'locus%s'%index+'\t0\n')
        else:
            for key in gb_seq.features:
                if key.type == "source":
                    index = index + 1
                    with open(resultdir+'/gb_location.txt','a') as w:
                        w.write(str(key.location).split(":")[0].split("[")[1].strip('<')+':'+str(key.location).split(":")[1].split("]")[0].strip('>')+'\t'+'locus%s'%index+'\t0\n')
        complete_file = resultdir+"/nucleotide.fasta"
        complete_file_obj = open(complete_file, "a")
        complete_file_obj.write(complete_fasta)
        break

def parse_fasta(name,nucleotide,resultdir):
    outfile = resultdir+'/nucleotide_all.fasta'
    outfile2 = resultdir+'/nucleotide.fasta'
    location = resultdir+'/gb_location.txt'
    all_len = 0
    with open(outfile,'a') as w:
        w.write('>%s\n'%name)
    for record in SeqIO.parse(nucleotide,"fasta"):
        id = record.id
        seq = str(record.seq)
        length = len(seq)
        start = all_len + 1
        end = length + all_len
        with open(location,'a') as w:
            w.write(str(start)+':'+str(end)+'\t'+id+'\t'+str(all_len)+'\n')
        with open(outfile,'a') as w:
            w.write(seq)
        all_len = all_len + length
    os.system('cp "'+nucleotide+'" '+outfile2)

def mkdir(resultdir):
    if not os.path.exists(resultdir):
        os.mkdir(resultdir)
    else:
        os.system("rm -rf %s/*"%resultdir)

def findResistanceGene(inFile,outdir,blastnPath,database,default_ident,default_cov):
    resistanceGenePosList = []
    cmd = blastnPath + ' -task blastn -dust no -evalue 1E-5 -culling_limit 1 -query ' + inFile + ' -db ' + database + ' -outfmt \"6 sseqid pident length  qstart qend slen sstrand qseqid gaps\" -out ' + outdir + '/blast.out'
    subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    blastResultDict = {}
    f = open(outdir+'/blast.out')
    lines = f.readlines()
    if len(lines) == 0:
        pass
    else:
        for line in lines:
            data = line.strip().split('\t')
            product = re.split('~~~',data[0])[3]
            acc = re.split('~~~',data[0])[2]
            ident = float(data[1])
            cov = float(data[2])*100/int(data[5])
            leftPos = data[3]
            rightPos = data[4]
            name = re.split('~~~',data[0])[1]+" "+leftPos+" "+rightPos
            strand = data[6]
            covlen = '1-'+data[2]+'/'+data[5]
            data = str(ident)+'|'+str(cov)+'|'+leftPos+'|'+rightPos+'|'+name+'|'+strand+'|'+data[-2].replace('|','~')+'|'+data[-1]+'|'+product+'|'+acc+'|'+covlen
            if name not in blastResultDict:
                blastResultDict[name] = []
                blastResultDict[name].append(data)
            else:
                blastResultDict[name].append(data)
    f.close()
    with open(outdir+'/AMR.possible.xls','w') as w:
        w.write('#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE\n')
    for name,dataList in blastResultDict.items():
        for data in dataList:
            ident = float(data.split('|')[0])
            cov = float(data.split('|')[1])
            leftPos = data.split('|')[2]
            rightPos = data.split('|')[3]
            name = data.split('|')[4].split(" ")[0]
            strand = data.split('|')[5]
            query = data.split('|')[6].replace('~','|')
            gaps = data.split('|')[7]
            product1 = data.split('|')[8]
            acc1 = data.split('|')[9]
            covlen1 = data.split('|')[10]
            if ident >= default_ident and cov >= default_cov:
                bestLeftPos = leftPos
                bestRightPos = rightPos
                pos = bestLeftPos + '-' + bestRightPos
                bestStand = strand
                bestIdent = format(ident, '0.2f')
                bestCov = format(cov, '0.2f')
                with open(outdir+'/AMR.possible.xls','a') as w:
                    w.write('nucleotide.fasta'+'\t'+query+'\t'+bestLeftPos+'\t'+bestRightPos+'\t'+bestStand+'\t'+name+'\t'+covlen1+'\t'+'.'+'\t'+gaps+'\t'+str(bestCov)+'\t'+str(bestIdent)+'\t'+'resDB'+'\t'+acc1+'\t'+product1+'\n')
                
    with open(outdir+'/AMR.result','w') as w:
        w.write('#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE\n')
    for name,dataList in blastResultDict.items():
        maxIdent = default_ident
        maxCov = default_cov
        pos = ''
        for data in dataList:
            ident = float(data.split('|')[0])
            cov = float(data.split('|')[1])
            leftPos = data.split('|')[2]
            rightPos = data.split('|')[3]
            name = data.split('|')[4].split(" ")[0]
            strand = data.split('|')[5]
            query = data.split('|')[6]
            gaps = data.split('|')[7]
            product1 = data.split('|')[8]
            acc1 = data.split('|')[9]
            covlen1 = data.split('|')[10]
            if ident >= maxIdent and cov >= maxCov:
                bestLeftPos = leftPos
                bestRightPos = rightPos
                pos = bestLeftPos + '-' + bestRightPos
                bestStand = strand
                bestIdent = format(ident, '0.2f')
                bestCov = format(cov, '0.2f')
        
        if pos != '':
            left1 = int(pos.split('-')[0])
            right1 = int(pos.split('-')[1])
            flag = True
            if resistanceGenePosList != []:
                for po in resistanceGenePosList:
                    left2 = int(po.split('-')[0])
                    right2 = int(po.split('-')[1])
                    if left1 < right2 and right1 > left2:
                        flag = False
                        break
                    else:
                        continue
            if flag:
                with open(outdir+'/AMR.result','a') as w:
                    w.write('nucleotide.fasta'+'\t'+query+'\t'+bestLeftPos+'\t'+bestRightPos+'\t'+bestStand+'\t'+name+'\t'+covlen1+'\t'+'.'+'\t'+gaps+'\t'+str(bestCov)+'\t'+str(bestIdent)+'\t'+'resDB'+'\t'+acc1+'\t'+product1+'\n')
                resistanceGenePosList.append(pos)
    os.remove(outdir + '/blast.out')


def AMR(nucleotide,resultdir,database,blastnPath,cov0,ident0):
    findResistanceGene(nucleotide,resultdir,blastnPath,database,ident0,cov0)
    infile = resultdir + "/AMR.result"
    outfile = resultdir + "/AMR.xls"
    tempfile = resultdir + "/AMR_temp.xls"
    f = open(infile)
    i = 0
    for lines in f:
        i = i + 1
        if i == 1:
            with open(tempfile,"a") as w:
                w.write(lines)
            continue
        lines = lines.strip()
        line = lines.split("\t")
        new_line = line[0].split("/")[-1] + '\t' + "\t".join(line[1:])
        with open(tempfile,"a") as w:
            w.write(new_line+'\n')
    f.close()
    f = open(infile)
    datas = f.readlines()
    f.close()
    with open(outfile,'w') as w:
        w.write(datas[0])
    newDatas = sorted(datas[1:],key = lambda data:(int(data.split('\t')[2])))
    for newData in newDatas:
        tmp_list = newData.split('\t')
        newData2 = tmp_list[0]+'\t'+tmp_list[1].replace('~','|')+'\t'+'\t'.join(tmp_list[2:])
        with open(outfile,'a') as w:
            w.write(newData2)
    os.remove(tempfile)
    os.remove(infile)
    out_file = resultdir + "/AMR.gff"
    flag = "AMR"
    format_genbank(outfile,out_file,flag,nucleotide,resultdir)
    print("AMR done")

def ISfinder(nucleotide,resultdir,blastn_path,cov1,ident1,database):
    cmd = blastn_path+" -num_threads 8 -task blastn -query "+nucleotide+" -db "+database+" -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq slen\" -out "+resultdir+"/ISfinder.xls"
    subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    with open(resultdir+"/ISfinder.xls", 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write('QUERY\tSUBJECT|FAMILY|GROUP|ORGANISM TYPE\t%IDENTITY\tLENGTH\tMISMATCH\tGAPOPEN\tQUERY START\tQUERY END\tSUBJECT START\tSUBJECT END\tEVALUE\tBITSCORE\tQUERY SEQ\tSUBJECT SEQ\tSUBJECT LENGTH\n'+content)
    blast_result = resultdir+'/ISfinder.xls'
    filter_result_temp = resultdir+'/ISfinder.filter.temp.xls'
    filter_result = resultdir+'/ISfinder.filter.xls'
    outfile = resultdir + "/IS.gff"
    flag = "IS"
    filter_blast_result(blast_result,filter_result_temp,cov1,ident1,flag)
    f = open(filter_result_temp)
    datas = f.readlines()
    f.close()
    with open(filter_result,'w') as w:
        w.write(datas[0])
    newDatas = sorted(datas[1:],key = lambda data:(int(data.split('\t')[6]),-float(data.split('\t')[14].strip()),-float(data.split('\t')[2])))
    for newData in newDatas:
        with open(filter_result,'a') as w:
            w.write(newData)
    os.remove(filter_result_temp)
    format_genbank(filter_result,outfile,flag,nucleotide,resultdir)
    print("ISfinder done")

def parse_IS(resultdir):
    family_list = []
    num_dict = {}
    file = resultdir+"/ISfinder.filter.xls"
    i = 0
    if os.path.exists(file):
        f = open(file)
        for lines in f:
            lines = lines.strip()
            if lines:
                i = i + 1
                if i == 1:
                    continue
                lin = lines.split("\t")[1]
                line = lin.split("|")
                name = line[0]
                family = line[1]
                group = line[2]
                origin = line[3]
                family_list.append(family)
        f.close()
    else:
        print("There is no ISfinder result,Please check run")
    family_set = set(family_list)
    for family in family_set:
        num = family_list.count(family)
        num_dict[family] = num 
    return num_dict

def plot(num_dict,resultdir):
    labels = []
    sizes = []
    num_set = sorted(num_dict.items(),key=lambda num_dict:num_dict[1],reverse=True)
    for val in num_set:
        labels.append(val[0])
        sizes.append(val[1])
    plt.rcdefaults()
    if len(labels) <= 3:
        fig, ax = plt.subplots(figsize=(15,1))
    elif len(labels) < 10 and len(labels) > 3:
        fig, ax = plt.subplots(figsize=(15,2))
    else:
        fig, ax = plt.subplots(figsize=(15,8))
    if max(sizes) >5:
        x_pos = np.arange(0,max(sizes)+1,5)
    else:
        x_pos = np.arange(0,max(sizes)+1)
    y_pos = np.arange(len(labels))
    bar = ax.barh(y_pos, sizes, align='center')
    ax.set_xticks(x_pos)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    plt.savefig(resultdir+"/IS.png")
    print("plot done")

def async_func(f):
    def wrapper(*args,**kwargs):
        thr = Thread(target=f,args=args,kwargs=kwargs)
        thr.start()
    return wrapper


@async_func
def integron_finder(current_dir,nucleotide,resultdir,blastn_path,database):
    merged_integron_path = "%s/nucleotide.integrons"%resultdir
    try:
        run_integron_finder(current_dir,resultdir,nucleotide)
        integron_pos_dict = parse_integron_result(resultdir,merged_integron_path)
        if integron_pos_dict:
            find_most_like_In(nucleotide,integron_pos_dict,resultdir,blastn_path,database)
        with open(resultdir+"/integron.finished",'w') as w:
            w.write("completed")
        print("integron_finder done")
    except:
        with open(resultdir+"/integron.finished",'w') as w:
            w.write("abnormal exit")
        print("integron_finder error")
    
    

def parse_integron_result(resultdir,merged_integron_path):
    integron_dict = {}
    if merged_integron_path:
        integron_file = merged_integron_path
        f = open(integron_file)
        lines = f.readlines()
        for line in lines[1:]:
            line = line.strip()
            if line:
                data = line.split("\t")
                name = data[1]+"|"+data[0]
                if name not in integron_dict:
                    integron_dict[name] = []
                    integron_dict[name].append(int(data[3]))
                    integron_dict[name].append(int(data[4]))
                else:
                    integron_dict[name].append(int(data[3]))
                    integron_dict[name].append(int(data[4]))
        f.close()
        integron_pos_dict = {}
        for k,v in integron_dict.items():
            v.sort()
            start = v[0]
            end = v[-1]
            integron_pos_dict[k] = str(start)+"|"+str(end)
        #shutil.move(integron_file,resultdir)
        return integron_pos_dict
    else:
        return {}

def find_most_like_In(nucleotide,integron_pos_dict,resultdir,blastn_path,database):
    integron_fasta = resultdir+"/integron.fa"
    seq_dict = {}
    for rec in SeqIO.parse(nucleotide,'fasta'):
        id = str(rec.id)
        seq = str(rec.seq)
        seq_dict[id] = seq
    for name,pos in integron_pos_dict.items():
        rec_id = name.split("|")[0]
        start = int(pos.split("|")[0])
        end = int(pos.split("|")[1])
        integron_seq = seq_dict[rec_id][start-1:end]
        with open(integron_fasta,'a') as w:
            w.write(">"+name.replace("|","~")+"\n"+integron_seq+"\n")
    blast_result = resultdir+'/integron.blast.out'
    cmd = blastn_path+" -num_threads 8 -query "+integron_fasta+" -db "+database+" -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq slen\" -out "+blast_result
    subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    
    size = os.path.getsize(blast_result)
    if size != 0:
        with open(blast_result, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write('QUERY\tINTEGRON TYPE|ACCESSION\t%IDENTITY\tLENGTH\tMISMATCH\tGAPOPEN\tQUERY START\tQUERY END\tSUBJECT START\tSUBJECT END\tEVALUE\tBITSCORE\tQUERY SEQ\tSUBJECT SEQ\tSUBJECT LENGTH\n'+content)
        filter_result = resultdir+'/integron.filter.xls'
        filter_result_temp = resultdir+'/integron.filter.temp.xls'
        outfile = resultdir + "/integron.gff"
        flag = "integron"
        cov2 = 0
        ident2 = 0
        flag = "integron"
        filter_blast_result(blast_result,filter_result_temp,cov2,ident2,flag)
        f = open(filter_result_temp)
        datas = f.readlines()
        f.close()
        with open(filter_result,'w') as w:
            w.write(datas[0])
        newDatas = sorted(datas[1:],key = lambda data:(int(data.split('\t')[6]),-float(data.split('\t')[14].strip()),-float(data.split('\t')[2])))
        name_tmp  = []
        for newData in newDatas:
            data_name = newData.split("\t")[0]
            if data_name not in name_tmp:
                name_tmp.append(data_name)
                with open(filter_result,'a') as w:
                    w.write(newData)
        #os.remove(blast_result)
        os.remove(filter_result_temp)
        format_genbank(filter_result,outfile,"integron",nucleotide,resultdir)
        #print("integron done")
    else:
        os.remove(blast_result)
        #print("no integron hits")
    




def transposon(nucleotide,resultdir,blastn_path,cov3,ident3,database):
    cmd = blastn_path+" -num_threads 8 -task blastn -query "+nucleotide+" -db "+database+" -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq slen\" -out "+resultdir+"/transposon.xls"
    subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    size = os.path.getsize(resultdir+"/transposon.xls")
    blast_result = resultdir+'/transposon.xls'
    if size != 0:
        with open(resultdir+"/transposon.xls", 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write('QUERY\tTRANSPOSON TYPE|ACCESSION\t%IDENTITY\tLENGTH\tMISMATCH\tGAPOPEN\tQUERY START\tQUERY END\tSUBJECT START\tSUBJECT END\tEVALUE\tBITSCORE\tQUERY SEQ\tSUBJECT SEQ\tSUBJECT LENGTH\n'+content)
        filter_result = resultdir+'/transposon.filter.xls'
        filter_result_temp = resultdir+'/transposon.filter.temp.xls'
        outfile = resultdir + "/transposon.gff"
        flag = "transposon"
        filter_blast_result(blast_result,filter_result_temp,cov3,ident3,flag)
        f = open(filter_result_temp)
        datas = f.readlines()
        f.close()
        with open(filter_result,'w') as w:
            w.write(datas[0])
        newDatas = sorted(datas[1:],key = lambda data:(int(data.split('\t')[6]),-float(data.split('\t')[14].strip()),-float(data.split('\t')[2])))
        for newData in newDatas:
            with open(filter_result,'a') as w:
                w.write(newData)
        os.remove(blast_result)
        os.remove(filter_result_temp)
        format_genbank(filter_result,outfile,flag,nucleotide,resultdir)
        print("transposon done")
    else:
        os.remove(blast_result)
        print("no transposon hits")

def filter_blast_result(blast_result,filter_result,coverage,identity,flag):
    if os.path.exists(filter_result):
        os.remove(filter_result)
    f = open(blast_result)
    i = 0
    for lines in f:
        i = i + 1
        if i == 1:
            if flag != "AMR":
                with open(filter_result,"a") as w:
                    title = '\t'.join(lines.split("\t")[0:14])+'\t'+'%COVERAGE'
                    w.write(title+'\n')
            else:
                with open(filter_result,"a") as w:
                    title = '\t'.join(lines.split("\t")[0:14])
                    w.write(title+'\n')
            continue
        lines = lines.strip()
        if lines:
            line = lines.split("\t")
            ident = float(line[2])
            slen = float(line[14])
            match_length = abs(int(line[9])-int(line[8])) + 1
            cov = round(match_length*100/slen,2)
            if flag != "AMR":
                if ident >= identity and cov >= coverage and match_length>= 50:
                    newline = '\t'.join(line[0:2]) + '\t' + str(round(ident,2)) + '\t' + '\t'.join(line[3:14]) + '\t' + str(cov)
                    with open(filter_result,"a") as w:
                        w.write(newline+'\n')
            else:
                if ident >= identity and match_length>= 50:
                    newline = '\t'.join(line[0:2]) + '\t' + str(round(ident,2)) + '\t' + '\t'.join(line[3:14])
                    with open(filter_result,"a") as w:
                        w.write(newline+'\n')
    f.close()

def format_genbank(infile,outfile,flag,fasta_file,resultdir):
    name_dict = {}
    query_dict = {}
    f = open(resultdir+"/gb_location.txt")
    for lines in f:
        lines = lines.strip()
        if lines:
            line = lines.split("\t")
            location = line[0]
            locus_tag = line[1]
            all_len = line[2]
            query_dict[locus_tag] = all_len
            name_dict[location] = locus_tag
    f.close()
    
    file = resultdir+'/nucleotide_all.fasta'
    if os.path.exists(file):
        for record in SeqIO.parse(file,"fasta"):
            filename = record.id
    f = open(infile)
    i = 0
    for lines in f:
        i = i + 1
        if i == 1:
            continue
        lines = lines.strip()
        if lines:
            line = lines.split("\t")
            #query
            if flag == "AMR":
                query = line[1]
                file = resultdir+'/nucleotide_all.fasta'
                if os.path.exists(file):
                    filename = filename
                else:
                    filename = query
            elif flag == "integron":
                query = "".join(line[0].split("~")[0:-1])
                file = resultdir+'/nucleotide_all.fasta'
                if os.path.exists(file):
                    filename = filename
                else:
                    filename = query
            else:
                query = line[0]
                file = resultdir+'/nucleotide_all.fasta'
                if os.path.exists(file):
                    filename = filename
                else:
                    filename = query
            
            #note
            if flag == "AMR":
                if line[-3] != "":
                    gene = line[5]
                    CDS = line[5]
                    db_xref = line[-3]
                    note = line[-1]
                else:
                    gene = line[5]
                    CDS = line[5]
                    note = line[-1]
                    db_xref = ""
            elif flag == "IS":
                text = line[1].split("|")
                gene = text[0]
                CDS = text[0]
                note = "family:"+text[1]+" group:"+text[2]+" organism:"+text[3]
                db_xref = ""
            elif flag == "integron":
                text = line[1].split("|")
                note = ""
                if len(text) > 1:
                    gene = text[0]
                    CDS = text[0]
                    db_xref = text[1]
                    
                else:
                    gene = text[0]
                    CDS = text[0]
                    db_xref = ""
                    
            elif flag == "transposon":
                text = line[1].split("|")
                note = ""
                if len(text) > 1:
                    gene = text[0]
                    CDS = text[0]
                    db_xref = text[1]
                else:
                    gene = text[0]
                    CDS = text[0]
                    db_xref = ""
            #start end
            if flag == "AMR":
                file = resultdir+'/nucleotide_all.fasta'
                if os.path.exists(file):
                    all_len = int(query_dict[query])
                    start = str(int(line[2]) + all_len)
                    end = str(int(line[3]) + all_len)
                else:
                    start = line[2]
                    end = line[3]
                #locus
                name = start + ":" + end
                locus = get_locus_tag(name_dict,name)
            else:
                file = resultdir+'/nucleotide_all.fasta'
                if os.path.exists(file):
                    all_len = int(query_dict[query.split("~")[0]])
                    start = str(int(line[6]) + all_len)
                    end = str(int(line[7]) + all_len)
                else:
                    start = line[6]
                    end = line[7]
                #locus
                name = start + ":" + end
                locus = get_locus_tag(name_dict,name)
            #source
            if flag == "AMR":
                source = "resistance"
            elif flag == "IS":
                source = "IS"
            elif flag == "integron":
                source = "integron"
            elif flag == "transposon":
                source = "transposon"
            #strand
            if flag == "AMR":
                strand = line[4]
            else:
                if int(line[8]) < int(line[9]):
                    strand = "+"
                else:
                    strand = "-"
            #feature
            features = ["gene","CDS","mobile_element","transposon"]
            if flag == "AMR":
                if db_xref == "":
                    with open(outfile,"a") as w:
                        w.write(filename+'\t'+source+'\t'+features[0]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';gene='+gene+';Note='+note+'\n')
                        w.write(filename+'\t'+source+'\t'+features[1]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';CDS='+CDS+';Note='+note+'\n')
                else:
                    with open(outfile,"a") as w:
                        w.write(filename+'\t'+source+'\t'+features[0]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';gene='+gene+';db_xref='+db_xref+';Note='+note+'\n')
                        
                        w.write(filename+'\t'+source+'\t'+features[1]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';CDS='+CDS+';db_xref='+db_xref+';Note='+note+'\n')
                        
            if flag == "IS":

                if db_xref == "":
                    with open(outfile,"a") as w:
                        w.write(filename+'\t'+source+'\t'+features[2]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';gene='+CDS+';Note='+note+'\n')
                else:
                    with open(outfile,"a") as w:
                        w.write(filename+'\t'+source+'\t'+features[1]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';gene='+CDS+';db_xref='+db_xref+';Note='+note+'\n')

            if flag == "integron":
                if db_xref == "":
                    with open(outfile,"a") as w:
                        w.write(filename+'\t'+source+'\t'+features[2]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';gene='+CDS+'\n')
                else:
                    with open(outfile,"a") as w:
                        w.write(filename+'\t'+source+'\t'+features[2]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';gene='+CDS+';db_xref='+db_xref+'\n')
            
            if flag == "transposon":
                if db_xref == "":
                    with open(outfile,"a") as w:
                        w.write(filename+'\t'+source+'\t'+features[3]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';transposon='+CDS+';gene='+gene+'\n')
                else:
                    with open(outfile,"a") as w:
                        w.write(filename+'\t'+source+'\t'+features[3]+'\t'+start+'\t'+end+'\t.\t'+strand+'\t.\tlocus_tag='+locus+';transposon='+CDS+';db_xref='+db_xref+';gene='+gene+'\n')
            
    f.close()
    if flag == "AMR":
        gb_file = resultdir+'/AMR.gb'
    elif  flag == "IS":
        gb_file = resultdir+'/IS.gb'
    elif  flag == "integron":
        gb_file = resultdir+'/integron.gb'
    elif  flag == "transposon":
        gb_file = resultdir+'/transposon.gb'
    if i > 1:
        f = open(outfile)
        datas = f.readlines()
        f.close()
        with open(outfile+'3','w') as w:
            w.write(datas[0])
        newDatas = sorted(datas[1:],key = lambda data:(int(data.split('\t')[3].strip())))
        for newData in newDatas:
            with open(outfile+'3','a') as w:
                w.write(newData)
        file = resultdir+'/nucleotide_all.fasta'
        if os.path.exists(file):
            fasta_file = file
        else:
            fasta_file = fasta_file
        convert(outfile+'3', fasta_file,resultdir,gb_file)

def get_locus_tag(name_dict,name):
    qstart = int(name.split(":")[0])
    qend = int(name.split(":")[1])
    max_len = 0
    max_dict = {}
    for k,v in name_dict.items():
        start = int(k.split(":")[0])
        end = int(k.split(":")[1])
        if qstart >= start and qend <= end:
            length = qend - qstart + 1
        elif qstart < start and qend > end:
            length = end - start + 1
        elif qstart < start and qend > start and qend <= end:
            length = qend - start + 1
        elif qstart >= start and qstart < end and qend > end:
            length = end - qstart
        else:
            length = 0
        if length > max_len:
            max_dict = {}
            rstart = str(start)
            rend = str(end)
            max_len = length
            max_dict[max_len] = rstart+':'+rend
    if max_dict != {}:
        ret = name_dict[max_dict[max_len]]
    else:
        ret = ""
    return ret
     
def getFileFormat(infile):
    try:
        for rec in SeqIO.parse(infile,'fasta'):
            seq = rec.seq
            if seq:
                return 'fasta'
    except:
        pass
        
    try:
        for rec in SeqIO.parse(infile,'genbank'):
            seq = rec.seq
            if seq:
                return 'genbank'
    except:
        pass
    
    return "unknown"



def run():
    nucleotide,genbank,resultdir,databases,coverages,identities,indir = arg_parse()
    curent_dir = os.path.dirname(os.path.abspath(__file__))
    database = databases.split("|")
    resDB = database[0]
    ISDatabase = database[1]
    IntegronDB = database[2]
    TransposonDB = database[3]
    coverage = coverages.split("|")
    cov0 = coverage[0]
    if cov0:
        cov0 = float(cov0)
    cov1 = coverage[1]
    if cov1:
        cov1 = float(cov1)
    cov2 = coverage[2]
    if cov2:
        cov2 = float(cov2)
    cov3 = coverage[3]
    if cov3:
        cov3 = float(cov3)
    identity = identities.split("|")
    ident0 = identity[0]
    if ident0:
        ident0 = float(ident0)
    ident1 = identity[1]
    if ident1:
        ident1 = float(ident1)
    ident2 = identity[2]
    if ident2:
        ident2 = float(ident2)
    ident3 = identity[3]
    if ident3:
        ident3 = float(ident3)
    blastn_path,blastp_path = check_dependencies()
    if blastn_path and blastp_path:
        mkdir(resultdir)
        #dir input
        if indir:
            for file in os.listdir(indir):
                filename = file.split('.')[0]
                outdir = resultdir + '/' + time.strftime('%Y%m%d%H%M%S',time.localtime())+'_'+filename
                mkdir(outdir)
                infile = os.path.join(indir,file)
                fileFormat = getFileFormat(infile)
                if fileFormat=="genbank":
                    parse_genbank(infile,outdir)
                    nucleotide = outdir+"/nucleotide.fasta"
                    protein = outdir+"/protein.fasta"
                elif fileFormat=="fasta":
                    name = re.sub('[^\w-]','',infile.split("/")[-1].split(".")[0])
                    parse_fasta(name,infile,outdir)
                    nucleotide = outdir+"/nucleotide.fasta"
                else:
                    print("Please check your file format:",infile)
                    sys.exit(-1)
                print ("Analysing file",infile," at "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                #-----
                #input sequence length
                seq_len = 0
                for rec in SeqIO.parse(nucleotide,'fasta'):
                    length = len(str(rec.seq))
                    seq_len = seq_len + length
                with open(outdir+'/seqLen.txt','w') as w:
                    w.write(str(seq_len))
                #replicon
                replicon_search(nucleotide,outdir,blastn_path,curent_dir)
                #integron
                if IntegronDB:
                    database = curent_dir+'/IntegronDB/Integron'
                    integron_finder(curent_dir,nucleotide,outdir,blastn_path,database)
                #resistance
                if resDB:
                    database = curent_dir+'/resDB/sequences'
                    AMR(nucleotide,outdir,database,blastn_path,cov0,ident0)
                #ISfinder
                if ISDatabase:
                    database = curent_dir+'/ISDatabase/IS'
                    ISfinder(nucleotide,outdir,blastn_path,cov1,ident1,database)
                    num_dict = parse_IS(outdir)
                    if num_dict != {}:
                        plot(num_dict,outdir)
                #transposon
                if TransposonDB:
                    database = curent_dir+'/TransposonDB/Transposon'
                    transposon(nucleotide,outdir,blastn_path,cov3,ident3,database)
                
                if IntegronDB:
                    while 1:
                        if os.path.exists(outdir+"/integron.finished"):
                            os.remove(outdir+"/integron.finished")
                            break
                        else:
                            time.sleep(1) 
                    #print("integron done")
                final_result = outdir+'/annotation.gff'
                catFlag = False
                for file in os.listdir(outdir):
                    if 'gff' in file:
                        catFlag = True
                if catFlag:
                    fname = codecs.open(final_result, "w",'utf-8')
                    for file in os.listdir(outdir):
                        infile1 = os.path.join(outdir,file)
                        if 'gff' in file and 'gff3' not in file:
                            x = codecs.open(infile1, "r",'utf-8')
                            fname.write(x.read().encode('utf-8').decode('utf-8'))
                            x.close()
                    fname.close()
                if os.path.exists(final_result):
                    f = open(final_result)
                    datas = f.readlines()
                    f.close()
                    with open(final_result+'3','w') as w:
                        w.write(datas[0])
                    newDatas = sorted(datas[1:],key = lambda data:(int(data.split('\t')[3].strip())))
                    for newData in newDatas:
                        with open(final_result+'3','a') as w:
                            w.write(newData)
                    gb_file = outdir+"/annotation.gb"
                    file = outdir+'/nucleotide_all.fasta'
                    if os.path.exists(file):
                        convert(final_result+"3",file,outdir,gb_file)
                    else:
                        convert(final_result+"3",nucleotide,outdir,gb_file)
                if catFlag:
                    os.system("rm "+outdir+"/*.gff")
                    os.system("rm "+outdir+"/*.gff3")
                os.remove(outdir+'/gb_location.txt')
                if os.path.exists(outdir+"/nucleotide_all.fasta"):
                    os.remove(outdir+"/nucleotide_all.fasta")
                if os.path.exists(outdir+"/nucleotide.fasta"):
                    os.remove(outdir+"/nucleotide.fasta")
                if os.path.exists(outdir+"/integron.blast.out"):
                    os.remove(outdir+"/integron.blast.out")
                if os.path.exists(outdir+"/integron.fa"):
                    os.remove(outdir+"/integron.fa")
                if os.path.exists(outdir+"/seqLen.txt"):
                    os.remove(outdir+"/seqLen.txt")
                if os.path.exists(outdir+"/nucleotide.integrons"):
                    os.rename(outdir+"/nucleotide.integrons",outdir+"/integrons.detail.xls")
                for file in os.listdir(outdir):
                    if '.xls' in file:
                        file1 = os.path.join(outdir,file)
                        file2 = file1.replace('.xls','.tsv')
                        os.system("mv %s %s"%(file1,file2))
                #-----
        else:
            #genbank input
            if genbank:
                fileFormat = getFileFormat(genbank)
                if fileFormat == "genbank":
                    pass
                else:
                    print("Please check your file format:",genbank)
                    sys.exit(-1)
                print ("Analysing file",genbank," at "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                parse_genbank(genbank,resultdir)
                nucleotide = resultdir+"/nucleotide.fasta"
                protein = resultdir+"/protein.fasta"
                
            #fasta input
            else:
                fileFormat = getFileFormat(nucleotide)
                if fileFormat == "fasta":
                    pass
                else:
                    print("Please check your file format:",nucleotide)
                    sys.exit(-1)
                print ("Analysing file",nucleotide," at "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                name = re.sub('[^\w-]','',nucleotide.split("/")[-1].split(".")[0])
                parse_fasta(name,nucleotide,resultdir)
                nucleotide = resultdir+"/nucleotide.fasta"
                
            
            #-----
            #input sequence length
            seq_len = 0
            for rec in SeqIO.parse(nucleotide,'fasta'):
                length = len(str(rec.seq))
                seq_len = seq_len + length
            with open(resultdir+'/seqLen.txt','w') as w:
                w.write(str(seq_len))
            #replicon
            replicon_search(nucleotide,resultdir,blastn_path,curent_dir)
            #integron
            if IntegronDB:
                database = curent_dir+'/IntegronDB/Integron'
                integron_finder(curent_dir,nucleotide,resultdir,blastn_path,database)
            #resistance
            if resDB:
                database = curent_dir+'/resDB/sequences'
                AMR(nucleotide,resultdir,database,blastn_path,cov0,ident0)
            #ISfinder
            if ISDatabase:
                database = curent_dir+'/ISDatabase/IS'
                ISfinder(nucleotide,resultdir,blastn_path,cov1,ident1,database)
                num_dict = parse_IS(resultdir)
                if num_dict != {}:
                    plot(num_dict,resultdir)
            #transposon
            if TransposonDB:
                database = curent_dir+'/TransposonDB/Transposon'
                transposon(nucleotide,resultdir,blastn_path,cov3,ident3,database)
            if IntegronDB:
                while 1:
                    if os.path.exists(resultdir+"/integron.finished"):
                        os.remove(resultdir+"/integron.finished")
                        break
                    else:
                        time.sleep(1)            
                #print("integron done")
            final_result = resultdir+'/annotation.gff'
            catFlag = False
            for file in os.listdir(resultdir):
                if 'gff' in file:
                    catFlag = True
            if catFlag:
                fname = codecs.open(final_result, "w",'utf-8')
                for file in os.listdir(resultdir):
                    infile1 = os.path.join(resultdir,file)
                    if 'gff' in file and 'gff3' not in file:
                        x = codecs.open(infile1, "r",'utf-8')
                        fname.write(x.read().encode('utf-8').decode('utf-8'))
                        x.close()
                fname.close()
            if os.path.exists(final_result):
                f = open(final_result)
                datas = f.readlines()
                f.close()
                with open(final_result+'3','w') as w:
                    w.write(datas[0])
                newDatas = sorted(datas[1:],key = lambda data:(int(data.split('\t')[3].strip())))
                for newData in newDatas:
                    with open(final_result+'3','a') as w:
                        w.write(newData)
                gb_file = resultdir+"/annotation.gb"
                file = resultdir+'/nucleotide_all.fasta'
                if os.path.exists(file):
                    convert(final_result+"3",file,resultdir,gb_file)
                else:
                    convert(final_result+"3",nucleotide,resultdir,gb_file)
            if catFlag:
                os.system("rm "+resultdir+"/*.gff")
                os.system("rm "+resultdir+"/*.gff3")
            os.remove(resultdir+'/gb_location.txt')
            if os.path.exists(resultdir+"/nucleotide_all.fasta"):
                os.remove(resultdir+"/nucleotide_all.fasta")
            if os.path.exists(resultdir+"/nucleotide.fasta"):
                os.remove(resultdir+"/nucleotide.fasta")
            if os.path.exists(resultdir+"/integron.blast.out"):
                os.remove(resultdir+"/integron.blast.out")
            if os.path.exists(resultdir+"/integron.fa"):
                os.remove(resultdir+"/integron.fa")
            if os.path.exists(resultdir+"/seqLen.txt"):
                os.remove(resultdir+"/seqLen.txt")
            if os.path.exists(resultdir+"/nucleotide.integrons"):
                os.rename(resultdir+"/nucleotide.integrons",resultdir+"/integrons.detail.xls")
            for file in os.listdir(resultdir):
                if '.xls' in file:
                    file1 = os.path.join(resultdir,file)
                    file2 = file1.replace('.xls','.tsv')
                    os.system("mv %s %s"%(file1,file2))
            #-----
        if os.path.exists(resultdir+"/Results_Integron_Finder_nucleotide"):
            os.system("rm -r %s/Results_Integron_Finder_nucleotide"%resultdir)
        print("All Done"," at "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    else:
        print("program terminated because lack of dependencies")
        sys.exit(-1)


def replicon_search(nucleotide,resultdir,blastn_path,curent_dir,cov=60,ident=90):
    database = curent_dir+'/repliconDB/replicon.fasta'
    cmd = blastn_path + ' -task blastn -dust no -perc_identity '+str(ident)+' -query ' + nucleotide + ' -db ' + database + ' -outfmt \"6 qseqid sseqid pident length  qstart qend slen sstrand  gaps\" -out ' + resultdir + '/replicon.blast.out'
    subprocess.run(cmd,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    blast_out = resultdir + '/replicon.blast.out'
    blast_tsv = resultdir+'/replicon.tsv'
    with open(blast_tsv,'w') as w:
        w.write('Plasmid\tIdentity\tQuery/Template\tContig\tPosition\tAccession\n')
    if os.path.exists(blast_out):
        if os.path.getsize(blast_out) != 0:
            f = open(blast_out)
            for lines in f:
                lines = lines.strip()
                if lines:
                    data = lines.split("\t")
                    qseqid = data[0]
                    sseqid = data[1]
                    pident = float(data[2])
                    length = int(data[3])
                    pos = data[4]+'..'+data[5]
                    slen = int(data[6])
                    if (length/slen)*100 >= cov:
                        tmp = []
                        plasmid = sseqid.split('_')[0]
                        acc = sseqid.split('_')[-1]
                        tmp.append(plasmid)
                        tmp.append(str(pident))
                        tmp.append(str(length)+'/'+str(slen))
                        tmp.append(qseqid)
                        tmp.append(pos)
                        tmp.append(acc)
                        with open(blast_tsv,'a') as w:
                            w.write('\t'.join(tmp)+'\n')
            f.close()
        os.system("rm "+blast_out)
    print("replicon search done")

if __name__ == "__main__":
    run()