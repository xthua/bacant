#!/bin/bash
if [ $# -lt 2 ];then
        echo "sh update_resdb.sh old_res.fasta new_res.fasta"
exit 0;
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
old_res_file=$1
new_res_file=$2
res_dir=`dirname $new_res_file`
mkdir -p $res_dir

#download table
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/latest/ReferenceGeneCatalog.txt -O amrfinder.xls

#download fasta
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS -O amrfinder.fasta

#add new fasta
python3 $DIR/get_update_data.py amrfinder.xls amrfinder.fasta $old_res_file $new_res_file

#build blast database
makeblastdb -in $new_res_file -out $res_dir/Res -dbtype nucl
