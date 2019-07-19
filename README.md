#Step_1: get patients_id
download data from [TCGA](https://portal.gdc.cancer.gov) and prepare the patients_id.txt
```bash
mkdir mutect
cat TCGALUADclinical.cases_selection.2019-03-13/clinical.tsv | cut -f 1 > ./mutect/patients_id.txt
```
#Step_2: get missense mutation from raw MAF file
```bash
cat TCGA.KIRC.mutect.2a8f2c83-8b5e-4987-8dbf-01f7ee24dc26.DR-10.0.somatic.maf | grep -v "^#" | grep Missense_Mutation | grep PASS|grep SNP| cut -f 1,5-13,55,56,68,111,116 > Missense_Mutation.PASS.txt
cp ../mutect/patients_id.txt  ./
```
#Step_3: get missense mutation for each patients
```bash
mkdir patients_missense_mutation
awk 'NR>1 {print "grep "$0" Missense_Mutation.PASS.txt > patients_missense_mutation/"$0".txt &"}' patients_id.txt > get_id.sh
chmod a+x get_id.sh
# the id looks like 745c699b-8ecf-4653-b41f-6620f11bcf39.txt 
./get_id.sh > get_id.log 2>&1 & 
# count 0 here 264M   #589 patients id txt in patients_missense_mutation
cd patients_missense_mutation 
ls -l | cut -d ":" -f 1 | cut -d " " -f 5- | sed "s/ //g" | grep "^0" | wc -l  
mkdir all wt mt 
#rm -rf all wt mt
cd ../
#download  Homo_sapiens.GRCh38.pep   data
wget http://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
cp Homo_sapiens.GRCh38.pep.all.fa  Homo_sapiens.GRCh38.pep.rename.fa
```
#Step_4: get peptite fragments  *.fa files，before run step5 
```bash
#need the right python envrioment for each *py script
source ~/conda.bashrc
source activate py2

cat >get_fragments.sh <<EOF
#!/bin/bash
##SBATCH -D /home/tangbo/scratch/work/bioin/KIDNEY_KIRC
#SBATCH -D /home/tangbo/scratch/work/bioin/skcm/patients_missense_mutationl
#SBATCH -J py_samtool
#SBATCH -o get_fragmets4_%A_%a.out
#SBATCH --partition=Lewis
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
##SBATCH --time=2:00:00   #2 hours maybe not enough use 2days
#SBATCH -t 2-00:00
#have to add 3 tools
module add samtools/samtools-1.7 	
module add bcftools/bcftools-1.7
module add htslib/htslib-1.7
srun hostname -s | sort -u >slurm.hosts
#echo "ok check this test"
python get_fragments.py 
EOF
sbatch get_fragments.sh
#check empty *fa exist
ls *.fa -l | cut -d ":" -f 1 | cut -d " " -f 5- | sed "s/ //g" | grep "^0" | wc -l  
```
#Step_5: Predict BA using netMHCPan4.0 
```bash 
#wget http://www.cbs.dtu.dk/services/NetMHCpan-4.0/data.Linux.tar.gz
#creat the linker of download data
#in netMHCpan-4.0 folder mkdir tmp as the install intruduction
export NMHOME=/home/tangbo/apps/netMHCpan-4.0
scp -r -PASS 22 wangdu@plsci1.rnet.missouri.edu:/scratch/BREAST_CANCER/HLA_list ./
#HLA_list split to A B C,as the next prediction will occupay huge size space
cat  HLA_list|grep  HLA-A>HLAa.list
cat  HLA_list|grep  HLA-B>HLAb.list
cat  HLA_list|grep  HLA-C>HLAc.list
cat  HLA_list|grep  HLA-A|wc -l   
cat  HLA_list|grep  HLA-B|wc -l   
cat  HLA_list|grep  HLA-C|wc -l  

#check the submit used py with right env
awk '{print "sbatch submit_single.sh "$0" >>all_HLA_job_ids"}' HLAa2.list > submit_HLAa_jobs.sh
chmod a+x submit_HLAa_jobs.sh
./submit_HLAa_jobs.sh  > step5a.log 2>&1 & 
#although the xls file is small but the fils is too many to over 100g
awk '{print "sbatch submit_single.sh "$0" >>all_HLA_job_ids"}' HLAb3.list > submit_HLAb_jobs.sh
chmod a+x submit_HLAb_jobs.sh
./submit_HLAb_jobs.sh  > step5b.log 2>&1 & 
awk '{print "sbatch submit_single.sh "$0" >>all_HLA_job_ids"}' HLAc2.list > submit_HLAc_jobs.sh
chmod a+x submit_HLAc_jobs.sh
./submit_HLAc_jobs.sh  > step5c.log 2>&1 & 

head -n 443 HLAa.list>HLAa1.list
awk 'FNR>=444 && FNR<=886' HLAa.list  >HLAa2.list
head -n 443 HLAb.list>HLAb1.list
awk 'FNR>=444 && FNR<=813' HLAb.list  >HLAb2.list
awk 'FNR>=814 && FNR<=1212' HLAb.list  >HLAb3.list
awk 'FNR>=1213 && FNR<=1412' HLAb.list  >HLAb4.list
head HLAc.list -n 400 >HLAc1.list 
awk 'FNR>=401 && FNR<=617' HLAc.list  >HLAc2.list

rm HLAb.wt.count HLAb.mt.count
#check run finished  HLA-C01:02 in the out  check by grep HLA-C01:02 *MHC.xls |wc -l
awk '{print "ls -l patients_missense_mutation/wt/"$0 "/* | wc -l >>  HLAb.wt.count"}' HLAb4.list  > count.wt.sh
chmod 775 count.wt.sh
./count.wt.sh
cat HLAb.wt.count|wc -l

awk '{print "ls -l patients_missense_mutation/mt/"$0 "/* | wc -l >>  HLAb.mt.count"}' HLAb4.list  > count.mt.sh
chmod 775 count.mt.sh
./count.mt.sh
cat HLAb.mt.count|wc -l        
#                             use the number == ls *MHC.xls |wc -l
paste HLAc1.list HLAc.wt.count | grep 168  | cut -f 1 > HLAc.donelist.wt
paste HLAc1.list HLAc.mt.count | grep 168  | cut -f 1 > HLAc.donelist.mt
cd patients_missense_mutation
mkdir wt_done mt_done
awk '{print "mv  wt/" $0 " wt_done/" }' ../HLAc.donelist.wt  > wt_mv.sh
chmod 775 wt_mv.sh
./wt_mv.sh
awk '{print "mv  mt/" $0 " mt_done/" }' ../HLAc.donelist.mt  > mt_mv.sh
chmod 775 mt_mv.sh
./mt_mv.sh
#store files on the remotely and reduce footprint in cluster 
rsync -avz --remove-source-files -e ssh  mt_done tangbo@plsci1.rnet.missouri.edu:/scratch/tangbo/ly/skcm/
rsync -avz --remove-source-files -e ssh  wt_done tangbo@plsci1.rnet.missouri.edu:/scratch/tangbo/ly/skcm/
#repeat step_5 step many times until all HLAs done
```
