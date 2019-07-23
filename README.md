# title(add some thing)

need a simple background introduction 
### Prerequisites

Some tools are required,

```
python v2.7.14
samtools v1.7
cftools v1.7
htslib v1.7
netMHCpan v4.0
```
## Example 
use the supported files in example to run 

* **example/skcm/patients_missense_mutationl/TCGA.SKCM.mutect.4b7a5729-b83e-4837-9b61-a6002dce1c0a.DR-10.0.somatic.maf** - *download data from TCGA.*
* **example/skcm/patients_missense_mutationl/patients_id.txt** - *the records of patients' ID.*
* **example/skcm/patients_missense_mutationl/Homo_sapiens.GRCh38.pep.rename.fa** - *procced sequnce file.*
* **example/skcm/patients_missense_mutationl/HLA_list** - *the list of HLA-A,HLA-B,HLA-C unit.*
* **example/skcm/patients_missense_mutationl/submit_single.sh** - *change the header to fit your cluster.*
* **example/skcm/patients_missense_mutationl/get_fragments.py** - *to get the fragment sequnce.*
* **example/skcm/patients_missense_mutationl/get_amp_single.py** - *to make prediction.*

```bash
cd example/skcm/patients_missense_mutationl
unzip TCGA.SKCM.mutect.4b7a5729-b83e-4837-9b61-a6002dce1c0a.DR-10.0.somatic.zip
cat TCGA.SKCM.mutect.4b7a5729-b83e-4837-9b61-a6002dce1c0a.DR-10.0.somatic.maf  | grep -v "^#" | grep Missense_Mutation | grep PASS|grep SNP| cut -f 1,5-13,55,56,68,111,116 > Missense_Mutation.PASS.txt
mkdir patients_missense_mutation
awk 'NR>1 {print "grep "$0" Missense_Mutation.PASS.txt > patients_missense_mutation/"$0".txt &"}' patients_id.txt > get_id.sh
chmod a+x get_id.sh
#get the file named by pation ID, the file look like 745c699b-8ecf-4653-b41f-6620f11bcf39.txt
./get_id.sh > get_id.log 2>&1 & 

cd patients_missense_mutation 
mkdir all wt mt 
#fragment sequnces will be genertsted in wt/ mt/ all/ resprectly
python get_fragments.py 
```

```bash
#prediction step
awk '{print "sbatch submit_single.sh "$0" >>all_HLA_job_ids"}' HLA_list > submit_HLA_jobs.sh
chmod a+x submit_HLA_jobs.sh
./submit_HLA_jobs.sh  > log  2>&1 & 
```
*Notice:*  if the space is not enough in running prediction step, please split HLA_list file to some small files, such as  "HLAa.list, HLAb.list, HLAc.list", and replace above "HLA_list" with them one by one.


## Authors

* **Yang Liu** - *University of Missouri, Columbia MO, USA*
* **Duolin Wang** - *University of Missouri, Columbia MO, USA*
* **Bowen Tang** - *University of Xiamen, Xiamen Fujian, China*
* **Email** - *ylmk2@mail.missouri.edu* 
* **Email** - *wangdu@mail.missouri.edu* 
* **Email** - *tangbo@mail.missouri.edu* 


## License
GNU v2.0


