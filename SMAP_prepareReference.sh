#!/bin/sh

# Program to prepare a reference genome combining human and mouse genomes (and transcritpomes)

# Usage: SMAP



usage="$(basename "$0") [-h] -t HumanGenome.fasta -u HumanTranscriptome.gtf -m MouseGenome.fasta -s MouseTranscriptome.gtf -o outputDirectory
Program to prepare a reference genome/transcriptome combining human and mouse genomes/transcriptomes

where:
-h  show this help text
-t is followed by the files defining the genome of the tumor compartment (usually human)
-u is followed by the files defining the transcriptome of the tumor compartment (usually human)
-m is followed by the files defining the genome of the stromal compartment (usually mouse)
-s is followed by the files defining the  transcriptome of the stromal compartment (usually mouse)
-o defines the directory that will contain all the combined genome information
"


tg=''
tt=''
sg=''
st=''
output=''

while getopts ':ht:u:m:s:o:' option; do
case "$option" in
h)
echo "$usage"
exit
;;

t)

tg="$OPTARG"
;;

u)

tt="$OPTARG"
;;

m)

sg="$OPTARG"
;;

s)

st="$OPTARG"
;;

o)

output=$OPTARG
;;

\?)
printf "invalid option: -%s\n" "$OPTARG" >&2
echo "$usage" >&2
exit 1
;;
esac
done
shift $((OPTIND - 1))






# Testing file existence
if [ "$tg" = '' ]; then printf "Option  -t is missing\n\n $usage" ;exit;fi
if [ "$tt" = '' ]; then printf "Option  -u is missing\n\n $usage";exit ;fi
if [ "$sg" = '' ]; then printf "Option  -m is missing\n\n $usage" ;exit ;fi
if [ "$st" = '' ]; then printf "Option  -s is missing\n\n $usage";exit ;fi
if [ "$output" = '' ]; then printf "Option  -o is missing\n\n $usage";exit ;fi



if [ ! -e "$tg" ]; then printf "Tumor genome file does not exists  (Option -t) \n" "$usage";exit ;fi
if [ ! -e "$tt" ]; then printf "Tumor transcriptome file does not exists  (Option -u) \n" "$usage";exit ;fi
if [ ! -e "$sg" ]; then printf "Stromal genome file does not exists  (Option -m) \n" "$usage";exit ;fi

if [ ! -e "$st" ]; then printf "Stromal transcriptome file does not exists  (Option -s) \n" "$usage";exit ;fi


if [[ ! -e $output ]]; then
mkdir $output
elif [[ ! -d $output ]]; then
echo "$dir already exists but is not a directory" 1>&2
exit
fi



# creating new combined reference genomes and transcriptomes

grep -v "^#" $tt | sed "s/^/hs/g" >$output"/combined.gtf"
grep -v "^#" $st | sed "s/^/mmu/g" >>$output"/combined.gtf"

sed "s/^>/>hs/g" $tg >$output"/combined.fasta"
sed "s/^>/>mmu/g" $sg >>$output"/combined.fasta"

printf "Combined genomes and transcriptomes are in directory:\n$output\n\n"
printf "Combined Genome file:\n$output/combined.fasta\n\n"
printf "Combined Transcriptome file:\n$output/combined.gtf\n\n"

printf "For further analysis run the STAR genome indexing :\n"

printf "mkdir $output/StarIndexed_combined\n"
printf " STAR --runMode genomeGenerate --genomeDir $output/StarIndexed_combined   \ \n"
printf " --genomeFastaFiles $output/combined.fasta   \ \n"
printf " --sjdbGTFfile $output/combined.gtf   \ \n"
printf " --sjdbOverhang 100   # set to the length of base pair sequencing \n"
printf "\n\n Optionally set --runThreadN option to the number of threads to use.\n\n"




