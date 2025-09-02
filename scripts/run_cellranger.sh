mkdir -p logs

inputDir="/groups/tanaka/download"

echo ${inputDir}

#samples=`echo 22861{1..6} 310697 310698`

samples=`echo 310697 310698`


echo $samples

dataDir='raw_data'
mkdir -p $dataDir

outDir="cellranger_out"
mkdir -p $outDir

for f in ${samples};
do
    echo "---"
    echo "sample -- $f"

    mkdir -p ${dataDir}/${f}
    mkdir -p ${outDir}/${f}
    
    R1=`find ${inputDir} -name "*fastq.gz"|grep "_R1_"|grep ${f}`
    R2=`find ${inputDir} -name "*fastq.gz"|grep "_R2_"|grep ${f}`

    echo $R1
    echo "-----"
    echo $R2

    name_R1=`basename ${R1}`
    name_R1=`echo ${name_R1}|tr '_' '\t'|cut -f2,4,5,6,7 | tr '\t' '_'`
    echo 'R1 name --'
    echo $name_R1

    name_R2=`basename ${R2}`
    name_R2=`echo ${name_R2}|tr '_' '\t'|cut -f2,4,5,6,7 | tr '\t' '_'`
    echo 'R2 name --'
    echo $name_R2
        
    script="cellranger_count_${f}.sh"
    
    cat <<EOF > $script
#!/usr/bin/bash	
#SBATCH --export=ALL	
#SBATCH --qos=short
#SBATCH --time=0-05:00:00
#SBATCH --partition=c
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1 
#SBATCH -o logs/cellranger_count_${f}.out
#SBATCH -e logs/cellranger_count_${f}.err
#SBATCH --job-name cellranger_${f}


ml load build-env/f2022
ml load cellranger/8.0.1

ln -s $R1 ${dataDir}/${f}/${name_R1}
ln -s $R2 ${dataDir}/${f}/${name_R2}

cellranger count \
--id=${f} \
--fastqs=${dataDir}/${f} \
--sample=${f} \
--transcriptome=/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/ens_mm10_10x \
--include-introns false \
--create-bam=false \
--output-dir ${outDir}/${f} \
--jobmode local \
--localcores=32 \
--localmem=120

    
EOF
    cat $script
    sbatch $script
    
    #break
    
done


