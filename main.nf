#!/usr/bin/env nextflow

// process illumina run and do basic bioinforamtics and qc

params.run = ""
params.plates = ""
Channel.fromPath(params.plates).set { plates_ch }
params.adapter = ""
params.outdir = ""
params.projectLocations = ""
params.runloc = ""

adapter_file = file(params.adapter)
project_locations = file(params.projectLocations)
run_locations = file(params.runloc)

process combine_plate_meta {
  container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/trim_reads'
  label 'combine_reads'

  input:
  file '*' from plates_ch.collect()

  output:
  file "plate_metadata.txt" into plate_meta

  """
  cat * | sed 's/10411pl/10411_pl/'g | sed 's/10447pl/10447_pl/'g | sed 's/10290pl/10290_pl/'g| sed 's/10292pl/10292_pl/'g | sed -e 's/\t/,/'g  > plate_metadata.txt
  
  """
}

process combine_reseqs {
  container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/trim_reads'
  label 'general_software'


  input:
  set row, well, barcode, samplename, reps, species, project, blank, gc, mb, cov from plate_meta.splitCsv()
  file project_locs from project_locations

  output:
  set project, barcode, file("${barcode}_*_1.fastq.gz"), file("${barcode}_*_2.fastq.gz"), 'project_url.txt', 'clean_samp.txt' into conc_reads, conc_reads2
  set barcode, project, species, gc, mb, cov, reps into sample_meta
  file 'project_url.txt' into project_path_ch
  file 'clean_samp.txt' into trim_ch

  """
  clean_samp=\$(echo "$samplename" | tr -cd '[:alnum:]\n ')
  echo ${barcode}_\$clean_samp > clean_samp.txt

  aws s3 sync $run_locations read1 --exclude "*" --include "*/${barcode}*_R1_001.fastq.gz"
  aws s3 sync $run_locations read1 --exclude "*" --include "climb/${barcode}*_1.fastq.gz"
  cat read1/*/*fastq.gz > ${barcode}_\${clean_samp}_1.fastq.gz

  aws s3 sync $run_locations read2 --exclude "*" --include "climb/${barcode}*_2.fastq.gz"
  aws s3 sync $run_locations read2 --exclude "*" --include "*/${barcode}*_R2_001.fastq.gz"
  cat read2/*/*fastq.gz > ${barcode}_\${clean_samp}_2.fastq.gz

  rm -rf read*

  grep $project $project_locs > project_url.txt

  """
}


conc_reads.map { it -> tuple(it[0], it[1], it[2] , it[3] , file(it[4]).text.trim(), file(it[5]).text.trim()) }.set {conc_reads_trim}


conc_reads_merge = conc_reads_trim.groupTuple(by:0).map { key, barcode ,r1,r2, project_path, samplename -> tuple( groupKey(key, barcode.size()), barcode, r1, r2, project_path, samplename) }
conc_reads_trans = conc_reads_merge.transpose()


process trim_reads {
  label 'general_software'
  container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/trim_reads'
  cpus 4

  publishDir "${params.outdir}/${project_path}/reads"

  input:
  set project, pair_id, file(read1), file(read2), project_path, samplename from conc_reads_trans
  file adapter from adapter_file

  output:
  set project, pair_id, samplename, file("${samplename}_1_trimmed.fastq.gz"), file("${samplename}_2_trimmed.fastq.gz"), file("${samplename}_U1_trimmed.fastq.gz"), file("${samplename}_U2_trimmed.fastq.gz"), project_path into read_upload, read_upload2, assemble, kraken

  """
  java -jar /tmp/software/trimmomatic-0.39.jar PE -phred33 -threads 4 ${read1} ${read2} ${samplename}_1_trimmed.fastq.gz ${samplename}_U1_trimmed.fastq.gz ${samplename}_2_trimmed.fastq.gz ${samplename}_U2_trimmed.fastq.gz ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  
  echo "trim_reads $read1 $samplename $project_path $project" > ${samplename}_trimming_check.txt
  echo "java -jar /tmp/software/trimmomatic-0.39.jar PE -phred33 -threads 4 ${read1} ${read2} ${samplename}_1_trimmed.fastq.gz ${samplename}_U1_trimmed.fastq.gz ${samplename}_2_trimmed.fastq.gz ${samplename}_U2_trimmed.fastq.gz ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" >> ${samplename}_trimming_check.txt

  aws s3 cp ${samplename}_trimming_check.txt "${params.outdir}/${project_path}/qc/"
  aws s3 cp $read1 "${params.outdir}/${project_path}/" --acl public-read
  aws s3 cp $read2 "${params.outdir}/${project_path}/" --acl public-read
  """

}


upload_reads = read_upload.groupTuple()


process assemble_samples {
  label 'assembly'
  publishDir "${params.outdir}/${project_path}/contigs", pattern: "${pair_id}.fast*"
  container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/assembly'

  input:
  set project, barcode, pair_id, file(read1_t), file(read2_t), file(a), file(b), project_path from assemble

  output:
  set project, barcode, pair_id, project_path, file(read1_t), file(read2_t), file("${pair_id}.fasta") into qc_samples
  set project, file("${pair_id}.fasta"), pair_id, project_path into prokka
  file("${pair_id}.fastg") into graph

  """
  spades.py -t 4 --careful --pe1-1 $read1_t --pe1-2 $read2_t -o ${pair_id}_spades && mv ${pair_id}_spades/contigs.fasta ${pair_id}.fasta && mv ${pair_id}_spades/assembly_graph.fastg ${pair_id}.fastg || { 
      touch ${pair_id}.fasta
      touch ${pair_id}.fastg
  }


  """
}

process annotate_samples {

   label 'general_software'
   publishDir "${params.outdir}/${project_path}/contigs"
   container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/annotation'
   cpus 4

   input:
   set project, file(assembly), name, project_path from prokka

   output:
   set project, name, file("${name}_annotation/*"), project_path into annotations, annotations2
   //set project, name, file("${name}_annotation/*.fna"), project_path into project_qc


   """
   sed -e 's/_length.*//' $assembly > cleanheaders.fasta
   /tmp/prokka/bin/prokka --locustag $name --outdir ${name}_annotation --prefix $name --rawproduct --cpus 4 --force cleanheaders.fasta && \
   rm ${name}_annotation/${name}.log || {
       mkdir ${name}_annotation
       touch ${name}_annotation/${name}.fna
   }

   """
}

// Combine metadata and assembly output data for sample by sample QC

// qc_and_meta = qc_samples.combine(sample_meta, by: 0)

process sample_qc {
   container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/sample_qc'
   label 'sina'

   input:
   set project, barcode, name, proj_loc, file(R1), file(R2), file(fasta) from qc_samples

   output:
   set barcode, project, name, proj_loc, file(R1), file(R2), file(fasta), file("16S_${name}.fasta"), file("16S_silva_${name}.fasta"), file("16S_silva_${name}.csv"), file("16S_silva_${name}.results"), file("cov.txt"), file("cov_noZeros.txt"), file("median_inserts.txt"), file("num_inserts_300.txt"), file("read_stats.txt"), file("size_gc.txt") into qcFiles_ch, test

   """
   if [ -s $fasta ]
   then

   	bwa index $fasta

	bwa mem -t 12 $fasta $R1 $R2 | samtools view -@ 12 -uS - | samtools sort -@ 12 - -o ${name}.self.sorted.bam

   	/tmp/bedtools2/bin/genomeCoverageBed -d -ibam ${name}.self.sorted.bam > ${name}_self_coverage.txt

   	awk '{sum+=\$3} END {print sum/NR}' ${name}_self_coverage.txt > cov.txt

   	awk '{if(\$3>0) print}' ${name}_self_coverage.txt | awk '{sum+=\$3} END {print sum/NR}' > cov_noZeros.txt

   	samtools view ${name}.self.sorted.bam | awk '{if(\$9>0 && \$9<1000) print \$9}' > all_inserts.txt

   	cat all_inserts.txt | datamash median 1 > median_inserts.txt

   	awk -v name=${name} '{if(\$1>300) print \$1}' all_inserts.txt | wc -l > num_inserts_300.txt

   	export PATH=$PATH:/tmp/barrnap/binaries/linux:/tmp/bedtools2/bin

   	/tmp/barrnap/bin/barrnap $fasta --outseq 16S_${name}.fasta

   	/tmp/sina-1.6.0-linux/bin/sina -i 16S_${name}.fasta -o 16S_silva_${name}.fasta --ptdb /tmp/SSURef_NR99_128_SILVA_07_09_16_opt.arb --search --search-db  /tmp/SSURef_NR99_128_SILVA_07_09_16_opt.arb --lca-fields tax_slv --meta-fmt csv >&2

   	cut -d',' -f8 16S_silva_${name}.csv | tr ' ' '\n' | sed '/^\$/'d | while read line ; do grep \${line%%.*} /tmp/SILVA_128_SSURef_Nr99_tax_silva.fasta | awk -F';' -v name=${name} '{print name"\t"\$(NF-1)"\t"\$NF"\t"line}' ; done | sort -u > 16S_silva_${name}.results

   	echo -e ">Seq1\n\$(sed '/>.*/d' $fasta)" | seqkit fx2tab --gc --length | awk '{print \$3"\t"\$4}' > size_gc.txt
   else
   	echo "0" > cov.txt
        echo "0" > cov_noZeros.txt
        echo "0" > all_inserts.txt 
        echo "0" > median_inserts.txt
        echo "0" > num_inserts_300.txt
        echo "-\t-" > 16S_silva_${name}.results
   	echo "0\t0" > size_gc.txt 
   	touch 16S_${name}.fasta 
   	touch 16S_silva_${name}.fasta
   	touch 16S_silva_${name}.csv
   	touch 16S_silva_${name}.results
   fi

   seqkit stats $R1 | tail -n+2 > read_stats.txt
  
   """ 

}


//combine qcoutputs with metadata

qc_and_meta = qcFiles_ch.combine(sample_meta, by: 0).map { it -> tuple(it[0], it[1], it[2], it[3] , file(it[4]), file(it[5]) , file(it[6]), file(it[7]), file(it[8]), file(it[9]), file(it[10]), file(it[11]).text.trim(), file(it[12]).text.trim(), file(it[13]).text.trim(), file(it[14]).text.trim(), file(it[15]).text.trim(), file(it[16]).text.trim(), it[17], it[18].tokenize(' '), it[19], it[20], it[21], it[22]) }.set { qc_and_meta_ch }



process qc_report {
   container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/project_qc'
   label 'general_software'
   publishDir "${params.outdir}/${project_loc}/qc"

   input:
   set barcode, project, name, project_loc, file(R1) , file(R2), file(fasta) , file(fasta16S), file(silva16Sfasta), file(silva16SCSV), file(silva16SResults), cov, covNoZeros, median_ins, num_ins300, read_stats, size_gc, proj2, species, gc, mb, coverage, reps from qc_and_meta_ch

   output:
   set project, barcode, name, project_loc, file("${name}_read_stats.txt"), file("${name}_coverage.txt"), file("${name}_gc.txt"), file("${name}_size.txt"), file("${name}_16S.txt"), file("$silva16SCSV") into group_qc_stats

   """
   array_read_stats=($read_stats)
   array_size_gc=($size_gc)   


   # for organisms that don't have an expected size make them 5mb
   megabases=$mb

   if [ -s \$megabases ]; then megabases=5 ; echo \$megabases ; fi

   printf "${name}\t$median_ins\t$cov\t$covNoZeros\t\${array_read_stats[3]//,/}\t$num_ins300\n" > ${name}_read_stats.txt
   
   predicted_coverage=\$(printf "\${array_read_stats[4]//,/}\t\$megabases" | awk -v megabases=\$megabases '{print (\$1 * 2) / (megabases * 1000000) }')

   echo | awk -v name=${name} -v pred_cov=\$predicted_coverage 'BEGIN { OFS="\t";} { if ($cov == "-" && $reps == 1 || $cov <= $coverage && $reps == 1 ) 
		 	print name,$cov,$coverage,pred_cov,$reps,"repeat"
	  	 else if ($reps >=2 && $cov == "-" || $cov <= $coverage && $reps >= 2)
		       	print name,$cov,$coverage,pred_cov,$reps,"fail"
		 else 
			print name,$cov,$coverage,pred_cov,$reps,"pass"
                }' > ${name}_coverage.txt
  

   exp_gc=$gc

   if [ -s \$exp_gc ]; then exp_gc=50 ; fi
 
   echo | awk -v exp_gc=\$exp_gc -v name=${name} -v pred_gc=\${array_size_gc[1]} 'BEGIN { OFS="\t";} { if (pred_gc > exp_gc+2 || pred_gc < exp_gc-2) 
							print name,pred_gc,exp_gc,"warning"
						 else
							print name,pred_gc,exp_gc,"pass"
						}' > "${name}"_gc.txt

   echo | awk -v megabases=\$megabases -v name=${name} -v pred_size=\${array_size_gc[0]} 'BEGIN { OFS="\t";} { if (\$megabases*1000000>pred_size+500000 || \$megabases*1000000<pred_size-500000 || pred_size == 0)
                                                        print name,pred_size,megabases*1000000,"warning"
                                                 else
                                                        print name,pred_size,megabases*1000000,"pass"
                                                }' > ${name}_size.txt

   while read line ; 
   do echo -e "\$line" |
	awk 'BEGIN { OFS="\t";} { if (\$3 == "${species[0]}" && \$4 == "${species[1]}") 
			print \$1,\$2,\$3" "\$4,"${species[0]} ${species[1]}","pass-species-match"
	       	else if (\$3 == "${species[0]}" && \$4 != "${species[1]}")
			print \$1,\$2,\$3" "\$4,"${species[0]} ${species[1]}","genus-match"
		else
			print \$1,\$2,\$3" "\$4,"${species[0]} ${species[1]}","warning"
	     }'
   done  < $silva16SResults > ${name}_16S.txt 

   """

}


group_qc_stats_ch = group_qc_stats.groupTuple()

process group_qc_stats {
   container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/project_qc'
   label 'general_software'

   input: 
   set project, barcode, name, project_loc, file(read_stats), file(coverage), file(gc), file(size), file(taxon_16S), file(taxon_16S_csv) from group_qc_stats_ch

   output:
   file "qc_report.txt" into qc_reports

   """
   aws s3 sync ${params.outdir}/${project_loc[0]}/qc/ qc
   
   echo -e "[General Stats]\n" >> qc_report.txt
   echo -e "SampleID\tMedianInsertSize\tMeanCoverage\tCoverageNoZeros\tNumberReads\tNumberReadsIns300" | cat - qc/*read_stats.txt >> qc_report.txt
   echo "\n[Coverage]\n" >> qc_report.txt
   echo -e "SampleID\tCoverage\tDesiredCoverage\tEstimatedCoverage\tRepetitions\tStatus" | cat - qc/*coverage.txt >> qc_report.txt
   echo "\n[GC]\n" >> qc_report.txt
   echo -e "SampleID\tPredictedGC\tExpectedSize\tStatus" | cat - qc/*gc.txt >> qc_report.txt
   echo "\n[Assembly Size]\n" >> qc_report.txt
   echo -e "SampleID\tAssembledSize\tExpectedSize\tStatus" | cat - qc/*size.txt >> qc_report.txt
   echo "\n[16S]\n" >> qc_report.txt 
   echo -e "SampleID\tPredictedGenus\tPredictedSpecies\tExpectedSpecies\tStatus" | cat - qc/*_16S.txt | uniq >> qc_report.txt

   aws s3 cp qc_report.txt "${params.outdir}/${project_loc[0]}/" --acl public-read
   
   cat $read_stats | jq --raw-input --slurp  'split("\n")| map(split("\t")) | .[0:-1]'| jq -s '{data:.[]}' > json_reads.txt 
   aws s3 cp json_reads.txt "${params.outdir}/${project_loc[0]}/" --acl public-read

   """

}


// Combine in the trimmed reads and the contigs folder (prokka) for group QC and upload

upload_contigs = annotations.groupTuple()


process zip_upload_reads {
   container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/project_qc'
   label 'general_software'

   input:
   set project, barcode, name, file(read1), file(read2), file(readU1), file(readU2), project_path from upload_reads

   """
   yum -y install zip
   aws s3 sync ${params.outdir}/${project_path[0]}/reads reads
   zip -r reads reads
   aws s3 cp reads.zip "${params.outdir}/${project_path[0]}/" --acl public-read

   """
}

process zip_upload_contigs {
   container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/project_qc'
   label 'general_software'

   input:
   set project, name, file(annotations), project_path from upload_contigs


   output:
   set project, name, file('contigs/*fasta'), project_path into group_qc
   """
   yum -y install zip
   aws s3 sync ${params.outdir}/${project_path[0]}/contigs contigs --exclude "cov.txt" --exclude "inserts.txt" --exclude "cov_noZeros.txt"
   zip -r contigs contigs
   aws s3 cp contigs.zip "${params.outdir}/${project_path[0]}/" --acl public-read

   """

}

process project_qc {
   container '281328365596.dkr.ecr.eu-west-1.amazonaws.com/project_qc'
   label 'general_software'

   input:
   set project, name, file(annotations), project_path from group_qc

   """
   if [ `cat $annotations | wc -l | cut -d" " -f1` > 1 ]
   then
   	quast.py $annotations -o quast_results

        aws s3 cp quast_results/report.html "${params.outdir}/${project_path[0]}/" --acl public-read

   	tail -n+2 quast_results/transposed_report.tsv | awk -F'\t' '{print \$1"\t"\$2"\t"\$3"\t"\$8"\t"\$9"\t"\$14"\t"\$15"\t"\$16"\t"\$17"\t"\$18"\t"\$19"\t"\$20"\t"\$21"\t"\$22}' | jq --raw-input --slurp  'split("\n")| map(split("\t")) | .[0:-1]'| jq -s '{data:.[]}' > assembly_json.txt
   else
   	
   	printf "assemblies_failed\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n" | jq --raw-input --slurp  'split("\n")| map(split("\t")) | .[0:-1]'| jq -s '{data:.[]}' > assembly_json.txt
   fi

   aws s3 cp assembly_json.txt "${params.outdir}/${project_path[0]}/" --acl public-read

   # make kraken file until proper kraken file is ready
   echo "'{data:.[]}'" > json_cont.txt

   

   # make empty kraken output until krakne part is ready.

   aws s3 cp json_cont.txt "${params.outdir}/${project_path[0]}/" --acl public-read
   

   """

}


