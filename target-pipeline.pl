use strict;
use warnings;
use feature 'say';
use Term::ANSIColor;
use Parallel::ForkManager;
use Getopt::Long qw/GetOptions/;

my ($sample,$threads,$r1,$r2,$panel,$lane,$pbs_id,$random,$email,$gender);
GetOptions(
	's=s' => \$sample,
	'r1=s' => \$r1,
	'r2=s' => \$r2,
	'c=s' => \$panel,
	'l=s' => \$lane,
	'r=s' => \$random,
	'e=s' => \$email,
	'g=s' => \$gender,
	't=s' => \$threads,
);
my $usage =<<EOF;

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

程序简介：      遗传病分析流程
使用方法：      perl target_pipeline_web.pl arg1 arg2 ...
更新日期：      2018年3月20日

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EOF
unless($sample && $r1 && $r2 && $lane && $random && $email){
        say colored($usage, "yellow on_magenta");
        exit 0;
}

#$threads = 7;

&report("创建目录");
unless(-e "/$tans/${lane}"){`mkdir -p "/$tans/${lane}"`}
unless(-e "/$tans/${lane}/BAM"){`mkdir -p "/$tans/${lane}/BAM"`}
unless(-e "/$tans/${lane}/Regioncount"){`mkdir -p "/$tans/${lane}/Regioncount"`}
unless(-e "/$tans/${lane}/Siteplot"){`mkdir -p "/$tans/${lane}/Siteplot"`}
unless(-e "/$tans/${lane}/SNP_INDEL"){`mkdir -p "/$tans/${lane}/SNP_INDEL"`}
unless(-e "/$tans/${lane}/Statistics"){`mkdir -p "/$tans/${lane}/Statistics"`}
unless(-e "/$tans/${lane}/SV"){`mkdir -p "/$tans/${lane}/SV"`}
unless(-e "/$tans/${lane}/VCF"){`mkdir -p "/$tans/${lane}/VCF"`}
unless(-e "/$tans/${lane}/mtDNA_HotSpot"){`mkdir -p "/$tans/${lane}/mtDNA_HotSpot"`}
unless(-e "/$tans/${lane}/CNV"){`mkdir -p "/$tans/${lane}/CNV"`}
unless(-e "/$tans/${lane}/GFF"){`mkdir -p "/$tans/${lane}/GFF"`}
unless(-e "/$tans/${lane}/RawDataSelectedGene"){`mkdir -p "/$tans/${lane}/RawDataSelectedGene"`}

if($panel ne ''){
	&normal;
}
&mysql_mail_manager;

sub normal{
	my $cutadapt = '/bin/cutadapt';
	my $ref = '/local_disk/DB/hg19/Sequence/BWAIndex/genome.fa';
	my $ref_fa = '/local_disk/DB/hg19/Sequence/WholeGenomeFasta/genome.fa';
	my $picard = '/disk1/software/picard-tools-2.2.3/picard.jar';
	my $bamtools = '/disk1/software/bin/bamtools';
	my $bwa = '/disk1/software/bin/bwa';
	my $samtools = '/disk1/software/bin/samtools';
	my $bamToBed = '/disk1/software/bin/bamToBed';
	my $coverageBed = '/disk1/software/bin/coverageBed';
	my $GATK = '/disk1/software/GATK-3.7/GenomeAnalysisTK.jar';
	my $dbsnp = '/local_disk/DB/dbsnp/dbsnp_147.hg19.vcf';
	my $varscan  = '/disk1/software/VarScan.v2.3.7.jar';
	my $beddir = '/disk1/bed/';
	my $local_disk = "/local_disk/tmp/web/batchmendelian/$random/$lane/$sample";	
	my $disk = "/ssd1/batchmendelian/$random/$lane/$sample";
	my $java = "java -Djava.io.tmpdir=$local_disk -Xmx15g -Xms5g -jar";
	unless(-e $disk){`mkdir -p $disk`;}
	unless(-e $local_disk){`mkdir -p $local_disk`;}
	chdir $disk;
	&report("拷贝fastq文件");
	if ($r1=~/,/ && $r2=~/,/){
		my @reads1=split ",", $r1;
		my @reads2=split ",", $r2;
		for(my $i=0; $i<@reads1; $i++){
			&parallel(
				"cat $reads1[$i] >> ${sample}_R1.fq.gz",
				"cat $reads2[$i] >> ${sample}_R2.fq.gz"
			);
		}
	} else {
		&parallel(
			"cp $r1 ${sample}_R1.fq.gz",
			"cp $r2 ${sample}_R2.fq.gz"
		);
	}
	&report("fastq文件质控");
	chomp(my $hardware = `zcat ${sample}_R1.fq.gz | head -1 | cut -d: -f1`);
	if($hardware =~ /NS500/){
		`$cutadapt -m 80 -q10,10 --nextseq-trim=10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}else{
		`$cutadapt -m 80 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
		#`$cutadapt -m 20 -q10,10 -o $sample.clean_R1.fq.gz -p $sample.clean_R2.fq.gz -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ${sample}_R1.fq.gz ${sample}_R2.fq.gz`;
	}
	&report("将fastq比对到人类参考基因组上");
	`$bwa mem -M -t $threads $ref $sample.clean_R1.fq.gz $sample.clean_R2.fq.gz | $samtools fixmate -O bam - - | $samtools sort -\@ $threads -m 1G - $sample.sort`;
	&check_ssd1_TNAS2;
	`$samtools index $sample.sort.bam > $sample.sort.bam.bai`;
	&report("添加BAM文件头信息");
	`$java $picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=$sample.sort.bam OUTPUT=$sample.sort.header.bam RGLB=$ref RGPL=ILLUMINA RGSM=GP1 RGPU=GRP1 SORT_ORDER=coordinate CREATE_INDEX=true`;
	&report("过滤BAM文件");
	`$bamtools filter -isMapped true -isPaired true -isProperPair true -in $sample.sort.header.bam -out $sample.sort.flt.bam`;
	`rm -f $sample.sort.header.bam $sample.sort.header.bai`;
	if ($panel =~ /DMD|Clinican|Canmute|SLC26A4|Deaf_V2|lung8|LAMA2|PKU|LBS/){
		my @parallel_command=(
			"$java $picard MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.rmdup.sorted.bam M=duplication-report2.txt",
		);
		if($panel =~ /PKU/){
			push @parallel_command,"perl /disk1/software/crest/extractSClip.pl -i $sample.sort.bam --ref_genome $ref -r chr6 > $sample.extractSClip.log 2>&1";
		}else{
			for my $chr (1 .. 22,'X','Y'){
				push @parallel_command,"perl /disk1/software/crest/extractSClip.pl -i $sample.sort.bam --ref_genome $ref -r chr$chr > $sample.extractSClip.log 2>&1";
			}
		}
		&report("去除PCR冗余序列/提取soft-clip reads");
		&parallel(@parallel_command);
		&report("合并$sample.sort.bam.*.cover/$sample.sort.bam.*.sclip.txt");
		`cat $sample.sort.bam.*.cover > $sample.sort.bam.cover`;
		`cat $sample.sort.bam.*.sclip.txt > $sample.sort.bam.sclip.txt`;
	}else{
		&report("去除PCR冗余序列");
		`$java $picard MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.rmdup.sorted.bam M=duplication-report2.txt`;
	}
	&report("标记PCR冗余序列/拷贝BAM文件到公共盘/bamtools stats/bamToBed");
	&parallel(
		"$java $picard MarkDuplicates VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=$sample.sort.flt.bam OUTPUT=$sample.sort.flt.mark.bam M=duplication-report.txt",
		"cp $sample.rmdup.sorted.bam /$tans/${lane}-web/BAM/",
		"cp $sample.rmdup.sorted.bai /$tans/${lane}-web/BAM/$sample.rmdup.sorted.bam.bai",
		"$bamtools stats -in $sample.sort.bam > $sample.sort.bam.stat",
		"$bamToBed -i $sample.rmdup.sorted.bam > $sample.rmdup.sorted.bamtobed",
	);
	`rm -f $sample.sort.flt.ba?`;
	if ($panel =~ /DMD|Clinican|Canmute|SLC26A4|Deaf_V2|lung8|LAMA2|PKU|LBS/){
		my @parallel_command=(
			"$java $GATK -T BaseRecalibrator -R $ref_fa -I $sample.sort.flt.mark.bam -knownSites $dbsnp -o recal.table -nct $threads $intervals"
		);
		if ($panel =~ /PKU/ || $panel =~ /LAMA2/){
			push @parallel_command,"perl /disk1/software/crest/CREST.pl -l 150 -f $sample.sort.bam.cover -d $sample.sort.bam --ref_genome $ref -t /local_disk/DB/hg19/Sequence/blat/genome.2bit -r chr6 -p $sample.sort.bam.chr6 > $sample.chr6.CREST.log 2>&1";
		}else{
			for my $chr (1 .. 22,'X','Y'){
				push @parallel_command,"perl /disk1/software/crest/CREST.pl -l 150 -f $sample.sort.bam.cover -d $sample.sort.bam --ref_genome $ref -t /local_disk/DB/hg19/Sequence/blat/genome.2bit -r chr$chr -p $sample.sort.bam.chr$chr > $sample.chr$chr.CREST.log 2>&1";
			}
		}			
		if($panel =~ /DMD/){
			push @parallel_command,"$samtools depth -aa -r chrX:31135345-33359192 -d 100000 $sample.rmdup.sorted.bam > $sample.DMD.depth";
		}elsif($panel =~ /SLC26A4/){
		    push @parallel_command,"$samtools depth -aa -r chr7:107301080-107358252 -d 100000 $sample.rmdup.sorted.bam > $sample.SLC26A4.depth";
		}elsif($panel =~ /PKU/ || $panel =~ /LAMA2/){
			push @parallel_command,"$samtools depth -aa -r chr6:129204100-129838000 -d 100000 $sample.rmdup.sorted.bam > $sample.LAMA2.depth";
		}
		&report("并行:校正平台误差/检测融合基因");
		&parallel(@parallel_command);
		`cat *.chr*.predSV.txt > $sample.sort.bam.predSV.txt`;
		`perl /disk2/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt`;
		if($panel =~ /DMD/){
			&parallel(
				"perl /disk2/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt",
				"perl /disk2/script/batchmendelian/siteplot_dmd.pl $sample.DMD.depth"
			);
		}elsif($panel =~ /SLC26A4/){
			&parallel(
				"perl /disk2/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt",
				"perl /disk2/script/batchmendelian/siteplot_SLC26A4.pl $sample.SLC26A4.depth SLC26A4 SLC26A4 chr7:107301080-107358252"
			);
		}elsif($panel =~ /PKU/ || $panel =~ /LAMA2/){
			&parallel(
                "perl /disk2/script/batchmendelian/extract_crest_v2.pl $sample.sort.bam.predSV.txt",
                "perl /disk2/script/batchmendelian/siteplot_LAMA2.pl $sample.LAMA2.depth"
            );
		}
		
		if($panel =~ /DMD/){
			&parallel(
				"cp *.pdf /$tans/${lane}-web/Siteplot/",
				"cp *.png /$tans/${lane}-web/Siteplot/",
				"cp $sample.smoothy_depth.txt /$tans/${lane}-web/Siteplot/$sample.DMD.depth.txt"
			);
		}elsif($panel =~ /SLC26A4/){
			&parallel(
				"cp *.pdf /$tans/${lane}-web/Siteplot/",
				"cp *.png /$tans/${lane}-web/Siteplot/",
				"cp $sample.SLC26A4.depth /$tans/${lane}-web/Siteplot/$sample.SLC26A4.depth.txt"
			);
		}elsif($panel =~ /PKU/ || $panel =~ /LAMA2/){
			&parallel(
                "cp *.pdf /$tans/${lane}-web/Siteplot/",
                "cp *.png /$tans/${lane}-web/Siteplot/",
                "cp $sample.LAMA2.depth /$tans/${lane}-web/Siteplot/$sample.LAMA2.depth.txt"
            );
		}
		if (glob "${sample}*crest.xls"){
			`cp ${sample}*crest.xls "/$tans/${lane}-web/SV"`;
		}
	}else{
		&report("校正平台误差");
		`$java $GATK -T BaseRecalibrator -R $ref_fa -I $sample.sort.flt.mark.bam -knownSites $dbsnp -o recal.table -nct $threads $intervals`;
	}
	&report("输出高质量比对序列");
	`$java $GATK -T PrintReads -R $ref_fa -I $sample.sort.flt.mark.bam -BQSR recal.table -o $sample.recal.bam -nct $threads`;
	`rm -f $sample.sort.flt.mark.ba?`;
	`rm -f $sample.rmdup.sorted.ba?`;
	for my $sub_chip ((split /,/,$panel)){
		next if $sub_chip eq '';
		if($sub_chip =~ /Wscancer|Canmute|Clinican|BloodV2|Immu_All|ImmuV2|MetFull|MetAdd|geneis|Can06|lung8/){
			#MPIleup+GATK
			&report("并行:mpileup低频突变/HaplotypeCaller常规突变/coverageBed");
			&parallel(
				"$samtools mpileup -s -f $ref -l $beddir/$sub_chip.bed -d 10000 -L 10000 $sample.recal.bam | $java $varscan mpileup2cns --min-reads2 1 --min-coverage 1 --strand-filter 1 --output-vcf 1 --variants 1 --min-var-freq 0.0001 --p-value 1 --min-avg-qual 20 > $sample.varscan.$sub_chip.vcf",
				"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/$sub_chip.bed -ip 100 -D $dbsnp -o $sample.gatk.$sub_chip.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
				"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.coverage",
				"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.target.coverage"
			);
			&report("捕获效率统计");
			`perl /disk2/script/batchmendelian/depth_count.pl $sample.sort.$sub_chip.coverage $sample.sort.$sub_chip.target.coverage > $sample.$sub_chip.sort.bam.depth.statistics.xls`;
			`perl /disk2/script/batchmendelian/read_count.pl $disk/$sample.$sub_chip.sort.bam.depth.statistics.xls $sub_chip`;
			`perl /disk2/script/batchmendelian/vcf2vcf.pl $sample.gatk.$sub_chip.vcf > $sample.gatk.$sub_chip.tmp.vcf`;
			`perl /disk2/script/batchmendelian/merge_vcf_for_gatk_varscan.pl $sample.gatk.$sub_chip.tmp.vcf $sample.varscan.$sub_chip.vcf > $sample.gatk.$sub_chip.split.vcf`;
			`mv $sample.varscan.$sub_chip.vcf /$tans/${lane}-web/VCF/$sample.varscan.$sub_chip.vcf`;
		}else{
			#GATK
			if ($sub_chip eq 'Deaf_V2'){
				&report("并行:HaplotypeCaller(常规区域和miRNA区域)/coverageBed");
				&parallel(
					"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/$sub_chip.bed -ip 100 -D $dbsnp -o $sample.gatk.$sub_chip.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
					"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/${sub_chip}_miRNA_region.bed -ip 100 -D $dbsnp -o $sample.gatk.${sub_chip}.miRNA.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
					"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.coverage",
					"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.target.coverage"
				);
				&parallel(
					"mv $sample.gatk.${sub_chip}.miRNA.vcf /$tans/${lane}-web/VCF/$sample.${sub_chip}.miRNA.vcf",
					"perl /disk2/script/batchmendelian/depth_count.pl $sample.sort.$sub_chip.coverage $sample.sort.$sub_chip.target.coverage > $sample.$sub_chip.sort.bam.depth.statistics.xls"
				);
				`perl /disk2/script/batchmendelian/read_count.pl $disk/$sample.$sub_chip.sort.bam.depth.statistics.xls $sub_chip`;
			}else{
				&report("并行:HaplotypeCaller/coverageBed");
				&parallel(
					"$java $GATK -T HaplotypeCaller -R $ref_fa -I $sample.recal.bam -L $beddir/$sub_chip.bed -ip 100 -D $dbsnp -o $sample.gatk.$sub_chip.vcf -A Coverage -stand_call_conf 30.0 -nct $threads",
					"$coverageBed -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.coverage",
					"$coverageBed -d -a $sample.rmdup.sorted.bamtobed -b $beddir/$sub_chip.bed > $sample.sort.$sub_chip.target.coverage"
				);
				`perl /disk2/script/batchmendelian/depth_count.pl $sample.sort.$sub_chip.coverage $sample.sort.$sub_chip.target.coverage > $sample.$sub_chip.sort.bam.depth.statistics.xls`;
				`perl /disk2/script/batchmendelian/read_count.pl $disk/$sample.$sub_chip.sort.bam.depth.statistics.xls $sub_chip`;
			}
			#不需要mpileup的通用
			&report("split vcf");
			`perl /disk2/script/batchmendelian/vcf2vcf.pl $sample.gatk.$sub_chip.vcf > $sample.gatk.$sub_chip.split.vcf`;
		}
		#所有Panel通用
		&report("并行:ANNOVAR功能注释/拷贝统计文件");
		&parallel(
			"perl /disk1/software/annovar/table_annovar.pl $sample.gatk.$sub_chip.split.vcf /local_disk/DB/annovar_db/humandb/ -buildver hg19 -protocol refGene,cytoBand,dbNsfpInterPro,Inhouse,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,cosmic74,dbnsfp30a,clinvar_20160302,spidex,mcap10,revel,dbscsnv11,dbnsfp31a_interpro,gnomad_exome,gnomad_genome,InterVar -operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -argument \"-splicing_threshold 5 --hgvs,,,,,,,,,,,,,,,,,,\" -vcfinput -nastring . --thread $threads --maxgenethread $threads",
			"perl /disk2/script/batchmendelian/depth_count3.pl $sample.sort.$sub_chip.target.coverage $sample $sub_chip",
			"mv $sample.$sub_chip.$sub_chip\_readcount.xls /$tans/${lane}-web/Regioncount/$sample.$sub_chip\_readcount.xls",
		);
		&report("并行:vcf2gff/拷贝捕获效率统计文件");
		&parallel(
			"perl /disk2/script/batchmendelian/ParseVCFHgvs.pl -v $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf -f $sample.gatk.$sub_chip.split.vcf.refGene.exonic_variant_function -g $gender > $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff",
			"mv $sample.${sub_chip}.count.xls /$tans/${lane}-web/Statistics/$sample.${sub_chip}.Statistics.xls"
		);
		&report("关联HGMD数据库");
		`perl /disk2/script/batchmendelian/HGMD.pl $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff $sample $sub_chip`;
		&report("拷贝SNP/INDEL文件");
		if($sub_chip eq 'MR_MUT-MR_CNV'){
			if(-e "/$tans/${lane}-web/Statistics/$sample.MR_MUT-MR_CNV.Statistics.xls"){
				`rm -f /$tans/${lane}-web/Statistics/$sample.MR_MUT-MR_CNV.Statistics.xls`;
			}
			if(-e "/$tans/${lane}-web/SNP_INDEL/${sample}.MR_MUT-MR_CNV.SNP_INDEL.xls"){
				`rm -f /$tans/${lane}-web/SNP_INDEL/${sample}.MR_MUT-MR_CNV.SNP_INDEL.xls`;
			}
			goto NEXT;
		}
		&report("并行:拷贝经过注释和过滤的xls/拷贝原始未经过滤的gff/备份vcf和gff文件");
		&parallel(
			"mv $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff.hgmd.xls /$tans/${lane}-web/SNP_INDEL/${sample}.${sub_chip}.SNP_INDEL.xls",
			"mv $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf.gff /$tans/${lane}-web/GFF/$sample.$sub_chip.gff",
			"mv $sample.gatk.$sub_chip.split.vcf.hg19_multianno.vcf /$tans/${lane}-web/VCF/$sample.$sub_chip.vcf"
		);
		&report("将捕获效率统计信息追加到Statistics.xls");
		`ls /$tans/${lane}-web/Statistics/*.xls | xargs -i sed 1d {} >> Statistics.xls`;
		`sed -i '1iSample\tRaw_data_bases(Mb)\tClean_data_bases(Mb)\tAligned_bases(Mb)\tAligned\tInitial bases on target\tBase covered on target\tCoverage of target region\tEffective bases on target\tFraction of effective bases on target\tAverage sequencing depth on target\tFraction of target covered with at least 4X\tFraction of target covered with at least 10X\tFraction of target covered with at least 20X\tduplication rate\tchip' Statistics.xls`;
		`mv Statistics.xls /$tans/${lane}-web/Statistics.xls`;
NEXT:
	}
	&report("清理空间");
	`rm -rf /local_disk/tmp/web/batchmendelian/$random/$lane/$sample`;
	`rm -rf *`;
	&report("该任务已运行完成");
}


sub report{
	my $comment = shift @_;
	print  STDERR $comment;
}



sub sent_email{
	`echo -e "恭喜！您提交的遗传病任务已经分析完成 | /usr/bin/mailx -s "遗传病任务分析完成" $email`;
}

sub parallel{
	my $pm = new Parallel::ForkManager(7);
	for (@_){
		my $pid = $pm->start and next;
		&do($_);
		$pm->finish;
	}
	$pm->wait_all_children;
	sub do{
		`$_[0]`;
	}
}

sub intervals{
	my $initial_intervals = $panel;
	$initial_intervals =~ s/MitoChip,//g;
	$initial_intervals =~ s/,MitoChip//g;
	$initial_intervals =~ s/,CNV//g;
	$initial_intervals =~ s/^CNV,//g;
	my $intervals = join " ",map {$_ = "-L /disk1/bed/$_.bed"} grep {$_ ne ""} (split /,/,$initial_intervals);
	return $intervals;
}
