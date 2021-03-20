$HOSTNAME = ""
params.outdir = 'results'  


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.inputparam){params.inputparam = ""} 

Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into{g_2_reads_g15_3;g_2_reads_g15_18}

Channel.value(params.mate).into{g_3_mate_g_13;g_3_mate_g15_3;g_3_mate_g15_11;g_3_mate_g15_16;g_3_mate_g15_18;g_3_mate_g15_19;g_3_mate_g15_20;g_3_mate_g15_21}
Channel.value(params.inputparam).into{g_19_mate_g24_26;g_19_mate_g24_30;g_19_mate_g24_38}

//* params.run_Adapter_Removal =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Adapter_Removal"
//* @style @multicolumn:{seed_mismatches, palindrome_clip_threshold, simple_clip_threshold} @condition:{Tool_for_Adapter_Removal="trimmomatic", seed_mismatches, palindrome_clip_threshold, simple_clip_threshold}, {Tool_for_Adapter_Removal="fastx_clipper", discard_non_clipped}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 1000
    $CPU  = 1
    $MEMORY = 24
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal)){
g_2_reads_g15_18.into{g15_18_reads_g15_19}
g15_18_log_file_g15_11 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Adapter_Removal {

input:
 set val(name), file(reads) from g_2_reads_g15_18
 val mate from g_3_mate_g15_18

output:
 set val(name), file("reads/*.fastq")  into g15_18_reads_g15_19
 file "*.{fastx,trimmomatic}.log"  into g15_18_log_file_g15_11

errorStrategy 'retry'

when:
(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
phred = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.phred
Tool_for_Adapter_Removal = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Tool_for_Adapter_Removal
Adapter_Sequence = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
min_length = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length
seed_mismatches = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches
palindrome_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold
simple_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
discard_non_clipped = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped
    
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.remove_previous_reads
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use Cwd qw();
 
runCmd("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);

system("!{runGzip}");
my $quality="!{phred}";
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic PE -threads 1 -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.2.fastq.unpaired ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper is not suitable for paired reads.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic SE -threads 1  -phred${quality} !{file1} reads/!{name}.fastq ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        runCmd("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir = join '/', @paths;
    $inputsdir .= "/work";
    print "INFO: inputs reads will be removed if they are located in the $workdir $inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}
}


//* params.run_Trimmer =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Trimmer"
//* @style @multicolumn:{trim_length_5prime,trim_length_3prime}, {trim_length_5prime_R1,trim_length_3prime_R1}, {trim_length_5prime_R2,trim_length_3prime_R2} @condition:{single_or_paired_end_reads="single", trim_length_5prime,trim_length_3prime}, {single_or_paired_end_reads="pair", trim_length_5prime_R1,trim_length_3prime_R1,trim_length_5prime_R2,trim_length_3prime_R2}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 500
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "long"
}
//* platform
//* autofill
if (!((params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer)){
g15_18_reads_g15_19.into{g15_19_reads_g15_20}
g15_19_log_file_g15_21 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Trimmer {

input:
 set val(name), file(reads) from g15_18_reads_g15_19
 val mate from g_3_mate_g15_19

output:
 set val(name), file("reads/*q")  into g15_19_reads_g15_20
 file "*.log" optional true  into g15_19_log_file_g15_21

errorStrategy 'retry'

when:
(params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer

shell:
phred = params.Adapter_Trimmer_Quality_Module_Trimmer.phred
single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads
trim_length_5prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime
trim_length_3prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime
trim_length_5prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1
trim_length_3prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1
trim_length_5prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2
trim_length_3prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads");
system("!{runGzip}");
my $file1 = "";
my $file2 = "";
if ("!{mate}" eq "pair") {
    $file1 = "!{file1}";
    $file2 = "!{file2}";
    my $trim1 = "!{trim_length_5prime_R1}:!{trim_length_3prime_R1}";
    my $trim2 = "!{trim_length_5prime_R2}:!{trim_length_3prime_R2}";
    my $len=getLength($file1);
    print "length of $file1: $len\\n";
    trimFiles($file1, $trim1, $len);
    my $len=getLength($file2);
    print "INFO: length of $file2: $len\\n";
    trimFiles($file2, $trim2, $len);
} else {
    $file1 = "!{file1}";
    my $trim1 = "!{trim_length_5prime}:!{trim_length_3prime}";
    my $len=getLength($file1);
    print "INFO: length of file1: $len\\n";
    trimFiles($file1, $trim1, $len);
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}



sub trimFiles
{
  my ($file, $trim, $len)=@_;
    my @nts=split(/[,:\\s\\t]+/,$trim);
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="";
    my $quality="-Q!{phred}";

    if (scalar(@nts)==2)
    {
      $param = "-f ".($nts[0]+1) if (exists($nts[0]) && $nts[0] >= 0 );
      $param .= " -l ".($len-$nts[1]) if (exists($nts[0]) && $nts[1] > 0 );
      $outfile="reads/$file";  
      $com="fastx_trimmer $quality -v $param -o $outfile -i $file > !{name}.fastx_trimmer.log" if ((exists($nts[0]) && $nts[0] > 0) || (exists($nts[0]) && $nts[1] > 0 ));
      print "INFO: $com\\n";
      if ($com eq ""){
          print "INFO: Trimmer skipped for $file \\n";
          system("mv $file reads/.");
      } else {
          runCmd("$com");
          print "INFO: Trimmer executed for $file \\n";
      }
    }

    
}


sub getLength
{
   my ($filename)=@_;
   open (IN, $filename);
   my $j=1;
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $len=length($line);
     }
     $j++;
   }
   close(IN);
   return $len;
}

sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}

'''

}
}



process Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary {

input:
 file logfile from g15_19_log_file_g15_21.collect()
 val mate from g_3_mate_g15_21

output:
 file "trimmer_summary.tsv"  into g15_21_outputFileTSV_g_16

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_trimmer\\.log/){
        $file =~ /(.*)\\.fastx_trimmer\\.log/;
        my $mapper   = "fastx_trimmer";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Trimmer" ];
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "trimmer_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* params.run_Quality_Filtering =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Quality_Filtering"
//* @style @multicolumn:{window_size,required_quality}, {leading,trailing,minlen}, {minQuality,minPercent} @condition:{tool="trimmomatic", minlen, trailing, leading, required_quality_for_window_trimming, window_size}, {tool="fastx", minQuality, minPercent}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill
if (!((params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering)){
g15_19_reads_g15_20.into{g15_20_reads_g_13}
g15_20_log_file_g15_16 = Channel.empty()
} else {


process Adapter_Trimmer_Quality_Module_Quality_Filtering {

input:
 set val(name), file(reads) from g15_19_reads_g15_20
 val mate from g_3_mate_g15_20

output:
 set val(name), file("reads/*q")  into g15_20_reads_g_13
 file "*.{fastx,trimmomatic}_quality.log" optional true  into g15_20_log_file_g15_16

errorStrategy 'retry'

when:
(params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering    

shell:
tool = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.tool
phred = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.phred
window_size = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size
required_quality_for_window_trimming = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming
leading = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.leading
trailing = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing
minlen = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen


// fastx parameters
minQuality = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality
minPercent = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent

remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 ="";
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads unpaired");
system("!{runGzip}");
my $param = "SLIDINGWINDOW:"."!{window_size}".":"."!{required_quality_for_window_trimming}";
$param.=" LEADING:"."!{leading}";
$param.=" TRAILING:"."!{trailing}";
$param.=" MINLEN:"."!{minlen}";

my $quality="!{phred}";

print "INFO: fastq quality: $quality\\n";
     
if ("!{tool}" eq "trimmomatic") {
    if ("!{mate}" eq "pair") {
        runCmd("trimmomatic PE -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.1.fastq.unpaired $param 2> !{name}.trimmomatic_quality.log");
    } else {
        runCmd("trimmomatic SE -phred${quality} !{file1} reads/!{name}.fastq $param 2> !{name}.trimmomatic_quality.log");
    }
} elsif ("!{tool}" eq "fastx") {
    if ("!{mate}" eq "pair") {
        print("WARNING: Fastx option is not suitable for paired reads. This step will be skipped.");
        system("mv !{file1} !{file2} reads/.");
    } else {
        runCmd("fastq_quality_filter  -Q $quality -q !{minQuality} -p !{minPercent} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx_quality.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}

##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}


'''

}
}



process merge_pairs {

input:
 set val(name),file(reads) from g15_20_reads_g_13
 val mate from g_3_mate_g_13

output:
 set val(name), file("reads/*.fastq")  into g_13_reads_g_14

errorStrategy 'retry'
maxRetries 2

script:
run_merge = params.run_Merge_Pairs
nameAll = reads.toString()
nameArray = nameAll.split(' ')
def file2;

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}

"""
#!/bin/bash
mkdir reads
if [ "${mate}" == "pair" ]; then
    if [ "${run_merge}" == "yes" ]; then
        pandaseq -F -N -T 40 -l 50 -f ${file1} -r ${file2} -w reads/${name}.fastq
    else
        cat ${file1} ${file2} > reads/${name}.fastq
    fi
else
    mv $file1 ${name}.fastq 2>/dev/null
    mv ${name}.fastq reads/.
fi

"""
 
}


if (!((params.run_Extract_Barcode && (params.run_Extract_Barcode == "yes")))){
g_13_reads_g_14.into{g_14_reads_g24_38}
} else {

process Extract_Barcode {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /ext\/.*.fastq$/) "aavbarcodes/$filename"
}

input:
 set val(name),file(reads) from g_13_reads_g_14

output:
 set val(name),file("ext/*.fastq")  into g_14_reads_g24_38

errorStrategy 'retry'
maxRetries 3

when:
(params.run_Extract_Barcode && (params.run_Extract_Barcode == "yes"))

script:
five_prime = params.Extract_Barcode.five_prime
three_prime = params.Extract_Barcode.three_prime
min_len = params.Extract_Barcode.min_len
max_len = params.Extract_Barcode.max_len
mismatch = params.Extract_Barcode.mismatch

"""
#!/bin/sh 
mkdir -p ext
cutadapt -a ${three_prime} -n ${mismatch} -o ext/${name}.tmp ${reads}
cutadapt -g ${five_prime} -m ${min_len} -M $max_len -n ${mismatch} -o ext/${name}.fastq ext/${name}.tmp

"""

}
}


//* params.run_Variant_Sequential_Mapping =   "yes"   //* @dropdown @options:"yes","no" @show_settings:"Sequential_Mapping" @description:"Filters out or quantify given sequence sets."

//both bowtie and bowtie2 indexes located in same path
bowtieIndexes = []
                 

//_nucleicAcidType="dna" should be defined in the autofill section of pipeline header in case dna is used.
_select_sequence = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping._select_sequence
index_directory = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping.index_directory
name_of_the_index_file = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping.name_of_the_index_file
_aligner = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping._aligner
aligner_Parameters = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping.aligner_Parameters
description = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping.description
filter_Out = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping.filter_Out
sense_antisense = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping.sense_antisense

desc_all=[]
description.eachWithIndex() {param,i -> 
    if (param.isEmpty()){
        desc_all[i] = name_of_the_index_file[i]
    }  else {
        desc_all[i] = param.replaceAll("[ |.|;]", "_")
    }
}
custom_index=[]
index_directory.eachWithIndex() {param,i -> 
    if (_select_sequence[i] == "custom"){
        custom_index[i] = param+"/"+name_of_the_index_file[i]
    }else {
        custom_index[i] = bowtieIndexes[_select_sequence[i]]
    }
}

mapList = []
paramList = []
alignerList = []
filterList = []
indexList = []
senseList = []

//concat default mapping and custom mapping
mapList = (desc_all) 
paramList = (aligner_Parameters)
alignerList = (_aligner)
filterList = (filter_Out)
indexList = (custom_index)
senseList = (sense_antisense)

mappingList = mapList.join(" ") // convert into space separated format in order to use in bash for loop
paramsList = paramList.join(",") // convert into comma separated format in order to use in as array in bash
alignersList = alignerList.join(",") 
filtersList = filterList.join(",") 
indexesList = indexList.join(",") 
senseList = senseList.join(",")

//* @style @condition:{remove_duplicates="yes",remove_duplicates_based_on_UMI_after_mapping},{remove_duplicates="no"},{_select_sequence="custom", index_directory,name_of_the_index_file,description,_aligner,aligner_Parameters,filter_Out,sense_antisense},{_select_sequence=("rRNA","ercc","miRNA","tRNA","piRNA","snRNA","rmsk","genome"),_aligner,aligner_Parameters,filter_Out,sense_antisense}  @array:{_select_sequence,_select_sequence, index_directory,name_of_the_index_file,_aligner,aligner_Parameters,filter_Out,sense_antisense,description} @multicolumn:{_select_sequence,_select_sequence,index_directory,name_of_the_index_file,_aligner,aligner_Parameters,filter_Out, sense_antisense, description},{remove_duplicates,remove_duplicates_based_on_UMI_after_mapping}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 4
    $MEMORY = 20
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 2000
    $CPU  = 4
    $MEMORY = 20
    $QUEUE = "long"
}
//* platform
//* autofill
if (!(params.run_Variant_Sequential_Mapping == "yes")){
g_14_reads_g24_38.into{g24_38_reads}
g24_38_bowfiles_g24_26 = Channel.empty()
g24_38_bam_file_g24_35 = Channel.empty()
g24_38_bam_file_g24_36 = Channel.empty()
g24_38_bam_index_g24_35 = Channel.empty()
g24_38_bam_index_g24_36 = Channel.empty()
g24_38_filter_g24_26 = Channel.empty()
g24_38_log_file_g24_30 = Channel.empty()
} else {


process Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping {

input:
 set val(name), file(reads) from g_14_reads_g24_38
 val mate from g_19_mate_g24_38

output:
 set val(name), file("final_reads/*q")  into g24_38_reads
 set val(name), file("bowfiles/?*") optional true  into g24_38_bowfiles_g24_26
 file "*/*_sorted.bam" optional true  into g24_38_bam_file_g24_35
 file "*/*_sorted.bam.bai" optional true  into g24_38_bam_index_g24_35
 val filtersList  into g24_38_filter_g24_26
 file "*/*_sorted.dedup.bam" optional true  into g24_38_bam_file_g24_36
 file "*/*_sorted.dedup.bam.bai" optional true  into g24_38_bam_index_g24_36
 file "*/*_duplicates_stats.log" optional true  into g24_38_log_file_g24_30

errorStrategy 'retry'
maxRetries 2

when:
params.run_Variant_Sequential_Mapping == "yes"

script:
nameAll = reads.toString()
nameArray = nameAll.split(' ')
def file2;

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}

remove_duplicates = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping.remove_duplicates
remove_duplicates_based_on_UMI_after_mapping = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping.remove_duplicates_based_on_UMI_after_mapping
remove_previous_reads = params.Sequential_Variant__Mapping_Module_Variant_Sequential_Mapping.remove_previous_reads

"""
#!/bin/bash
mkdir reads final_reads bowfiles
workflowWorkDir=\$(cd ../../ && pwd)
if [ -n "${mappingList}" ]; then
    $runGzip
    #rename files to standart format
    if [ "${mate}" == "pair" ]; then
        mv $file1 ${name}.1.fastq 2>/dev/null
        mv $file2 ${name}.2.fastq 2>/dev/null
        mv ${name}.1.fastq ${name}.2.fastq reads/.
    else
        mv $file1 ${name}.fastq 2>/dev/null
        mv ${name}.fastq reads/.
    fi
    #sequential mapping
    k=0
    prev="reads"
    IFS=',' read -r -a paramsListAr <<< "${paramsList}" #create comma separated array 
    IFS=',' read -r -a filtersListAr <<< "${filtersList}"
    IFS=',' read -r -a indexesListAr <<< "${indexesList}"
    IFS=',' read -r -a alignersListAr <<< "${alignersList}"
    IFS=',' read -r -a senseListAr <<< "${senseList}"
    wrkDir=\$(pwd)
    for rna_set in ${mappingList}
    do
        ((k++))
        printf -v k2 "%02d" "\$k" #turn into two digit format
        mkdir -p \${rna_set}/unmapped
        cd \$rna_set
        ## create link of the target file to prevent "too many symlinks error"
        for r in \${wrkDir}/\${prev}/*; do
            targetRead=\$(readlink -e \$r)
            rname=\$(basename \$r)
            echo "INFO: ln -s \$targetRead \$rname"
            ln -s \$targetRead \$rname
        done
        genomeDir=`dirname "\${indexesListAr[\$k-1]}"`
        echo "INFO: genomeDir: \$genomeDir"
        if [ -e "\${indexesListAr[\$k-1]}.1.bt2" -o  -e "\${indexesListAr[\$k-1]}.fa"  -o  -e "\${indexesListAr[\$k-1]}.fasta"  -o  -e "\$genomeDir/SAindex" ]; then
            if [ -e "\${indexesListAr[\$k-1]}.fa" ] ; then
                fasta=\${indexesListAr[\$k-1]}.fa
            elif [ -e "\${indexesListAr[\$k-1]}.fasta" ] ; then
                fasta=\${indexesListAr[\$k-1]}.fasta
            fi
            echo "INFO: fasta: \$fasta"
            if [ -e "\${indexesListAr[\$k-1]}.1.bt2" -a "\${alignersListAr[\$k-1]}" == "bowtie2" ] ; then
                echo "INFO: \${indexesListAr[\$k-1]}.1.bt2 Bowtie2 index found."
            elif [ -e "\${indexesListAr[\$k-1]}.1.ebwt" -a "\${alignersListAr[\$k-1]}" == "bowtie" ] ; then
                echo "INFO: \${indexesListAr[\$k-1]}.1.ebwt Bowtie index found."
            elif [ -e "\$genomeDir/SAindex" -a "\${alignersListAr[\$k-1]}" == "STAR" ] ; then
                echo "INFO: \$genomeDir/SAindex STAR index found."
            elif [ -e "\${indexesListAr[\$k-1]}.fa" -o  -e "\${indexesListAr[\$k-1]}.fasta" ] ; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2-build \$fasta \${indexesListAr[\$k-1]}
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    if [ -e "\${indexesListAr[\$k-1]}.gtf" ]; then
                        STAR --runMode genomeGenerate --genomeDir \$genomeDir --genomeFastaFiles \$fasta --sjdbGTFfile \${indexesListAr[\$k-1]}.gtf --genomeSAindexNbases 5
                    else
                        echo "WARNING: \${indexesListAr[\$k-1]}.gtf not found. STAR index is not generated."
                    fi
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie-build \$fasta \${indexesListAr[\$k-1]}
                fi
            fi
                
            if [ "${mate}" == "pair" ]; then
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${indexesListAr[\$k-1]} --no-unal --un-conc unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq --al-conc ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.1.fastq ${name}.2.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.1.fastq
                    mv ${name}.starUnmapped.out.mate2 unmapped/${name}.unmapped.2.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}   \${indexesListAr[\$k-1]}  --un  unmapped/${name}.unmapped.fastq -1 ${name}.1.fastq -2 ${name}.2.fastq -S  \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow1_\${rna_set}  
                    mv unmapped/${name}.unmapped_1.fastq unmapped/${name}.unmapped.1.fastq
                    mv unmapped/${name}.unmapped_2.fastq unmapped/${name}.unmapped.2.fastq
                fi
            else
                if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                    bowtie2 \${paramsListAr[\$k-1]} -x \${indexesListAr[\$k-1]} --no-unal --un  unmapped/${name}.unmapped.fastq -U ${name}.fastq --al ${name}.fq.mapped -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow_\${rna_set}  
                elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                    STAR \${paramsListAr[\$k-1]}  --genomeDir \$genomeDir --readFilesIn ${name}.fastq --outSAMtype SAM  --outFileNamePrefix ${name}.star --outReadsUnmapped Fastx
                    mv ${name}.starAligned.out.sam \${rna_set}_${name}_alignment.sam
                    mv ${name}.starUnmapped.out.mate1 unmapped/${name}.unmapped.fastq
                    mv ${name}.starLog.final.out \${k2}_${name}.star_\${rna_set}
                elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                    bowtie \${paramsListAr[\$k-1]}  \${indexesListAr[\$k-1]}  --un  unmapped/${name}.unmapped.fastq  ${name}.fastq  -S \${rna_set}_${name}_alignment.sam 2>&1 | tee \${k2}_${name}.bow1_\${rna_set}  
                    
                fi
            fi
            echo "INFO: samtools view -bT \${fasta} \${rna_set}_${name}_alignment.sam > \${rna_set}_${name}_alignment.bam"
            samtools view -bT \${fasta} \${rna_set}_${name}_alignment.sam > \${rna_set}_${name}_alignment.bam
            rm -f \${rna_set}_${name}_alignment.sam
            if [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                mv \${rna_set}_${name}_alignment.bam \${rna_set}_${name}_tmp0.bam
                echo "INFO: samtools view -F 0x04 -b \${rna_set}_${name}_tmp0.bam > \${rna_set}_${name}_alignment.bam"
                samtools view -F 0x04 -b \${rna_set}_${name}_tmp0.bam > \${rna_set}_${name}_alignment.bam  # Remove unmapped reads
                if [ "${mate}" == "pair" ]; then
                    echo "# unique mapped reads: \$(samtools view -f 0x40 -F 0x4 -q 255 \${rna_set}_${name}_alignment.bam | cut -f 1 | sort -T '.' | uniq | wc -l)" >> \${k2}_${name}.bow1_\${rna_set}
                else
                    echo "# unique mapped reads: \$(samtools view -F 0x40 -q 255 \${rna_set}_${name}_alignment.bam | cut -f 1 | sort -T '.' | uniq | wc -l)" >> \${k2}_${name}.bow1_\${rna_set}
                fi
            fi
            if [ "${mate}" == "pair" ]; then
                mv \${rna_set}_${name}_alignment.bam \${rna_set}_${name}_alignment.tmp1.bam
                echo "INFO: samtools sort -n -o \${rna_set}_${name}_alignment.tmp2 \${rna_set}_${name}_alignment.tmp1.bam"
                samtools sort -n -o \${rna_set}_${name}_alignment.tmp2.bam \${rna_set}_${name}_alignment.tmp1.bam 
                echo "INFO: samtools view -bf 0x02 \${rna_set}_${name}_alignment.tmp2.bam >\${rna_set}_${name}_alignment.bam"
                samtools view -bf 0x02 \${rna_set}_${name}_alignment.tmp2.bam >\${rna_set}_${name}_alignment.bam
                rm \${rna_set}_${name}_alignment.tmp1.bam \${rna_set}_${name}_alignment.tmp2.bam
            fi
            echo "INFO: samtools sort -o \${rna_set}@${name}_sorted.bam \${rna_set}_${name}_alignment.bam"
            samtools sort -o \${rna_set}@${name}_sorted.bam \${rna_set}_${name}_alignment.bam 
            echo "INFO: samtools index \${rna_set}@${name}_sorted.bam"
            samtools index \${rna_set}@${name}_sorted.bam
            
             # split sense and antisense bam files. 
            if [ "\${senseListAr[\$k-1]}" == "Yes" ]; then
                if [ "${mate}" == "pair" ]; then
                    echo "INFO: paired end sense antisense separation"
                	samtools view -f 65 -b \${rna_set}@${name}_sorted.bam >\${rna_set}@${name}_forward_sorted.bam
	                samtools index \${rna_set}@${name}_forward_sorted.bam
	                samtools view -F 16 -b \${rna_set}@${name}_forward_sorted.bam >\${rna_set}@${name}_sense_sorted.bam
	                samtools index \${rna_set}@${name}_sense_sorted.bam
	                samtools view -f 16 -b \${rna_set}@${name}_forward_sorted.bam >\${rna_set}@${name}_antisense_sorted.bam
	                samtools index \${rna_set}@${name}_antisense_sorted.bam
                else
	                echo "INFO: single end sense antisense separation"
	                samtools view -F 16 -b \${rna_set}@${name}_sorted.bam >\${rna_set}@${name}_sense_sorted.bam
	                samtools index \${rna_set}@${name}_sense_sorted.bam
	                samtools view -f 16 -b \${rna_set}@${name}_sorted.bam >\${rna_set}@${name}_antisense_sorted.bam
	                samtools index \${rna_set}@${name}_antisense_sorted.bam
                fi
            fi
            
            if [ "${remove_duplicates}" == "yes" ]; then
                ## check read header whether they have UMI tags which are separated with underscore.(eg. NS5HGY:2:11_GTATAACCTT)
                umiCheck=\$(samtools view \${rna_set}@${name}_sorted.bam |head -n 1 | awk 'BEGIN {FS="\\t"}; {print \$1}' | awk 'BEGIN {FS=":"}; \$NF ~ /_/ {print \$NF}')
                
                # based on remove_duplicates_based_on_UMI_after_mapping
                if [ "${remove_duplicates_based_on_UMI_after_mapping}" == "yes" -a ! -z "\$umiCheck" ]; then
                    echo "INFO: umi_mark_duplicates.py will be executed for removing duplicates from bam file"
                    echo "python umi_mark_duplicates.py -f \${rna_set}@${name}_sorted.bam -p 4"
                    python umi_mark_duplicates.py -f \${rna_set}@${name}_sorted.bam -p 4
                else
                    echo "INFO: Picard MarkDuplicates will be executed for removing duplicates from bam file"
                    if [ "${remove_duplicates_based_on_UMI_after_mapping}" == "yes"  ]; then
                        echo "WARNING: Read header have no UMI tags which are separated with underscore. Picard MarkDuplicates will be executed to remove duplicates from alignment file (bam) instead of remove_duplicates_based_on_UMI_after_mapping."
                    fi
                    echo "INFO: picard MarkDuplicates OUTPUT=\${rna_set}@${name}_sorted.deumi.sorted.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=\${rna_set}@${name}_sorted.bam"
                    picard MarkDuplicates OUTPUT=\${rna_set}@${name}_sorted.deumi.sorted.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=\${rna_set}@${name}_sorted.bam 
                fi
                #get duplicates stats (read the sam flags)
                samtools flagstat \${rna_set}@${name}_sorted.deumi.sorted.bam > \${k2}@\${rna_set}@${name}_duplicates_stats.log
                #remove alignments marked as duplicates
                samtools view -b -F 0x400 \${rna_set}@${name}_sorted.deumi.sorted.bam > \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup
                #sort deduplicated files by chrom pos
                echo "INFO: samtools sort -o \${rna_set}@${name}_sorted.dedup.bam \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup"
                samtools sort -o \${rna_set}@${name}_sorted.dedup.bam \${rna_set}@${name}_sorted.deumi.sorted.bam.x_dup 
                samtools index \${rna_set}@${name}_sorted.dedup.bam
                #get flagstat after dedup
                echo "##After Deduplication##" >> \${k2}@\${rna_set}@${name}_duplicates_stats.log
                samtools flagstat \${rna_set}@${name}_sorted.dedup.bam >> \${k2}@\${rna_set}@${name}_duplicates_stats.log
                # split sense and antisense bam files. 
	            if [ "\${senseListAr[\$k-1]}" == "Yes" ]; then
	                if [ "${mate}" == "pair" ]; then
	                    echo "INFO: paired end sense antisense separation"
	                	samtools view -f 65 -b \${rna_set}@${name}_sorted.dedup.bam >\${rna_set}@${name}_forward_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_forward_sorted.dedup.bam
		                samtools view -F 16 -b \${rna_set}@${name}_forward_sorted.dedup.bam>\${rna_set}@${name}_sense_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_sense_sorted.dedup.bam
		                samtools view -f 16 -b \${rna_set}@${name}_forward_sorted.dedup.bam >\${rna_set}@${name}_antisense_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_antisense_sorted.dedup.bam
	                else
		                echo "INFO: single end sense antisense separation"
		                samtools view -F 16 -b \${rna_set}@${name}_sorted.dedup.bam >\${rna_set}@${name}_sense_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_sense_sorted.dedup.bam
		                samtools view -f 16 -b \${rna_set}@${name}_sorted.dedup.bam >\${rna_set}@${name}_antisense_sorted.dedup.bam
		                samtools index \${rna_set}@${name}_antisense_sorted.dedup.bam
	                fi
	            fi
            fi
            
        
            for file in unmapped/*; do mv \$file \${file/.unmapped/}; done ##remove .unmapped from filename
            if [ "\${alignersListAr[\$k-1]}" == "bowtie2" ]; then
                grep -v Warning \${k2}_${name}.bow_\${rna_set} > ${name}.tmp
                mv ${name}.tmp \${k2}_${name}.bow_\${rna_set}
                cp \${k2}_${name}.bow_\${rna_set} ./../bowfiles/.
            elif [ "\${alignersListAr[\$k-1]}" == "bowtie" ]; then
                cp \${k2}_${name}.bow1_\${rna_set} ./../bowfiles/.
            elif [ "\${alignersListAr[\$k-1]}" == "STAR" ]; then
                cp \${k2}_${name}.star_\${rna_set} ./../bowfiles/.
            fi
            cd ..
            # if filter is on, remove previously created unmapped fastq. 
            if [ "\${filtersListAr[\$k-1]}" == "Yes" ]; then
                if [ "\${prev}" != "reads" ]; then
                    echo "INFO: remove prev: \${prev}/*"
                    rm -rf \${prev}/*
                elif  [ "${remove_previous_reads}" == "true" ]; then
                    echo "INFO: inputs reads will be removed if they are located in the workdir"
                    for f in \${prev}/*; do
                        targetFile=\$(readlink -e \$f)
                        echo "INFO: targetFile: \$targetFile"
                        if [[ \$targetFile == *"\${workflowWorkDir}"* ]]; then
                            rm -f \$targetFile
                            echo "INFO: \$targetFile located in workdir and deleted."
                        fi
                    done
                fi
            # if filter is off remove current unmapped fastq
            else
                echo "INFO: remove \${rna_set}/unmapped/*"
                rm -rf \${rna_set}/unmapped/*
            fi
        else
            echo "WARNING: \${indexesListAr[\$k-1]} Mapping skipped. File not found."
            cd unmapped 
            ln -s \${wrkDir}/\${rna_set}/*fastq .
            cd ..
            cd ..
        fi
        
        if [ "\${filtersListAr[\$k-1]}" == "Yes" ]; then
            prev=\${rna_set}/unmapped
        fi
    done
    cd final_reads && ln -s \${wrkDir}/\${prev}/* .
else 
    mv ${reads} final_reads/.
fi
"""

}
}


mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process Sequential_Variant__Mapping_Module_Sequential_Mapping_Bam_dedup_count {

input:
 file bam from g24_38_bam_file_g24_36.collect()
 file index from g24_38_bam_index_g24_36.collect()

output:
 file "*.counts.tsv"  into g24_36_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my @header_antisense;
my @header_sense;
my %all_files;
my %sense_files;
my %antisense_files;

my @mappingList = (!{mappingListQuoteSep});
my @rawIndexList = (!{rawIndexList});
my %indexHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;

chomp(my $contents = `ls *.bam`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        $file=~/(.*)@(.*)_sorted(.*)\\.bam/;
        my $mapper = $1; 
        my $name = $2; ##header
        #print $3;
        if ($3 eq ".dedup"){
            $dedup = ".dedup";
        }
        if ($name=~/_antisense$/){
        	push(@header_antisense, $name) unless grep{$_ eq $name} @header_antisense; #mapped element header
        	$antisense_files{$mapper} .= $file." ";
        }
        elsif ($name=~/_sense$/){
        	push(@header_sense, $name) unless grep{$_ eq $name} @header_sense; #mapped element header
        	$sense_files{$mapper} .= $file." ";
        }
        else{
			push(@header, $name) unless grep{$_ eq $name} @header; #mapped element header
	        $all_files{$mapper} .= $file." ";
        }
}

runCov(\\%all_files, \\@header, \\@indexHash, "", $dedup);
runCov(\\%sense_files, \\@header_sense, \\@indexHash, "sense", $dedup);
runCov(\\%antisense_files, \\@header_antisense, \\@indexHash, "antisense", $dedup);

sub runCov {
	my ( \$files, \$header, \$indexHash, \$sense_antisense, \$dedup) = @_;
	open OUT, ">header".\$sense_antisense.".tsv";
	print OUT join ("\\t", "id","len",@{\$header}),"\\n";
	close OUT;
	my $par = "";
	if ($sense_antisense=~/^sense\$/){
      $par = "-s";
    }elsif($sense_antisense=~/^antisense\$/){
      $par = "-S";
    }
	
	foreach my $key (sort keys %{\$files}) {  
	   my $bamFiles = ${\$files}{$key};
		unless (-e ${indexHash}{$key}.".bed") {
            print "2: bed not found run makeBed\n";
                if (-e ${indexHash}{$key}.".fa") {
                    makeBed(${indexHash}{$key}.".fa", $key, ${indexHash}{$key}.".bed");
                } elsif(-e ${indexHash}{$key}.".fasta"){
                    makeBed(${indexHash}{$key}.".fasta", $key, ${indexHash}{$key}.".bed");
                }
        }
	    
		my $com =  "bedtools multicov $par -bams $bamFiles -bed ".${indexHash}{$key}.".bed >$key${dedup}${sense_antisense}.counts.tmp\n";
        print $com;
        `$com`;
        my $iniResColumn = int(countColumn(${indexHash}{$key}.".bed")) + 1;
	    `awk -F \\"\\\\t\\" \\'{a=\\"\\";for (i=$iniResColumn;i<=NF;i++){a=a\\"\\\\t\\"\\$i;} print \\$4\\"\\\\t\\"(\\$3-\\$2)\\"\\"a}\\' $key${dedup}${sense_antisense}.counts.tmp> $key${dedup}${sense_antisense}.counts.tsv`;
	    `sort -k3,3nr $key${dedup}${sense_antisense}.counts.tsv>$key${dedup}${sense_antisense}.sorted.tsv`;
        `cat header${sense_antisense}.tsv $key${dedup}${sense_antisense}.sorted.tsv> $key${dedup}${sense_antisense}.counts.tsv`;
	}
}

sub countColumn {
    my ( \$file) = @_;
    open(IN, \$file);
    my $line=<IN>;
    chomp($line);
    my @cols = split('\\t', $line);
    my $n = @cols;
    close OUT;
    return $n;
}

sub makeBed {
    my ( \$fasta, \$type, \$bed) = @_;
    print "makeBed $fasta\\n";
    print "makeBed $bed\\n";
    open OUT, ">$bed";
    open(IN, \$fasta);
    my $name="";
    my $seq="";
    my $i=0;
    while(my $line=<IN>){
        chomp($line);
        if($line=~/^>(.*)/){
            $i++ if (length($seq)>0);
            print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
            $name="$1";
            $seq="";
        } elsif($line=~/[ACGTNacgtn]+/){
            $seq.=$line;
        }
    }
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''
}


process Sequential_Variant__Mapping_Module_Deduplication_Summary {

input:
 file flagstat from g24_38_log_file_g24_30.collect()
 val mate from g_19_mate_g24_30

output:
 file "deduplication_summary.tsv"  into g24_30_outputFileTSV

errorStrategy 'retry'
maxRetries 2

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i=0;
chomp(my $contents = `ls *_duplicates_stats.log`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $i++;
    $file=~/(.*)@(.*)@(.*)_duplicates_stats\\.log/;
    my $mapOrder = int($1); 
    my $mapper = $2; #mapped element 
    my $name = $3; ##sample name
    push(@header, $mapper) unless grep{$_ eq $mapper} @header; 
        
    # my $duplicates;
    my $aligned;
    my $dedup; #aligned reads after dedup
    my $percent=0;
    if ("!{mate}" eq "pair" ){
        #first flagstat belongs to first bam file
        chomp($aligned = `cat $file | grep 'properly paired (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        #second flagstat belongs to dedup bam file
        chomp($dedup = `cat $file | grep 'properly paired (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    } else {
        chomp($aligned = `cat $file | grep 'mapped (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        chomp($dedup = `cat $file | grep 'mapped (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    }
    # chomp($duplicates = `cat $file | grep 'duplicates' | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    # $dedup = int($aligned) - int($duplicates);
    if ("!{mate}" eq "pair" ){
       $dedup = int($dedup/2);
       $aligned = int($aligned/2);
    } 
    $percent = "0.00";
    if (int($aligned)  > 0 ){
       $percent = sprintf("%.2f", ($aligned-$dedup)/$aligned*100); 
    } 
    $tsv{$name}{$mapper}=[$aligned,$dedup,"$percent%"];
    $headerHash{$mapOrder}=$mapper;
    $headerText{$mapOrder}=["$mapper (Before Dedup)", "$mapper (After Dedup)", "$mapper (Duplication Ratio %)"];
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary = "deduplication_summary.tsv";
open(OUT, ">$summary");
print OUT "Sample\\t";
my @headArr = ();
for my $mapOrder (@sortedOrderArray) {
    push (@headArr, @{$headerText{$mapOrder}});
}
my $headArrAll = join("\\t", @headArr);
print OUT "$headArrAll\\n";

foreach my $name (keys %tsv){
    my @rowArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push (@rowArr, @{$tsv{$name}{$headerHash{$mapOrder}}});
    }
    my $rowArrAll = join("\\t", @rowArr);
    print OUT "$name\\t$rowArrAll\\n";
}
close(OUT);
'''
}

mappingListQuoteSep = mapList.collect{ '"' + it + '"'}.join(",") 
rawIndexList = indexList.collect{ '"' + it + '"'}.join(",") 
process Sequential_Variant__Mapping_Module_Sequential_Mapping_Bam_count {

input:
 file bam from g24_38_bam_file_g24_35.collect()
 file index from g24_38_bam_index_g24_35.collect()

output:
 file "*.counts.tsv"  into g24_35_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my @header_antisense;
my @header_sense;
my %all_files;
my %sense_files;
my %antisense_files;

my @mappingList = (!{mappingListQuoteSep});
my @rawIndexList = (!{rawIndexList});
my %indexHash;
my $dedup = "";
@indexHash{@mappingList} = @rawIndexList;

chomp(my $contents = `ls *.bam`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        $file=~/(.*)@(.*)_sorted(.*)\\.bam/;
        my $mapper = $1; 
        my $name = $2; ##header
        #print $3;
        if ($3 eq ".dedup"){
            $dedup = ".dedup";
        }
        if ($name=~/_antisense$/){
        	push(@header_antisense, $name) unless grep{$_ eq $name} @header_antisense; #mapped element header
        	$antisense_files{$mapper} .= $file." ";
        }
        elsif ($name=~/_sense$/){
        	push(@header_sense, $name) unless grep{$_ eq $name} @header_sense; #mapped element header
        	$sense_files{$mapper} .= $file." ";
        }
        else{
			push(@header, $name) unless grep{$_ eq $name} @header; #mapped element header
	        $all_files{$mapper} .= $file." ";
        }
}

runCov(\\%all_files, \\@header, \\@indexHash, "", $dedup);
runCov(\\%sense_files, \\@header_sense, \\@indexHash, "sense", $dedup);
runCov(\\%antisense_files, \\@header_antisense, \\@indexHash, "antisense", $dedup);

sub runCov {
	my ( \$files, \$header, \$indexHash, \$sense_antisense, \$dedup) = @_;
	open OUT, ">header".\$sense_antisense.".tsv";
	print OUT join ("\\t", "id","len",@{\$header}),"\\n";
	close OUT;
	my $par = "";
	if ($sense_antisense=~/^sense\$/){
      $par = "-s";
    }elsif($sense_antisense=~/^antisense\$/){
      $par = "-S";
    }
	
	foreach my $key (sort keys %{\$files}) {  
	   my $bamFiles = ${\$files}{$key};
		unless (-e ${indexHash}{$key}.".bed") {
            print "2: bed not found run makeBed\n";
                if (-e ${indexHash}{$key}.".fa") {
                    makeBed(${indexHash}{$key}.".fa", $key, ${indexHash}{$key}.".bed");
                } elsif(-e ${indexHash}{$key}.".fasta"){
                    makeBed(${indexHash}{$key}.".fasta", $key, ${indexHash}{$key}.".bed");
                }
        }
	    
		my $com =  "bedtools multicov $par -bams $bamFiles -bed ".${indexHash}{$key}.".bed >$key${dedup}${sense_antisense}.counts.tmp\n";
        print $com;
        `$com`;
        my $iniResColumn = int(countColumn(${indexHash}{$key}.".bed")) + 1;
	    `awk -F \\"\\\\t\\" \\'{a=\\"\\";for (i=$iniResColumn;i<=NF;i++){a=a\\"\\\\t\\"\\$i;} print \\$4\\"\\\\t\\"(\\$3-\\$2)\\"\\"a}\\' $key${dedup}${sense_antisense}.counts.tmp> $key${dedup}${sense_antisense}.counts.tsv`;
	    `sort -k3,3nr $key${dedup}${sense_antisense}.counts.tsv>$key${dedup}${sense_antisense}.sorted.tsv`;
        `cat header${sense_antisense}.tsv $key${dedup}${sense_antisense}.sorted.tsv> $key${dedup}${sense_antisense}.counts.tsv`;
	}
}

sub countColumn {
    my ( \$file) = @_;
    open(IN, \$file);
    my $line=<IN>;
    chomp($line);
    my @cols = split('\\t', $line);
    my $n = @cols;
    close OUT;
    return $n;
}

sub makeBed {
    my ( \$fasta, \$type, \$bed) = @_;
    print "makeBed $fasta\\n";
    print "makeBed $bed\\n";
    open OUT, ">$bed";
    open(IN, \$fasta);
    my $name="";
    my $seq="";
    my $i=0;
    while(my $line=<IN>){
        chomp($line);
        if($line=~/^>(.*)/){
            $i++ if (length($seq)>0);
            print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
            $name="$1";
            $seq="";
        } elsif($line=~/[ACGTNacgtn]+/){
            $seq.=$line;
        }
    }
    print OUT "$name\\t1\\t".length($seq)."\\t$name\\t0\\t+\\n" if (length($seq)>0); 
    close OUT;
}

'''
}


process Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary {

input:
 file logfile from g15_20_log_file_g15_16.collect()
 val mate from g_3_mate_g15_16

output:
 file "quality_filter_summary.tsv"  into g15_16_outputFileTSV_g_16

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i = 0;
chomp( my $contents = `ls *.log` );
my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapper   = "";
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_quality\\.log/){
        $mapper   = "fastx";
        $file =~ /(.*)\\.fastx_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    } elsif ($file =~ /(.*)\\.trimmomatic_quality\\.log/){
        $mapper   = "trimmomatic";
        $file =~ /(.*)\\.trimmomatic_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "quality_filter_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}


process Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary {

input:
 file logfile from g15_18_log_file_g15_11.collect()
 val mate from g_3_mate_g15_11

output:
 file "adapter_removal_summary.tsv"  into g15_11_outputFileTSV_g_16
 file "adapter_removal_detailed_summary.tsv" optional true  into g15_11_outputFile

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx\\.log/){
        $file =~ /(.*)\\.fastx\\.log/;
        my $mapper   = "fastx";
        my $name = $1;    ##sample name
        push( @header, $mapper );

        my $in;
        my $out;
        my $tooshort;
        my $adapteronly;
        my $noncliped;
        my $Nreads;

        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $tooshort =`cat $file | grep 'too-short reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $adapteronly =`cat $file | grep 'adapter-only reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $noncliped =`cat $file | grep 'non-clipped reads.' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $Nreads =`cat $file | grep 'N reads.' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        $tsvDetail{$name}{$mapper} = [ $in, $tooshort, $adapteronly, $noncliped, $Nreads, $out ];
        $headerTextDetail{$mapOrder} = ["Total Reads","Too-short reads","Adapter-only reads","Non-clipped reads","N reads","Reads After Adapter Removal"];
    } elsif ($file =~ /(.*)\\.trimmomatic\\.log/){
        $file =~ /(.*)\\.trimmomatic\\.log/;
        my $mapper   = "trimmomatic";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        
        my $in;
        my $out;

        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        


        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "adapter_removal_summary.tsv";
my $detailed_summary = "adapter_removal_detailed_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );
if (%headerTextDetail){
    writeFile( $detailed_summary, \\%headerTextDetail, \\%tsvDetail );  
}

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

g15_11_outputFileTSV_g_16= g15_11_outputFileTSV_g_16.ifEmpty([""]) 
g15_21_outputFileTSV_g_16= g15_21_outputFileTSV_g_16.ifEmpty([""]) 
g15_16_outputFileTSV_g_16= g15_16_outputFileTSV_g_16.ifEmpty([""]) 

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 30
    $CPU  = 1
    $MEMORY = 10
    $QUEUE = "short"
}
//* platform
//* autofill

process Overall_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /overall_summary.tsv$/) "summary/$filename"
}

input:
 file adapterSum from g15_11_outputFileTSV_g_16
 file trimmerSum from g15_21_outputFileTSV_g_16
 file qualitySum from g15_16_outputFileTSV_g_16

output:
 file "overall_summary.tsv"  into g_16_outputFileTSV

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @rawFiles = split(/[\\n]+/, $contents);
my @files = ();
# order must be in this order for chipseq pipeline: bowtie->dedup
# rsem bam pipeline: dedup->rsem, star->dedup
# riboseq ncRNA_removal->star
my @order = ("adapter_removal","trimmer","quality","extractUMI","extractValid","sequential_mapping","ncRNA_removal","bowtie","star","hisat2","tophat2", "dedup","rsem","kallisto","esat","count");
for ( my $k = 0 ; $k <= $#order ; $k++ ) {
    for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
        if ( $rawFiles[$i] =~ /$order[$k]/ ) {
            push @files, $rawFiles[$i];
        }
    }
}

print Dumper \\@files;
##add rest of the files
for ( my $i = 0 ; $i <= $#rawFiles ; $i++ ) {
    push(@files, $rawFiles[$i]) unless grep{$_ == $rawFiles[$i]} @files;
}
print Dumper \\@files;

##Merge each file according to array order

foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @header) = ( split("\\t", $line1) );
        push @seen_cols, @header;

        while (my $line=<IN>) {
        chomp($line);
        my ( $ID, @fields ) = ( split("\\t", $line) ); 
        my %this_row;
        @this_row{@header} = @fields;

        #print Dumper \\%this_row;

        foreach my $column (@header) {
            if (! exists $all_rows{$ID}{$column}) {
                $all_rows{$ID}{$column} = $this_row{$column}; 
            }
        }   
    }
    close IN;
}

#print for debugging
#print Dumper \\%all_rows;
#print Dumper \\%seen_cols;

#grab list of column headings we've seen, and order them. 
my @cols_to_print = uniq(@seen_cols);
my $summary = "overall_summary.tsv";
open OUT, ">$summary";
print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
foreach my $key ( keys %all_rows ) { 
    #map iterates all the columns, and gives the value or an empty string. if it's undefined. (prevents errors)
    print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
}
close OUT;

sub uniq {
    my %seen;
    grep ! $seen{$_}++, @_;
}

'''


}

//* autofill
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 8
    $QUEUE = "short"
}
//* platform
//* autofill

process Sequential_Variant__Mapping_Module_Sequential_Mapping_Summary {

input:
 set val(name), file(bowfile) from g24_38_bowfiles_g24_26
 val mate from g_19_mate_g24_26
 val filtersList from g24_38_filter_g24_26

output:
 file '*.tsv'  into g24_26_outputFileTSV_g24_13
 val "sequential_mapping_sum"  into g24_26_name_g24_13

errorStrategy 'retry'
maxRetries 2

shell:
'''
#!/usr/bin/env perl
open(my \$fh, '>', "!{name}.tsv");
print $fh "Sample\\tGroup\\tTotal Reads\\tReads After Sequential Mapping\\tUniquely Mapped\\tMultimapped\\tMapped\\n";
my @bowArray = split(' ', "!{bowfile}");
my $group= "\\t";
my @filterArray = (!{filtersList});
foreach my $bowitem(@bowArray) {
    # get mapping id
    my @bowAr = $bowitem.split("_");
    $bowCount = $bowAr[0] + -1;
    # if bowfiles ends with underscore (eg. bow_rRNA), parse rRNA as a group.
    my ($RDS_In, $RDS_After, $RDS_Uniq, $RDS_Multi, $ALGN_T, $a, $b, $aPer, $bPer)=(0, 0, 0, 0, 0, 0, 0, 0, 0);
    if ($bowitem =~ m/bow_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN, $bowitem);
        my $i = 0;
        while(my $line=<IN>){
            chomp($line);
            $line=~s/^ +//;
            my @arr=split(/ /, $line);
            $RDS_In=$arr[0] if ($i=~/^1$/);
            # Reads After Filtering column depends on filtering type
            if ($i == 2){
                if ($filterArray[$bowCount] eq "Yes"){
                    $RDS_After=$arr[0];
                } else {
                    $RDS_After=$RDS_In;
                }
            }
            if ($i == 3){
                $a=$arr[0];
                $aPer=$arr[1];
                $aPer=~ s/([()])//g;
                $RDS_Uniq=$arr[0];
            }
            if ($i == 4){
                $b=$arr[0];
                $bPer=$arr[1];
                $bPer=~ s/([()])//g;
                $RDS_Multi=$arr[0];
            }
            $ALGN_T=($a+$b);
            $i++;
        }
        close(IN);
    } elsif ($bowitem =~ m/star_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN2, $bowitem);
        my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $bowitem | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($uniqAligned = `cat $bowitem | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $bowitem | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		## Here we exclude "Number of reads mapped to too many loci" from multimapped reads since in bam file it called as unmapped.
		## Besides, these "too many loci" reads exported as unmapped reads from STAR.
		$RDS_In = int($inputCount);
		$RDS_Multi = int($multimapped);
        $RDS_Uniq = int($uniqAligned);
        $ALGN_T = $RDS_Uniq+$RDS_Multi;
		if ($filterArray[$bowCount] eq "Yes"){
            $RDS_After=$RDS_In-$ALGN_T;
        } else {
            $RDS_After=$RDS_In;
        }
    } elsif ($bowitem =~ m/bow1_([^\\.]+)$/){
        $group = "$1\\t";
        open(IN2, $bowitem);
        my $multimapped;
		my $aligned;
		my $inputCount;
		my $uniqAligned;
		chomp($inputCount = `cat $bowitem | grep '# reads processed:' | awk '{sum+=\\$4} END {print sum}'`);
		chomp($aligned = `cat $bowitem | grep '# reads with at least one reported alignment:' | awk '{sum+=\\$9} END {print sum}'`);
		chomp($uniqAligned = `cat $bowitem | grep '# unique mapped reads:' | awk '{sum+=\\$5} END {print sum}'`);
		## Here we exclude "Number of reads mapped to too many loci" from multimapped reads since in bam file it called as unmapped.
		## Besides, these "too many loci" reads exported as unmapped reads from STAR.
		$RDS_In = int($inputCount);
		$RDS_Multi = int($aligned) -int($uniqAligned);
		if ($RDS_Multi < 0 ){
		    $RDS_Multi = 0;
		}
        $RDS_Uniq = int($uniqAligned);
        $ALGN_T = int($aligned);
		if ($filterArray[$bowCount] eq "Yes"){
            $RDS_After=$RDS_In-$ALGN_T;
        } else {
            $RDS_After=$RDS_In;
        }
    }
    
    print $fh "!{name}\\t$group$RDS_In\\t$RDS_After\\t$RDS_Uniq\\t$RDS_Multi\\t$ALGN_T\\n";
}
close($fh);



'''

}


process Sequential_Variant__Mapping_Module_Merge_TSV_Files {

input:
 file tsv from g24_26_outputFileTSV_g24_13.collect()
 val outputFileName from g24_26_name_g24_13.collect()

output:
 file "${name}.tsv"  into g24_13_outputFileTSV_g24_14

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


process Sequential_Variant__Mapping_Module_Sequential_Mapping_Short_Summary {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /sequential_mapping_short_sum.tsv$/) "short_counts/$filename"
	else if (filename =~ /sequential_mapping_detailed_sum.tsv$/) "detailed_counts/$filename"
}

input:
 file mainSum from g24_13_outputFileTSV_g24_14

output:
 file "sequential_mapping_short_sum.tsv"  into g24_14_outputFileTSV
 file "sequential_mapping_detailed_sum.tsv"  into g24_14_outputFile

errorStrategy 'retry'
maxRetries 2

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_rows;
my @seen_cols_short;
my @seen_cols_detailed;
my $ID_header;

chomp(my $contents = `ls *.tsv`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
        open IN,"$file";
        my $line1 = <IN>;
        chomp($line1);
        ( $ID_header, my @h) = ( split("\\t", $line1) );
        my $totalHeader = $h[1];
        my $afterFilteringHeader = $h[2];
        my $uniqueHeader = $h[3];
        my $multiHeader = $h[4];
        my $mappedHeader = $h[5];
        push(@seen_cols_short, $totalHeader) unless grep{$_ eq $totalHeader} @seen_cols_short; #Total reads Header
        push(@seen_cols_detailed, $totalHeader) unless grep{$_ eq $totalHeader} @seen_cols_detailed; #Total reads Header

        my $n=0;
        while (my $line=<IN>) {
                
                chomp($line);
                my ( $ID, @fields ) = ( split("\\t", $line) ); 
                #SHORT
                push(@seen_cols_short, $fields[0]) unless grep{$_ eq $fields[0]} @seen_cols_short; #mapped element header
                $all_rows{$ID}{$fields[0]} = $fields[5];#Mapped Reads
                #Grep first line $fields[1] as total reads.
                if (!exists $all_rows{$ID}{$totalHeader}){    
                        $all_rows{$ID}{$totalHeader} = $fields[1];
                } 
                $all_rows{$ID}{$afterFilteringHeader} = $fields[2]; #only use last entry
                #DETAILED
                $uniqueHeadEach = "$fields[0] (${uniqueHeader})";
                $multiHeadEach = "$fields[0] (${multiHeader})";
                $mappedHeadEach = "$fields[0] (${mappedHeader})";
                push(@seen_cols_detailed, $mappedHeadEach) unless grep{$_ eq $mappedHeadEach} @seen_cols_detailed;
                push(@seen_cols_detailed, $uniqueHeadEach) unless grep{$_ eq $uniqueHeadEach} @seen_cols_detailed;
                push(@seen_cols_detailed, $multiHeadEach) unless grep{$_ eq $multiHeadEach} @seen_cols_detailed;
                $all_rows{$ID}{$mappedHeadEach} = $fields[5];
                $all_rows{$ID}{$uniqueHeadEach} = $fields[3];
                $all_rows{$ID}{$multiHeadEach} = $fields[4];
    }
    close IN;
    push(@seen_cols_short, $afterFilteringHeader) unless grep{$_ eq $afterFilteringHeader} @seen_cols_short; #After filtering Header
}


#print Dumper \\%all_rows;
#print Dumper \\%seen_cols_short;

printFiles("sequential_mapping_short_sum.tsv",@seen_cols_short,);
printFiles("sequential_mapping_detailed_sum.tsv",@seen_cols_detailed);


sub printFiles {
    my($summary, @cols_to_print) = @_;
    
    open OUT, ">$summary";
    print OUT join ("\\t", $ID_header,@cols_to_print),"\\n";
    foreach my $key ( keys %all_rows ) { 
        print OUT join ("\\t", $key, (map { $all_rows{$key}{$_} // '' } @cols_to_print)),"\\n";
        }
        close OUT;
}

'''


}

//* params.run_FastQC =  "no"  //* @dropdown @options:"yes","no" @description:"FastQC provides quality control checks on raw sequence data."



process Adapter_Trimmer_Quality_Module_FastQC {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.(html|zip)$/) "fastQC/$filename"
}

input:
 val mate from g_3_mate_g15_3
 set val(name), file(reads) from g_2_reads_g15_3

output:
 file '*.{html,zip}'  into g15_3_FastQCout

errorStrategy 'retry'
maxRetries 3

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
${runGzip}
fastqc ${file} 
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
