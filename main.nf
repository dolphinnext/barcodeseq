$HOSTNAME = ""
params.outdir = 'results'  

//* params.count_reverse_complement =  "no"   //* @dropdown @options:"yes","no" @show_settings:"Count_Variants"

if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.mergedmate){params.mergedmate = ""} 
if (!params.variants_file){params.variants_file = ""} 

Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.into{g_2_reads_g15_3;g_2_reads_g15_18}

Channel.value(params.mate).into{g_3_mate_g_13;g_3_mate_g15_3;g_3_mate_g15_11;g_3_mate_g15_16;g_3_mate_g15_18;g_3_mate_g15_19;g_3_mate_g15_20;g_3_mate_g15_21}
Channel.value(params.mergedmate).set{g_19_mate_g_45}
Channel.value(params.variants_file).set{g_35_filePath_g_45}

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
g_13_reads_g_14.into{g_14_reads_g_45}
} else {

process Extract_Barcode {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /ext\/.*.fastq$/) "aavbarcodes/$filename"
}

input:
 set val(name),file(reads) from g_13_reads_g_14

output:
 set val(name),file("ext/*.fastq")  into g_14_reads_g_45

errorStrategy 'retry'
maxRetries 3

when:
(params.run_Extract_Barcode && (params.run_Extract_Barcode == "yes"))

shell:
count_reverse_complement = params.count_reverse_complement
five_prime = params.Extract_Barcode.five_prime
three_prime = params.Extract_Barcode.three_prime
min_len = params.Extract_Barcode.min_len
max_len = params.Extract_Barcode.max_len
mismatch = params.Extract_Barcode.mismatch

'''
#!/usr/bin/env perl

my $count_reverse_complement = "!{count_reverse_complement}";
my $five_prime = "!{five_prime}";
my $three_prime = "!{three_prime}";
my $min_len = "!{min_len}";
my $max_len = "!{max_len}";
my $mismatch = "!{mismatch}";

runCmd("mkdir -p ext");

if ($count_reverse_complement eq "no"){
	runCmd("cutadapt -a $three_prime -n $mismatch -o ext/!{name}.tmp !{reads}");
	runCmd("cutadapt -g $five_prime -m $min_len -M $max_len -n $mismatch -o ext/!{name}.fastq ext/!{name}.tmp");
}
else {
	my $rev_threeprime=$three_prime;
	my $rev_fiveprime=$five_prime;
	$rev_threeprime=~tr/ACGT/TGCA/;
	$rev_fiveprime=~tr/ACGT/TGCA/;
	$rev_threeprime = reverse($rev_threeprime);
	$rev_fiveprime = reverse($rev_fiveprime);
	runCmd("cutadapt -a $three_prime -n $mismatch -o ext/!{name}.tmp.1 !{reads}");
	runCmd("cutadapt -g $five_prime -m $min_len -M $max_len -n $mismatch -o ext/!{name}.tmp.2 ext/!{name}.tmp.1");
	runCmd("cutadapt -a $rev_fiveprime -n $mismatch -o ext/!{name}.tmp.3 !{reads}");
	runCmd("cutadapt -g $rev_threeprime -m $min_len -M $max_len -n $mismatch -o ext/!{name}.tmp.4 ext/!{name}.tmp.3");
	runCmd("cat ext/!{name}.tmp.2 ext/!{name}.tmp.4 > ext/!{name}.fastq");
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


//* params.variants_file =  "/project/umw_biocore/kucukura/data/akuous/variants/variants.txt"  //* @input  @description:"Variant sequences Format: VariantSeq\tVariantname" @show_settings:"Count_Variants" 


process count_variants {

input:
 set val(name),file(reads) from g_14_reads_g_45
 val mate from g_19_mate_g_45
 val variants_file from g_35_filePath_g_45

output:
 set "*.tsv"  into g_45_outputFileTSV_g_44
 file "*.fa"  into g_45_fastaFile_g_50

errorStrategy 'retry'
maxRetries 2

shell:
count_reverse_complement = params.count_reverse_complement

'''
#!/usr/bin/env perl
use POSIX;

open(VARIANTS, "!{variants_file}");
open(READS, "!{reads}");
my $count_revcomp = "!{count_reverse_complement}";
my $name = "!{name}";
my $count_revcomp = "no";
my %lib=();
my %libtotal=();
my %variants=();
my %variantCounts=();
my %variantNames=();
my %variantlib=();
my %unknownvariants=();
my %unknownvariantCounts=();
my $variantlen=0;
my $i=0;
while($line=<VARIANTS>)
{
  chomp($line);
  if ($line=~/([ACGT]*)[\\t\\s](.*)/) {
    if (length($2)>0 && length($1)>0) {
        $variants{$1}=0;
        $variantlib{$1} = $2;
        my @libarr=split(/_/, $2);
        $lib{$libarr[0]}=0;
        $libtotal{$libarr[0]}=0;
        $variantCounts{$libarr[0]}=();
        $variantNames{$libarr[0]}=();
    }
  }
}
while($line=<READS>)
{
  chomp($line);
  $i++;
  if ($i % 100000 == 0){
     print "$i reads processed!\\n"; 
  }
  if ($line=~/^([ACGT]*)$/) {
  	if ($count_revcomp eq "yes"){
  		my $revcom = $line;
        $revcom=~tr/ACGT/TGCA/;
        $revcom = reverse($revcom);
        $line=$revcom if (exists($variants{$revcom}));
    }
    if (exists($variants{$line})) {
        $variants{$line}++;
    }
    else {
         if (exists($unknownvariants{$line})) {
             $unknownvariants{$line}++;
         }else{
             $unknownvariants{$line}=1;
         }
      }
   }
}
my $totalunknown=0;
my $maxunknown=0;
foreach my $key (keys %unknownvariants)
{
   push(@unknownvariantCounts, $unknownvariants{$key});
   if ($maxunknown < $unknownvariants{$key}) { $maxunknown =  $unknownvariants{$key} };
   $totalunknown += $unknownvariants{$key};
}

open(OUTKNOWN, ">", $name."_known.tsv");
print OUTKNOWN "lib\\tseq\\tcount\\n";
$i = 0;
foreach my $key (keys %variants)
{
   if ($variants{$key}>0){
      $variantlen = length($key);
      print OUTKNOWN $variantlib{$key}."\\t".$key."\\t".$variants{$key}."\\n";
      my @libarr=split(/_/,$variantlib{$key});
      $lib{$libarr[0]}++;
      $libtotal{$libarr[0]}+=$variants{$key};
      push(@{$variantCounts{$libarr[0]}}, $variants{$key});
      push(@{$variantNames{$libarr[0]}}, $variantlib{$key});
   }
   $i++;
}
close(OUTKNOWN);

open(OUTKNOWNSUM, ">", $name."_knownsum.tsv");
print OUTKNOWNSUM "lib\\tcount\\n";

my $identified =  0;
my $totalreads = 0;
my $avg = 0;
my $std = 0;
my $cv = 0;
my $variant = "";
my $highestvariant="";
my $lowestvariant="";

foreach my $key (keys %lib)
{
   print OUTKNOWNSUM $key."\\t".$lib{$key}."\\n";
   if ($identified < $lib{$key}){
      $identified = $lib{$key};
      $totalreads = $libtotal{$key};
      if ($identified > 0) {
         $variant = $key;
         $avg =  average(\\@{$variantCounts{$key}}); 
         $std =  stdev(\\@{$variantCounts{$key}});
	 $cv = $std*100/$avg;
         ($highestvariant, $lowestvariant) = highlow(\\@{$variantCounts{$key}}, \\@{$variantNames{$key}});
      }
   }
}
close(OUTKNOWNSUM);

open(OUTKNOWNSTAT, ">", $name."_knownstats.tsv");
print OUTKNOWNSTAT "Sample\\t$variant\\nidentified\\t$identified\\ntotalreads\\t$totalreads\\naverage\\t".round($avg,2)."\\nstd\\t".round($std,2)."\\n%CV\\t".round($cv,2)."\\nHighest_Read(variant)\\t$highestvariant\\nLowest_Read(variant)\\t$lowestvariant\\nsize(nt)\\t$variantlen\\n";
close(OUTKNOWNSTAT);

open(OUTUNKNOWN, ">", $name."_unknown.tsv");
open(OUTUNKNOWNFA, ">", $name."_unknown.fa");
print OUTUNKNOWN "seq\\tcount\\n";

foreach my $key (keys %unknownvariants)
{
   print OUTUNKNOWN $key."\\t".$unknownvariants{$key}."\\n";
   print OUTUNKNOWNFA ">".$key.":".$unknownvariants{$key}."\\n".$key."\\n";
}
close(OUTUNKNOWN);
close(OUTUNKNOWNFA);
open(OUTUNKNOWNSTAT, ">", $name."_unknownstats.tsv");
print OUTUNKNOWNSTAT "Sample\\t$name\\nTotal_Reads\\t$totalunknown\\nAverage\\t".round(average(\\@unknownvariantCounts),2)."\\nHighest_Read\\t$maxunknown\\n";

close(OUTUNKNOWNSTAT);

sub round{
     my ($num, $dec) = @_;
     return floor($num*10**$dec)/10**$dec;
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub highlow{
        my($data, $names) = @_;
        if (@data == 1){
            return 0;
        }
        my $lownum = 100000;
        my $highnum = 0;
        my $lowind = "";
        my $highind = "";
        my $i=0;
        foreach(@$data) {
           my $name = ${$names}[$i];
           if ($_ >= $highnum) { 
                $highind = ($highnum == $_ ? "$highind & $name":$name);
                $highnum = $_; 
           }
           if ($_ <= $lownum) { 
                $lowind = ($lownum == $_ ? "$lowind & $name":$name);
                $lownum = $_; 
           }
           $i++;
        }
        return ("$lownum ($lowind)", "$highnum ($highind)");
}
'''
}


process Blastn {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "blast_out/$filename"
}

input:
 set val(name),file(query) from g_45_fastaFile_g_50

output:
 file "*.tsv"  into g_50_outputFile

errorStrategy 'retry'
maxRetries 3

when:
(params.run_Blastn && (params.run_Blastn == "yes"))

shell:
db = params.Blastn.db
outfmt = params.Blastn.outfmt
max_target_seqs = params.Blastn.max_target_seqs
num_threads = params.Blastn.num_threads
task = params.Blastn.task
evalue = params.Blastn.evalue

blastn -task ${task} -query ${query} -evalue ${evalue} -db ${db} -num_threads ${num_threads} -max_target_seqs ${max_target_seqs} -out ${name}_blastout.tsv -outfmt ${outfmt}




}


process mergeCounts {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /out\/.*.tsv$/) "variant_counts/$filename"
}

input:
 file tsv from g_45_outputFileTSV_g_44.collect()

output:
 file "out/*.tsv"  into g_44_outputFile

errorStrategy 'retry'
maxRetries 1
	
script:
'''
#!/usr/bin/env perl

use strict;

#################### VARIABLES ######################
 my %tf = (
	"known" => 3,
	"unknown" => 2,
	"knownsum" => 2,
	"knownstats" => 2,
	"unknownstats" => 2
 );
my $indir=$ENV{'PWD'};
my $outdir=$ENV{'PWD'}."/out";
################### PARAMETER PARSING ####################
`mkdir -p $outdir`;
foreach my $type (keys %tf) {
	my $out = "out/".${type}.".tsv";
	opendir I, $indir or die "Could not open ./\\n";
	my @files = grep(/_${type}.tsv\$/, readdir(I));
	closedir I; 
	print "type ".$type." ... col: ".$tf{$type}."\\n";
	my @a=();
	my %b=();
	my %c=();
	my $i=0;
	foreach my $f (@files) { 
	  my $ind=$indir."/".$f;
	  my $libname=$f;
	  $libname=~s/_${type}.tsv//;
	  print "$libname\\n";
	  $i++;
	  $a[$i]=$libname;
	  print "${ind}\\n";
	  open IN, $ind or die "Could not open ".$ind."\\n";
	  my $lineNum = 0;
	
	while(<IN>)
	  {
	      next if ($lineNum++ == 0);
	      my @v=split; 
	      #print "Using $type $tf{$type}\\n";
	      my $count = $v[$tf{$type}-1];
	      #print $count."\\t";
	      $b{$v[0]}{$i}= $count;
	      $c{$v[0]}=$v[1] if($tf{$type}>2);
	  }
	  print "\\n";
	  close IN;
	}
	
	open OUT, ">".${out} or die "Could not open ".${out}." file for writting\\n";
	
	if ($tf{$type} == 2) {
	  print OUT "Seq";
	} else {
	  print OUT "Id\\tSeq";
	}
	
	for(my $j=1;$j<=$i;$j++)
	{
	 print OUT "\\t$a[$j]";
	}
	print OUT "\\n";
	
	foreach my $key (keys %b)
	{
	 if ($tf{$type} == 2) {
	   print OUT $key; #only print out the first column;
	 }
	 else
	 {
	   print OUT $key."\\t".$c{$key}; #only print out the first two columns
	 }
	 for(my $j=1;$j<=$i;$j++)
	 {
	   if (exists($b{$key}) && exists($b{$key}{$j})){
	     print OUT "\\t".$b{$key}{$j};
	   } else {
	     print OUT "\\t0";
	   }
	 }
	 print OUT "\\n";
	}
	
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
