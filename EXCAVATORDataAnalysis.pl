#!/usr/bin/perl

use warnings;
use Pod::Usage;
use Getopt::Long;
use strict;
use File::Path;
use Cwd 'abs_path';

######################################################################
#
#	Variables initialization with default values
#
######################################################################

my ($record,@words,@unlinkfiles,@SamplesVect);

my ($Program_Folder_Path,$Main_Output_Folder_Path,$Sample_Out_Folder_Path,$Sample_Label,$Analysis_Label,$Assembly,$Target_Name,$Target_File_Path,$RC_Folder_Path,$RCNorm_Folder_Path,$Images_Folder_Path,$R_Target_Folder,$R_Target_Path,$R_Norm_Path,$Input_File_Path,$Mode,$Sample_Res_Folder_Path,$Sample_Plot_Folder_Path,$Sample_Data_Folder_Path,$processors);

my ($verbose,$help,$man);


######################################################################
#
#  Defining system variables
#
######################################################################

my ($myscript,$workingfolder,$L1,$L2);

$myscript = abs_path($0);
$L1=length($myscript);
$L2=length($0);
$workingfolder=substr $myscript, 0, ($L1 - $L2 - 1);


$Program_Folder_Path="$workingfolder";


######################################################################
#
#	Reading user's options
#
######################################################################

GetOptions('processors=s'=>\$processors,'assembly|a=s'=>\$Assembly,'target|t=s'=>\$Target_Name,'output|o=s'=>\$Main_Output_Folder_Path,'mode|e:s'=>\$Mode,'verbose|v'=>\$verbose,'help|h'=>\$help,'man|m'=>\$man,) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error: the number of arguments found at command line is incorrect.");

$Input_File_Path=$ARGV[0];




die ("User ERROR: Please set option --mode with a valid string ('pooling' or 'paired').\n") if( $Mode ne "pooling" && $Mode ne "paired" );


#die ("User ERROR: the target specified with --target wasn't found for assembly $Assembly. \nPlease check --target and --assembly. \nError found") if (! -d "$Program_Folder_Path/data/targets/$Assembly/$Target_Name" );

######################################################################
#
#  Checking system folders
#
######################################################################

print "Checking output folders...\n";

if(-e $Main_Output_Folder_Path){ 
  print "'$Main_Output_Folder_Path' folder ready!\n";
}
else{
  mkpath($Main_Output_Folder_Path);
  if(-e $Main_Output_Folder_Path){
    print "'$Main_Output_Folder_Path' folder created!\n";
  }  
}

print "Checking output subfolders...\n";


if(-e "$Main_Output_Folder_Path/Results"){
  print "'$Main_Output_Folder_Path/Results' folder ready!\n";
}
else{
  mkdir ("$Main_Output_Folder_Path/Results");
  print "Creating '$Main_Output_Folder_Path/Results' folder...\n";
}

if(-e "$Main_Output_Folder_Path/Plots"){
  print "'$Main_Output_Folder_Path/Plots' folder ready!\n";
}
else{
  mkdir ("$Main_Output_Folder_Path/Plots");
  print "Creating '$Main_Output_Folder_Path/Plots' folder...\n";
}

if(-e "$Main_Output_Folder_Path/.tmp"){
  print "'$Main_Output_Folder_Path/.tmp' folder ready!\n";
}
else{
  mkdir ("$Main_Output_Folder_Path/.tmp");
  print "Creating '$Main_Output_Folder_Path/.tmp' folder...\n";
}


######################################################################
#
#  When mode is pooling generates the Read Count for Controls
#
######################################################################


if($Mode eq "pooling"){ 
  print "Creating Pooling Control!\n";
  system qq(R --slave --args $Program_Folder_Path,$Main_Output_Folder_Path,$Program_Folder_Path/data/targets/$Assembly/$Target_Name,$Input_File_Path,$Mode,$Target_Name < $Program_Folder_Path/lib/R/PoolingCreateControl.R);
  print "Pooling Control Created!\n";
}




######################################################################
#
#  Creating Multi Processor Analysis
#
######################################################################


print "Preparing Multiprocessor Analysis...\n";
system qq(R --slave --args $Program_Folder_Path,$Input_File_Path,$Main_Output_Folder_Path,$Target_Name,$Assembly,$processors,$Mode < $Program_Folder_Path/lib/R/DataAnalysisParallel.R);
print "Multiprocessor Analysis Complete\n";
print "Starting Multiprocessor Analysis!\n";


######################################################################
#
#  Running Multi Processor Analysis
#
######################################################################

my $Input_File_Parallel="$Main_Output_Folder_Path/.tmp/ParallelAnalyzePerla.sh";


open(CHECKBOOK,"$Input_File_Parallel") || die "Couldn't open the input file $Input_File_Parallel.";
my @pids;
while($record=<CHECKBOOK>){
  my $childpid = fork() or exec($record);
  push(@pids,$childpid);
}

print "My Children: ", join(' ',@pids), "\n";
waitpid($_,0) for @pids;

print "Multiprocessor Analysis Complete!\n";




######################################################################
#
#  Documentation
#
######################################################################

=head1 SYNOPSIS 

 perl EXCAVATORDataAnalysis.pl [arguments] [options]
 
  Options:

       -h, --help                   Print help message.
       -m, --man                    Print complete documentation.
       -v, --verbose                Use verbose output.

 Function:
 
 EXCAVATORDataAnalysis.pl performs segmentation of the WMRC and classify each segmented region as one of 5 possible discrete states (2-copy deletion, 1-copy deletion, normal, 1-copy duplication and N-copy amplification).

 Example: 
 
 EXCAVATOR2> perl EXCAVATORDataAnalysis.pl ExperimentalFileAnalysis.w50K.txt --processors 6 --target MyTarget_w50K --assembly hg19 --output /.../OutEXCAVATOR2/Results_MyProject_w50K --mode pooling/paired

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief usage message and detailed explanation of options.

=item B<--man>

Print the complete manual of the script.

=item B<--verbose>

Use verbose output.

=item B<--processors>

The number of thread to use for the analysis.

=item B<--output>

The output folder for resulting files.

=item B<--assembly>

The assembly exploited for read mapping and target initialization.

=item B<--target>

The "target name" used for target initialization with TargetPerla.pl.

=item B<--mode>

The experimental design mode to use. The possible options are "pooling" or "paired".

=back

=head1 DESCRIPTION

EXCAVATORDataAnalysis.pl perform the segmentation of the WMRC by means of the Shifting Level Model algorithm and exploits FastCall algorithm to classify each segmented region as one of the five possible discrete states (2-copy deletion, 1-copy deletion, normal, 1-copy duplication and N-copy amplification). The FastCall calling procedure takes into account sample heterogeneity and exploits the Expectation Maximization algorithm to estimate the parameters of a five gaussian mixture model and to provide the probability that each segment belongs to a specific copy number state.


EXCAVATOR2 is freely available to the community for non-commercial use. For questions or comments, please contact "romina.daurizio@gmail.com".
=cut
