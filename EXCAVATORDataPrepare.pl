#!/usr/bin/perl

use warnings;
use Pod::Usage;
use Getopt::Long;
use strict;
use File::Path;
use Cwd 'abs_path';
#use Forks::Super;

my ($target,$assembly,$MAPQ,$Program_Folder_Path,$Input_File_Path);
my ($record,@words,@unlinkfiles,@SamplesVect);
my ($verbose,$help,$man);
my ($Bam_File_Path,$Sample_Out_Folder_Path,$Sample_Label,$RC_Folder_Path);
my ($RCNorm_Folder_Path,$RCImages_Folder_Path,$processors);

$MAPQ = 0;

######################################################################
#
#  Reading user's options
#
######################################################################

GetOptions('processors=s'=>\$processors,'target=s'=>\$target,'assembly=s'=>\$assembly,'mapq=i'=>\$MAPQ,'verbose|v'=>\$verbose,'help|h'=>\$help,'man|m'=>\$man,) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error: too many arguments were presented at command line.");


######################################################################
#
#	Defining system variables
#
######################################################################

my ($myscript,$workingfolder,$L1,$L2);

$myscript = abs_path($0);
$L1=length($myscript);
$L2=length($0);
$workingfolder=substr $myscript, 0, ($L1 - $L2 - 1);


$Program_Folder_Path="$workingfolder";


$Input_File_Path=$ARGV[0];

my $range=99999;
my $ID=int(rand($range));
my $labeltemp=".tmp.$ID";


### Creating Temporary Folder for MultiProcessor analysis ###
my $Temp_Folder="$Program_Folder_Path/$labeltemp";
mkpath ("$Temp_Folder");


my $Input_File_Parallel="$Temp_Folder/ParallelReadPerla.sh";


print "Program folder is: $Program_Folder_Path\n";


######################################################################
#
#  Creating Multi Processor Analysis
#
######################################################################


  print "Preparing Multiprocessor Analysis...\n";
  system qq(R --slave --args $Program_Folder_Path,$Input_File_Path,$labeltemp,$target,$assembly,$processors < $Program_Folder_Path/lib/R/DataPrepareParallel.R);
  print "Starting Multiprocessor Analysis!\n";

open(CHECKBOOK,"$Input_File_Parallel") || die "Couldn't open the input file $Input_File_Parallel.";
my @pids;
while($record=<CHECKBOOK>){
my $childpid = fork() or exec($record);
push(@pids,$childpid);
}

#print "My Children: ", join(' ',@pids), "\n";
waitpid($_,0) for @pids;

print "Multiprocessor Analysis Complete!\n";

######################################################################
#
#  Documentation
#
######################################################################

=head1 SYNOPSIS 

 perl EXCAVATORDataPrepare.pl [arguments] [options]

 Options:

       -h, --help                   Print help message.
       -m, --man                    Print complete documentation.
       -v, --verbose                Use verbose output.
           --mapq <integer>         Select mapping quality for .bam file filtering; if omitted default value is 0.
           
 Function:
 
 perl EXCAVATORDataPrepare.pl performs RC calculations, data normalization and data analysis on multiple .bam files.

 Example: 
 
 EXCAVATOR2> perl EXCAVATORDataPrepare.pl ExperimentalFilePrepare.w50000.txt --processors 6 --target MyTarget_w50000 --assembly hg19 


=head1 OPTIONS

=over 8

=item B<--help>

Print a brief usage message and detailed explanation of options.

=item B<--man>

Print the complete manual of the script.

=item B<--verbose>

Use verbose output.

=item B<--mapq>

Sets the numeric value of the mapping quality for .bam file filtering; must be an integer number. If omitted default value is 0.

=item B<--assembly>

The assembly exploited for read mapping and target initialization (hg19 or hg38).

=item B<--target>

The "target name" used for target initialization with TargetPerla.pl.


=back

=head1 DESCRIPTION

ReadPerla.pl is a Perl script which is part of the EXCAVATOR2 package. It performs RC calculations for In-target and Off-target regions, data normalization and data analysis on multiple .bam files.

The mapping quality value which is used by SAMtools can be set by means of the option --mapq when running ReadPerla.pl. If omitted default value is 0.

EXCAVATOR is freely available to the community for non-commercial use. For questions or comments, please contact "romina.daurizio@gmail.com".
=cut


