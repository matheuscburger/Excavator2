#!/usr/bin/perl

use warnings;
use Pod::Usage;
use Getopt::Long;
use strict;
use File::Path;
use Cwd 'abs_path';

my ($target,$assembly,$MAPQ,$Program_Folder_Path,$Input_File_Path);
my ($record,@words,@unlinkfiles,@SamplesVect);
my ($verbose,$help,$man);
my ($Bam_File_Path,$Sample_Out_Folder_Path,$Sample_Label,$RC_Folder_Path);
my ($RCNorm_Folder_Path,$RCImages_Folder_Path);

$MAPQ = 20;

######################################################################
#
#	Reading user's options
#
######################################################################

GetOptions('target=s'=>\$target,'assembly=s'=>\$assembly,'mapq=i'=>\$MAPQ,'verbose|v'=>\$verbose,'help|h'=>\$help,'man|m'=>\$man,) or pod2usage ();

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
$workingfolder=substr $myscript, 0, ($L1 - $L2 - 22);

$Program_Folder_Path="$workingfolder";

#print "Program folder is: $Program_Folder_Path\n";


$Input_File_Path=$ARGV[0];


open(CHECKBOOK,"$Input_File_Path") || die "Couldn't open the input file $Input_File_Path.";

while($record=<CHECKBOOK>){
  chomp($record);
  @words=split(' ',$record);
  
  $Bam_File_Path=$words[0];
  $Sample_Out_Folder_Path=$words[1];
  $Sample_Label=$words[2];
  
  $RC_Folder_Path="$Sample_Out_Folder_Path/RC";
  $RCNorm_Folder_Path="$Sample_Out_Folder_Path/RCNorm";
  $RCImages_Folder_Path="$Sample_Out_Folder_Path/Images";
  
  if(-e "$RC_Folder_Path"){
    print "'$RC_Folder_Path' folder ready!\n";
  }
  else{
    mkpath ("$RC_Folder_Path");
    print "Creating '$RC_Folder_Path' folder...\n";
  }
  
  
  if(-e "$RCNorm_Folder_Path"){
    print "'$RCNorm_Folder_Path' folder ready!\n";
  }
  else{
    mkpath ("$RCNorm_Folder_Path");
    print "Creating '$RCNorm_Folder_Path' folder...\n";
  }
  
  
  if(-e "$RCImages_Folder_Path"){
    print "'$RCImages_Folder_Path' folder ready!\n";
  }
  else{
    mkpath ("$RCImages_Folder_Path");
    print "Creating '$RCImages_Folder_Path' folder...\n";
  }
  
  print "Working on sample $Sample_Label.\n";
  
  push(@SamplesVect,"$Sample_Label");
  
  ######################################################################
  #
  #  Read count
  #
  ######################################################################
  
  print "Creating Read Count Data...\n";
  system qq($Program_Folder_Path/lib/bash/./FiltBam.sh $Bam_File_Path $MAPQ $Sample_Out_Folder_Path $Program_Folder_Path $Sample_Label $target $assembly);
  print "Removing temporary files...\n";
  rmtree("$Sample_Out_Folder_Path/.tmp");
  print "Read Count done!\n";
  print "Normalizing Read Count Data...\n";
  system qq(R --slave --args $Program_Folder_Path,$Sample_Out_Folder_Path,$Sample_Label,$target,$assembly < $Program_Folder_Path/lib/R/EXCAVATORNormalizationExome.R);
  print "Normalization done!\n";
}

print "Performed RC calculations and Normalization on samples: \n";
print "@SamplesVect";
@unlinkfiles=("$Input_File_Path");
unlink @unlinkfiles;


