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

my ($Program_Folder_Path,$Main_Output_Folder_Path,$Sample_Out_Folder_Path,$Sample_Label,$Analysis_Label,$Assembly,$Target_Name,$Target_File_Path,$RC_Folder_Path,$RCNorm_Folder_Path,$Images_Folder_Path,$R_Target_Folder,$R_Target_Path,$R_Norm_Path,$Input_File_Path,$Mode,$Sample_Res_Folder_Path,$Sample_Plot_Folder_Path,$Sample_Data_Folder_Path);

my ($verbose,$help,$man);


######################################################################
#
#	Reading user's options
#
######################################################################

GetOptions('assembly|a=s'=>\$Assembly,'target|t=s'=>\$Target_Name,'output|o=s'=>\$Main_Output_Folder_Path,'mode|e:s'=>\$Mode,'verbose|v'=>\$verbose,'help|h'=>\$help,'man|m'=>\$man,) or pod2usage ();

@ARGV == 1 or pod2usage ("Syntax error: the number of arguments found at command line is incorrect.");

$Input_File_Path=$ARGV[0];

#die ("User ERROR: Please set option --mode with a valid string ('pooling' or 'somatic' or 'nocontrol').\n") if( $Mode ne "pooling" && $Mode ne "somatic" && $Mode ne "nocontrol" );

#die ("User ERROR: Please set option --assembly with a valid string ('hg18' or 'hg19').\n") if( $Assembly ne "hg18" && $Assembly ne "hg19" );

#print "Output folder is: $Main_Output_Folder_Path\n";

######################################################################
#
#	Defining system variables
#
######################################################################

my ($myscript,$workingfolder,$L1,$L2);

$myscript = abs_path($0);
$L1=length($myscript);
$L2=length($0);
$workingfolder=substr $myscript, 0, ($L1 - $L2 - 25);

$Program_Folder_Path="$workingfolder";

#print "Program folder is: $Program_Folder_Path\n";

die ("User ERROR: the target specified with --target wasn't found for assembly $Assembly. \nPlease check --target and --assembly. \nError found") if (! -d "$Program_Folder_Path/data/targets/$Assembly/$Target_Name" );

$R_Target_Folder="$Program_Folder_Path/data/targets/$Assembly/$Target_Name";

######################################################################
#
#	Reading input file
#
######################################################################

open(CHECKBOOK,"$Input_File_Path") || die "Couldn't open the system input file.";

while($record=<CHECKBOOK>){
    chomp($record);
    @words=split(' ',$record);

    $Analysis_Label=$words[0];
    $RC_Folder_Path=$words[1];
    $Sample_Label=$words[2];


    $Sample_Res_Folder_Path="$Main_Output_Folder_Path/Results/$Sample_Label";
    $Sample_Plot_Folder_Path="$Main_Output_Folder_Path/Plots/$Sample_Label";

    if($Analysis_Label =~ m/T/){ 
     mkpath($Sample_Res_Folder_Path);
     mkpath($Sample_Plot_Folder_Path);
   }

}
# 
# ######################################################################
# #
# #	Data calculations
# #
# ######################################################################
# 
# 
# 
print "Starting Segmentation and Calling...\n";
 system qq(R --slave --args $Main_Output_Folder_Path,$R_Target_Folder,$Input_File_Path,$Mode,$Target_Name,$Program_Folder_Path,$Assembly < $Program_Folder_Path/lib/R/EXCAVATORInferenceExome.R);
print "Segmentation and Calling Complete\n";



  system qq(R --slave --args $Main_Output_Folder_Path,$Input_File_Path < $Program_Folder_Path/lib/R/EXCAVATORPlotsExome.R);
# 
# #system qq(cp $Program_Folder_Path/ParameterFile.txt $Main_Output_Folder_Path/ParameterFile.$ID.txt);
# 
# #system qq(cp $Program_Folder_Path/.Ainput.$ID $Main_Output_Folder_Path/Input.$ID.Copy.txt);
# 
# print "Performed CNVs analysis with samples: \n";
# print "@SamplesVect";
# 
# close(CHECKBOOK);
# @unlinkfiles=("$Program_Folder_Path/.Ainput.$ID");
# unlink @unlinkfiles;


