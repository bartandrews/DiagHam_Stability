#!/usr/bin/perl -w
#
# script for plotting output of FCIHofstadterModel using gnuplot
#
use strict 'vars';


my $ExtendSymmetries=0;
while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-s/ )
      {
	$ExtendSymmetries=1;
      }
    shift(@ARGV);
  }



if (!defined($ARGV[0]))
  {
    print("usage PlotHofstadterSpectrum.pl [-s] [Directory|single file]\n");
    exit(1);
  }


my @ListFiles;
my $ReadMore = 1;
if ( defined($ARGV[0]) )
  {
    if ( -d$ARGV[0] )
      {
	chdir($ARGV[0]);
      }
    else # single file:
      {
	my $TmpFile2;
	foreach $TmpFile2 (@ARGV)
	  {
	    if (-e $TmpFile2)
	      {
		push (@ListFiles, $TmpFile2);
	      }
	  }
	$ReadMore=0;
      }
  }
my $TmpFile;
if ($ReadMore == 1) # array still empty?
  {
    foreach $TmpFile (<*>)
      {
	if ($TmpFile =~ /bosons_hofstadter_X_\d*_Y_\d*_q_\d*_n_\d*_x_\d*_y_\d*.*_gx_0_gy_0.dat/)
	  {
	    push (@ListFiles, $TmpFile);
	  }
      }
  }


foreach $TmpFile (@ListFiles)
  {
    my $Max;
    my $Min;
    print ("Processing $TmpFile\n");
    if ($ExtendSymmetries)
      {
	my $ExtendedFile = $TmpFile;
	$ExtendedFile =~ s/.dat/\_ext.dat/;
	print ("Creating $ExtendedFile\n");
	&ExtendSymmetry($TmpFile, $ExtendedFile);
	&CreatePostScript($ExtendedFile);
      }
    else
      {
	&CreatePostScript($TmpFile);
      }  
    print ("\n\n");
#    &FindMinMax($TmpFile, 1, \$Min, \$Max, 0, 0, 15);
#    print ($TmpFile." ".$Min." ".$Max."\n");
  }


# find minimum and maximum values in a file
#
# $_[0] = file name
# $_[1] = column where to search
# $_[2] = reference on min value
# $_[3] = reference on max value
# $_[4] = additonnal constraint on another column
# $_[5] = max value for the constraint
# $_[6] = min value for the constraint

sub FindMinMax
  {
    my $FileName = $_[0];
    my $Column = $_[1];
    my $Min = $_[2];
    my $Max = $_[3];
    if (defined($_[4]))
      {
	my $ColumnConstraint = $_[4];
	my $MinConstraint = $_[5];
	my $MaxConstraint = $_[6];
	my $Flag = 0;
	open (INFILE, $FileName);
	my $TmpLine;
	my @TmpArray;
#	$TmpLine = <INFILE>;
#	$TmpLine = <INFILE>;
	foreach $TmpLine (<INFILE>)
	  {
	    chomp ($TmpLine);
	    if ($TmpLine =~ m/^\s*#+/)
	      {
		# print("ignoring comment line ".$TmpLine."\n");
	      }
	    else
	      {
		@TmpArray = split (/ /, $TmpLine);
		if (($TmpArray[$ColumnConstraint] <= $MaxConstraint) && ($TmpArray[$ColumnConstraint] >= $MinConstraint)) {
		  if ($Flag == 0) {
		    $$Min = $TmpArray[$Column];
		    $$Max = $$Min;
		    $Flag = 1;
		  } else {
		    if ($TmpArray[$Column] < $$Min) {
		      $$Min = $TmpArray[$Column];
		    }
		    if ($TmpArray[$Column] > $$Max) {
		      $$Max = $TmpArray[$Column];
		    }
		  }
		}
	      }
	  }
      }
    else
      {
	open (INFILE, $FileName);
	my $TmpLine;
	my @TmpArray;
#	$TmpLine = <INFILE>;
#	$TmpLine = <INFILE>;
	$TmpLine = <INFILE>;
	@TmpArray = split (/ /, $TmpLine);
	$$Min = $TmpArray[$Column];
	$$Max = $$Min;
	foreach $TmpLine (<INFILE>)
	  {
	    chomp ($TmpLine);
	    @TmpArray = split (/ /, $TmpLine);
	    if ($TmpArray[$Column] < $$Min)
	      {
		$$Min = $TmpArray[$Column];
	      }
	    if ($TmpArray[$Column] > $$Max)
	      {
		$$Max = $TmpArray[$Column];
	      }
	  }
	close (INFILE);
      }
}


# create postscript graph from data file
#
# $_[0] = input file name
# $_[1] = output file name
sub ExtendSymmetry
    {
       my $TmpFile=$_[0];
       my $ExtendedFile=$_[1];
       
       $TmpFile =~ /X\_(\d+)\_Y\_(\d+)\_q\_(\d+)\_n\_(\d+)\_x\_(\d+)\_y\_(\d*)\_/;
       my $FileX = $1;
       my $FileY = $2;
       my $Fileq = $3;
       my $Filen = $4;
       my $FileLx = $5;
       my $FileLy = $6;

       print("Trying to open outfile $ExtendedFile\n");

       open (OUTFILE, ">$ExtendedFile") or die ("Cannot open outfile $ExtendedFile\n");
       open (INFILE, $TmpFile);
       
       foreach my $TmpLine (<INFILE>)
	 {
	   if ($TmpLine =~ m/^\s*#+/)
	     {
	       print OUTFILE $TmpLine;
	     }
	   else
	     {
	       chomp $TmpLine;
	       my @TmpArray = split (/ /,  $TmpLine);	   	   
	       my $kx = $TmpArray[0];
	       my $ky = $TmpArray[1];
	       my $E =  $TmpArray[2];
	       
	       if ( ($kx==0) || (( $FileLx%2 == 0) && ($kx== $FileLx/2)))
		 {
		   if ( ($ky==0) || (( $FileLy%2 == 0) && ($ky== $FileLy/2)))
		     {
		       print OUTFILE $TmpLine."\n";
		     }
		   else
		     {
		       print OUTFILE ($TmpLine."\n"."$kx ".($FileLy-$ky)." $E\n");
		     }
		 }
	       else
		 {
		   if ( ($ky==0) || (( $FileLy%2 == 0) && ($ky== $FileLy/2)))
		     {
		       print OUTFILE ($TmpLine."\n".($FileLx-$kx)." $ky $E\n");
		     }
		   else
		     {
		       print OUTFILE ($TmpLine."\n".($FileLx-$kx)." $ky $E\n");
		       print OUTFILE ("$kx ".($FileLy-$ky)." $E\n");
		       print OUTFILE (($FileLx-$kx)." ".($FileLy-$ky)." $E\n");
		       
		     }
		 }
	       
	     }
	   
	 }
       close (INFILE);
       close (OUTFILE);
     }
    
# create postscript graph from data file
#
# $_[0] = file name
# /deprecated option: $_[1] = print flag (1 if true)

sub CreatePostScript
  {
    my $FileName = $_[0];
    my $PrintFlag = 0;
    my $Max;
    my $Min;
    &FindMinMax($FileName, 2, \$Min, \$Max, 0, 0, 14);
    my $Delta = ($Max - $Min) / 20.0;
    $Max += $Delta;
    $Min -= $Delta;
    my $PlotFileName =  $FileName; # previously random: "tmp".time().".p";
    $PlotFileName =~ s/\.dat/\.gp/;
    #my @TmpArray = split (/_/,  $OutputFile);
  
    $FileName =~ /X\_(\d+)\_Y\_(\d+)\_q\_(\d+)\_n\_(\d+)\_x\_(\d+)\_y\_(\d*)\_/;
    my $FileX = $1;
    my $FileY = $2;
    my $Fileq = $3;
    my $Filen = $4;
    my $Filex = $5;
    my $Filey = $6;

    #my $Title = "N = ".$Filen."  $Filex x $Filey MUC with ($Fileq; $FileX x $FileY)";
    my $Title = "C=, {/Symbol n}=:  N = ".$Filen."  L_x=$Filex, L_y=$Filey, MUC: p=$Fileq, q=$FileX x $FileY";

    my $MaxX = $Filex * $Filey;

    my $Delta = ($Max-$Min)/10.0;

    my $OutputFile = $FileName;
    $OutputFile =~ s/\.dat/\.ps/;
    open (OUTFILE, ">$PlotFileName");
    print OUTFILE ("set xrange [-1:".($MaxX+1)."]
set yrange [".$Min.":".($Max+$Delta)."]
set xlabel \"k_x*L_y + k_y\"
set ylabel \"E(k)\"
set size 1.0, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
plot \"".$FileName."\" using (\$1*".$Filey."+\$2):3 title \"".$Title."\" with points lc rgbcolor \"red\" lt 6 pt 1 ps 2
");
    close (OUTFILE);
    # print plotfile on screen
    open (INFILE, $PlotFileName);
    my $TmpLine;
    foreach $TmpLine (<INFILE>)
      {
	print $TmpLine;
      }
    close (INFILE);

    `gnuplot $PlotFileName`;
    if ($PrintFlag == 1)
      {
	`lpr $OutputFile`;
      }
    #`rm -f $PlotFileName`;
  }


