#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;

#if (!(defined($ARGV[4])))
#  {
#    die "usage: FindFermionLaplacianNeutralGap.pl ";
#  }

my $FlagStrict = 0;
my $Result = GetOptions ("strict+" => \$FlagStrict); 

my @IncompressibleState2S;
my @IncompressibleStateN;
&FindIncompressibleState(\@IncompressibleStateN, \@IncompressibleState2S);

my $MaxQ = 20;
my $QValue = 2;
while ($QValue <= 20)
  {
    my $PValue = 1;
    while ($PValue < $QValue)
      {
	my $GCDPQ = &FindGCD($PValue, $QValue);
	my $ReducedP = $PValue / $GCDPQ;
	my $ReducedQ = $QValue / $GCDPQ;
	if ($ReducedQ >= $QValue)
	  {
	    my $Pos1 = 0;
	    while ($Pos1 < $#IncompressibleStateN)
	      {
		my $Pos2 = $Pos1 + 1;
		my @TmpFraction;
		while ($Pos2 <= $#IncompressibleStateN)
		  {
		    if (($PValue * ($IncompressibleState2S[$Pos1] - $IncompressibleState2S[$Pos2])) == 
			($QValue * ($IncompressibleStateN[$Pos1] - $IncompressibleStateN[$Pos2])))
		      {
			push (@TmpFraction, $Pos2);
		      }
		    $Pos2++;
		  }
		if ((($#TmpFraction >= 1) && ($FlagStrict == 1)) ||
		    (($#TmpFraction >= 0) && ($FlagStrict == 0)))
		  {
		    $Pos2 = 0;
		    my $Flag = 0;
		    while (($Flag == 0) && ($Pos2 < $Pos1))
		      {
			if (($PValue * ($IncompressibleState2S[$Pos1] - $IncompressibleState2S[$Pos2])) == 
			    ($QValue * ($IncompressibleStateN[$Pos1] - $IncompressibleStateN[$Pos2])))
			  {
			    $Flag = 1;
			  }
			$Pos2++
		      }
		    if ($Flag == 0)
		      {
			push (@TmpFraction, $Pos1);
			print ("possible candidate for fraction ".$PValue."/".$QValue.": ");
			foreach $Pos2 (@TmpFraction)
			  {
			    print ("(".$IncompressibleStateN[$Pos2].",".$IncompressibleState2S[$Pos2].") ");
			  }
			print (" with shift ");
			$Pos2 = ($QValue * $IncompressibleStateN[$TmpFraction[0]]) - ($PValue * $IncompressibleState2S[$TmpFraction[0]]);
			my $GCDPQ = &FindGCD($Pos2, $PValue);
			$Pos2 /= $GCDPQ;
			print $Pos2;
			if  (($PValue / $GCDPQ)!= 1)
			  {
			    print ("/".($PValue / $GCDPQ));
			  }
			print ("\n");
		      }
		  }
		$Pos1++;
	      }
	  }
	$PValue++;
      }
    $QValue++;
  }

# find all possible candidates for incompressible state
#
# $_[0] =  reference on the array containing N value of the incompressible states
# $_[1] =  reference on the array containing 2S value of the incompressible states

sub FindIncompressibleState
  {
    my $IncompressibleStateN = $_[0];
    my $IncompressibleState2S = $_[1];
    my $NbrFermions = 1000;
    my $MaxN = 0;
    my $TmpFileName;
    
    foreach $TmpFileName (<*>)
      {
	if ((-d $TmpFileName) && ($TmpFileName =~ /n\_\d*/))
	  {
	    $TmpFileName =~ s/n\_//;
	    if ($TmpFileName < $NbrFermions)
	      {
		$NbrFermions = $TmpFileName;
	      }
	    else
	      {
		if ($TmpFileName > $MaxN)
		  {
		    $MaxN = $TmpFileName;
		  }
	      }
	  }
      }
    
    while ($NbrFermions <= $MaxN)
      {
	$TmpFileName = "n_".$NbrFermions;
	if (-e $TmpFileName)
	  {
	    chdir ("n_".$NbrFermions);
	    foreach $TmpFileName (<*>)
	      {
		if ((-f $TmpFileName) && ($TmpFileName =~ /^fermions\_.*\_n\_\d*\_2s\_\d*\_l\.dat/))
		  {
		    my $SValue = $TmpFileName;
		    $SValue =~ s/^.*\_2s\_(\d*)\_.*/$1/;
		    open (INFILE, $TmpFileName);
		    my $TmpLine2 = <INFILE>;
		    chomp ($TmpLine2);
		    my @TmpValue = split (/ /, $TmpLine2);		
		    my $GroundEnergy = $TmpValue[1];
		    my $GroundMomentum = $TmpValue[0];
		    while (defined($TmpLine2 = <INFILE>))
		      {
			chomp ($TmpLine2);
			@TmpValue = split (/ /, $TmpLine2);		
			if (($TmpValue[1] < $GroundEnergy) && 
			    ((abs($TmpValue[1] - $GroundEnergy) > (abs($TmpValue[1]) * 1e-10)) 
			     || (abs($TmpValue[1] - $GroundEnergy) > 1e-10)))
			  {
			    $GroundEnergy = $TmpValue[1];
			    $GroundMomentum = $TmpValue[0];			
			  }
		      }
		    close(INFILE);
		    if ($GroundMomentum == 0)
		      {
			push (@$IncompressibleStateN, $NbrFermions);
			push (@$IncompressibleState2S, $SValue);
		      }
		  }
	      }
	    chdir ("..");
	  }
	++$NbrFermions
      }
  }


sub FindGCD
  {
    my $TmpA = $_[0];
    my $TmpB = $_[1];
    if ($TmpA > $TmpB)
      {
	return &FindGCD($TmpB, $TmpA);
      }
    my $Tmp;
    while (($TmpB != $TmpA) && ($TmpA != 0))
      {
	$Tmp = $TmpB % $TmpA;
	$TmpB = $TmpA;
	$TmpA = $Tmp;	
      }
    return $TmpB;
  }
