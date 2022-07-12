#!/usr/bin/perl -w

use strict 'vars';

my $TmpFile;
foreach $TmpFile (<*>)
  {
    if ($TmpFile =~ /bosons\_delta\_exc.*\.dat/ )
      {
	my $NewName;
	my @TmpArray = split (/_/, $TmpFile);
	if ($TmpArray[4] eq "n")
	    {
	      my $N = $TmpArray[5];
	      my $S = 0;
	      if ($TmpArray[3] =~ /^m/)
		{
		  $TmpArray[3] =~ s/^m//;
		  $S -= $TmpArray[3];
		}
	      else
		{
		  $S += $TmpArray[3];
		}
	      if ($TmpArray[6] eq "nu")
		{
		  if ($TmpArray[7] == 1)
		    {
		      $S += $TmpArray[5] - 1;
		    }
		  else
		    {
		      if ($TmpArray[7] == 5)
			{
			  $S += 2 * $TmpArray[5] - 2;
			}
		    }
		  $NewName = "bosons_delta_n_".$N."_2s_".$S."_".$TmpArray[8];
		}
	      else
		{ 
		  $S += 2 * $TmpArray[5] - 2;
		  $NewName = "bosons_delta_n_".$N."_2s_".$S."_".$TmpArray[6];
		}
	      print ($TmpFile." ".$S." ".$N." ".$NewName."\n");
	      if (!(-e "n_$N"))
		{
		  mkdir("n_$N")
		}
	      `mv $TmpFile n_$N/$NewName`;
	    }
      }
    else
      {
	if ($TmpFile =~ /bosons\_delta.*\_2s\_.*\..*/ )
	  {
	    my @TmpArray = split (/_/, $TmpFile);
	    if ($TmpArray[2] eq "n")
	      {
		my $N = $TmpArray[3];
		if (!(-e "n_$N"))
		  {
		    mkdir("n_$N")
		  }
		print ($TmpFile."\n");
		`mv $TmpFile n_$N/$TmpFile`;
	      }
	  }
	else
	  {
	    if ($TmpFile =~ /bosons\_exc.*\.dat/ )
	      {
		my $NewName;
		my @TmpArray = split (/_/, $TmpFile);
		if ($TmpArray[3] eq "n")
		  {
		    my $N = $TmpArray[4];
		    my $S = 0;
		    if ($TmpArray[2] =~ /^m/)
		      {
			$TmpArray[2] =~ s/^m//;
			$S -= $TmpArray[2];
		      }
		    else
		      {
			$S += $TmpArray[2];
		      }
		    if ($TmpArray[5] eq "nu")
		      {
			if ($TmpArray[6] == 1)
			  {
			    $S += $TmpArray[4] - 1;
			  }
			else
			  {
			    if ($TmpArray[6] == 5)
			      {
				$S += 2 * $TmpArray[4] - 2;
			      }
			  }
			$NewName = "bosons_delta_n_".$N."_2s_".$S."_".$TmpArray[7];
		      }
		    else
		      { 
			$S += 2 * $TmpArray[4] - 2;
			$NewName = "bosons_delta_n_".$N."_2s_".$S."_".$TmpArray[5];
		      }
		    print ($TmpFile." ".$S." ".$N." ".$NewName."\n");
		    if (!(-e "n_$N"))
		      {
			mkdir("n_$N")
		      }
		    `mv $TmpFile n_$N/$NewName`;
		  }
	      }
	    else
	      {
		if ($TmpFile =~ /bosons\_delta\_n\_.*\..*/ )
		  {
		    my @TmpArray = split (/_/, $TmpFile);
		    my $N = $TmpArray[3];
		    if (!(-e "n_$N"))
		      {
			mkdir("n_$N")
		      }
		    my $S = 2 * ($N - 1);
		    my $NewName = "bosons_delta_n_".$N."_2s_".$S."_".$TmpArray[4];
		    print ($TmpFile." ".$S." ".$N." ".$NewName."\n");
		    `mv $TmpFile n_$N/$NewName`;
		  }
	      }
	  }
      }
  }

