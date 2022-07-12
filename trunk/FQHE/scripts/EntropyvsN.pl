#!/usr/bin/perl -w

use strict 'vars';

my $TmpFile;
my $TmpMaxBlockSize = 0;
my $MaxN = 0;
my $InteractionName;
foreach $TmpFile (<*>)
  {
     if ($TmpFile =~ /^fermions\_(.+)\_n\_(\d+)\_2s\_(\d+)\_lz\_0\.0\.ent$/)
       {
	 if ($2 > $MaxN)
	   {
	     $MaxN = $2;
	   }
	 if ($3 > $TmpMaxBlockSize)
	   {
	     $TmpMaxBlockSize = $3;
	   }
	 $InteractionName = $1;
       }
  }


my @Entropies;
foreach $TmpFile (<*>)
  {
     if ($TmpFile =~ /^fermions\_$InteractionName\_n\_(\d+)\_2s\_(\d+)\_lz\_0\.0\.ent$/)
       {
	 my $TmpN = $1;
	 my $Tmp2S = $2;
	 unless (open (INFILE, $TmpFile))
	   {
	     die ("can't open ".$TmpFile."\n");
	   }
	 my $TmpLine;
	 while (defined($TmpLine = <INFILE>))
	   {
	     chomp ($TmpLine);
	     my @TmpArray = split (/\s+/, $TmpLine);
	     if (defined($Entropies[$TmpArray[0]]))
	       {
		 my $TmpValue = $Entropies[$TmpArray[0]];
		 $$TmpValue{$TmpN} = $TmpArray[1];
	       }
	     else
	       {
		 my %TmpValue;
		 $TmpValue{$TmpN} = $TmpArray[1];
		 $Entropies[$TmpArray[0]] = \%TmpValue;
	       }
	     $TmpArray[0] = $Tmp2S + 1 -$TmpArray[0];
	     if (defined($Entropies[$TmpArray[0]]))
	       {
		 my $TmpValue = $Entropies[$TmpArray[0]];
		 $$TmpValue{$TmpN} = $TmpArray[1];
	       }
	     else
	       {
		 my %TmpValue;
		 $TmpValue{$TmpN} = $TmpArray[1];
		 $Entropies[$TmpArray[0]] = \%TmpValue;
	       }
	   }
	 close (INFILE);
       }
   }

my $TmpLa = 1;
my $SLaResults = "";
my $LatexFile = "";
my $LatexParity = 1;
my $NbrLatexFile = 1;
while ($TmpLa <= $TmpMaxBlockSize)
  {
    my $TmpValue = $Entropies[$TmpLa];
    my @NValues = sort {$a <=> $b} (keys(%$TmpValue));
    if ($#NValues > 3)
      {
	my $OutputFile = "fermions_".$InteractionName."_la_".$TmpLa.".ent";
	unless (open(OUTFILE, ">".$OutputFile))
	  {
	    die ("can't create ".$OutputFile."\n");
	  }
	my $TmpN;
	foreach $TmpN (@NValues)
	  {
	    if ($TmpN > 6)
	      {
		print OUTFILE $TmpN." ".(1.0 / $TmpN)." ".$$TmpValue{$TmpN}."\n";
	      }
	  }
	close (OUTFILE);
	unless (open(OUTFILE, ">".$OutputFile.".p"))
	  {
	    die ("can't create ".$OutputFile.".p\n");
	  }
	print OUTFILE "set xrange [0:0.15]
set terminal postscript landscape enhanced eps \"Helvetica\" 22
set output \"".$OutputFile.".eps\"
set xlabel \"1/N\" font \"Helvetica,26\"
set ylabel \"S_la\" font \"Helvetica,26\"
set size 1.0, 1.0
f1(x)= a0 +  a1 * x ** 2.5
fit f1(x) \"".$OutputFile."\" using 2:3 via a0,a1
plot \"".$OutputFile."\" using 2:3 notitle with points 31, f1(x) notitle with lines lt 1 lw 2
";
	close (OUTFILE);
	unlink("fit.log");
	`gnuplot $OutputFile.p`;
	my $TmpLine;
	unless (open(INFILE, "fit.log"))
	  {
	    die ("can't open fit.log\n");
	  }
	$LatexFile .= "\\begin{figure}
\\includegraphics[width=12cm,angle=0]{".$OutputFile.".eps}
\\caption{\$l_a=".$TmpLa."\$. fit : ";
	foreach $TmpLine (<INFILE>)
	  {
	    chomp ($TmpLine);
	    if ($TmpLine =~ /^a0\s+\=\s+(\-?\d+\.?\d*e?\-?\d*)\s+\+\/\-\s+(\-?\d+\.?\d*e?\-?\d*)/)
	      {
		 $SLaResults .= $TmpLa." ".(sqrt($TmpLa))." ".$1." ".$2."\n";
		 $LatexFile .= "\$a_0=".$1."\\pm".$2."\$ and ";
	      }
	    if ($TmpLine =~ /^a1\s+\=\s+(\-?\d+\.?\d*e?\-?\d*)\s+\+\/\-\s+(\-?\d+\.?\d*e?\-?\d*)/)
	      {
		 $LatexFile .= "\$a_1=".$1."\\pm".$2."\$";
	      }
	  }
	$LatexFile .= "}\n\\end{figure}\n";
	if ($LatexParity == 2)
	  {
	    &WriteLatexFile($LatexFile, $NbrLatexFile);
	    $LatexFile = "";
	    ++$NbrLatexFile;
	    $LatexParity = 1;
	  }
	else
	  {
	    $LatexParity = 2;
	  }
	close (INFILE);
	unlink("fit.log");
      }
    $TmpLa++;
  }

unless (open(OUTFILE, ">fermions_".$InteractionName.".ent"))
  {
    die ("can't create fermions_".$InteractionName.".ent\n");
  }
print OUTFILE $SLaResults;
close (OUTFILE);

unless (open(OUTFILE, ">fermions_".$InteractionName.".ent.p"))
  {
    die ("can't create fermions_".$InteractionName.".ent.p\n");
  }
print OUTFILE "set xrange [0:]
set terminal postscript landscape enhanced eps \"Helvetica\" 22
set output \"fermions_".$InteractionName.".ent.eps\"
set xlabel \"{l_a}^{1/2}\" font \"Helvetica,26\"
set ylabel \"S_{l_a}\" font \"Helvetica,26\"
set size 1.0, 1.0
f1(x)= a0 +  a1 * x
fit f1(x) \"fermions_".$InteractionName.".ent\" using 2:3 via a0,a1
plot \"fermions_".$InteractionName.".ent\" using 2:3 notitle with points 31, \"fermions_".$InteractionName.".ent\" using 2:3:4 notitle with errorbars, f1(x) notitle with lines lt 1 lw 2
";
close (OUTFILE);
`gnuplot fermions_$InteractionName.ent.p`;


unless (open(INFILE, "fit.log"))
  {
    die ("can't open fit.log\n");
  }

if ($LatexParity == 2)
  {
    &WriteLatexFile($LatexFile, $NbrLatexFile);
    $NbrLatexFile++;
  }

$LatexFile = "\\begin{figure}
\\includegraphics[width=12cm,angle=0]{fermions_".$InteractionName.".ent.eps}
\\caption{fit : ";

my $TmpLine;
foreach $TmpLine (<INFILE>)
  {
    chomp ($TmpLine);
    if ($TmpLine =~ /^a0\s+\=\s+(\-?\d+\.?\d*e?\-?\d*)\s+\+\/\-\s+(\-?\d+\.?\d*e?\-?\d*)/)
      {
	$LatexFile .= "\$a_0=".$1."\\pm".$2."\$ and ";
      }
    if ($TmpLine =~ /^a1\s+\=\s+(\-?\d+\.?\d*e?\-?\d*)\s+\+\/\-\s+(\-?\d+\.?\d*e?\-?\d*)/)
      {
	$LatexFile .= "\$a_1=".$1."\\pm".$2."\$";
      }
  }
$LatexFile .= "}\n\\end{figure}\n";

&WriteLatexFile($LatexFile, 0);


sub WriteLatexFile()
  {
    my $LatexFile = $_[0];
    my $FileIndex = $_[1];
    my $LatexFileName = "tmp".$FileIndex.".tex";
    unless (open (OUTFILE, ">".$LatexFileName))
      {
	die ("can't create ".$LatexFileName."\n");
      }
    
    print OUTFILE "\\documentclass[prb,aps,epsfig]{revtex4}

\\usepackage{graphicx,epsf}
\\usepackage{amsmath}
\\usepackage{bm}
\\usepackage{pstricks}
\\usepackage{psfrag}

\\begin{document}\n\n";

    print OUTFILE $LatexFile;
    print OUTFILE "
\\end{document}\n";

    close (OUTFILE);

    `latex $LatexFileName`;
    `latex $LatexFileName`;
    $LatexFileName = "tmp".$FileIndex.".dvi";
    my $PsFileName = "tmp".$FileIndex.".ps";
    `dvips $LatexFileName -o $PsFileName`;
    unlink ($LatexFileName);
  }
