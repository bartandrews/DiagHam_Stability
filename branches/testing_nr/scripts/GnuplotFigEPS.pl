#!/usr/bin/perl -w

########################################
#  Copyright NGUYEN Duc Phuong (2005)  #
########################################

#  Requirement: the gnuplot script file must contain "set term fig"
#  Required packages: gnuplot 4.0 (other versions do not work), transfig, tetex
#  $ARGV[0] = name of the gnuplot script file as input
#  $ARGV[1] = name of the output file saved as eps format

$GnuplotFile = shift (@ARGV);
$EPSFile = shift (@ARGV);

# open the Gnuplot file and replace \ by \\, save in a temporary file
$GnuplotTmpFile = "/tmp/" . $GnuplotFile;
open (INFILE, $GnuplotFile)
  or die "Cannot open the gnuplot script file: $GnuplotFile \n"; 
open (GNUPLOT, ">$GnuplotTmpFile")
  or die "Cannot open the temporary gnuplot script file: $GnuplotTmpFile \n";
while (<INFILE>)
  {
    chomp;
    $_ =~ s/\s*$//; # to remove all redundant space at the end of each line
    $_ =~ s/\\/\\\\/g; # to replace \ at the end of the line by \\ so that gnuplot can recognize Latex symbols
    $_ =~ s/\\\\$/\\/; # to replace \\ by \ at the end of each line, in the case we have \ as break line flag in gnuplot script
    $_ =~ s/\\\\n/\\n/g; # to replace \\n by \n so Gnuplot can recognize a break line flag! Pay attention to Latex symbols such as \nu -> does not work
    if ($_ =~ /set out/)
      {
	@List = split (/\s+/, $_);
	$FigFile = $List[$#List];
	$FigFile =~ s/\"//g;
      }
    print GNUPLOT ($_ . "\n");
  }
close (INFILE); close (GNUPLOT);
if (!$FigFile)
  {
    print "The name for fig file is not defined in the gnuplot script. Will use default.fig instead\n";
    $FigFile = "default.fig";
  }

# run gnuplot from a script file to obtain a fig file
`gnuplot $GnuplotTmpFile > $FigFile`;

# open the Fig file and replace the 4 by 6 (special text flag), save in a temporary file
$FigTmpFile = "/tmp/" . $FigFile . ".tmp";
open (INFILE, $FigFile)
  or die "Cannot open the fig file: $FigFile";
open (FIG, ">$FigTmpFile")
  or die "Cannot open the temporary fig file: ;$FigTmpFile";
while ($Line = <INFILE>)
  {    
    if (($Line =~ /\\\\/) || ($Line =~ /\$/))
      {	
	chomp ($Line);	
	@List = split (/\s/, $Line);	
	$List[9] = 6;	 
	$Line = join (" ", @List) . "\n";
      }
    print FIG ($Line);
  }

close (INFILE); close (FIG);

# create pstex and pstex_t files
$PstexFile = "/tmp/" . $EPSFile . ".pstex";
$Pstex_tFile = "/tmp/" . $EPSFile . ".pstex_t";
`fig2dev -L pstex $FigTmpFile $PstexFile`;
`fig2dev -L pstex_t -p $PstexFile $FigTmpFile $Pstex_tFile`;

# create the Latex file
$LatexFile = "/tmp/" . $EPSFile . ".tex";
open (LATEX, ">$LatexFile") 
  or die "Cannot open the Latex file: $LatexFile";
print LATEX (
"\\documentclass{article}
\\usepackage{epsfig}
\\usepackage{color}
\\setlength{\\textwidth}{100cm}
\\setlength{\\textheight}{100cm}
\\begin{document}
\\pagestyle{empty}
\\input{$Pstex_tFile}
\\end{document}"
);
close (LATEX);

# compile the Latex file and create the corresponding eps file
`cd /tmp && latex $LatexFile`;
`cd /tmp && latex $LatexFile`;
$DVIFile = $EPSFile . ".dvi";
`cd /tmp && dvips -E -o $EPSFile $DVIFile`;
`mv /tmp/$EPSFile .`;

#unlink $GnuplotTmpFile;
#unlink $FigFile;
#unlink $FigTmpFile;
unlink $PstexFile;
unlink $Pstex_tFile;
unlink $LatexFile;



