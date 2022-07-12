#!/usr/bin/perl -w

use strict 'vars';

use XML::Simple;
use Getopt::Long;

my $HeaderFile = "";
my $FooterFile = "";
my $XMLFile = "";
my $PDFFile = "";
my $PathToXMLDoc2Tex = "";
my $Result = GetOptions ("header:s" => \$HeaderFile, "footer:s" => \$FooterFile, "xml=s" => \$XMLFile, "pdf:s" => \$PDFFile, "path:s" => \$PathToXMLDoc2Tex); 

unless (-e $XMLFile)
  {
    die ("can't find $XMLFile\n");
  }


if ($HeaderFile eq "")
  {
    $HeaderFile = "docs/built_in_programs/header.tex";
  }
if ($FooterFile eq "")
  {
    $FooterFile = "docs/built_in_programs/footer.tex";
  }
if ($PathToXMLDoc2Tex eq "")
  {
    $PathToXMLDoc2Tex = "scripts/docs";
  }
my $TruncatedXMLFileName = $XMLFile;
$TruncatedXMLFileName =~ s/.*\///g;
if ($PDFFile eq "")
  {    
    $PDFFile = $TruncatedXMLFileName;
    $PDFFile =~ s/\.xml/\.pdf/i;
  }

my $TmpLatexFileName = $TruncatedXMLFileName;
$TmpLatexFileName =~ s/\.xml$//i;
$TmpLatexFileName .= "_tmp_".time().".tex";
my $TmpLatexFile = `$PathToXMLDoc2Tex/XMLDoc2Tex.pl --xml=$XMLFile --header=$HeaderFile --footer=$FooterFile`;
open (OUTFILE, ">$TmpLatexFileName");
print OUTFILE $TmpLatexFile;
close (OUTFILE);
`pdflatex $TmpLatexFileName`;
$TmpLatexFileName =~ s/\.tex$/\.pdf/i;
`mv -f $TmpLatexFileName $PDFFile`;
$TmpLatexFileName =~ s/\.pdf$/\.\*/i;
`rm -f $TmpLatexFileName`;
