#!/usr/bin/perl -w

use strict 'vars';

use XML::Simple;
use Data::Dumper;
use Getopt::Long;

my $HeaderFile = "";
my $FooterFile = "";
my $XMLFile = "";
my $Result = GetOptions ("header:s" => \$HeaderFile, "footer:s" => \$FooterFile, "xml=s" => \$XMLFile); 

if ($HeaderFile eq "")
  {
    $HeaderFile = "docs/built_in_programs/header.php";
  }
if ($FooterFile eq "")
  {
    $FooterFile = "docs/built_in_programs/footer.php";
  }

unless (open (INFILE, $HeaderFile))
  {
    die ("can't open ".$HeaderFile."\n");
  }
my $TmpLine;
foreach $TmpLine (<INFILE>)
  {
    print $TmpLine;
  }
close (INFILE);

print "echo \"";

my $XMLContent = XMLin($XMLFile);

if (defined($$XMLContent{'name'}))
  {
    print "<div class=\\\"programdoctitle\\\">".$$XMLContent{'name'}."<br /><br /></div>\n";
  }
if (defined($$XMLContent{'shortdesc'}))
  {
    print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Summary : </span>".&ConvertLatex2HTML($$XMLContent{'shortdesc'})."<br /><br /></div>\n";
  }
if (defined($$XMLContent{'location'}))
  {
    print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Location :</span> <span class=\\\"programdocdirectoryname\\\">".$$XMLContent{'location'}."</span><br /><br /></div>\n";
  }

if (defined($$XMLContent{'name'}))
  {
    print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Usage : </span> <span class=\\\"programdoccommand\\\">".$$XMLContent{'name'}." [options]</span><br /><br />\n";
  }

if (defined($$XMLContent{'optiongroup'}))
  {
    print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Options :</span><br />

<ul>\n";
    my $Key;
    my $OptionGroupRef = $$XMLContent{'optiongroup'};
    my @OptionGroupNames;
    if (defined($$XMLContent{'optiongroupsort'}))
      {
	@OptionGroupNames = split (/\s*\,\s*/, $$XMLContent{'optiongroupsort'});
      }
    else
      {
	@OptionGroupNames = sort (keys(%$OptionGroupRef));
      }
    foreach $Key (@OptionGroupNames)
      {
	my $Value = $$OptionGroupRef{$Key};
	print "<li>".&ConvertLatex2HTML($Key)." :
<ul>\n";
	my $OptionRef = $$Value{'option'};
	if (defined($$OptionRef{'name'}))
	  {
	    print PrintOptionInfo($OptionRef);
	  }
	else
	  {
	    my @OptionSortArray;
	    if (defined($$Value{'sort'}))
	      {
		@OptionSortArray = split (/\s*\,\s*/, $$Value{'sort'});
	      }
	    else
	      {
		@OptionSortArray = sort (keys (%$OptionRef));
	      }
	    my $OptionName;
	    foreach $OptionName (@OptionSortArray)
	      {
		my $OptionValue = $$OptionRef{$OptionName};
		$$OptionValue{'name'} = $OptionName;
		print PrintOptionInfo($OptionValue);
	      }
	    print "</li>\n";
	  }
	print "</ul></li>\n";
      }
    print "</ul>\n";
  }
if (defined($$XMLContent{'longdesc'}))
  {
    print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Description :</span><br />";
    print &ConvertLatex2HTML($$XMLContent{'longdesc'})."<br /><br /></div>\n";
  }
if (defined($$XMLContent{'accuracy'}))
  {
    print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Accuracy :</span><br />";
    print &ConvertLatex2HTML($$XMLContent{'accuracy'})."<br /><br /></div>\n";
  }
if (defined($$XMLContent{'remarks'}))
  {
    print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Remarks :</span><br />";
    print &ConvertLatex2HTML($$XMLContent{'remarks'})."<br /><br /></div>\n";
  }
if (defined($$XMLContent{'relatedprog'}))
  {
    print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Related programs and scripts :</span><br />

<ul>\n";
    my $RelatedProgram = $$XMLContent{'relatedprog'};
    if (defined($$RelatedProgram{'name'}))
      {
	print PrintRelatedProgramInfo($RelatedProgram);
      }
    else
      {
	my @RelatedProgramNames;
	if (defined($$XMLContent{'relatedprogsort'}))
	  {
	    @RelatedProgramNames = split (/\s*\,\s*/, $$XMLContent{'relatedprogsort'});
	  }
	else
	  {
	    @RelatedProgramNames = sort (keys(%$RelatedProgram));
	  }
	my $Key;
	foreach $Key (@RelatedProgramNames)
	  {
	    my $Value = $$RelatedProgram{$Key};
	    $$Value{'name'} = $Key;
	    print PrintRelatedProgramInfo($Value);
	  }
      }
    print "</ul></div>\n";
  }

if (defined($$XMLContent{'author'}))
  {
    my $AuthorInfo = $$XMLContent{'author'};
    if (defined($$AuthorInfo{'name'}))
      {
	print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Author : </span>".&PrintAuthorInfo($AuthorInfo)."<br /><br /></div>\n";
      }
    else
      {
	print "<div class=\\\"programdocsection\\\"><span class=\\\"programdocsectiontitle\\\">Authors : </span>\n";
	print "<ul>\n";
	my $Key;
	my $Value;
	while (($Key, $Value) = each (%$AuthorInfo))
	  {
	    $$Value{'name'} = $Key;
	    print "<li> ".PrintAuthorInfo($Value)."</li>\n";
	  }
	print "</ul><br /></div>\n";
      }
  }

print "\";";


unless (open (INFILE, $FooterFile))
  {
    die ("can't open ".$FooterFile."\n");
  }
foreach $TmpLine (<INFILE>)
  {
    print $TmpLine;
  }
close (INFILE);

# print author information in a string
#
# $_[0] = reference on the hash table containing author information
# return value = string corresponding to the author information

sub PrintAuthorInfo
  {
    my $AuthorInfo = $_[0];
    my $TmpString = "";
    if (defined($$AuthorInfo{'homepage'}))
      {
	$TmpString .= "<a href=\\\"".$$AuthorInfo{'homepage'}."\\\">".$$AuthorInfo{'name'}."</a>";	
      }
    else
      {
	$TmpString .= $$AuthorInfo{'name'};
      }
    return $TmpString;
  }

# print option information in a string
#
# $_[0] = reference on the hash table containing option information
# return value = string corresponding to the option information

sub PrintOptionInfo
  {
    my $OptionValue = $_[0];
    my $TmpString = "<li>";
    if (defined($$OptionValue{'short'}))
      {
	$TmpString .= "-".$$OptionValue{'short'}.", ";
      }
    $TmpString .= "--".$$OptionValue{'name'}." : ".&ConvertLatex2HTML($$OptionValue{'description'});
    if (defined($$OptionValue{'default'}))
      {
	$TmpString .= " (default value is set to ".$$OptionValue{'default'}.")";
      }
    $TmpString .= "</li>\n";
    return $TmpString;
  }

# print related program information in a string
#
# $_[0] = reference on the hash table containing related program information
# return value = string corresponding to the related program information

sub PrintRelatedProgramInfo
  {
    my $Value = $_[0];
    return "<li><span class=\\\"programdocrelatedprogram\\\">".$$Value{'name'}."</span><br />\n <span class=\\\"programdocrelatedprogramsection\\\">location</span> : <span class=\\\"programdocdirectoryname\\\">".$$Value{'location'}."</span><br />\n  <span class=\\\"programdocrelatedprogramsection\\\">usage</span> : ".&ConvertLatex2HTML($$Value{'usage'})."<br />\n  <span class=\\\"programdocrelatedprogramsection\\\">description</span> : ".&ConvertLatex2HTML($$Value{'description'})."<br /><br /></li>\n";
  }

# convert a Latex part of text to HTML
#
# $_[0] = string containing the text to convert
# return value = converted text

sub ConvertLatex2HTML
  {
    my $TmpString = $_[0];    
    $TmpString =~ s/\\\\//g;
    $TmpString =~ s/\n/\<br \/\>/g;
    $TmpString =~ s/\\\_/\_/g;
    $TmpString =~ s/\\\'(.)/\&$1acute\;/g;
    $TmpString =~ s/\\\`(.)/\&$1grave\;/g;
    $TmpString =~ s/\\\^(.)/\&$1circ\;/g;

    $TmpString =~ s/\\href\{([^\}]*)\}\{([^\}]*)\}/\<a href\=\\\"$2\\\"\>$1\<\/a\>/g;
    $TmpString =~ s/\\program\{([^\}]*)\}/\<span class\=\\\"programdocprogram\\\">$1\<\/span\>/g;
    $TmpString =~ s/\\option\{([^\}]*)\}/\<span class\=\\\"programdocoption\\\">\-$1\<\/span\>/g;
    $TmpString =~ s/\\longoption\{([^\}]*)\}/\<span class\=\\\"programdoclongoption\\\">\-\-$1\<\/span\>/g;
    $TmpString =~ s/\\filename\{([^\}]*)\}/\<span class\=\\\"programdocfilename\\\">$1\<\/span\>/g;
    $TmpString =~ s/\\directoryname\{([^\}]*)\}/\<span class\=\\\"programdocdirectoryname\\\">$1\<\/span\>/g;
    $TmpString =~ s/\\paper\{([^\}]*)\}/\<span class\=\\\"programdocpaper\\\">$1\<\/span\>/g;
    $TmpString =~ s/\\command\{([^\}]*)\}/\<span class\=\\\"programdoccommand\\\">$1\<\/span\>/g;
    $TmpString =~ s/\\commandarray\{([^\}]*)\}/\<br \/\>\<div class\=\\\"programdoccommandarray\\\">$1\<\/div\>/g;

    my $StartPos = 0;
    my $EndPos = 0;
    while ($EndPos >= 0)
      {	
	$StartPos = index ($TmpString, "\$", 0);
	$EndPos = index ($TmpString, "\$", $StartPos + 1);	
	if ($EndPos > $StartPos)
	  {
	    substr($TmpString, $StartPos, $EndPos - $StartPos + 1, "<span class=\\\"programdocformula\\\">".&ConvertLatexSmallFormula2HTML(substr($TmpString, $StartPos + 1, $EndPos - $StartPos - 1))."</span>");
	  }
      }

    $EndPos = 0;
    while ($EndPos >= 0)
      {	
	$StartPos = index ($TmpString, "\\begin{eqnarray}", 0);
	$EndPos = index ($TmpString, "\\end{eqnarray}", $StartPos + 16);	
	if ($EndPos > $StartPos)
	  {
	    substr($TmpString, $StartPos, $EndPos - $StartPos + 14, &ConvertLatexLongFormula2HTML(substr($TmpString, $StartPos + 16, $EndPos - $StartPos - 16)))
	  }
      }

    $TmpString =~ s/\<br \/\>/\<br \/\>\n/g;
    $TmpString =~ s/\\programwithlink\{([^\}]*)\}\{([^\}]*)\}/\".ShowProgramWithLink\(\$PathToRoot, \"$1\"\, \"$2\"\).\"/g;
    $TmpString =~ s/\\paperwithlink\{([^\}]*)\}\{([^\}]*)\}/\".ShowPaperLinkFromArxiv\(\$PathToRoot, \"$2\"\, \"$1\"\).\"/g;

    return $TmpString;
  }


# convert a Latex formula embedded in a $ $ to HTML
#
# $_[0] = string containing the formula to convert
# return value = converted formula

sub ConvertLatexSmallFormula2HTML
  {
    my $TmpString = $_[0];
    $TmpString =~ s/\\alpha/\&alpha\;/g;
    $TmpString =~ s/\\Alpha/\&Alpha\;/g;
    $TmpString =~ s/\\beta/\&beta\;/g;
    $TmpString =~ s/\\Beta/\&Beta\;/g;
    $TmpString =~ s/\\gamma/\&gamma\;/g;
    $TmpString =~ s/\\Gamma/\&Gamma\;/g;
    $TmpString =~ s/\\delta/\&delta\;/g;
    $TmpString =~ s/\\Delta/\&Delta\;/g;
    $TmpString =~ s/\\epsilon/\&epsilon\;/g;
    $TmpString =~ s/\\Epsilon/\&Epsilon\;/g;
    $TmpString =~ s/\\zeta/\&zeta\;/g;
    $TmpString =~ s/\\Zeta/\&Zeta\;/g;
    $TmpString =~ s/\\eta/\&eta\;/g;
    $TmpString =~ s/\\Eta/\&Eta\;/g;
    $TmpString =~ s/\\theta/\&theta\;/g;
    $TmpString =~ s/\\Theta/\&Theta\;/g;
    $TmpString =~ s/\\iota/\&iota\;/g;
    $TmpString =~ s/\\Iota/\&Iota\;/g;
    $TmpString =~ s/\\kappa/\&kappa\;/g;
    $TmpString =~ s/\\Kappa/\&Kappa\;/g;
    $TmpString =~ s/\\lambda/\&lambda\;/g;
    $TmpString =~ s/\\Lambda/\&Lambda\;/g;
    $TmpString =~ s/\\mu/\&mu\;/g;
    $TmpString =~ s/\\Mu/\&Mu\;/g;
    $TmpString =~ s/\\nu/\&nu\;/g;
    $TmpString =~ s/\\Nu/\&nu\;/g;
    $TmpString =~ s/\\omicron/\&omicron\;/g;
    $TmpString =~ s/\\Omicron/\&Omicron\;/g;
    $TmpString =~ s/\\xi/\&xi\;/g;
    $TmpString =~ s/\\Xi/\&Xi\;/g;
    $TmpString =~ s/\\pi/\&pi\;/g;
    $TmpString =~ s/\\Pi/\&Pi\;/g;
    $TmpString =~ s/\\rho/\&rho\;/g;
    $TmpString =~ s/\\Rho/\&Rho\;/g;
    $TmpString =~ s/\\sigma/\&sigma\;/g;
    $TmpString =~ s/\\Sigma/\&Sigma\;/g;
    $TmpString =~ s/\\tau/\&Tau\;/g;
    $TmpString =~ s/\\upsilon/\&upsilon\;/g;
    $TmpString =~ s/\\Upsilon/\&\Upsilon;/g;
    $TmpString =~ s/\\chi/\&chi\;/g;
    $TmpString =~ s/\\Chi/\&Chi\;/g;
    $TmpString =~ s/\\psi/\&psi\;/g;
    $TmpString =~ s/\\Psi/\&Psi\;/g;
    $TmpString =~ s/\\phi/\&phi\;/g;
    $TmpString =~ s/\\Phi/\&Phi\;/g;
    $TmpString =~ s/\\omega/\&omega\;/g;
    $TmpString =~ s/\\Omega/\&Omega\;/g;

    $TmpString =~ s/\\pm/\&plusmn\;/g;

    $TmpString =~ s/\\\;/ /g;
    $TmpString =~ s/\\nonumber//g;

    $TmpString =~ s/\\frac\{([^\}]*)\}\{([^\}]*)\}/\($1\)\/\($2\)/g;
    $TmpString =~ s/\{([^\}]*)\}\^\{([^\}]*)\}/\($1\)\<sup\>$2\<\/sup\>/g;
    $TmpString =~ s/\^\{([^\}]*)\}/\<sup\>$1\<\/sup\>/g;
    $TmpString =~ s/\^(.)/\<sup\>$1\<\/sup\>/g;
    $TmpString =~ s/\_\{([^\}]*)\}/\<sub\>$1\<\/sub\>/g;
    $TmpString =~ s/\_(.)/\<sub\>$1\<\/sub\>/g;

    return $TmpString;
  }

# convert a Latex formula embedded in a equation array to HTML
#
# $_[0] = string containing the formula to convert
# return value = converted formula

sub ConvertLatexLongFormula2HTML
  {
    $_[0] =~ s/^[\s\n]*//g;
    $_[0] =~ s/[\s\n]*$//g;
    my @Lines = split (/\\\\/, $_[0]);
    my $TmpLine;
    my $TmpString = "<br /> <div class=\\\"programdocformulaarray\\\"><table style=\\\"margin-left: auto;margin-right: auto;\\\">\n";
    foreach $TmpLine (@Lines)
      {
	$TmpString .= "<tr>\n";
	$TmpLine =~ s/^[\s\n]*//g;
	$TmpLine =~ s/[\s\n]*$//g;
	my @Columns = split (/\&/, $TmpLine);
	my $TmpColumn;
	foreach $TmpColumn (@Columns)
	  {
	    $TmpString .= "<td>";	    
	    $TmpColumn =~ s/^[\s]*//g;
	    $TmpColumn =~ s/[\s]*$//g;
	    $TmpString .= &ConvertLatexSmallFormula2HTML($TmpColumn)."</td>"	    
	    
	  }
	$TmpString .= "</tr>\n";
      }
    $TmpString .= "</table></div>\n";
    return $TmpString;
  }
