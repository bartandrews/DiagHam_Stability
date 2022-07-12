#!/usr/bin/perl -w

use strict 'vars';

use XML::Simple;
use Getopt::Long;

my $HeaderFile = "";
my $FooterFile = "";
my $XMLFile = "";
my $Result = GetOptions ("header:s" => \$HeaderFile, "footer:s" => \$FooterFile, "xml=s" => \$XMLFile); 

if ($HeaderFile eq "")
  {
    $HeaderFile = "docs/built_in_programs/header.tex";
  }
if ($FooterFile eq "")
  {
    $FooterFile = "docs/built_in_programs/footer.tex";
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


my $XMLContent = XMLin($XMLFile);
my $Key;
my $Value;

if (defined($$XMLContent{'name'}))
  {
    print "\\begin{center}{\\Large{\\bf{\\underline{".$$XMLContent{'name'}."}:}}}\\end{center}\n\\vspace{0.5cm}\n";
  }
if (defined($$XMLContent{'shortdesc'}))
  {
    print "\\begin{flushleft}{\\bf{\\underline{Summary} : }}".$$XMLContent{'shortdesc'}."\\end{flushleft}\n\n";
  }
if (defined($$XMLContent{'location'}))
  {
    print "\\begin{flushleft}{\\bf{\\underline{Location} : }} \\directoryname{".$$XMLContent{'location'}."}\\end{flushleft}\n\n";
  }

if (defined($$XMLContent{'name'}))
  {
    print "\\begin{flushleft}{\\bf{\\underline{Usage} : }}\\command{".$$XMLContent{'name'}." [options]}\\end{flushleft}\n\n";
  }

if (defined($$XMLContent{'optiongroup'}))
  {
    print "\\begin{flushleft}{\\bf{\\underline{Options} : }}\\end{flushleft}

\\begin{description}\n";
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
	print "\\item[\\textbullet]{\\bf{".$Key." :}}\n
\\begin{description}\n";
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
	  }
	print "\\end{description}\n\n";
      }
    print "\\end{description}\n";
  }
if (defined($$XMLContent{'longdesc'}))
  {
    print "\\begin{flushleft}{\\bf{\\underline{Description} : }}\\end{flushleft}\n\n";
    print $$XMLContent{'longdesc'}."\\\\\n\n";
  }
if (defined($$XMLContent{'accuracy'}))
  {
    print "\\begin{flushleft}{\\bf{\\underline{Accuracy} : }}\\end{flushleft}\n\n";
    print $$XMLContent{'accuracy'}."\\\\\n\n";
  }
if (defined($$XMLContent{'remarks'}))
  {
    print "\\begin{flushleft}{\\bf{\\underline{Remarks} : }}\\end{flushleft}\n\n";
    print $$XMLContent{'remarks'}."\\\\\n";
  }
if (defined($$XMLContent{'relatedprog'}))
  {
    print "\\begin{flushleft}{\\bf{\\underline{Related programs and scripts} : }}\\end{flushleft}\n\n

\\begin{description}\n";
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
    print "\\end{description}\n";
  }

if (defined($$XMLContent{'author'}))
  {
    my $AuthorInfo = $$XMLContent{'author'};
    if (defined($$AuthorInfo{'name'}))
      {
	print "\\begin{flushleft}{\\bf{\\underline{Author} : }}".&PrintAuthorInfo($AuthorInfo)."\\end{flushleft}\n\n";
      }
    else
      {
	print "\\begin{flushleft}{\\bf{\\underline{Authors} : }}\\end{flushleft}\n";
	print "\\begin{description}\n";
	my $Key;
	my $Value;
	while (($Key, $Value) = each (%$AuthorInfo))
	  {
	    $$Value{'name'} = $Key;
	    print "\\item[\\textbullet] ".PrintAuthorInfo($Value)."\n";
	  }
	print "\\end{description}\n\n";
      }
  }


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
    my $TmpString = $$AuthorInfo{'name'};
    if ((defined($$AuthorInfo{'email'})) || (defined($$AuthorInfo{'homepage'})))
      {
	$TmpString .= " (";
	if (defined($$AuthorInfo{'email'}))
	  {
	    $TmpString .= $$AuthorInfo{'email'};
	    if (defined($$AuthorInfo{'homepage'}))
	      {
		$TmpString .= ", ".$$AuthorInfo{'homepage'};
	      }
	  }
	else
	  {
	    $TmpString .= $$AuthorInfo{'homepage'};
	  }
	 $TmpString .= ")";
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
    my $TmpString = "\\item[";
    if (defined($$OptionValue{'short'}))
      {
	$TmpString .= "-".$$OptionValue{'short'}.", ";
      }
    $TmpString .= "- -".$$OptionValue{'name'}." : ] ".$$OptionValue{'description'};
    if (defined($$OptionValue{'default'}))
      {
	$TmpString .= " (default value is set to ".$$OptionValue{'default'}.")";
      }
    $TmpString .= "\n";
    return $TmpString;
  }

# print related program information in a string
#
# $_[0] = reference on the hash table containing related program information
# return value = string corresponding to the related program information

sub PrintRelatedProgramInfo
  {
    my $Value = $_[0];
    return "\\item[\\textbullet]{\\bf{".$$Value{'name'}."}}\\\\\n  {\\underline{location}} : \\directoryname{".$$Value{'location'}."}\\\\\n  {\\underline{usage}} : {\\it {".$$Value{'usage'}."}}\\\\\n  {\\underline{description}} : ".$$Value{'description'}."\\\\\n";
  }
