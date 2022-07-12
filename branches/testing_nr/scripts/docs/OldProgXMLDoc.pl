#!/usr/bin/perl -w

use strict 'vars';

if (!defined($ARGV[0]))
  {
    die "usage: ProgXMLDoc path/program\n";
  }

if (!(-e $ARGV[0]))
  {
    die $ARGV[0]." does not exist\n";
  }
if (!(-x $ARGV[0]))
  {
    die $ARGV[0]." is not a program\n";
  }

my $Path = $ARGV[0];
my $ProgramName = $ARGV[0];
$ProgramName =~ s/^.*\/([^\/]*)$/$1/;
$Path =~ s/^\.\///;
$Path =~ s/\/[^\/]*$//;
my $Options = `$ARGV[0] --help`;
my @TmpLines = split (/\n/, $Options);
my $TmpLine;
print "<?xml version=\"1.0\" ?>
<program name=\"".$ProgramName."\">
  <location>".$Path."</location>
  <authors name=\"\">
    <email></email>
    <homepage></homepage>
  </authors>
  <optiongroupsort></optiongroupsort>
  <optiongroup name=\"misc options\">
  <optiongroup name=\"\">
";
foreach $TmpLine (@TmpLines)
  {
    chomp ($TmpLine);
    $TmpLine =~ s/^\s*//;
    if ($TmpLine =~ /^\-/)
      {
	my $OptionLong = $TmpLine;
	$OptionLong =~ s/ \: .*//;
	my $OptionShort = $OptionLong;
	$OptionShort =~ s/\-\-.*//;
	$OptionShort =~ s/^\-([a-z,A-Z])\,.*/$1/;
	$OptionLong =~ s/^\-[a-z,A-Z]\, //;
	$OptionLong =~ s/^\-\-//;
	$TmpLine =~ s/^.* \: //;
	my $DefaultValue = $TmpLine;
	$DefaultValue =~ s/.*(\(default value.*\)).*/$1/;
	$TmpLine  =~ s/\(default value.*\)//;
	print "  <option name=\"".$OptionLong."\">\n";
	if ($OptionShort ne "")
	  {
	    print "    <short>".$OptionShort."</short>\n";
	  }
	$TmpLine =~ s/\</\\&lt;/m;
	$TmpLine =~ s/\>/\\&gt;/m;
	$TmpLine =~ s/\s*$//;
	print "    <description>".$TmpLine."</description>\n";
	if ($DefaultValue ne $TmpLine)
	  {
	    $DefaultValue  =~ s/\(default value = (.*)\)/$1/;
	    print "    <default>".$DefaultValue."</default>\n";
	  }
	print "  </option>\n"
      }
  }
print "  </optiongroup>

  <shortdescription></shortdescription>

  <longdescription></longdescription>

  <accuracy></accuracy>

  <relatedprog name=\"\" type=\"\">
    <location></location>
    <usage></usage>
    <description></description>
  </relatedprog>

  <remarks></remarks>
</program>\n";
