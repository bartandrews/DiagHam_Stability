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
my $OptionText = `$ARGV[0] --help`;
my @OptionGroups;

print "<?xml version=\"1.0\"   encoding=\"ISO-8859-1\"?>
<program name=\"".$ProgramName."\">
  <location>".$Path."</location>
  <authors name=\"\">
    <email></email>
    <homepage></homepage>
  </authors>
  <optiongroupsort>";

&ParseHelp($OptionText, \@OptionGroups);

my $Pos = 0;
while ($Pos < $#OptionGroups)
  {
    print $OptionGroups[$Pos]->GetName().", ";
    $Pos++;
  }
print $OptionGroups[$Pos]->GetName()."</optiongroupsort>\n";
my $TmpOptionGroup;
foreach $TmpOptionGroup (@OptionGroups)
  {
    print $TmpOptionGroup->ShowInformationXML();
  }

print "
  <shortdescription></shortdescription>

  <longdescription></longdescription>

  <accuracy></accuracy>

  <relatedprogsort></relatedprogsort>

  <relatedprog name=\"\" type=\"\">
    <location></location>
    <usage></usage>
    <description></description>
  </relatedprog>

  <remarks></remarks>
</program>\n";



# parse help returned by a program
#
# $_[0] = help message displayed by the program
# $_[1] = reference on an ordered array where option groups will be stored
# $_[2] = reference on a hash table where an option array will ber stored for each option group

sub ParseHelp
  {
    my @TmpLines = split (/\n/, $_[0]);
    my $OptionGroups = $_[1];
    my $TmpLine;
    my $Pos = 0;
    while (($Pos <= $#TmpLines) && !($TmpLines[$Pos] =~ /^Options\:/))
      {
	$Pos++;
      }
    $Pos++;
    while ($Pos <= $#TmpLines)
      {
	chomp ($TmpLines[$Pos]);
	$TmpLines[$Pos] =~ s/\:\s*$//;
	my $TmpOptionGroupName = $TmpLines[$Pos];
	if (($Pos <= $#TmpLines) && ($TmpOptionGroupName ne ""))
	  {
	    $Pos++;	    
	    my $TmpOptionGroup =  OptionGroup->new();
	    $TmpOptionGroup->SetName($TmpOptionGroupName);
	    chomp ($TmpLines[$Pos]);
	    $TmpLines[$Pos] =~ s/^\s//;
	    $TmpLines[$Pos] =~ s/\s*//;
	    while (($Pos <= $#TmpLines) && ($TmpLines[$Pos] ne ""))
	      {
		my $TmpOption =  Option->new();
		$TmpOption->ParseFromString($TmpLines[$Pos]);
		$TmpOptionGroup->AddOption($TmpOption);
		$Pos++;
		if ($Pos <= $#TmpLines)
		  {
		    chomp ($TmpLines[$Pos]);
		    $TmpLines[$Pos] =~ s/^\s//;
		    $TmpLines[$Pos] =~ s/\s*//;
		  }
	      }
	    push (@$OptionGroups, $TmpOptionGroup);
	  }
	else
	  {
	    $Pos++;
	  }
      }
  }


{package OptionGroup;

 # constructor
 #

 sub new
   {
     my $Class = shift;
     my $Self  = {};
     $Self->{NAME} = undef;
     $Self->{OPTIONS} = [];
     bless ($Self, $Class);
     return $Self;
   }
 
 # get an ordered comma-separated list of options
 #
 # return value = string containing ordered comma-separated list of options
 
 sub GetOptionOrderedList
   {
     my $Self = shift;
     my $String = "";
     if ($#{@{$Self->{OPTIONS}}} < 0)
       {
	 return $String;
       }
     my $Pos = 0;
     while ($Pos < $#{@{$Self->{OPTIONS}}})
       {
	 $String .= ${@{$Self->{OPTIONS}}}[$Pos]->GetName().", ";
	 $Pos++;
       }
     $String .= ${@{$Self->{OPTIONS}}}[$Pos]->GetName();
     return $String;
   }
 
 # display option group information using XML format
 #
 # return value = string containing option group information

 sub ShowInformationXML
   {
     my $Self = shift;
     my $String = "";
     $String .= "  <optiongroup name=\"".$Self->{NAME}."\">\n";
     $String .= "    <optiongroupsort>".$Self->GetOptionOrderedList()."</optiongroupsort>\n";
     my $TmpOption;
     foreach $TmpOption (@{$Self->{OPTIONS}})
       {
	 $String .= $TmpOption->ShowInformationXML();
       }
     $String .= "  </optiongroup>\n";
     return $String;
   }
 
 # add an option to the option group
 #
 
 sub AddOption
   {
     my $Self = shift;
     push (@{$Self->{OPTIONS}}, $_[0]);
   }

 # get name of the option group
 #
 # return value = option group name

 sub GetName
   {
     my $Self = shift;
     return $Self->{NAME};
   }

 # set name of the option group
 #
 # $_[0] = option group name

 sub SetName
   {
     my $Self = shift;
     $Self->{NAME} = $_[0];
   }

}

{package Option;

 # constructor
 #
 
 sub new
   {
     my $Class = shift;
     my $Self  = {};
     $Self->{NAME} = undef;
     $Self->{SHORTCUT} = undef;
     $Self->{DESCRIPTION} = undef;
     $Self->{DEFAULTVALUE} = undef;
     bless ($Self, $Class);
     return $Self;
   }
 

 # parse option parameter from a string
 #
 
 sub ParseFromString
   {
     my $Self = shift;
     my $TmpLine = $_[0];
     chomp ($TmpLine);
     $TmpLine =~ s/^\s*//;
     if ($TmpLine =~ /^\-/)
       {
	 $Self->{NAME} = $TmpLine;
	 $Self->{NAME} =~ s/ \: .*//;
	 $Self->{SHORTCUT} = $Self->{NAME};
	 $Self->{SHORTCUT} =~ s/\-\-.*//;
	 $Self->{SHORTCUT} =~ s/^\-([a-z,A-Z])\,.*/$1/;
	 $Self->{NAME} =~ s/^\-[a-z,A-Z]\, //;
	 $Self->{NAME} =~ s/^\-\-//;
	 $TmpLine =~ s/^.* \: //;
	 $Self->{DEFAULTVALUE} = $TmpLine;
	 $Self->{DEFAULTVALUE} =~ s/.*\(.*default value \= (.*)\).*/$1/;
	 $TmpLine  =~ s/\(default value \=.*\)//;
	 $TmpLine  =~ s/\, default value \=.*\)/\)/;
	 $TmpLine =~ s/\</\\&lt;/m;
	 $TmpLine =~ s/\>/\\&gt;/m;
	 $TmpLine =~ s/\s*$//;
	 $Self->{DESCRIPTION} = $TmpLine;
	 if (($Self->{DEFAULTVALUE} eq "") || ($Self->{DEFAULTVALUE} eq $Self->{DESCRIPTION}))
	   {
	     $Self->{DEFAULTVALUE}  = undef;
	   }
       }
   }

 # display option information using XML format
 #
 
 sub ShowInformationXML
   {
     my $Self = shift;
     my $String = "";
     $String .= "    <option name=\"".$Self->{NAME}."\">\n";
     if ((defined($Self->{SHORTCUT})) && ($Self->{SHORTCUT} ne ""))
       {
	 $String .= "    <short>".$Self->{SHORTCUT}."</short>\n";
       }
     $String .= "    <description>".$Self->{DESCRIPTION}."</description>\n";
     if ((defined($Self->{DEFAULTVALUE})) && ($Self->{DEFAULTVALUE} ne ""))
       {
	 $String .= "    <default>".$Self->{DEFAULTVALUE}."</default>\n";
       }
     $String .= "    </option>\n";
     return $String;
   }

 # get name of the option group
 #
 # return value = option group name

 sub GetName
   {
     my $Self = shift;
     return $Self->{NAME};
   }

}
