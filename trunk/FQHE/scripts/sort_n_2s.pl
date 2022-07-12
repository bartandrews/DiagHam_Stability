#!/usr/bin/perl -w
use strict 'vars';

my $TmpFileName;
foreach $TmpFileName (<*>)
{
    if ($TmpFileName =~ /\_n\_\d*\_2s\_\d*\_/)
    {
	my $NValue = $TmpFileName;
	my $SValue = $TmpFileName;
	$NValue =~ s/.*\_n\_(\d*)\_.*/$1/;
	$SValue =~ s/.*\_2s\_(\d*)\_.*/$1/;
	if (!((-e "n_$NValue") && (-d "n_$NValue")))
	{
	    `mkdir n_$NValue`;
	}
	if (!((-e "n_$NValue/2s_$SValue") && (-d "n_$NValue/2s_$SValue")))
	{
	    `mkdir n_$NValue/2s_$SValue`;
	}
	`mv -f $TmpFileName n_$NValue/2s_$SValue`;
    }
}


