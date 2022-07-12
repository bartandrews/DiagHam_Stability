#!/usr/bin/perl -w

use strict 'vars';


if ((!(defined($ARGV[1]))) || (defined($ARGV[2])))
  {
    die ("usage: SvnTarball https://my-svn-depot project\n");
  }

my $SVNDepot = $ARGV[0];
my $ProjectName = $ARGV[1];

my $TmpDir = &MakeTmpDir();
chdir ($TmpDir);
if (&SVNCheckout($SVNDepot, $ProjectName) != 0)
  {
    print "error while retrieving project ".$ProjectName." at ".$SVNDepot."\n";
    rmdir ($TmpDir);
    die;
  }

&SVNClean();
my $TarballName = &CreateTarball($ProjectName);

`mv $TarballName ../`;
chdir ("..");
`rm -rf $TmpDir`;



# do a svn checkout for a given project
#
# $_[0] = depot location
# $_[1] = project name
# return value = 0 if no error occured, else 1

sub SVNCheckout()
  {
    my $SVNDepot = $_[0];
    my $ProjectName = $_[1];
    my $Command = "svn checkout ".$SVNDepot." ".$ProjectName;
    `$Command`;
    return 0;
  }


# clean any svn related files from a directory and its subdirectory
#

sub SVNClean()
  {
    my $TmpFile;
    foreach $TmpFile (<*>)
      {
	if ((-d $TmpFile) && ($TmpFile ne "..") && ($TmpFile ne "."))
	  {
	    chdir ($TmpFile);
	    &SVNClean();
	    chdir ("..");
	  }
      }
    if (-e ".svn")
      {
	`rm -rf .svn`;
      }
  }

# create a tarball of given directory
#
# $_[0] = directory to be tarballed
# return name = tarball name

sub CreateTarball()
  {
    my $DirectoryName = $_[0];
    my $Day;
    my $Month;
    my $Year;
    ($Day,$Month,$Year) = (gmtime(time))[3,4,5];
    if ($Day < 10)
      {
	$Day = "0".$Day;
      }
    $Year += 1900;
    $Month++;
    if ($Month < 10)
      {
	$Month = "0".$Month;
      }
    my $TarballName = $DirectoryName.$Year.$Month.$Day.".tar.gz";
    my $Command = "tar -czf ".$TarballName." ".$DirectoryName;
    `$Command`;    
    return $TarballName;
  }


# create a temporary directory
#

sub MakeTmpDir()
  {
    my $tmpDirName = "tmp".time();
    mkdir ($tmpDirName);
    return $tmpDirName;

  }
