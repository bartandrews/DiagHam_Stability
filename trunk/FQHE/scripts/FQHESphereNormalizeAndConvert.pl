#!/usr/bin/perl -w

my $InputFile = $ARGV[0];
my $UnsqueezedBasis = "/home/regnault/development/Physics/DiagHam/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereConvertHaldaneBasis";
my $NormalizeProgram = "/home/regnault/development/Physics/DiagHam/build/FQHE/src/Programs/FQHEOnSphere/FQHESphereUnnormalizeState";

$InputFile =~ /^([^\_]+)\_unnormalized\_haldane\_(.*)\_n\_(\d+)\_2s\_(\d+)(\_lz\_\d+\.\d+\.vec)$/;
my $Statistics = $1;
my $InteractionName = $2;
my $NbrParticles = $3;
my $LzMax = $4;
my $Suffix = $5;
if ((-e $InputFile) && (defined($Statistics)))
  {
    if (!(-e $InteractionName."_n_".$NbrParticles."_2s_".$LzMax.".dat"))
      {
	die "can't find root description file ".$InteractionName."_n_".$NbrParticles."_2s_".$LzMax.".dat\n";
      }
    my $Output = $Statistics."_haldane_".$InteractionName."_n_".$NbrParticles."_2s_".$LzMax.$Suffix;
    my $Command = $NormalizeProgram." -i ".$InputFile." --haldane --reference-file ".$InteractionName."_n_".$NbrParticles."_2s_".$LzMax.".dat --normalize -o ".$Output;
    if (-e $InteractionName."_n_".$NbrParticles.".hil")
      {
	$Command .= " --load-hilbert ".$InteractionName."_n_".$NbrParticles.".hil";
      }
    system ($Command);
    if (!(-e $Output))
      {
	die "error while normalizing state ".$InputFile."\n";
      }
    my $Output2 = $Statistics."_".$InteractionName."_n_".$NbrParticles."_2s_".$LzMax.$Suffix;
    $Command = $UnsqueezedBasis." ".$Output." --reference-file ".$InteractionName."_n_".$NbrParticles."_2s_".$LzMax.".dat -o".$Output2;
    system ($Command);
    if (!(-e $Output2))
      {
	die "error while converting state ".$Output."\n";
      }    
  }
