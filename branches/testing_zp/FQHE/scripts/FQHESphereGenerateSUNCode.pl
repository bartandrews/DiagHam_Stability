#!/usr/bin/perl -w

use strict 'vars';

my $NValue = $ARGV[0];


my $FunctionDeclaration = "// evaluate Hilbert space dimension for fermions with SU(".$NValue.") spin
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value\n";
my $Index = 1;
while ($Index <= $NValue)
  {
    $FunctionDeclaration .= "// nbrN".$Index." = number of particles of type ".$Index."\n";
    $Index++;
  }
$FunctionDeclaration .= "// return value = Hilbert space dimension
long FermionSU".$NValue."ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz";
$Index = 1;
while ($Index <= $NValue)
  {
    $FunctionDeclaration .= ", int nbrN".$Index."\n";
    $Index++;
  }
$FunctionDeclaration .= ");";


my $FunctionBody = "{
  if ((nbrFermions < 0) || (totalLz < 0) || (lzMax < 0) ";
$Index = 1;
while ($Index <= $NValue)
  {
    $FunctionBody .= "|| (nbrN".$Index." < 0) || ((nbrN".$Index." - 1)> lzMax)";
    $Index++;
  }
$FunctionBody .= " || (lzMax < 0) || ((nbrFermions * lzMax - ((";
$Index = 1;
while ($Index <= $NValue)
  {
    $FunctionBody .= "(nbrN".$Index." * nbrN".$Index.")";
    if ($Index != $NValue)
      {
	 $FunctionBody .= " + ";
      }
    $Index++;
  }
$FunctionBody .= " - nbrFermions) >> 1)) < totalLz))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if (nbrFermions == 1) 
    if (lzMax >= totalLz)
      return 1l;
    else
      return 0l;
  unsigned long Tmp = 0l;
  if (nbrFermions > ".$NValue.")
    {
      Tmp += FermionSU".$NValue."ShiftedEvaluateHilbertSpaceDimension(nbrFermions - ".$NValue.", lzMax - 1, totalLz - (lzMax * ".$NValue.")";
$Index = 1;
while ($Index <= $NValue)
  {
    $FunctionBody .= ", nbrN".$Index." - 1";
    $Index++;
  }
$FunctionBody .= ");
    }\n";
my $Step = $NValue - 1;
while ($Step >= 2)
  {
    $FunctionBody .= "  if (nbrFermions >= ".$Step.")
    {
";
    $FunctionBody .= "      Tmp += FermionSU".$NValue."ShiftedEvaluateHilbertSpaceDimension(nbrFermions - ".$Step.", lzMax - 1, totalLz - (lzMax * ".$Step.")";
    $Index = 1;
    while ($Index <= $NValue)
      {
	$FunctionBody .= ", NbrN".$Index;
	$Index++;
      }
    $FunctionBody .= ");\n";
    $FunctionBody .= "    }
";
    $Step--;
  }
$FunctionBody .= "  return  (Tmp \n";
$Index = 1;
while ($Index <= $NValue)
  {
    $FunctionBody .= "	     + FermionSU".$NValue."ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax";
    for (my $Index2 = 1; $Index2 < $Index; $Index2++)
      {
	$FunctionBody .= ", nbrN".$Index2
      }
    $FunctionBody .= ", nbrN".$Index." - 1";
    for (my $Index2 = $Index + 1; $Index2 <= $NValue; $Index2++)
      {
	$FunctionBody .= ", nbrN".$Index2
      }
    $FunctionBody .= ")\n";
    $Index++;
  }
$FunctionBody .= "     + FermionSU".$NValue."ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, totalLz";
$Index = 1;
while ($Index <= $NValue)
  {
    $FunctionBody .= ", NbrN".$Index;
    $Index++;
  }
$FunctionBody .= "));
}
";

print $FunctionDeclaration;
print $FunctionBody;
