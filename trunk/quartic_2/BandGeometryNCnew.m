(* ::Package:: *)

2rt


(* ::Section:: *)
(*Public function declarations -- finish descriptions*)


BeginPackage["BandGeometryNCnew`"]


VersionString[]= "BandGeometry-nocompile v2016-6-29";
If[$Notebooks,
	Print[VersionString[]]
];
ModelList[]={
	"Haldane","Hofstadter","HofPyLandau",
	"HofPySymmetric","HofSymmetric","HofGeneral","p2x4"
};


bzMinValue::usage = "bzMinValue[function,otherParameters] finds the minimum value that \[OpenCurlyQuote]function\[CloseCurlyQuote] takes over the 2D Brillouin zone. The first two arguments of \[OpenCurlyQuote]function\[CloseCurlyQuote] are assummed to be the x and y components of momentum, normalized to the range [0...2Pi]. Any remaining parameters can be passed to \[OpenCurlyQuote]function\[CloseCurlyQuote] by including them as extra arguments.";
bzMaxValue::usage ="bzMaxValue[function,otherParameters] finds the maximum value that \[OpenCurlyQuote]function\[CloseCurlyQuote] takes over the 2D Brillouin zone. The first two arguments of \[OpenCurlyQuote]function\[CloseCurlyQuote] are assummed to be the x and y components of momentum, normalized to the range [0...2Pi]. Any remaining parameters can be passed to \[OpenCurlyQuote]function\[CloseCurlyQuote] by including them as extra arguments.";
bzIntegrate::usage="";
bzRMSKnownMean::usage="bzRMSKnownMean[function,otherParameters,Mean] computes the RMS fluctuations of \[OpenCurlyQuote]function\[CloseCurlyQuote] over the 2D Brillouin zone when the function\[CloseCurlyQuote]s average is known and supplied as the third argument. The first two arguments of \[OpenCurlyQuote]function\[CloseCurlyQuote] are assummed to be the x and y components of momentum, normalized to the range [0...2Pi]. Any remaining parameters can be passed to \[OpenCurlyQuote]function\[CloseCurlyQuote] by including them as extra arguments.";
bzMeanAndRMS::usage="bzRMS[function,otherParameters] computes the RMS fluctuations of \[OpenCurlyQuote]function\[CloseCurlyQuote] over the 2D Brillouin zone. The first two arguments of \[OpenCurlyQuote]function\[CloseCurlyQuote] are assummed to be the x and y components of momentum, normalized to the range [0...2Pi]. Any remaining parameters can be passed to \[OpenCurlyQuote]function\[CloseCurlyQuote] by including them as extra arguments.";


latticeDim::usage = "";
aa::usage = "";
bb::usage = "";
cc::usage = "";
bzXRange::usage = "";
bzYRange::usage = "";
bzRegionFunction::usage = "";
bzAspectRatio::usage = "";
latticeFourierBasis::usage="";
identifyWallpaper::usage = "";
intJacobian::usage = "";
intSymFactor::usage = "";


blochN::usage="";
blochNConst::usage="";
blochDN::usage="";
blochN2::usage="";
NblochN::usage="";
NblochDN::usage="";
NblochN2::usage="";
NblochDN2::usage="";
nBands::usage = "nBands[\\[CloseCurlyDoubleQuote]Model Name\\[CloseCurlyDoubleQuote]] returns the number of bands in the corresponding tight-binding model. Currently implemented models are given by the command modelList.";
nParams::usage = "nBands[\\[CloseCurlyDoubleQuote]Model Name\\[CloseCurlyDoubleQuote]] returns the number of bands in the corresponding tight-binding model. Currently implemented models are given by the command modelList.";
lattice::usage = "lattice[\\[CloseCurlyDoubleQuote]Model Name\\[CloseCurlyDoubleQuote]] returns the Bravais lattice used by the corresponding tight-binding model. Currently implemented models are given by the command modelList.";
paramShortNames::usage="";
unitCellVol::usage="";
H::usage ="[\\[CloseCurlyDoubleQuote]Model Name\\[CloseCurlyDoubleQuote], {kx, ky}, {parameters...}] gives the band Hamiltonian at the corresponding point in momentum space.";
dH::usage ="[\\[CloseCurlyDoubleQuote]Model Name\\[CloseCurlyDoubleQuote], dir, {kx, ky}, {parameters...}] gives the gradient in the dir\[CloseCurlyQuote]th momentum coordinate of the band Hamiltonian at the corresponding point in momentum space.";
NH::usage ="[\\[CloseCurlyDoubleQuote]Model Name\\[CloseCurlyDoubleQuote], {kx, ky}, {parameters...}] gives the band Hamiltonian at the corresponding point in momentum space.";
NdH::usage ="[\\[CloseCurlyDoubleQuote]Model Name\\[CloseCurlyDoubleQuote], dir, {kx, ky}, {parameters...}] gives the gradient in the dir\[CloseCurlyQuote]th momentum coordinate of the band Hamiltonian at the corresponding point in momentum space.";
termsInH::usage="";
dTermsInH::usage="";
energyUpperBound::usage = "[\\[CloseCurlyDoubleQuote]Model Name\\[CloseCurlyDoubleQuote], {\[Theta]_, \[Theta]2_, r2_}]";
bandEnergy::usage = "[kx_?NumericQ, ky_?NumericQ, name_String, params_, band_: 1]";
bandGap::usage = "[name_, params_, band_: 1]";
bandwidth::usage = "[name_, params_, band_: 1]";
gapOverWidth::usage = "[name_, params_, band_: 1]";
bandTouchingRatio::usage="";
gsVectors::usage = "[name_String, k_, params_, nOccBands_, \[Delta]prec_: 0]";


chernNumber::usage = "[name_, params_, nGrid_: 4]";
restaBilinears::usage="restaFactors[kx_?NumericQ,ky_?NumericQ,name_,params_,nOccBands_Integer,\[Delta]prec_:0]";
curvature::usage = "[kx_?NumericQ, ky_?NumericQ, name_, params_, nBands_: 1, \[Delta]prec_: 0]";
metric::usage = "[kx_?NumericQ, ky_?NumericQ, name_, params_, nBands_: 1, \[Delta]prec_: 0]";
metric11::usage="";
metric12::usage="";
metric22::usage="";
trG::usage = "[kx_?NumericQ, ky_?NumericQ, name_, params_, nBands_: 1, \[Delta]prec_: 0]";
detG::usage = "[kx_?NumericQ, ky_?NumericQ, name_, params_, nBands_: 1, \[Delta]prec_: 0]";
trGInequality::usage = "[kx_?NumericQ, ky_?NumericQ, name_, params_, nBands_: 1, \[Delta]prec_: 0]";
detGInequality::usage = "[kx_?NumericQ, ky_?NumericQ, name_, params_, nBands_: 1, \[Delta]prec_: 0]";
normTr::usage = "";


rmsCurvature::usage = "[name_, params_]";
trRMSMetric::usage="";
bTimesTrG::usage="";
bTimesDetG::usage="";
trGTimesDetG::usage="";
covarianceBTrG::usage="";
covarianceBDetG::usage="";
covarianceTrGDetG::usage="";


paramGradientBrms::usage="";
paramGradientBrms2::usage="";
paramHessianBrms::usage="";
paramGradientTrGrms::usage="";


sl2Ranisotropy::usage="";
metricEccentricity::usage="";
gradCurvature::usage="";
covDerivCurvature::usage="";
metric11Fourier::usage="";
metric12Fourier::usage="";
metric22Fourier::usage="";
sqrtDetGFourier::usage="";
curvatureFourier::usgae="";


lineOfData::usage="";
parallelLineOfData::usage="";
curvatureFourier10::usage="";
curvatureFourier11::usage="";
metricFourier10::usage="";
sigmaBAndIneqs::usage="";
intMetricComponents::usage="";
covariancesBTrGDetG::usage="";


functionShortName::usage="";
BZIntegratedFunctions[]={
	"metric11","metric12","metric22",
	"trG","trGInequality","detG","detGInequality",
	"metricEccentricity","sl2Ranisotropy",
	"gradCurvature","covDerivCurvature"
};
Do[
	With[{
		fn=Symbol[h],
		fnName=functionShortName[Symbol[h]],
		intFn=Symbol["int"<>ToUpperCase[StringTake[h,1]]<>StringDrop[h,1]],
		rmsFn=Symbol["rms"<>ToUpperCase[StringTake[h,1]]<>StringDrop[h,1]],
		avgWrmsFn=Symbol["avgWrms"<>ToUpperCase[StringTake[h,1]]<>StringDrop[h,1]]
	},
	intFn::usage="";
	rmsFn::usage="";
	avgWrmsFn::usage="";
	],
	{h,BZIntegratedFunctions[]}
];


(* ::Text:: *)
(*Select from opts only the options applicable to the functions listed in heads. Optionally, replace a subset of the options in opts with those in override. The latter option lets user-specified options override custom defaults for, eg, bzIntegrate, which in turn override the built-in defaults for NIntegrate.*)


filterOptions[opts_List,heads_List,override_List:{}]:=
Sequence@@FilterRules[
	DeleteDuplicates[Join[override,opts]],
	Flatten[Options/@heads]
];


(* ::Text:: *)
(*We use the following symmetrization commands when defining the larger Hamiltonians below -- we then only need to type in the entries below the diagonal, which reduces the possiblity of typos. Note that when PadRight is applied to a list of list of different lengths, it appends sufficiently many zeros to each list in order to produce a rectangular array.*)


matrixSym[M_] := (#1 + Transpose[#1] & )[PadRight[M]]; 
matrixAntiSym[M_] := (#1 - Transpose[#1] & )[PadRight[M]]; 
hermitianSym[M_] := (#1 + ConjugateTranspose[#1] & )[PadRight[M]]; 
hermitianSkewSym[M_] := (#1 - ConjugateTranspose[#1] & )[PadRight[M]]; 


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Lattice data*)


(* ::Text:: *)
(*So far, we don't really make use of any facts about the basis for any of the crystal systems used, so for now just provide information on*)


(* ::Text:: *)
(*Convention: a are unit vectors spanning real-space unit cell; then b define reciprocal lattice unit cell.*)
(*Also, functions defining a Brillouin zone, etc., are maps from the vector {k_x, k_y, ...} to a vector angular variables on the torus, in [0,2\[Pi])^n (which is what the Hamiltonian function is expecting).*)


(* ::Subsection:: *)
(*Definitions applicable to all lattices*)


(* ::Text:: *)
(*reference for the following trick: tutorial/FunctionsThatRememberValuesTheyHaveFound*)


aa[lattice_String,j_Integer]:= aa[lattice,j]= aa[lattice][[j]];
bb[lattice_String,j_Integer]:= bb[lattice,j]= bb[lattice][[j]];

bb[lattice_String]:=
2Pi*Inverse[
	Transpose[aa[lattice]]
];

unitCellVol[lattice_String]:=Abs[Det[bb[lattice]]];


bzRegionFunction[{lattice_String,"UnitCell"},x_,y_,f_]:=
And@@Table[
	(0 <= aa[lattice,i].{x,y} <= 2Pi),
	{i, latticeDim[lattice]}
];

bzXRange[{lattice_String,"UnitCell"}]:= bzXRange[{lattice,"UnitCell"}]=
({Min[#],Max[#]}&)[
	First/@{
				{0,0},
				bb[lattice,1],
				bb[lattice,2],
				bb[lattice,1]+bb[lattice,2]
			}
];

bzYRange[{lattice_String,"UnitCell"}]:= bzYRange[{lattice,"UnitCell"}]=
({Min[#],Max[#]}&)[
	Last/@{
				{0,0},
				bb[lattice,1],
				bb[lattice,2],
				bb[lattice,1]+bb[lattice,2]
			}
];

bzAspectRatio[{lattice_String,"UnitCell"}]:= bzAspectRatio[{lattice,"UnitCell"}]= 
With[{
	bx= bzXRange[{lattice,"UnitCell"}],
	by= bzYRange[{lattice,"UnitCell"}]
},
Abs[by[[2]]-by[[1]]]/Abs[bx[[2]]-bx[[1]]]
];


(* ::Text:: *)
(*The Zoology gauge choice breaks the 6=fold rotation symmetry of the triangular lattice down to reflection over the k_x axis. In order to have that symmetry in our basis of Fourier modes, we need to take a2 and a3 as the lattice basis vectors (instead of a1 and a2). *)


latticeFourierBasis[pt_,latt_String,{m_,n_}]:=
With[{
	aa1=If[latt==="Triangular",aa[latt,3],aa[latt,1]]
},
Exp[I*pt.(m*aa1+n*aa[latt,2])]
];


(* ::Subsection:: *)
(*Square lattice*)


latticeDim["Square"]=2;

aa["Square"]=
(2Pi)*{
	{1, 0},
	{0, 1}
};


bzRegionFunction[{"Square","BZ"},x_,y_,f_]:=
And@@Table[
	(-1/2 <= bb["Square",i].{x,y} <= 1/2),
	{i,2}
];

bzXRange[{"Square","BZ"}]= {-1/2, 1/2};
bzYRange[{"Square","BZ"}]= {-1/2, 1/2};
bzAspectRatio[{"Square","BZ"}]= 1;


(* ::Subsection:: *)
(*Triangular Lattice*)


latticeDim["Triangular"]=2;

aa["Triangular"]=
(4Pi)/Sqrt[3]*{
	{1, 0},
	(1/2)*{-1, Sqrt[3]}
};

aa["Triangular",3]= -aa["Triangular",2]-aa["Triangular",1];
bb["Triangular",3]= -bb["Triangular",2]+bb["Triangular",1];


bzRegionFunction[{"Triangular","BZ"},x_,y_,f_]:=
And@@Table[
	(-1/2 <= bb["Triangular",i].{x,y} <= 1/2),
	{i,3}
];

bzXRange[{"Triangular","BZ"}]= {-1/Sqrt[3], 1/Sqrt[3]};
bzYRange[{"Triangular","BZ"}]= {-1/2, 1/2};
bzAspectRatio[{"Triangular","BZ"}]= Sqrt[3]/2;


(* ::Section::Closed:: *)
(*Wallpaper group data (used by integration)*)


(* ::Text:: *)
(*For the most part, conventions have been adapted from the diagrams in*)
(*http://commons.wikimedia.org/wiki/Wallpaper_group_diagrams*)
(*http://en.wikipedia.org/wiki/List_of_planar_symmetry_groups#Wallpaper_groups*)


(* ::Subsection:: *)
(*Definitions applicable for all symmetries*)


cc[latticeAndGroup_List,j_Integer]:= cc[latticeAndGroup,j]= cc[latticeAndGroup][[j]];

linearBasisChangeQ[latticeAndGroup_List]:=
Equal@@Dimensions[
	cc[latticeAndGroup]
];

intJacobian[latticeAndGroup_List]:=intJacobian[latticeAndGroup]=
intSymFactor[latticeAndGroup]*
If[
	linearBasisChangeQ[latticeAndGroup],
(*then*)
	Det[cc[latticeAndGroup]],
(*else*)
	{
	Det[
		cc[latticeAndGroup][[{1,2}]]
	],
	Det[
		cc[latticeAndGroup][[{1,3}]]
	],
	Det[
		cc[latticeAndGroup][[{3,2}]]
	]
	}
];


intSymFactor[{_String,"p1"}]= 1;
cc[{lattice_String,"p1"}]:=bb[lattice];


(* ::Subsection:: *)
(*Square lattice*)


intSymFactor[{"Square","p2"}]= 2;
cc[{"Square","p2"}]:=
DiagonalMatrix[{1,1/2}];

intSymFactor[{"Square","cmm"}]= 4;
cc[{"Square","cmm"}]:= 
{
	{1/2,0},
	{1/2,1/2},
	{0,-1/2}
};

intSymFactor[{"Square","p4"}] = 4;
cc[{"Square","p4"}]= (1/2)IdentityMatrix[2];

intSymFactor[{"Square","pmm"}] = intSymFactor[{"Square","p4"}];
cc[{"Square","pmm"}]= cc[{"Square","p4"}];

intSymFactor[{"Square","p4m"}] = 8;
cc[{"Square","p4m"}]=
{
	{1/2, 0},
	{1/4, 1/4},
	{-1/4,1/4}
};

intSymFactor[{"Square","p4g"}] = 8;
cc[{"Square","p4g"}]=
{
	{1/2, 0},
	{0, 1/2},
	{-1/4, -1/4}
};


(* ::Subsection:: *)
(*Triangular lattice*)


intSymFactor[{"Triangular","p2"}]= 2;
cc[{"Triangular","p2"}]:=
{
	{Sqrt[3]/2,0},
	{Sqrt[3]/2,1/2},
	{0,-1/2}
};

intSymFactor[{"Triangular","cm"}]= 2;
cc[{"Triangular","cm"}]=cc[{"Triangular","p2"}];

intSymFactor[{"Triangular","cmm"}]= 4;
cc[{"Triangular","cmm"}]:= 
{
	{Sqrt[3]/4,0},
	{Sqrt[3]/2,1/2},
	{-Sqrt[3]/4,-1/2}
};

intSymFactor[{"Triangular","p3"}] = 3;
cc[{"Triangular","p3"}]=
{
	{1/Sqrt[3],0},
	{ 1/(2 Sqrt[3]), 1/2}
};

intSymFactor[{"Triangular","p3m1"}] = 6;
cc[{"Triangular","p3m1"}]=
{
	{1/Sqrt[3],0},
	{1/(2Sqrt[3]),1/2},
	{-Sqrt[3]/4,-1/4}
};

intSymFactor[{"Triangular","p31m"}] = 6;
cc[{"Triangular","p31m"}]=
{
	{Sqrt[3]/4,1/4},
	{1/(2Sqrt[3]),1/2},
	{1/(4Sqrt[3]),-1/4}
};

intSymFactor[{"Triangular","p6"}]=intSymFactor[{"Triangular","p31m"}];
cc[{"Triangular","p6"}]=cc[{"Triangular","p31m"}];

intSymFactor[{"Triangular","p6m"}] = 12;
cc[{"Triangular","p6m"}]=
{
	{Sqrt[3]/4, 1/4},
	{1/(4Sqrt[3]),1/4},
	{-1/(2Sqrt[3]),0}
};


(* ::Subsection:: *)
(*Symmetry testing*)


approxEqual[values_List]:=
	Max[Abs/@Differences[values]]<=1.0*^-8;
approxEqual[{}]=True;

rotocenterTest1[k_,f_,{pt_,f0_}]:=
approxEqual[
	Append[
		f/@Rest[
			NestList[Dot[RotationMatrix[2Pi/k],#]&,pt,k-1]
		],
		f0
	]
];

rotocenterTest[k_,f_,ptsFs_]:=
And@@(rotocenterTest1[k,f,#]&/@ptsFs);

mirrorTest1[axis_,f_,{pt_,f0_}]:=
approxEqual[
	{f0, f[ReflectionMatrix[axis].pt]}
];

mirrorTest[axis_,f_,ptsFs_]:=
And@@(mirrorTest1[axis,f,#]&/@ptsFs);

identifyWallpaper::sym="Inconsistent symmetry behavior for function `1`.";


(* ::Text:: *)
(*NB we don't test for pgg properly*)


identifyWallpaper[f_,"Square",nSamples_:3]:=
With[{
	b=bb["Square"],
	samples={#,f[#]}&/@RandomReal[{-1,1},{nSamples,2}]
},
	Which[
		rotocenterTest[4,f,samples],
			Switch[{
					mirrorTest[b[[1]],f,samples],
					mirrorTest[b[[1]]+b[[2]],f,samples]
				},
				{True,True},"p4m",
				{True,False},Message[identifyWallpaper::sym,f]; "p1",
				{False,True},"p4g",
				_,"p4"
			],
		rotocenterTest[2,f,samples],
			Switch[{
					mirrorTest[b[[1]],f,samples],
					mirrorTest[b[[1]]+b[[2]],f,samples]
				},
				{True,True},Message[identifyWallpaper::sym,f]; "p1",
				{True,False},"pmm",
				{False,True},"cmm",
				_,"p2"
			],
		True,"p1"
	]
];


(* ::Text:: *)
(*Again, use a basis consisting of -b3, b1 because the gauge choice in the Zoology paper breaks sym. down to reflection over the x-axis*)


identifyWallpaper[f_,"Triangular",nSamples_:3]:=
With[{
	b={bb["Triangular",3],bb["Triangular",1]},
	samples={#,f[#]}&/@RandomReal[{-1,1},{nSamples,2}]
},
	Which[
		rotocenterTest[2,f,samples],
			Switch[{
					rotocenterTest[3,f,samples],
					mirrorTest[b[[2]]+b[[1]],f,samples]
				},
				{True,True},"p6m",
				{True,False},"p6",
				{False,True},"cmm",
				_,"p2"
			],
		rotocenterTest[3,f,samples],
			Switch[{
					mirrorTest[b[[2]]+b[[1]],f,samples],
					mirrorTest[b[[2]]-b[[1]],f,samples]
				},
				{True,True},Message[identifyWallpaper::sym,f]; "p1",
				{True,False},"p31m",
				{False,True},"p3m1",
				_,"p3"
			],
		mirrorTest[b[[2]]-b[[1]],f,samples],"cm",
		True,"p1"
	]
];


(* ::Section::Closed:: *)
(*Integration and minimization over the BZ*)


(* ::Subsection:: *)
(*Minimization and maximization*)


Options[bzMinMaxValue]=
{
	"Symmetry"->True,
	AccuracyGoal->4, 
	PrecisionGoal->4
};
bzMinMaxValue[minOrMax_Symbol,F_,modelParams_List,otherParams___,
	opts:OptionsPattern[{bzMinMaxValue,NMinValue}]]:=
With[{
	wallpaper=
	{
		lattice[modelParams],
		If[
			OptionValue["Symmetry"],
		(*then*)
			identifyWallpaper[
				(F[#,modelParams,otherParams]&),
				lattice[modelParams]
			],	
		(*else*)
			"p1"]
	}
},
If[
	linearBasisChangeQ[wallpaper],
(*then*)
	minOrMax[{
		Re@F[{x,y}.cc[wallpaper],modelParams,otherParams],
		0<=x<=1&&0<=y<=1
	},
	{x,y},
	filterOptions[{opts},{NMinValue}]
	],
(*else*)
	minOrMax[{
		Re@F[{x,y,x*y}.cc[wallpaper],modelParams,otherParams],
		0<=x<=1&&0<=y<=1
	},
	{x,y},
	filterOptions[{opts},{NMinValue}]
	]
]];

bzMinValue[F_Symbol,modelParams_List,otherParams___,opts:OptionsPattern[bzMinMaxValue]]:=
	bzMinMaxValue[NMinValue,F,modelParams,otherParams,opts];

bzMaxValue[F_Symbol,modelParams_List,otherParams___,opts:OptionsPattern[bzMinMaxValue]]:=
	bzMinMaxValue[NMaxValue,F,modelParams,otherParams,opts];


(* ::Subsection:: *)
(*Integration*)


(* ::Text:: *)
(*alternatively: Method->{"MultiDimensionalRule","Generators"->5,"SymbolicProcessing"->0},*)
(**)
(*Limit Max Recursion by default in order to limit time-consuming runaway evaluations... previous error messages mentioned failure to converge after 18 recursive bisections.'*)
(**)
(*We use the With[ ] statements in order to ensure that the integrand is as evaluated/simple as possible before calling NIntegrate, since we turned off the symbolic simplification that it does.*)


Options[bzIntegrate]=
{
	"BZ Area"->1,
	"Symmetry"->True,
	AccuracyGoal->6,
	PrecisionGoal->5,
	MaxRecursion->9,
	Method->{"GlobalAdaptive","SymbolicProcessing"->0}
};

bzIntegrate[F_,modelParams_List,FParams___,
	opts:OptionsPattern[{bzIntegrate,NIntegrate}]]:=
Module[{
	wallpaper=
	{
		lattice[modelParams],
		If[OptionValue["Symmetry"],
		(*then*)
			identifyWallpaper[
				(F[#,modelParams,FParams]&),
				lattice[modelParams]
			],	
		(*else*) "p1"
		]
	},x,y
},
If[
	NumericQ[OptionValue["BZ Area"]],
	(*then*) OptionValue["BZ Area"]/unitCellVol[modelParams],
	(*else*) 1
]*If[
	linearBasisChangeQ[wallpaper],
(*then*)
	With[{
		XY={x,y}.cc[wallpaper],
		OPTS=filterOptions[Options[bzIntegrate],{NIntegrate},{opts}]
	},
	intJacobian[wallpaper]*
	NIntegrate[
		F[XY,modelParams,FParams],
		{x,0,1},{y,0,1},
		OPTS
	]],
(*else*)
	With[{
		JAC=intJacobian[wallpaper].{1,x,y},
		XY={x,y,x*y}.cc[wallpaper],
		OPTS=filterOptions[Options[bzIntegrate],{NIntegrate},{opts}]
	},
	NIntegrate[
		JAC*F[XY,modelParams,FParams],
		{x,0,1},{y,0,1},
		OPTS
	]]
]];


bzRMSKnownMean[F_,modelParams_List,FParams___,meanOfF_Real,
	opts:OptionsPattern[{bzIntegrate,NIntegrate}]]:=
Sqrt[
	bzIntegrate[
		((F[#1,#2]-meanOfF)^2&),
		modelParams,FParams,"BZ Area"->1,opts
	]
];

bzMeanAndRMS[F_,modelParams_List,FParams___,
	opts:OptionsPattern[{bzIntegrate,NIntegrate}]]:=
With[{
	intF=bzIntegrate[F,modelParams,FParams,"BZ Area"->1,opts]
},
{intF, bzRMSKnownMean[F,modelParams,FParams,intF,opts]}
];


(* ::Section:: *)
(*Band Hamiltonians*)


(* ::Text:: *)
(*Formatting convention: since, in principle, we want to study different lattice models, we identify them with a string. This saves us the trouble of having to explicitly pass each of these functions and their gradients, etc. to the routines above. *)
(*A Hamiltonian is specified by a function of the form*)
(*H[ "Name", {Subscript[\[Theta], 1], Subscript[\[Theta], 2], Subscript[\[Theta], 3], ...}, {parameters}]*)
(*which should return a numerical matrix of size nBands["Name"].  Here "Name" is any identification string, {\[Theta]1...} is the coordinates of a point on the reciprocal space torus -- Note: our convention is that each component runs from 0 to 2\[Pi], regradless of the lattice or shape of the Brillouin zone, for compatibility with the "Zoology..." paper and in order to make doing BZ integrals a bit easier. Converting between these angular coordinates and the physical wavenumber vector is really only needed for making nice plots, and is done by the kToTorus function in the section defining lattices.*)
(*Finally, "params" is an arbitrary vector of any remaining couplings needed to specify the Hamiltonian.*)
(*The gradient of each Hamiltonian is computed automatically.*)


(* ::Subsection:: *)
(*Haldane model *)


(* ::Text:: *)
(*F. D. M. Haldane, Phys. Rev. Lett. 61, 2015 (1988). *)
(*REMARK: as defined above, our a["Triangular"] real space Bravais vectors are equal to Haldane's b[1], b[2]. The rescaled aa vectors used below are Haldane's a[1],a[2],a[3].*)


nBands["Haldane"]=2;
lattice["Haldane"]="Triangular";
paramShortNames["Haldane"]={"t1","t2","phi","m"};

blochN[k_List,{"Haldane",t1_,t2_,\[Phi]_,M_}]=
Module[{nn},
nn[1]= -(4Pi/3)bb["Triangular",2];
nn[2]= (4Pi/3)bb["Triangular",1];
nn[3]= -(4Pi/3)bb["Triangular",3];
{
	t1*Sum[Cos[k.nn[i]],{i,3}],
	t1*Sum[Sin[k.nn[i]],{i,3}],
	M-2t2*Sin[\[Phi]]Sum[Sin[k.aa["Triangular",i]],{i,3}]
}];

blochNConst[k_List,{"Haldane",t1_,t2_,\[Phi]_,M_}]:=
2t2*Cos[\[Phi]]Sum[Cos[k.aa["Triangular",i]],{i,3}];

blochN2[k_List,{"Haldane",t1_,t2_,\[Phi]_,M_}]:=
(M-2t2*Sin[\[Phi]]Sum[Sin[k.aa["Triangular",i]],{i,3}])^2+
t1^2*(3+2*Sum[Cos[k.aa["Triangular",i]],{i,3}]);


(* ::Subsection:: *)
(*Hofstadter Model*)


(* ::Text:: *)
(*q x1 magnetic unit cell, Landau gauge*)


nParams["Hofstadter"]=2; 
lattice["Hofstadter"]="Square";
nBands[{"Hofstadter",p_,q_}]:=q;


H[{kx_,ky_},{"Hofstadter",p_,q_Integer}]:=
SparseArray[{
		Band[{1,1}]->Table[
			2*Cos[2Pi*(ky+m*p/q)],
			{m,q}
		],
		Band[{2,1}]->Exp[2Pi*I*kx],
		Band[{1,2}]->Exp[-2Pi*I*kx],
		{q,1}->Exp[-2Pi*I*kx],
		{1,q}->Exp[2Pi*I*kx]
	},
	{q,q}
];


NdH[{kx_,ky_},{"Hofstadter",p_,q_Integer},derivs:{1..}]:=
With[{n=Length[derivs]},
SparseArray[{
		Band[{2,1}]->(2Pi*I)^n*Exp[2Pi*I*kx],
		Band[{1,2}]->(-2Pi*I)^n*Exp[-2Pi*I*kx],
		{q,1}->(-2Pi*I)^n*Exp[-2Pi*I*kx],
		{1,q}->(2Pi*I)^n*Exp[2Pi*I*kx]
	},
	{q,q}
]];

NdH[{kx_,ky_},{"Hofstadter",p_,q_Integer},derivs:{2..}]:=
With[{n=Length[derivs]},
SparseArray[{
		Band[{1,1}]->Table[
			(2Pi*I)^n*Exp[2Pi*I(ky+m*p/q)] + 
			(-2Pi*I)^n*Exp[-2Pi*I(ky+m*p/q)],
			{m,q}
		]
	},
	{q,q}
]];

(* All mixed derivatives vanish *)
NdH[{kx_,ky_},{"Hofstadter",p_,q_Integer},derivs_List]:=
SparseArray[{},{q,q}];


nParams["QHofstadter"]=3; 
lattice["QHofstadter"]="Square";
nBands[{"QHofstadter",p_,q_,t2_}]:=q;


H[{kx_,ky_},{"QHofstadter",p_,q_Integer,t2_}]:=

If[q>4,

	SparseArray[{
			Band[{1,1}]->Table[
				-2*Cos[2Pi*(ky+m*p/q)]-2*t2*Cos[4Pi*(ky+m*p/q)],
				{m,q}
			],
			Band[{2,1}]->-Exp[2Pi*I*kx],
			Band[{1,2}]->-Exp[-2Pi*I*kx],
			Band[{3,1}]->-t2*Exp[4Pi*I*kx],
			Band[{1,3}]->-t2*Exp[-4Pi*I*kx],
			{q,1}->-Exp[-2Pi*I*kx],
			{1,q}->-Exp[2Pi*I*kx],
			{q,2}->-t2*Exp[-4Pi*I*kx],
			{2,q}->-t2*Exp[4Pi*I*kx],
			{1,q-1}->-t2*Exp[4Pi*I*kx],
			{q-1,1}->-t2*Exp[-4Pi*I*kx]
		},
		{q,q}
	],
	If[q==4,
			SparseArray[{
				Band[{1,1}]->Table[
					-2*Cos[2Pi*(ky+m*p/q)]-2*t2*Cos[4Pi*(ky+m*p/q)],
					{m,q}
				],
				Band[{2,1}]->-Exp[2Pi*I*kx],
				Band[{1,2}]->-Exp[-2Pi*I*kx],
				Band[{3,1}]->-t2*Exp[4Pi*I*kx]-t2*Exp[-4Pi*I*kx],
				Band[{1,3}]->-t2*Exp[-4Pi*I*kx]-t2*Exp[4Pi*I*kx],
				{q,1}->-Exp[-2Pi*I*kx],
				{1,q}->-Exp[2Pi*I*kx]
			},
			{q,q}
		],
		If[q==3,
			SparseArray[{
					Band[{1,1}]->Table[
						-2*Cos[2Pi*(ky+m*p/q)]-2*t2*Cos[4Pi*(ky+m*p/q)],
						{m,q}
					],
					Band[{2,1}]->-Exp[2Pi*I*kx]-t2*Exp[-4Pi*I*kx],
					Band[{1,2}]->-Exp[-2Pi*I*kx]-t2*Exp[4Pi*I*kx],
					{q,1}->-Exp[-2Pi*I*kx]-t2*Exp[4Pi*I*kx],
					{1,q}->-Exp[2Pi*I*kx]-t2*Exp[-4Pi*I*kx]
				},
				{q,q}
			]
		]
	]		

];


NdH[{kx_,ky_},{"QHofstadter",p_,q_Integer,t2_},derivs:{1..}]:=
With[{n=Length[derivs]},
If[q>4,

	SparseArray[{
			Band[{2,1}]->-(2Pi*I)^n*Exp[2Pi*I*kx],
			Band[{1,2}]->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
			Band[{3,1}]->-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
			Band[{1,3}]->-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx],
			{q,1}->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
			{1,q}->-(2Pi*I)^n*Exp[2Pi*I*kx],
			{q,2}->-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx],
			{2,q}->-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
			{1,q-1}->-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
			{q-1,1}->-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx]
		},
		{q,q}
	],
	If[q==4,
		SparseArray[{
				Band[{2,1}]->-(2Pi*I)^n*Exp[2Pi*I*kx],
				Band[{1,2}]->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
				Band[{3,1}]->-t2*(4Pi*I)^n*Exp[4Pi*I*kx]-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx],
				Band[{1,3}]->-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx]-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
				{q,1}->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
				{1,q}->-(2Pi*I)^n*Exp[2Pi*I*kx]
			},
			{q,q}
		],
		If[q==3,
			SparseArray[{
				Band[{2,1}]->-(2Pi*I)^n*Exp[2Pi*I*kx]-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx],
				Band[{1,2}]->-(-2Pi*I)^n*Exp[-2Pi*I*kx]-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
				{q,1}->-(-2Pi*I)^n*Exp[-2Pi*I*kx]-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
				{1,q}->-(2Pi*I)^n*Exp[2Pi*I*kx]-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx]
			},
			{q,q}
			]
		]
	]		
]];

NdH[{kx_,ky_},{"QHofstadter",p_,q_Integer,t2_},derivs:{2..}]:=
With[{n=Length[derivs]},
SparseArray[{
		Band[{1,1}]-> Table[
			-(2Pi*I)^n*Exp[2Pi*I(ky+m*p/q)]- (-2Pi*I)^n*Exp[-2Pi*I(ky+m*p/q)]-t2*((4Pi*I)^n*Exp[4Pi*I(ky+m*p/q)] + (-4Pi*I)^n*Exp[-4Pi*I(ky+m*p/q)]),
			{m,q}]
	 },
	{q,q}
]];

(* All mixed derivatives vanish *)
NdH[{kx_,ky_},{"QHofstadter",p_,q_Integer,t2_},derivs_List]:=
SparseArray[{},{q,q}];


nParams["p2x4"]=3; 
lattice["p2x4"]="Square";
nBands[{"p2x4",p_,q_,t2_}]:=q;
(* *** *)
H[{kx_,ky_},{"p2x4",p_,q_Integer,t2_}]:=

If[q>4,

	SparseArray[{
			Band[{1,1}]->Table[
				-2*Cos[2Pi*(ky+m*p/q)]-2*t2*Cos[4Pi*(ky+m*p/q)],
				{m,q}
			],
			Band[{2,1}]->-Exp[2Pi*I*kx],
			Band[{1,2}]->-Exp[-2Pi*I*kx],
			{q,1}->-Exp[-2Pi*I*kx],
			{1,q}->-Exp[2Pi*I*kx]
		},
		{q,q}
	],
	If[q==4,
			SparseArray[{
				Band[{1,1}]->Table[
					-2*Cos[2Pi*(ky+m*p/q)]-2*t2*Cos[4Pi*(ky+m*p/q)],
					{m,q}
				],
				Band[{2,1}]->-Exp[2Pi*I*kx],
				Band[{1,2}]->-Exp[-2Pi*I*kx],
				{q,1}->-Exp[-2Pi*I*kx],
				{1,q}->-Exp[2Pi*I*kx]
			},
			{q,q}
		],
		If[q==3,
			SparseArray[{
					Band[{1,1}]->Table[
						-2*Cos[2Pi*(ky+m*p/q)]-2*t2*Cos[4Pi*(ky+m*p/q)],
						{m,q}
					],
					Band[{2,1}]->-Exp[2Pi*I*kx]-t2*Exp[-4Pi*I*kx],
					Band[{1,2}]->-Exp[-2Pi*I*kx]-t2*Exp[4Pi*I*kx],
					{q,1}->-Exp[-2Pi*I*kx]-t2*Exp[4Pi*I*kx],
					{1,q}->-Exp[2Pi*I*kx]-t2*Exp[-4Pi*I*kx]
				},
				{q,q}
			]
		]
	]		

];

(***)
NdH[{kx_,ky_},{"p2x4",p_,q_Integer,t2_},derivs:{1..}]:=
With[{n=Length[derivs]},
If[q>4,

	SparseArray[{
			Band[{2,1}]->-(2Pi*I)^n*Exp[2Pi*I*kx],
			Band[{1,2}]->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
			Band[{3,1}]->-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
			Band[{1,3}]->-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx],
			{q,1}->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
			{1,q}->-(2Pi*I)^n*Exp[2Pi*I*kx]

		},
		{q,q}
	],
	If[q==4,
		SparseArray[{
				Band[{2,1}]->-(2Pi*I)^n*Exp[2Pi*I*kx],
				Band[{1,2}]->-(-2Pi*I)^n*Exp[-2Pi*I*kx],

				{q,1}->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
				{1,q}->-(2Pi*I)^n*Exp[2Pi*I*kx]
			},
			{q,q}
		],
		If[q==3,
			SparseArray[{
				Band[{2,1}]->-(2Pi*I)^n*Exp[2Pi*I*kx]-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx],
				Band[{1,2}]->-(-2Pi*I)^n*Exp[-2Pi*I*kx]-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
				{q,1}->-(-2Pi*I)^n*Exp[-2Pi*I*kx]-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
				{1,q}->-(2Pi*I)^n*Exp[2Pi*I*kx]-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx]
			},
			{q,q}
			]
		]
	]		
]];
(****)
NdH[{kx_,ky_},{"p2x4",p_,q_Integer,t2_},derivs:{1..}]:=
With[{n=Length[derivs]},
If[q>4,

	SparseArray[{
			Band[{2,1}]->-(2Pi*I)^n*Exp[2Pi*I*kx],
			Band[{1,2}]->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
			Band[{3,1}]->-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
			Band[{1,3}]->-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx],
			{q,1}->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
			{1,q}->-(2Pi*I)^n*Exp[2Pi*I*kx]

		},
		{q,q}
	],
	If[q==4,
		SparseArray[{
				Band[{2,1}]->-(2Pi*I)^n*Exp[2Pi*I*kx],
				Band[{1,2}]->-(-2Pi*I)^n*Exp[-2Pi*I*kx],

				{q,1}->-(-2Pi*I)^n*Exp[-2Pi*I*kx],
				{1,q}->-(2Pi*I)^n*Exp[2Pi*I*kx]
			},
			{q,q}
		],
		If[q==3,
			SparseArray[{
				Band[{2,1}]->-(2Pi*I)^n*Exp[2Pi*I*kx]-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx],
				Band[{1,2}]->-(-2Pi*I)^n*Exp[-2Pi*I*kx]-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
				{q,1}->-(-2Pi*I)^n*Exp[-2Pi*I*kx]-t2*(4Pi*I)^n*Exp[4Pi*I*kx],
				{1,q}->-(2Pi*I)^n*Exp[2Pi*I*kx]-t2*(-4Pi*I)^n*Exp[-4Pi*I*kx]
			},
			{q,q}
			]
		]
	]		
]];


NdH[{kx_,ky_},{"p2x4",p_,q_Integer,t2_},derivs:{2..}]:=
With[{n=Length[derivs]},
SparseArray[{
		Band[{1,1}]-> Table[
			-(2Pi*I)^n*Exp[2Pi*I(ky+m*p/q)]- (-2Pi*I)^n*Exp[-2Pi*I(ky+m*p/q)]-t2*((4Pi*I)^n*Exp[4Pi*I(ky+m*p/q)] + (-4Pi*I)^n*Exp[-4Pi*I(ky+m*p/q)]),
			{m,q}]
	 },
	{q,q}
]];

(* All mixed derivatives vanish *)
NdH[{kx_,ky_},{"p2x4",p_,q_Integer,t2_},derivs_List]:=
SparseArray[{},{q,q}];


(* ::Subsection::Closed:: *)
(*Hofstadter Model -- from Python code*)


posIndex[XY_,UxUy_]:=
With[{
	xy=MapThread[Mod,{XY,UxUy}]
},
	xy[[1]]+(First@UxUy)*xy[[2]]
];


(* ::Text:: *)
(*Landau gauge*)


nParams["HofPyLandau"]=3; 
lattice["HofPyLandau"]="Square";
nBands[{"HofPyLandau",ux_,uy_,nFlux_}]:= ux*uy;


hoppingPhaseLandau[KxKy_,XY_,{ux_,uy_},nFlux_]:=
Module[{
	numTrans,xy
},
numTrans=-MapThread[Quotient,{XY,{ux,uy}}];
xy=MapThread[Mod,{XY,{ux,uy}}];
0*KxKy.numTrans-2Pi*xy[[2]]*numTrans[[1]]nFlux/uy
];

dHoppingPhaseLandau[XY_,{ux_,uy_}]:=
-0*MapThread[Quotient,{XY,{ux,uy}}];

H[{kx_,ky_},{"HofPyLandau",ux_Integer,uy_Integer,nFlux_}]:=
Module[{i,j,\[CapitalDelta],
	initIndex,finalIndex,phase,phaseFactors,
	lattice\[Delta]s={{1,0},{-1,0},{0,1},{0,-1}},
	tmpHamiltonian=ConstantArray[0.0,{ux*uy,ux*uy}]
},
Do[
	phaseFactors=Pi*{0,2i}*nFlux/(ux*uy);
	initIndex=posIndex[{i,j},{ux,uy}]+1;
	Do[
		finalIndex=posIndex[{i,j}+\[CapitalDelta],{ux,uy}]+1;
		phase=hoppingPhaseLandau[2Pi*{kx,ky},{i,j}+\[CapitalDelta],{ux,uy},nFlux];
		If[initIndex>=finalIndex,
			tmpHamiltonian[[initIndex,finalIndex]]+=
				Exp[I*(phase-\[CapitalDelta].(2*Pi*{kx,ky}+phaseFactors))];
			If[initIndex!=finalIndex,
				tmpHamiltonian[[finalIndex,initIndex]]+=
					Exp[-I*(phase-\[CapitalDelta].(2*Pi*{kx,ky}+phaseFactors))]
			];
		];,
		{\[CapitalDelta],lattice\[Delta]s}
	];,
	{j,0,uy-1},
	{i,0,ux-1}
];
Return[tmpHamiltonian];
];


NdH[{kx_,ky_},{"HofPyLandau",ux_Integer,uy_Integer,nFlux_},{1}]:=
Module[{i,j,\[CapitalDelta],
	initIndex,finalIndex,phase,phaseFactors,
	lattice\[Delta]s={{1,0},{-1,0},{0,1},{0,-1}},
	tmpHamiltonian=ConstantArray[0.0,{ux*uy,ux*uy}]
},
Do[
	phaseFactors=Pi*{0,2i}*nFlux/(ux*uy);
	initIndex=posIndex[{i,j},{ux,uy}]+1;
	Do[
		finalIndex=posIndex[{i,j}+\[CapitalDelta],{ux,uy}]+1;
		phase=hoppingPhaseLandau[2Pi*{kx,ky},{i,j}+\[CapitalDelta],{ux,uy},nFlux];
		If[initIndex>=finalIndex,
			tmpHamiltonian[[initIndex,finalIndex]]+=
				(2Pi*I (dHoppingPhaseLandau[{i,j}+\[CapitalDelta],{ux,uy}]-\[CapitalDelta]).{1,0})*
				Exp[I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))];
			If[initIndex!=finalIndex,
				tmpHamiltonian[[finalIndex,initIndex]]+=
					(-2Pi*I (dHoppingPhaseLandau[{i,j}+\[CapitalDelta],{ux,uy}]-\[CapitalDelta]).{1,0})*
					Exp[-I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))]
			];
		];,
		{\[CapitalDelta],lattice\[Delta]s}
	];,
	{j,0,uy-1},
	{i,0,ux-1}
];
Return[tmpHamiltonian];
];

NdH[{kx_,ky_},{"HofPyLandau",ux_Integer,uy_Integer,nFlux_},{2}]:=
Module[{i,j,\[CapitalDelta],
	initIndex,finalIndex,phase,phaseFactors,
	lattice\[Delta]s={{1,0},{-1,0},{0,1},{0,-1}},
	tmpHamiltonian=ConstantArray[0.0,{ux*uy,ux*uy}]
},
Do[
	phaseFactors=Pi*{0,2i}*nFlux/(ux*uy);
	initIndex=posIndex[{i,j},{ux,uy}]+1;
	Do[
		finalIndex=posIndex[{i,j}+\[CapitalDelta],{ux,uy}]+1;
		phase=hoppingPhaseLandau[2Pi*{kx,ky},{i,j}+\[CapitalDelta],{ux,uy},nFlux];
		If[initIndex>=finalIndex,
			tmpHamiltonian[[initIndex,finalIndex]]+=
				(2Pi*I (dHoppingPhaseLandau[{i,j}+\[CapitalDelta],{ux,uy}]-\[CapitalDelta]).{0,1})*
				Exp[I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))];
			If[initIndex!=finalIndex,
				tmpHamiltonian[[finalIndex,initIndex]]+=
					(-2Pi*I (dHoppingPhaseLandau[{i,j}+\[CapitalDelta],{ux,uy}]-\[CapitalDelta]).{0,1})*
					Exp[-I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))]
			];
		];,
		{\[CapitalDelta],lattice\[Delta]s}
	];,
	{j,0,uy-1},
	{i,0,ux-1}
];
Return[tmpHamiltonian];
];


(* ::Text:: *)
(*Symmetric gauge*)


nParams["HofPySymmetric"]=3; 
lattice["HofPySymmetric"]="Square";
nBands[{"HofPySymmetric",ux_,uy_,nFlux_}]:=ux*uy;


hoppingPhaseSymmetric[KxKy_,XY_,{ux_,uy_},nFlux_]:=
Module[{
	numTrans,xy
},
numTrans=-MapThread[Quotient,{XY,{ux,uy}}];
xy=MapThread[Mod,{XY,{ux,uy}}];
0*KxKy.numTrans+
Pi*xy[[1]]*numTrans[[2]]nFlux/uy-Pi*xy[[2]]*numTrans[[1]]nFlux/ux
];

dHoppingPhaseSymmetric[XY_,{ux_,uy_}]:=
-0*MapThread[Quotient,{XY,{ux,uy}}];

H[{kx_,ky_},{"HofPySymmetric",ux_Integer,uy_Integer,nFlux_}]:=
Module[{i,j,\[CapitalDelta],
	initIndex,finalIndex,phase,phaseFactors,
	lattice\[Delta]s={{1,0},{-1,0},{0,1},{0,-1}},
	tmpHamiltonian=ConstantArray[0.0,{ux*uy,ux*uy}]
},
Do[
	phaseFactors=Pi*{-j,i}*nFlux/(ux*uy);
	initIndex=posIndex[{i,j},{ux,uy}]+1;
	Do[
		finalIndex=posIndex[{i,j}+\[CapitalDelta],{ux,uy}]+1;
		phase=hoppingPhaseSymmetric[2Pi*{kx,ky},{i,j}+\[CapitalDelta],{ux,uy},nFlux];
		If[initIndex>=finalIndex,
			tmpHamiltonian[[initIndex,finalIndex]]+=
				Exp[I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))];
			If[initIndex!=finalIndex,
				tmpHamiltonian[[finalIndex,initIndex]]+=
					Exp[-I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))]
			];
		];,
		{\[CapitalDelta],lattice\[Delta]s}
	];,
	{j,0,uy-1},
	{i,0,ux-1}
];
Return[tmpHamiltonian];
];


NdH[{kx_,ky_},{"HofPySymmetric",ux_Integer,uy_Integer,nFlux_},{1}]:=
Module[{i,j,\[CapitalDelta],
	initIndex,finalIndex,phase,phaseFactors,
	lattice\[Delta]s={{1,0},{-1,0},{0,1},{0,-1}},
	tmpHamiltonian=ConstantArray[0.0,{ux*uy,ux*uy}]
},
Do[
	phaseFactors=Pi*{-j,i}*nFlux/(ux*uy);
	initIndex=posIndex[{i,j},{ux,uy}]+1;
	Do[
		finalIndex=posIndex[{i,j}+\[CapitalDelta],{ux,uy}]+1;
		phase=hoppingPhaseSymmetric[2Pi*{kx,ky},{i,j}+\[CapitalDelta],{ux,uy},nFlux];
		If[initIndex>=finalIndex,
			tmpHamiltonian[[initIndex,finalIndex]]+=
				(2Pi*I (dHoppingPhaseSymmetric[{i,j}+\[CapitalDelta],{ux,uy}]-\[CapitalDelta]).{1,0})*
				Exp[I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))];
			If[initIndex!=finalIndex,
				tmpHamiltonian[[finalIndex,initIndex]]+=
					(-2Pi*I (dHoppingPhaseSymmetric[{i,j}+\[CapitalDelta],{ux,uy}]-\[CapitalDelta]).{1,0})*
					Exp[-I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))]
			];
		];,
		{\[CapitalDelta],lattice\[Delta]s}
	];,
	{j,0,uy-1},
	{i,0,ux-1}
];
Return[tmpHamiltonian];
];

NdH[{kx_,ky_},{"HofPySymmetric",ux_Integer,uy_Integer,nFlux_},{2}]:=
Module[{i,j,\[CapitalDelta],
	initIndex,finalIndex,phase,phaseFactors,
	lattice\[Delta]s={{1,0},{-1,0},{0,1},{0,-1}},
	tmpHamiltonian=ConstantArray[0.0,{ux*uy,ux*uy}]
},
Do[
	phaseFactors=Pi*{-j,i}*nFlux/(ux*uy);
	initIndex=posIndex[{i,j},{ux,uy}]+1;
	Do[
		finalIndex=posIndex[{i,j}+\[CapitalDelta],{ux,uy}]+1;
		phase=hoppingPhaseSymmetric[2Pi*{kx,ky},{i,j}+\[CapitalDelta],{ux,uy},nFlux];
		If[initIndex>=finalIndex,
			tmpHamiltonian[[initIndex,finalIndex]]+=
				(2Pi*I (dHoppingPhaseSymmetric[{i,j}+\[CapitalDelta],{ux,uy}]-\[CapitalDelta]).{0,1})*
				Exp[I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))];
			If[initIndex!=finalIndex,
				tmpHamiltonian[[finalIndex,initIndex]]+=
					(-2Pi*I (dHoppingPhaseSymmetric[{i,j}+\[CapitalDelta],{ux,uy}]-\[CapitalDelta]).{0,1})*
					Exp[-I*(phase-\[CapitalDelta].(2Pi*{kx,ky}+phaseFactors))]
			];
		];,
		{\[CapitalDelta],lattice\[Delta]s}
	];,
	{j,0,uy-1},
	{i,0,ux-1}
];
Return[tmpHamiltonian];
];


(* ::Subsection:: *)
(*General gauges*)


(* ::Text:: *)
(*Symmetric gauge, faster -- see symmetric gauge hof.nb*)


nParams["HofSymmetric"]=3; 
lattice["HofSymmetric"]="Square";
nBands[{"HofSymmetric",ux_,uy_,nFlux_}]:=ux*uy;


hofSymmetricXBlock[kx_,j_,{m_,n_}]:=
Normal@ SparseArray[
	Which[
		m==1,{
			Band[{1,1}]->Exp[-2Pi*I(kx-j (n+1)/(2m*n))]+Exp[2Pi*I(kx-j (n+1)/(2m*n))]
		},
		m==2,{
			Band[{1,2}]->Exp[-2Pi*I(kx-j/(2m*n))]+Exp[2Pi*I(kx-j (n+1)/(2m*n))],
			Band[{2,1}]->Exp[-2Pi*I(kx-j (n+1)/(2m*n))]+Exp[2Pi*I(kx-j/(2m*n))]
		},
		True,{
			Band[{1,2}]->Exp[-2Pi*I(kx-j/(2m*n))],
			Band[{m,1}]->Exp[-2Pi*I(kx-j (n+1)/(2m*n))],
			Band[{2,1}]->Exp[2Pi*I(kx-j/(2m*n))],
			Band[{1,m}]->Exp[2Pi*I(kx-j (n+1)/(2m*n))]
		}
	],
	{m,m}
];

DHofSymmetricXBlock[kx_,j_,{m_,n_}]:=
2Pi*I*Normal@ SparseArray[
	Which[
		m==1,{
			Band[{1,1}]->-Exp[-2Pi*I(kx-j (n+1)/(2m*n))]+Exp[2Pi*I(kx-j (n+1)/(2m*n))]
		},
		m==2,{
			Band[{1,2}]->-Exp[-2Pi*I(kx-j/(2m*n))]+Exp[2Pi*I(kx-j (n+1)/(2m*n))],
			Band[{2,1}]->-Exp[-2Pi*I(kx-j (n+1)/(2m*n))]+Exp[2Pi*I(kx-j/(2m*n))]
		},
		True,{
			Band[{1,2}]->-Exp[-2Pi*I(kx-j/(2m*n))],
			Band[{m,1}]->-Exp[-2Pi*I(kx-j (n+1)/(2m*n))],
			Band[{2,1}]->Exp[2Pi*I(kx-j/(2m*n))],
			Band[{1,m}]->Exp[2Pi*I(kx-j (n+1)/(2m*n))]
		}
	],
	{m,m}
];

hofSymmetricYBlock[ky_,ij_,{m_,n_}]:=
DiagonalMatrix@ Table[
	Exp[2Pi*I(ky+k*ij/(2m*n))],
	{k,0,m-1} 
];

H[{kx_,ky_},{"HofSymmetric",m_Integer,n_Integer,nFlux_}]:=
ArrayFlatten@ Table[
	Which[
		n==1, 
			hofSymmetricXBlock[kx,j*nFlux,{m,n}]+
			hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}]+
			Conjugate@hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}],
		n==2&&i==0&&j==1,
			hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}]+
			Conjugate@hofSymmetricYBlock[ky,nFlux,{m,n}],
		n==2&&i==1&&j==0,
			Conjugate@hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}]+
			hofSymmetricYBlock[ky,nFlux,{m,n}],
		i==j, hofSymmetricXBlock[kx,j*nFlux,{m,n}],
		i==0&&j==n-1, hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}],
		i==j+1, hofSymmetricYBlock[ky,nFlux,{m,n}],
		i==n-1&&j==0, Conjugate@hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}],
		j==i+1, Conjugate@hofSymmetricYBlock[ky,nFlux,{m,n}],
		True, ConstantArray[0,{m,m}]
	],
	{i,0,n-1},
	{j,0,n-1}
];

NdH[{kx_,ky_},{"HofSymmetric",m_Integer,n_Integer,nFlux_},{1}]:=
ArrayFlatten@ Table[
	Which[
		i==j, DHofSymmetricXBlock[kx,j*nFlux,{m,n}],
		True, ConstantArray[0,{m,m}]
	],
	{i,0,n-1},
	{j,0,n-1}
];

NdH[{kx_,ky_},{"HofSymmetric",m_Integer,n_Integer,nFlux_},{2}]:=
2Pi*I*ArrayFlatten@ Table[
	Which[
		n==1, 
			hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}]+
			-Conjugate@hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}],
		n==2&&i==0&&j==1,
			hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}]+
			-Conjugate@hofSymmetricYBlock[ky,nFlux,{m,n}],
		n==2&&i==1&&j==0,
			-Conjugate@hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}]+
			hofSymmetricYBlock[ky,nFlux,{m,n}],
		i==0&&j==n-1, hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}],
		i==j+1, hofSymmetricYBlock[ky,nFlux,{m,n}],
		i==n-1&&j==0, -Conjugate@hofSymmetricYBlock[ky,(m+1)nFlux,{m,n}],
		j==i+1, -Conjugate@hofSymmetricYBlock[ky,nFlux,{m,n}],
		True, ConstantArray[0,{m,m}]
	],
	{i,0,n-1},
	{j,0,n-1}
];


(* ::Text:: *)
(*General gauge, intrpolating between Symmetric and Landau, for testing gauge invariance*)


nParams["HofGeneral"]=4; 
lattice["HofGeneral"]="Square";
nBands[{"HofGeneral",ux_,uy_,__}]:= ux*uy;

generalPhaseFn[1,{x_,y_},{P_,Q_},\[Gamma]_]:=
y (1-\[Gamma])/(P*Q)+ If[x == P-1, y*\[Gamma]/Q, 0];
	
generalPhaseFn[2,{x_,y_},{P_,Q_},\[Gamma]_]:=
-x*\[Gamma]/(P*Q)+ If[y == Q-1, -x(1-\[Gamma])/P, 0];	

hofXTransBlock[kx_,y_,{P_,Q_},xPhaseFn_]:=
Normal@ SparseArray[
	Which[
		P==1,{
			Band[{1,1}]->Exp[2Pi*I*kx]Exp[2Pi*I*-xPhaseFn[{P-1,y},{P,Q}]]
		},
		True,{
			Band[{2,1}]->Table[
					Exp[2Pi*I*kx]Exp[2Pi*I*-xPhaseFn[{x,y},{P,Q}]],
					{x,0,P-2}
				],
			Band[{1,P}]->Exp[2Pi*I*kx]Exp[2Pi*I*-xPhaseFn[{P-1,y},{P,Q}]]
		}
	],
	{P,P}
];

hofTrans[1,{kx_,ky_},{P_Integer,Q_Integer,nFlux_},phaseFn_]:=
ArrayFlatten@ Table[
	If[ i==j, 
		hofXTransBlock[kx,j,{P,Q},(phaseFn[1,#1,#2]&)],
		ConstantArray[0,{P,P}]
	],
	{i,0,Q-1},
	{j,0,Q-1}
];

hofYTransBlock[ky_,y_,{P_,Q_},yPhaseFn_]:=
DiagonalMatrix@ Table[
	Exp[2Pi*I*ky]Exp[2Pi*I*-yPhaseFn[{x,y},{P,Q}]],
	{x,0,P-1} 
];

hofTrans[2,{kx_,ky_},{P_Integer,Q_Integer,nFlux_},phaseFn_]:=
ArrayFlatten@ Table[
	If[ i == Mod[j+1,Q],
		hofYTransBlock[ky,j,{P,Q},(phaseFn[2,#1,#2]&)],
		ConstantArray[0,{P,P}]
	],
	{i,0,Q-1},
	{j,0,Q-1}
];

H[{kx_,ky_},{"HofGeneral2",P_Integer,Q_Integer,nFlux_,phaseFn_}]:=
(hofTrans[1,{kx,ky},{P,Q,nFlux},phaseFn]+
hofTrans[2,{kx,ky},{P,Q,nFlux},phaseFn]+
ConjugateTranspose@ hofTrans[1,{kx,ky},{P,Q,nFlux},phaseFn]+
ConjugateTranspose@ hofTrans[2,{kx,ky},{P,Q,nFlux},phaseFn])/.{Conjugate->Identity};

NdH[{kx_,ky_},{"HofGeneral2",m_Integer,n_Integer,nFlux_,\[Gamma]_},{1}]:=
ArrayFlatten@ Table[
	Which[
		i==j, DHofSymmetricXBlock2[kx,j*nFlux,{m,n},\[Gamma]],
		True, ConstantArray[0,{m,m}]
	],
	{i,0,n-1},
	{j,0,n-1}
];

NdH[{kx_,ky_},{"HofGeneral2",m_Integer,n_Integer,nFlux_,\[Gamma]_},{2}]:=
2Pi*I*ArrayFlatten@ Table[
	Which[
		n==1, 
			hofSymmetricYBlock2[ky,-(\[Gamma]+n(1-\[Gamma]))nFlux,{m,n}]+
			-ConjugateHofSymmetricYBlock2[ky,-(\[Gamma]+n(1-\[Gamma]))nFlux,{m,n}],
		n==2&&i==0&&j==1,
			hofSymmetricYBlock2[ky,-(\[Gamma]+n(1-\[Gamma]))nFlux,{m,n}]+
			-ConjugateHofSymmetricYBlock2[ky,-\[Gamma]*nFlux,{m,n}],
		n==2&&i==1&&j==0,
			-ConjugateHofSymmetricYBlock2[ky,-(\[Gamma]+n(1-\[Gamma]))nFlux,{m,n}]+
			hofSymmetricYBlock2[ky,-\[Gamma]*nFlux,{m,n}],
		i==0&&j==n-1, hofSymmetricYBlock2[ky,-(\[Gamma]+n(1-\[Gamma]))nFlux,{m,n}],
		i==j+1, hofSymmetricYBlock2[ky,-\[Gamma]*nFlux,{m,n}],
		i==n-1&&j==0, -ConjugateHofSymmetricYBlock2[ky,-(\[Gamma]+n(1-\[Gamma]))nFlux,{m,n}],
		j==i+1, -ConjugateHofSymmetricYBlock2[ky,-\[Gamma]*nFlux,{m,n}],
		True, ConstantArray[0,{m,m}]
	],
	{i,0,n-1},
	{j,0,n-1}
];


(* ::Subsection:: *)
(*Definitions for all models*)


modelName[params_List]:= First[params];
modelCouplings[params_List]:= Rest[params];

nParams[model_String]:=nParams[model]=Length[paramShortNames[model]];
nParams[params_List]:=nParams[modelName[params]];
(* nBands[params_List]:=nBands[modelName[params]]; *)
lattice[params_List]:=lattice[modelName[params]];
latticeDim[params_List]:=latticeDim[lattice[modelName[params]]];
unitCellVol[params_List]:=unitCellVol[lattice[modelName[params]]];


(* ::Text:: *)
(*We can (and do) write any two-band Hamiltonian as a 3-vector (blochN) dotted into the three Pauli sigma matrices.*)


H[k_List,params_List]:=
blochNConst[k,params]*IdentityMatrix[2]+
	blochN[k,params].Array[PauliMatrix,3]/;
		nBands[params]==2;


(* ::Text:: *)
(*We want to make dH/dk\[Mu] a predefined function rather than recomputing the derivative each time, but defining it "by hand" leaves the possibility for error if someone changes the definition of H but forgets to make the corresponding change in dH/dk\[Mu]. Instead, given an analytic function of momentum H[k,params], this routine calculates the first and second momentum derivatives and defines a new function dH[k, params, \[Mu]] (and dH[k, params, \[Mu],\[Nu]] = d^2 H/dk\[Mu] dk\[Nu]). *)
(**)
(*We're putting the derivative component labels \[Mu], \[Nu], \[Rho] etc. on the end of the list of arguments to dH, because our general convention for momentum-dependent functions is that the first parameter is always the momentum, the second is always the list of parameters defining the model, followed by any extra arguments or parameters the function itself requires. In particular, this is what bzIntegrate expects.*)


(* ::Text:: *)
(*We use Set (=) instead of SetDelayed (:=) here and afterwards in this section for defining the Hamiltonians, on the theory that the former will evaluate slightly faster...*)


Block[{dummy},
	NH[{kx_?NumericQ,ky_?NumericQ},params_List]:=
		N@H[{kx,ky},params];
	NdH[{kx_?NumericQ,ky_?NumericQ},params_List,{\[Mu]_Integer}]:=
	N@D[
		H[ReplacePart[{kx,ky},\[Mu]->dummy],params],
		dummy
		]/.{dummy->({kx,ky}[[\[Mu]]])}
	];
SetAttributes[NH,NumericFunction];
SetAttributes[NdH,NumericFunction];


(* ::Section::Closed:: *)
(*Analytic solution, curvature and metric for two-band Hamiltonians*)


(* ::Subsection::Closed:: *)
(*GS energy and eigenvector*)


bandEnergy2[k_,params_,sgn_]:=
With[{
	n0=blochNConst[k,params],
	n=blochN[k,params]
},
n0+sgn*Sqrt[n.n]
];

gsVectors2[k_,params_,sgn_]:=
With[{
	n=blochN[k,params],
	nn=(Norm[blochN[k,params]]/.Abs->Identity)
},
1/Sqrt[2nn(nn+sgn*n[[3]])]*
{
	n[[3]]+sgn*nn,
	n[[1]]+I*n[[2]]
}];


(* ::Subsection:: *)
(*Analytic curvature and metric*)


(* ::Text:: *)
(*Conventions: Curvature and metric aren't affected by a shift in energy so w/log we can parameterize any two-band hamiltonian by a 3-vector n:*)
(*H = n . \[Sigma]*)
(*We don't assume the n given to us is normalized, since computing that analytically is more trouble than it's worth, and instead normalize the curvature and metric explicitly.*)
(**)
(*The following functions expect n to be given by the function*)
(*N2[ "Model Name", momentum, parameters],*)
(*where "Model Name" is any string labeling the Hamiltonian we're investigating (in the hopes that this would simplify computing these things for lots of different models), momentum is a vector in the reciprocal space, and parameters is a list of all other coupling constants needed to specify the Hamiltonian and can be of any length.*)
(*We also need the momentum derivatives of the Hamiltonian, which are provided by another function*)
(*\[PartialD]H/\[PartialD]Subscript[k, \[Mu]]= dN2["Model Name", \[Mu], momentum, parameters].*)


curvature2[k_,params_]:=
With[{
	n=NblochN[k,params]
},
(1/2)(
	n.Cross[NblochDN[k,params,{1}],NblochDN[k,params,{2}]]
)/(n.n)^(3/2)
];

metricComponent2[k_,params_,{\[Mu]_,\[Nu]_},OptionsPattern[{"Unimodular"->False}]]:=
With[{
	n=NblochN[k,params],
	dn=Table[NblochDN[k,params,{\[Rho]}],{\[Rho],2}]
},(
Cross[n,dn[[\[Mu]]]].Cross[n,dn[[\[Nu]]]]
)/
If[OptionValue["Unimodular"],
	(*true*)  Norm[n]Max[1.*^-8,Abs[n.Cross[dn[[1]],dn[[2]]]]],
	(*false*) (2*n.n)^2
]
];

trMetric2[k_,params_,OptionsPattern[{"Unimodular"->False,"Inequality"->False}]]:=
With[{
	n=NblochN[k,params],
	dn=Table[NblochDN[k,params,{\[Mu]}],{\[Mu],2}]
},(
	Cross[n,dn[[1]]].Cross[n,dn[[1]]]+
	Cross[n,dn[[2]]].Cross[n,dn[[2]]]-
	If[OptionValue["Inequality"],
		(*true*)  2Norm[n]Abs[n.Cross[dn[[1]],dn[[2]]]],
		(*false*) 0
])/
If[OptionValue["Unimodular"],
	(*true*)  Norm[n]Max[1.*^-8,Abs[n.Cross[dn[[1]],dn[[2]]]]],
	(*false*) (2*n.n)^2
]
];

detMetric2[k_,params_,OptionsPattern[{"Inequality"->False}]]:=
With[{
	n=NblochN[k,params],
	dn=Table[NblochDN[k,params,{\[Mu]}],{\[Mu],2}]
},
(
	(Norm[Cross[n,dn[[1]]]]Norm[Cross[n,dn[[2]]]])^2-
	Norm[Cross[n,dn[[1]]].Cross[n,dn[[2]]]]^2-
	If[OptionValue["Inequality"],
		(*true*)  (Norm[n](n.Cross[dn[[1]],dn[[2]]]))^2,
		(*false*) 0
])/(2*n.n)^4
];


(* ::Section:: *)
(*Numerical diagonalization of general Hamiltonians*)


(* ::Subsection:: *)
(*Spectrum*)


(* ::Text:: *)
(*The first thing we want is a flat band with a large gap to excitations -- for more complicated models like the Ruby model the only way to keep tabs on this is numerically.*)
(*Unfortunately the code below has more duplication than I'd like, but this is due to the fact that any function one wants Mathematica to treat purely numerically must have the NumericQ test on its inputs, and there's no way to do this for a pure function (which could have then been passed to a general minimization/integration routine...).*)


(* ::Text:: *)
(*Convention: we've sorted the band energies in ascending order, so the ground state is band #1.*)


bandEnergy[{kx_?NumericQ,ky_?NumericQ},params_,band_Integer:1]:=
First@bandEnergy[{kx,ky},params,{band}];

bandEnergy[{kx_?NumericQ,ky_?NumericQ},params_,bands_List]:=
If[
	nBands[params]==2,
	(*then*)
		(bandEnergy2[{kx,ky},params,(-1)^#]&)/@bands,
	(*else*)
		(Chop@Eigenvalues[
				NH[{kx,ky},params],
				Max[bands]
			])[[bands]]
];

SetAttributes[bandEnergy,NumericFunction];


bandGap[params_,band_:1]:=
bzMinValue[bandEnergy,params,band]-
	bzMaxValue[bandEnergy,params,band+1];

bandwidth[params_,band_:1]:=
bzMaxValue[bandEnergy,params,band]-
	bzMinValue[bandEnergy,params,band];

gapOverWidth[params_,band_:1]:=
With[{
	unoccEMin=bzMaxValue[bandEnergy,params,band+1],
	occEMax=bzMinValue[bandEnergy,params,band],
	occEMin=bzMaxValue[bandEnergy,params,band]
},
	(unoccEMin-occEMax)/(occEMax-occEMin)
];

bandTouchingRatio[params_,band_Integer]:=
With[{
	pts=Flatten[
		Table[
			{i,j}.bb[lattice[params]]/#,
			{i,0,#-1},{j,0,#-1}
		]&/@ Switch[
			lattice[params],
			"Triangular",{2,3},
			"Square",{4},
			_,{2}
		],
	2]
},
(Min[#]/Mean[#]&)[
	Abs[#1-#2]&@@@(
		bandEnergy[#,params,{band,band+1}]&/@ pts
	)
]
];


(* ::Subsection:: *)
(*Numerical diagonalization *)


(* ::Text:: *)
(*partitionEigensystem returns two lists of pairs {\[Lambda]_i, v_i} (where \[Lambda]_i is an eigenvalue and v_i the associated eigenvector); the first list ("ground states") consists the pairs associated with the nBands smallest (real) eigenvalues, and the second ("excited states") containing the rest. Mathematica sorts eigenvalues by absolute value, even when the spectrum is real, so unless we know a priori that the band Hamiltonian is positive definite, we need to examine all the eigenvalues to identify the ground state.*)
(*I've changed this routine to use Eigensystem instead of singular value decomposition -- further testing showed that the former method was slightly more numerically accurate, with no discernable difference in speed, and we don't need the Hamiltonian to be positive definite.*)


(* ::Text:: *)
(*The important new addition is that we "gauge fix" this answer -- when given a purely numerical matrix, mathematica returns eigenvectors normalized to unity; in addition, it seems to choose the phase so that the final component of each vector is real. The sign of this component is positive or negative seemly at random, and this can throw off the discretized Chern number computation below (in regions where the GS eigenvector isn't changing much, the inner product on links of the grid becomes close to -1, which hits the branch cut in the Arg function). Although that algorithm is supposed to be gauge invariant, its numerical stability is much improved if we make sure the eigenvectors have a common phase so that this problem doesn't occur. *)


(* ::Text:: *)
(*gsVectors only returns the list of ground state eigenvectors. *)


numericalGaugeFix[vect_]:=
vect/(Sign@@(
	Select[Reverse[vect],RealExponent[#]>-8.0&,1]
));

numericalGaugeFix[val_,vect_List]:=
{
	val,
	vect/(Sign@@(
		Select[Reverse[vect],RealExponent[#]>-8.0&,1]
	))
};


partitionEigensystem[{kx_,ky_},params_,occBands_:1]:=
With[{
	occ=If[IntegerQ[occBands],{occBands},occBands],
	unocc=Complement[Range[nBands[params]],#]&[
		If[IntegerQ[occBands],{occBands},occBands]
	],
	eigensys=Transpose[
		Eigensystem[NH[{kx,ky},params]]
	]
},
Function[
	{sortOrder},
	{
		numericalGaugeFix@@@
			eigensys[[ sortOrder[[occ]]]],
		numericalGaugeFix@@@
			eigensys[[ sortOrder[[unocc]]]]
	}
][
	Ordering[eigensys[[All,1]]]
]
];

gsVectors[{kx_,ky_},params_,nOccBands_]:=
If[
	nBands[params]==2&&nOccBands==1,
(*then*)
	{gsVectors2[{kx,ky},params,-1]},
(*else*)
	partitionEigensystem[{kx,ky},params,nOccBands][[1,All,2]]
];


(* ::Section:: *)
(*Numerical Chern number, curvature and metric for general band Hamiltonians*)


(* ::Subsection:: *)
(*Discretized Chern number computation*)


(* ::Text:: *)
(*[1]	T. Fukui, Y. Hatsugai, and H. Suzuki, J. Phys. Soc. Jpn. 74, 1674 (2005), cond-mat/0503172*)


(* ::Text:: *)
(*After lots of failed attempts to numerically evalute the Chern number through a less expensive line integral of the Berry connection, we decided to go with the method used in the above paper. This computes the Chern number as a BZ integral of the Berry curvature, so it has the advantage of being independent of gauge, but the authors point out that since the Chern number is quantized, the "integral" may safely be evaluated on a coarse mesh, as the sum of the discretized curvature on each plaquette (which is evaluated as the product of the phases of wavefunction overlaps around the plaquette).*)


(* ::Text:: *)
(*vects is a 2d array of the ground state eigenvectors at each point on the sampling grid -- to improve stability, we "nudge" the grid by a small amount to ensure that we don't run into problems from evaluating the ground state eigenvectors at points of abnormally high symmetry. The argument of the "/@" in the next line makes a list of plaquettes whose entries are the GS eigenvector at each of the vertices. The statement defining links  chops up this 2D list into a list of the paths of the (directed) links that we need to evaluate the wavefunction overlaps on in order to obtain the plaquette curvatures. Note that we need to use FractionalPart on the total u(1) flux in each plaquette, so links on the interior of the BZ don't necessarily cancel out.*)
(**)
(*Remark: checked that we get correct normalization and sign by comparing against known results for the Haldane model.*)


(* ::Text:: *)
(*Addendum: In actual use, apparently the algorithm as written can break when we're too close to a Chern number transition, and return "Indeterminate". In case this happens, fall back to a numeric integration of the Berry curvature. *)


chernNumber/: functionShortName[chernNumber]="chern";
chernNumber[params_List,nOcc_:1]:=
unitCellVol[params]/(2Pi)*bzIntegrate[(curvature[#1,##2]&),params,nOcc];
SetAttributes[chernNumber,NumericFunction];


(*chernNumber/: functionShortName[chernNumber]="chern";
chernNumber[params_List,nOcc_:1,OptionsPattern[{"nGrid"->7}]]:=
Module[{
	nudge= 1.*^-2*RandomReal[{-1,1},2],
	vects,links,ans,
	n=OptionValue["nGrid"]
},
vects=Table[
	gsVectors[
		(({i,j}+nudge)/n).cc[{lattice[params],"p1"}],
		params,
		nOcc
	],
{i,0,n},{j,0,n}
];
links=Partition[
	#[[{1,2,4,3,1}]],
	2,1
]&/@
	Flatten[
		Partition[
			vects,
			{2,2},{1,1}
		],
		{{1,2},{3,4}}
	];
ans= Total[
	FractionalPart[Total[#]]&/@
		Apply[
			(Arg[Det[
				Outer[
					Function[{v1,v2},Conjugate[v1].v2],
					#1,#2,
					1
				]
			]]/Pi&),
			links,
			{2}
		]
];
If[!NumericQ[ans],
	(*then*)Round[unitCellVol[params]/(2Pi)*
				bzIntegrate[curvature,params,nOcc,"BZ Area"->1]
			],
	(*else*) ans
]];
SetAttributes[chernNumber,NumericFunction];*)


(* ::Subsection:: *)
(*Curvature/metric subroutines*)


(* ::Text:: *)
(*For lack of a better name, I'm referring to the quantity*)
(*Subscript[R, \[Mu]]^(i,j) = \[LeftAngleBracket] Subscript[e, i] | \!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]\ H\) | Subscript[g, j] \[RightAngleBracket]/(Subscript[E, i] - Subscript[E, j])*)
(*as a "Resta factor," since it's a useful intermediate step in computing the curvature and metric tensors. Here e_i and g_j refer to an excited and ground state, respectively. The value of this for one pair (i,j) is computed by restaFactorEntry; restaFactor calls this to compute the tensor for all pairs of ground and excited states. It turned out to be a bit cleaner to compute these for all values of \[Mu] in restaFactors, since we wind up needing all that information to compute the curvature, metric etc.*)


restaFactorEntry[{gsEnergy_,gsVect_},{esEnergy_,esVect_},dH_]:=
(Conjugate[esVect].dH.gsVect)/(esEnergy-gsEnergy);

restaFactor[{gsData_,esData_},dH_]:=
Outer[
	restaFactorEntry[#1,#2,dH]&,
	gsData,
	esData,
	1
];

restaFactors[{kx_?NumericQ,ky_?NumericQ},params_,nOccBands_:1]:=
With[{
	sortedEigenstates=
		partitionEigensystem[{kx,ky},params,nOccBands]
},
Table[
	restaFactor[
		sortedEigenstates,
		NdH[{kx,ky},params,{\[Mu]}]
	],	
{\[Mu],latticeDim[lattice[params]]}
]
];

restaBilinears[{kx_?NumericQ,ky_?NumericQ},params_,nOccBands_:1]:=
With[{
	r=restaFactors[{kx,ky},params,nOccBands],
	d=latticeDim[lattice[params]]
},
Table[
	Tr[ConjugateTranspose[r[[\[Mu]]]].r[[\[Nu]]]],
	{\[Mu],d}, {\[Nu],d}
]];

SetAttributes[restaFactorEntry,NumericFunction]
SetAttributes[restaFactor,NumericFunction]
SetAttributes[restaFactors,NumericFunction]
SetAttributes[restaBilinears,NumericFunction]


(* ::Text:: *)
(*NB "curvature" normalization is division by |B|/2, in order to be comparable with the unimodular normalization of Sqrt[Det[G]]*)


normDenominator[resta_,norm_]:=
Max[1.*^-8,
	Switch[norm,
		"Unimodular", Sqrt[Abs[Det[Re[resta]]]],
		"Curvature", Abs[Im[resta[[2,1]]]],
		_, 1
	]
];


(* ::Subsection:: *)
(*Components of the metric*)


(* ::Text:: *)
(*possible issue: bzIntegrate of metric[[1,2]] gives incorrect answer; integrating the component directly does, though.*)


curvature/: functionShortName[curvature]="B";
curvature[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1]:=
If[
	nBands[params]==2,
(*then*)
	curvature2[{kx,ky},params],
(*else*)
	2*Im[restaBilinears[{kx,ky},params,nOcc][[2,1]]]
];

metric[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1]:=
If[
	nBands[params]==2,
(*then*)
	Table[
		metricComponent2[{kx,ky},params,{\[Mu],\[Nu]}],
		{\[Mu],2},{\[Nu],2}
	],
(*else*)
	With[{
		r=restaBilinears[{kx,ky},params,nOcc]
	},
	Table[
		Re[r[[\[Mu],\[Nu]]]],
		{\[Mu],2},{\[Nu],2}
	]]
];

metric11/: functionShortName[metric11]="g11";
Options[metric11]={"Norm"->None};
metric11[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1,
	opts:OptionsPattern[]]:=
If[
	nBands[params]==2,
(*then*) 
	metricComponent2[{kx,ky},params,{1,1}]/
	If[OptionValue["Norm"]==="Unimodular",
		detMetric2[{kx,ky},params,"Unimodular"->False,"Inequality"->False],
		1
	],
(*else*) 
	With[{r=restaBilinears[{kx,ky},params,nOcc]},
		Re[r[[1,1]]]/normDenominator[r,OptionValue["Norm"]]
	]
];

metric12/: functionShortName[metric12]="g12";
Options[metric12]={"Norm"->None};
metric12[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1,
	opts:OptionsPattern[]]:=
If[
	nBands[params]==2,
(*then*) 
	metricComponent2[{kx,ky},params,{1,2}]/
	If[OptionValue["Norm"]==="Unimodular",
		detMetric2[{kx,ky},params,"Unimodular"->False,"Inequality"->False],
		1
	],
(*else*) 
	With[{r=restaBilinears[{kx,ky},params,nOcc]},
		Re[r[[1,2]]]/normDenominator[r,OptionValue["Norm"]]
	]
];

metric22/: functionShortName[metric22]="g22";
Options[metric22]={"Norm"->None};
metric22[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1,
	opts:OptionsPattern[]]:=
If[
	nBands[params]==2,
(*then*) 
	metricComponent2[{kx,ky},params,{2,2}]/
	If[OptionValue["Norm"]==="Unimodular",
		detMetric2[{kx,ky},params,"Unimodular"->False,"Inequality"->False],
		1
	],
(*else*) 
	With[{r=restaBilinears[{kx,ky},params,nOcc]},
		Re[r[[2,2]]]/normDenominator[r,OptionValue["Norm"]]
	]
];

SetAttributes[curvature,NumericFunction];
SetAttributes[metric,NumericFunction];


(* ::Subsection:: *)
(*Geometric inequalities*)


(* ::Text:: *)
(*Lower bounds on the trace and determinant of the quantum metric, proved by RR in the "Band geometry..." paper. Both of the following quantities are >=0; we know the determinant inequality is always saturated for any two-band Hamiltonian. *)


trG/: functionShortName[trG]="trG";
Options[trG]={"Norm"->None};
trG[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1,
	OptionsPattern[]]:=
If[
	nBands[params]==2,
(*then*)
	If[OptionValue["Norm"]===None,
		trMetric2[{kx,ky},params,"Unimodular"->False,"Inequality"->False],
		trMetric2[{kx,ky},params,"Unimodular"->True,"Inequality"->False]
	],
(*else*)
	With[{
		r=restaBilinears[{kx,ky},params,nOcc]
	},
	Tr[Re[r]]/normDenominator[r,OptionValue["Norm"]]
	]
];

detG/: functionShortName[detG]="detG";
Options[detG]={"Norm"->None};
detG[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1,
	OptionsPattern[]]:=
If[
	nBands[params]==2,
(*then*)
	If[OptionValue["Norm"]===None,
		detMetric2[{kx,ky},params,"Inequality"->False],
		1.
	],
(*else*)
	With[{
		r=restaBilinears[{kx,ky},params,nOcc]
	},
	Det[Re[r]]/(normDenominator[r,OptionValue["Norm"]])^2
	]
];

trGInequality/: functionShortName[trGInequality]="trIneq";
Options[trGInequality]={"Norm"->None};
trGInequality[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1,
	OptionsPattern[]]:=
If[
	nBands[params]==2,
(*then*)
	If[OptionValue["Norm"]===None,
		trMetric2[{kx,ky},params,"Unimodular"->False,"Inequality"->True],
		trMetric2[{kx,ky},params,"Unimodular"->True,"Inequality"->True]
	],
(*else*)
	With[{
		r=restaBilinears[{kx,ky},params,nOcc]
	},(
	Tr[Re[r]]-Abs[2*Im[r[[2,1]]]]
	)/normDenominator[r,OptionValue["Norm"]]
	]
];

detGInequality/: functionShortName[detGInequality]="detIneq";
Options[detGInequality]={"Norm"->None};
detGInequality[{kx_?NumericQ,ky_?NumericQ},params_List,nOccBands_:1,
	OptionsPattern[]]:=
If[
	nBands[params]==2,
(*then*)
	If[OptionValue["Norm"]===None,
		detMetric2[{kx,ky},params,"Inequality"->True],
		1.
	],
(*else*)
	With[{
		r=restaBilinears[{kx,ky},params,nOccBands]
	},(
	Det[Re[r]]-(2*Im[r[[2,1]]])^2/4
	)/(normDenominator[r,OptionValue["Norm"]])^2
	]
];


(* ::Subsection:: *)
(*Covariance matrix*)


(* ::Text:: *)
(*Think that reducing AccuracyGoal keeps NIntegrate from wasting time refining estimates in regions where the integrand is near zero (since if one measures zero to n significant digits, one still has zero digits of precision, since we still don't know the scale of the quantity (ie where the nonzero digits start) but one has n digits of accuracy.*)


rmsCurvature/: functionShortName[rmsCurvature]="B-RMS";
rmsCurvature[params_,nOcc_:1,
	opts:OptionsPattern[{bzIntegrate,NIntegrate}]]:=
unitCellVol[params]/(2Pi)*Last@bzMeanAndRMS[
	(curvature[#1,##2]/(2Pi)&),
	params,
	nOcc,
	opts
];
(* rmsCurvature[params_,nOcc_:1,
	opts:OptionsPattern[{bzIntegrate,NIntegrate}]]:=
unitCellVol[params]*bzRMSKnownMean[
	(curvature[#1,##2]/(2Pi)&),
	params,
	chernNumber[params,nOcc,"nGrid"->7],
	opts
];*)

trRMSMetric/: functionShortName[trRMSMetric]="trRMSg";
trRMSMetric[params_,nOcc_:1,
	opts:OptionsPattern[{bzIntegrate,NIntegrate}]]:=
With[{
	diag=rmsMetric11[params,nOcc,opts],
	offdiag=bzRMSKnownMean[
		metric12,
		params,
		nOcc,
		0.0,
		opts
	]
},
Sqrt[(1/2)(offdiag^2+diag^2)]
];


Options[bTimesTrG]={"Norm"->None};
bTimesTrG[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1,
	means_List:{0,0},OptionsPattern[]]:=
With[{
	r=restaBilinears[{kx,ky},params,nOcc]
},
(2Im[r[[2,1]]]-First[means])(
	Tr[Re[r]]/normDenominator[r,OptionValue["Norm"]]- 
	Last[means]
)];

Options[bTimesDetG]={"Norm"->None};
bTimesDetG[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1,
	means_List:{0,0},OptionsPattern[]]:=
With[{
	r=restaBilinears[{kx,ky},params,nOcc]
},
(2Im[r[[2,1]]]-First[means])(
	Det[Re[r]]/(normDenominator[r,OptionValue["Norm"]])^2- 
	Last[means]
)];

Options[trGTimesDetG]={"Norm"->None};
trGTimesDetG[{kx_?NumericQ,ky_?NumericQ},params_List,nOcc_:1,
	means_List:{0,0},OptionsPattern[]]:=
With[{
	r=restaBilinears[{kx,ky},params,nOcc]
},
(Tr[Re[r]]-First[means])(Det[Re[r]]-Last[means])/
	(normDenominator[r,OptionValue["Norm"]])^3 
];

SetAttributes[bTimesTrG,NumericFunction];
SetAttributes[bTimesDetG,NumericFunction];
SetAttributes[trGTimesDetG,NumericFunction];

covarianceBTrG/: functionShortName[covarianceBTrG]="covBTrG";
Options[covarianceBTrG]=
	Join[Options[bTimesTrG],Options[bzIntegrate]];
covarianceBTrG[params_,nOcc_:1,opts:OptionsPattern[]]:=
With[{
	meanB=2Pi*chernNumber[params,nOcc,"nGrid"->7]/unitCellVol[params],
	trg=bzIntegrate[
			(trG[#1,##2,FilterRules[{opts},Options@trG]]&),
			params,nOcc,
			FilterRules[{opts},Options@bzIntegrate]
		]
},
bzIntegrate[
	(bTimesTrG[#1,##2,{meanB,trg},
		FilterRules[{opts},Options@bTimesTrG]]&),
	params,nOcc,
	FilterRules[{opts},Options@bzIntegrate]
]
];

covarianceBDetG/: functionShortName[covarianceBDetG]="covBDetG";
Options[covarianceBDetG]=
	Join[Options[bTimesDetG],Options[bzIntegrate]];
covarianceBDetG[params_,nOcc_:1,opts:OptionsPattern[]]:=
With[{
	meanB=2Pi*chernNumber[params,nOcc,"nGrid"->7]/unitCellVol[params],
	detg=bzIntegrate[
			(detG[#1,##2,FilterRules[{opts},Options@detG]]&),
			params,nOcc,
			FilterRules[{opts},Options@bzIntegrate]
		]
},
bzIntegrate[
	(bTimesDetG[#1,##2,{meanB,detg},
		FilterRules[{opts},Options@bTimesDetG]]&),
	params,nOcc,
	FilterRules[{opts},Options@bzIntegrate]
]
];

covarianceTrGDetG/: functionShortName[covarianceTrGDetG]="covTrGDetG";
Options[covarianceTrGDetG]=
	Join[Options[trGTimesDetG],Options[bzIntegrate]];
covarianceTrGDetG[params_,nOcc_:1,opts:OptionsPattern[]]:=
With[{
	trg=bzIntegrate[
			(trG[#1,##2,FilterRules[{opts},Options@trG]]&),
			params,nOcc,
			FilterRules[{opts},Options@bzIntegrate]
		],
	detg=bzIntegrate[
			(detG[#1,##2,FilterRules[{opts},Options@detG]]&),
			params,nOcc,
			FilterRules[{opts},Options@bzIntegrate]
		]
},
bzIntegrate[
	(trGTimesDetG[#1,##2,{trg,detg},
		FilterRules[{opts},Options@trGTimesDetG]]&),
	params,nOcc,
	FilterRules[{opts},Options@bzIntegrate]
]
];


(* ::Section::Closed:: *)
(*misc. observables*)


sl2Ranisotropy/: functionShortName[sl2Ranisotropy]="sl2rA";
sl2Ranisotropy[{kx_?NumericQ,ky_?NumericQ},params_List]:=
With[{
	t=trG[{kx,ky},params,"Norm"->"Unimodular"]
},
(2+t-Sqrt[t(t+4)])/2
];

metricEccentricity/: functionShortName[metricEccentricity]="eccG";
metricEccentricity[{kx_?NumericQ,ky_?NumericQ},params_List]:=
With[{
	rR=Re[restaBilinears[{kx,ky},params,1]]
},
Function[{tr,sqrt},
	If[tr<1.*^-6&&sqrt<1.*^-6,0.,
	2/(tr+sqrt)*Sqrt[tr*sqrt]
]][
	rR[[1,1]]+rR[[2,2]],
	Sqrt[(rR[[1,1]]-rR[[2,2]])^2+4*rR[[1,2]]^2]
]
];


(* ::Subsection:: *)
(*Fourier components*)


curvatureFourier[{kx_?NumericQ,ky_?NumericQ},params_,modes_]:=
latticeFourierBasis[{kx,ky},lattice[params],-modes]curvature[{kx,ky},params];

metric11Fourier[{kx_?NumericQ,ky_?NumericQ},params_,modes_]:=
latticeFourierBasis[{kx,ky},lattice[params],-modes]metric11[{kx,ky},params];
metric12Fourier[{kx_?NumericQ,ky_?NumericQ},params_,modes_]:=
latticeFourierBasis[{kx,ky},lattice[params],-modes]metric12[{kx,ky},params];
metric22Fourier[{kx_?NumericQ,ky_?NumericQ},params_,modes_]:=
latticeFourierBasis[{kx,ky},lattice[params],-modes]metric22[{kx,ky},params];
sqrtDetGFourier[{kx_?NumericQ,ky_?NumericQ},params_,modes_]:=
latticeFourierBasis[{kx,ky},lattice[params],-modes]Sqrt[detG[{kx,ky},params]];

SetAttributes[curvatureFourier,NumericFunction];
SetAttributes[metric11Fourier,NumericFunction];
SetAttributes[metric12Fourier,NumericFunction];
SetAttributes[metric22Fourier,NumericFunction];
SetAttributes[sqrtDetGFourier,NumericFunction];


(* ::Subsection:: *)
(*Combined lines of data*)


(* ::Text:: *)
(*functionSpec is a list of either function names (as Symbols) or two-element lists of the name followed by the list of options to pass to that function.*)


lineOfData[params_,functionSpec_List]:=
Flatten@Table[
	If[ListQ[fn],
		(First[fn])[params,Sequence@@Last[fn]],
		fn[params]
	],
	{fn,functionSpec}
];

parallelLineOfData[params_,functionSpec_List]:=
Flatten@ParallelTable[
	If[ListQ[fn],
		(First[fn])[params,Sequence@@Last[fn]],
		fn[params]
	],
	{fn,functionSpec}
];


curvatureFourier10[params_]:=
Through[{Re,Im}[
bzIntegrate[ 
	(curvatureFourier[#1,#2,{1,0}]&),
	params
]
]];

curvatureFourier11[params_]:=
Flatten@Join[
	curvatureFourier10[params],
	Through[{Re,Im}[
		{
		bzIntegrate[ 
			(curvatureFourier[#1,#2,{1,1}]&),
			params
		],
		bzIntegrate[ 
			(curvatureFourier[#1,#2,{1,-1}]&),
			params
		]}
	]],
	{intGradCurvature2[params]}
];

metricFourier10[params_]:=
Flatten@Through[{Re,Im}[{
	bzIntegrate[ 
		(metric11Fourier[#1,#2,{1,0}]&),
		params
	],
	bzIntegrate[ 
		(metric12Fourier[#1,#2,{1,0}]&),
		params
	],
	bzIntegrate[ 
		(metric22Fourier[#1,#2,{1,0}]&),
		params
	],
	bzIntegrate[ 
		(sqrtDetGFourier[#1,#2,{1,0}]&),
		params
	]
}]];


sigmaBAndIneqs[params_,opts:OptionsPattern[{intTrGInequality}]]:=
Flatten@If[
	nBands[params]==2,
	{
		rmsCurvature[params],
		intGradCurvature2[params],
		intTrGInequality[params,opts]
	},{
		rmsCurvature[params],
		intTrGInequality[params,opts],
		intDetGInequality[params,opts]
	}
];

intMetricComponents[params_,opts:OptionsPattern[{intMetric11}]]:=
Flatten@{
	intMetric11[params,opts],
	intMetric12[params,opts],
	intMetric22[params,opts]
};

covariancesBTrGDetG[params_,opts:OptionsPattern[{covarianceBTrG}]]:=
With[{
	meanB=2Pi*chernNumber[params,"nGrid"->7]/unitCellVol[params],
	trg=bzIntegrate[
			(trG[#1,#2,FilterRules[{opts},Options@trG]]&),
			params,
			FilterRules[{opts},Options@bzIntegrate]
		],
	detg=bzIntegrate[
			(detG[#1,#2,FilterRules[{opts},Options@detG]]&),
			params,
			FilterRules[{opts},Options@bzIntegrate]
		]
},
{
	bzIntegrate[
		(bTimesTrG[#1,#2,{meanB,trg},
			FilterRules[{opts},Options@bTimesTrG]]&),
		params,
		FilterRules[{opts},Options@bzIntegrate]
	],
	bzIntegrate[
		(bTimesDetG[#1,#2,{meanB,detg},
			FilterRules[{opts},Options@bTimesDetG]]&),
		params,
		FilterRules[{opts},Options@bzIntegrate]
	],
	bzIntegrate[
		(trGTimesDetG[#1,#2,{trg,detg},
			FilterRules[{opts},Options@trGTimesDetG]]&),
		params,
		FilterRules[{opts},Options@bzIntegrate]
	]
}
];


(* ::Section:: *)
(*Parameter derivatives for steepest descent -- just momenta for now*)


(* ::Text:: *)
(*NB: The general expressions are derived from degenerate perturbation theory (ie, the physical H taken as input to the functions much be such that all the occupied states have the same energy) and will need to be havily modified in the case of general occupied bands.*)
(**)
(*The integrands are written out assuming two spatial dimensions in the name of efficiency only; the expressions are valid for any number of dimensions as long as one is willing to compute all the derivatives.*)
(**)
(*See "Note on perturbation theory" for the derivation of the expressions for the gradient and Hessian of RMS curavture with respect to model parameters.*)


(* ::Text:: *)
(*Also note that the list called "P" in the subroutines below consists of the numerical values of {P, R, R^2, R^3, ...} so that eg P[[3]] actually refers to R^2.*)


(* ::Text:: *)
(*partitionEigensystem returns two lists of pairs {\[Lambda]_i, v_i} (where \[Lambda]_i is an eigenvalue and v_i the associated eigenvector); the first list ("ground states") *)
(**)
(*{{{egs,vgs}},{{es1},{es2},...}}*)


(* ::Text:: *)
(*NB: Minus sign in the expression for ResolventTerm arises because we've inverted the spectrum of compiled Hams (which are, currently, the Hams for all models) in the interest of numerical efficiency.*)


resolventTerm[{gsEnergy_,gsVect_},{esEnergy_,esVect_}]:=
(-1)Outer[
	Times,
	esVect,
	Conjugate[esVect]
]/(gsEnergy-esEnergy);

projectorAndResolvent[{kx_?NumericQ,ky_?NumericQ},params_,
	rMax_Integer, band_Integer:1]:=
With[{
	gses=partitionEigensystem[{kx,ky},params,{band}]
},
Prepend[
	FoldList[Dot,#,ConstantArray[#,rMax-1]]&[
		Total[
			resolventTerm[gses[[1,1]],#]&/@ gses[[2]]
		]
	],
	Outer[
		Times,
		gses[[1,1,2]],
		Conjugate[gses[[1,1,2]]]
	]
]];

SetAttributes[resolventTerm,NumericFunction]
SetAttributes[projectorAndResolvent,NumericFunction]


(* ::Subsection:: *)
(*Projector derivatives*)


dxP[{kx_?NumericQ,ky_?NumericQ},params_]:=
With[{
	x=NdH[{kx,ky},params,{1}],
	P=projectorAndResolvent[{kx,ky},params,1]
},
P[[2]].x.P[[1]]+P[[1]].x.P[[2]]
]

dyP[{kx_?NumericQ,ky_?NumericQ},params_]:=
With[{
	y=NdH[{kx,ky},params,{2}],
	P=projectorAndResolvent[{kx,ky},params,1]
},
P[[2]].y.P[[1]]+P[[1]].y.P[[2]]
]


(* ::Subsection:: *)
(*Momentum gradients*)


gradCurvatureSub[P_,{x_,y_},{xx_,xy_,yy_}]:={
2*Total[Im@{
	Tr[P[[1]].y.P[[3]].xx],
	-Tr[P[[1]].x.P[[3]].xy],
	Tr[P[[1]].y.P[[2]].x.P[[3]].x],
	2Tr[P[[1]].y.P[[3]].x.P[[2]].x],
	-Tr[P[[1]].x.P[[3]].y.P[[2]].x],
	-3*Tr[P[[1]].x]Tr[P[[1]].y.P[[4]].x]
}],
2*Total[Im@{
	Tr[P[[1]].y.P[[3]].xy],
	Tr[P[[1]].yy.P[[3]].x],
	2Tr[P[[1]].y.P[[2]].y.P[[3]].x],
	Tr[P[[1]].y.P[[3]].y.P[[2]].x],
	Tr[P[[1]].y.P[[3]].x.P[[2]].y],
	-3*Tr[P[[1]].y]Tr[P[[1]].y.P[[4]].x]
}]
};

gradCurvature/: functionShortName[gradCurvature]="dB2";
gradCurvature[{kx_?NumericQ,ky_?NumericQ},params_]:=
With[{
	x=NdH[{kx,ky},params,{1}],
	y=NdH[{kx,ky},params,{2}],
	xx=NdH[{kx,ky},params,{1,1}],
	xy=NdH[{kx,ky},params,{1,2}],
	yy=NdH[{kx,ky},params,{2,2}],
	L=projectorAndResolvent[{kx,ky},params,4]
},
Norm@gradCurvatureSub[L,{x,y},{xx,xy,yy}]
];


gradDetGSub[P_,{x_,y_},{xx_,xy_,yy_}]:=
With[{
	x3x=Tr[P[[1]].x.P[[3]].x],
	y3x=Tr[P[[1]].y.P[[3]].x],
	y3y=Tr[P[[1]].y.P[[3]].y],
	x4x=Tr[P[[1]].x.P[[4]].x],
	y4x=Tr[P[[1]].y.P[[4]].x],
	y4y=Tr[P[[1]].y.P[[4]].y],
	dxE=Tr[P[[1]].x],
	dyE=Tr[P[[1]].x]
},{
2*Total[{
	-Re[y3x](
		Re@Tr[P[[1]].y.P[[3]].xx]+
		Re@Tr[P[[1]].x.P[[3]].xy]+
		Re@Tr[P[[1]].y.P[[2]].x.P[[3]].x]+
		2*Re@Tr[P[[1]].y.P[[3]].x.P[[2]].x]+
		Re@Tr[P[[1]].x.P[[3]].y.P[[2]].x]+
		-3dxE*Re[y4x]+
		-dyE*x4x
	),
	x3x(
		Re@Tr[P[[1]].y.P[[3]].xy]+
		Re@Tr[P[[1]].y.P[[3]].y.P[[2]].x]+
		Re@Tr[P[[1]].y.P[[3]].x.P[[2]].y]+
		-dyE*Re[y4x]+
		-dxE*y4y
	),
	y3y(
		Re@Tr[P[[1]].x.P[[3]].xx]+
		2*Re@Tr[P[[1]].x.P[[3]].x.P[[2]].x]+
		-2dxE*x4x
	)
}],
2*Total[{
	-Re[y3x](
		Re@Tr[P[[1]].y.P[[3]].xy]+
		Re@Tr[P[[1]].yy.P[[3]].x]+
		2*Re@Tr[P[[1]].y.P[[2]].y.P[[3]].x]+
		Re@Tr[P[[1]].y.P[[3]].y.P[[2]].x]+
		Re@Tr[P[[1]].x.P[[3]].x.P[[2]].y]+
		-3dyE*Re[y4x]+
		-dxE*y4y
	),
	y3y(
		Re@Tr[P[[1]].x.P[[3]].xy]+
		Re@Tr[P[[1]].y.P[[2]].x.P[[3]].x]+
		Re@Tr[P[[1]].x.P[[3]].y.P[[2]].x]+
		-dxE*Re[y4x]+
		-dyE*x4x
	),
	x3x(
		Re@Tr[P[[1]].y.P[[3]].yy]+
		2*Re@Tr[P[[1]].y.P[[3]].y.P[[2]].y]+
		-2dyE*y4y
	)
}]
}];


(* ::Text:: *)
(*NB: covDerivCurvature as defined includes Jacobian factor of sqrt det G for integration measure*)


covDerivCurvature/: functionShortName[covDerivCurvature]="covDB2";
covDerivCurvature[{kx_?NumericQ,ky_?NumericQ},params_]:=
With[{
	x=NdH[{kx,ky},params,{1}],
	y=NdH[{kx,ky},params,{2}],
	xx=NdH[{kx,ky},params,{1,1}],
	xy=NdH[{kx,ky},params,{1,2}],
	yy=NdH[{kx,ky},params,{2,2}],
	L=projectorAndResolvent[{kx,ky},params,4]
},
Function[{
	B,dB,detg,dDetG,Ginverse
},
Re[
(detg*dB-(B/2)dDetG).(Ginverse/(detg^(5/2))).(detg*dB-(B/2)dDetG)
]][
	2Im@Tr[L[[1]].y.L[[3]].x],
	gradCurvatureSub[L,{x,y},{xx,xy,yy}],
	-(Re@Tr[L[[1]].y.L[[3]].x])^2+
		Tr[L[[1]].x.L[[3]].x]*Tr[L[[1]].y.L[[3]].y],
	gradDetGSub[L,{x,y},{xx,xy,yy}],
	{{
		Tr[L[[1]].y.L[[3]].y],
		-Re@Tr[L[[1]].y.L[[3]].x]
	},{
		-Re@Tr[L[[1]].y.L[[3]].x],
		Tr[L[[1]].x.L[[3]].x]
	}}
]];


(* ::Section:: *)
(*Automatic definition of BZ-integrated quantities*)


Do[
	With[{
		fn=Symbol[h],
		fnName=functionShortName[Symbol[h]],
		intFn=Symbol["int"<>ToUpperCase[StringTake[h,1]]<>StringDrop[h,1]],
		rmsFn=Symbol["rms"<>ToUpperCase[StringTake[h,1]]<>StringDrop[h,1]],
		avgWrmsFn=Symbol["avgWrms"<>ToUpperCase[StringTake[h,1]]<>StringDrop[h,1]]
	},
	SetAttributes[fn,NumericFunction];

	intFn/: functionShortName[intFn]=fnName;
	Options[intFn]=Join[{"RMS"->False},Options[fn],Options[bzIntegrate]];
	intFn[params_List, nOcc_:1, opts:OptionsPattern[{intFn,NIntegrate}]]:=
	If[OptionValue["RMS"],bzMeanAndRMS,bzIntegrate][
		(fn[#1,##2,Sequence@@FilterRules[{opts},Options@fn]]&),
		params,nOcc,
		FilterRules[{opts},Flatten[Options/@{bzIntegrate,NIntegrate}]]
	];

	rmsFn/: functionShortName[rmsFn]=fnName<>"-RMS";
	Options[rmsFn]=Join[Options[fn],Options[bzIntegrate]];
	rmsFn[params_List, nOcc_:1, opts:OptionsPattern[{rmsFn,NIntegrate}]]:=
	Last@bzMeanAndRMS[
		(fn[#1,##2,Sequence@@FilterRules[{opts},Options@fn]]&),
		params,nOcc,
		FilterRules[{opts},Flatten[Options/@{bzIntegrate,NIntegrate}]]
	];

	avgWrmsFn/: functionShortName[avgWrmsFn]={fnName,fnName<>"-RMS"};
	Options[avgWrmsFn]=Join[Options[fn],Options[bzIntegrate]];
	avgWrmsFn[params_List, nOcc_:1, opts:OptionsPattern[{avgWrmsFn,NIntegrate}]]:=
	bzMeanAndRMS[
		(fn[#1,##2,Sequence@@FilterRules[{opts},Options@fn]]&),
		params,nOcc,
		FilterRules[{opts},Flatten[Options/@{bzIntegrate,NIntegrate}]]
	];
	],
{h,BZIntegratedFunctions[]}
];


(* ::Section:: *)
(*End*)


End[]


EndPackage[]



