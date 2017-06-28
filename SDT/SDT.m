(* ::Package:: *)

(* ::Section:: *)
(*SDT*)


(* :Title: Signal Detection Tools *)
(* :Context: SDT` *)
(* :Author: Flip Phillips, OSU, Skidmore College *)
(* :Summary: 
	This package provides various signal detection theory function to Mathematica.
*)
(* :Package Version: $Revision: 28 $ *)
(* :Mathematica Version: 8.0+ *)
(* :Copyright: Copyright 1999-2017, Flip Phillips, All Rights Reserved.  *)
(* :History: *)
(* :Keywords: *)
(* :Limitations: 
	some of the h/f code still needs the 0/1 corrections *)
(* :Discussion:  *)
(* :References: 
	This package is designed for doing signal detection theory calculations, based on equations from:
 	N. A. Macmillan & C. D. Creelman (1991).Detection Theory: A User's Guide, Cambridge University Press.

	mAFC stuff from:
 	Smith, J. E. Keith (1982). Simple Algorithms for M-Alternative Forced Choice, P&P.
*)


(* ::Text:: *)
(*LR+ = sensitivity / (1-specificity) = (a/(a+c)) / (b/(b+d))*)
(*LR- = (1-sensitivity) / specificity = (c/(a+c)) / (d/(b+d))*)
(*Post-test odds = pre-test odds * LR*)
(*Pre-test odds = pre-test probability / (1-pre-test probability)*)
(*Post-test probability = post-test odds / (post test odds+1)*)


BeginPackage["SDT`"]


DPrime::usage=
"DPrime[h,f,opts] computes sensitivity as \!\(\*SuperscriptBox[\(d\), \(\[Prime]\)]\) for a hit rate H and false-alarm rate F.
DPrime[pc,opts] computes sensitivity as \!\(\*SuperscriptBox[\(d\), \(\[Prime]\)]\) from p(c).
DPrime[pc,m,Design->mAFC] computes sensitivity as \!\(\*SuperscriptBox[\(d\), \(\[Prime]\)]\) for a mAFC design where m = number of alternatives
DPrime[mtx,opts] computes sensitivity as \!\(\*SuperscriptBox[\(d\), \(\[Prime]\)]\) for a response matrix.";


LogAlpha::usage=
"LogAlpha[h,f,opts] computes sensitivity as log(\[Alpha]) for a hit rate H and false-alarm rate F
LogAlpha[pc,opts] computes sensitivity as log(\[Alpha]) for a p(c).";


Criterion::usage="Criterion[h,f,opts] Computes C for a hit rate H and false-alarm rate F. For SameDifferent designs you need to call it with its matrix form - ala Criterion[m,opts].";


ROC::usage="ROC[dprime,f,opts] Computer the ROC curve location for a given \!\(\*SuperscriptBox[\(d\), \(\[Prime]\)]\) and false-alarm rate.";
(* PlotROC::usage="PlotROC[dPrime,opts] creates a ROC plot for the given \!\(\*SuperscriptBox[\(d\), \(\[Prime]\)]\)"; *)


Make2x2::usage="Make2x2[srmatrix,responses] creates a traditional 2x2 SDT matrix from a list of stimuli and responses (srmatrix) in the form {{stim,resp},{stim,resp}...}. The responses variable defaults to {\"yes\",\"no\"} and should contain the list of 2 possible responses if they are different from this. The resulting matrix contains rows for each stimulus class and columns for each response. The matrix is left unnormalized (counts only) by default and in the form {{hits,miss},{false alarm,correct rejections}}.";
Normalize2x2::usage="Normalize is an option to Make2x2 that specifies that the results should be expressed in percentages rather than in raw results counts.";


Design::usage="Option to SDT routines describing the experimental design. One of: YesNo, SameDifferent, ABX, TwoAFC, ThreeAFC, mAFC, Oddity.";
HitsAndFalseAlarms::usage="HitsAndFalseAlarms[{{s,r},...}] calculates {h,f} for a given list of signal present / response {{s,r},...}.
These can be specified by {1,0} or {True,False}";
(* Options: ZeroCorrection\[Rule]True treats {0,0} as a 'hit' e.g., for non-true-SDT trials. Defaults to False";*)


Sensitivity::usage="Sensitivity[h,f,opts], Sensitivity[matrix,opts] calculates sensitivity (q) for a given SDT matrix or h/f rate.";
PercentCorrect::usage="PercentCorrect[h,f,opts], PercentCorrect[matrix,opts] calculates sensitivity in the form p(c) for a given SDT matrix or h/f rate.";


(* Z::usage="Z[p] returns the quantile of p for a Gaussian distribution with \[Mu]=0, \[Sigma]=1."; *)

(* ZListPlot::usage="ZPlot[list,opts] transformed plot, {x,z[y]}. Linear when y is the CDF of a Gaussian with \[Mu]=0, \[Sigma]=1.";
ZZListPlot::usage="ZZPlot[list,opts] transformed plot, {z[x],z[y]}. Linear when x and y are the CDF of a Gaussian with \[Mu]=0, \[Sigma]=1.";
*)


YesNo::usage="Design of SDT task where the response is with respect to the stimulus: was it present?";
SameDifferent::usage="Design of SDT task where the response is Same or Different stimulus after initial stimulus presentation.";
TwoAFC::usage="Two alternative forced choice design.";
ThreeAFC::usage="Three alternative forced choice design.";
mAFC::usage="multi-alternative forced choice design.";
ABX::usage="SDT design where A and B are presented, then X is presented and response is whether it is A or B.";
Oddity::usage="SDT design that is odd.";
ZeroCorrection::usage="Correct for {0,1}";


Begin["`Private`"]


(* internal stuff *)
z[p_]:= Quantile[NormalDistribution[0,1],Re[p]];
phi[p_]:= PDF[NormalDistribution[0,1],z[p]];
p[z_]:=(1+Erf[z/Sqrt[2]])/2;

toprop[m_?MatrixQ,opts___?OptionQ]:=Module[{n},
	Map[(n=N[Plus @@ #];
		If[!FreeQ[#,{0,_}],{1/(2 n),1-(1/(2 n))},
			If[!FreeQ[#,{_,0}],{1-(1/(2 n)),1/(2 n)},#/n],#/n])&,m]]


(* Error messages *)


SDT::unity="Hit and False-Alarm rates must be between 0 and 1";
SDT::unity2="p(c) must be between 0 and 1";
SDT::unity3="p must be between 0 and 1";
SDT::baddsgn="Invalid or Unimplimented Design type \"`1`\"";


Make2x2::badresp="The SR matrix does not contain responses from the list `1`.";


(* Hits/fas *)
Options[HitsAndFalseAlarms]={ZeroCorrection->False};

(* HitsAndFalseAlarms[m_?MatrixQ,opts___?OptionQ]:=Module[{pp,corr},
	corr=ZeroCorrection/.{opts}/.Options[HitsAndFalseAlarms];
	pp=If[corr,toprop[m,opts],Map[(#/N[Plus @@ #])&,m]];
	{pp[[1,1]],pp[[2,1]]}]
*)

HitsAndFalseAlarms[m_?MatrixQ,opts___?OptionQ]:=Module[{bm,hits,fas,nS},
	bm=m/.{0->False,1->True};
	hits = Count[bm, {True, True}];
	fas = Count[bm, {False, True}];
	nS = Count[bm, {True, _}];
	{hits/nS,fas/nS}]


(* p(c) *)
Options[PercentCorrect]={Design->YesNo};

PercentCorrect[dPrime_Real,opts___?OptionQ]:=Module[{design},
	design=Design/.{opts}/.Options[PercentCorrect];
	Switch[design,
		YesNo,		p[dPrime/2],
		SameDifferent, 2(p[dPrime/2]^2)-2p[dPrime/2]+1,
		_,Message[SDT::baddsgn,design]
	]]

PercentCorrect[h_Real,f_Real,opts___?OptionQ]:=Module[{design},
	design=Design/.{opts}/.Options[PercentCorrect];
	Switch[design,
		YesNo,		PercentCorrect[DPrime[h,f]],
		SameDifferent,p[(z[h]-z[f])/2] ,
		TwoAFC,    	p[(z[h]-z[f])/2] ,
		_,Message[SDT::baddsgn,design]
	]]

PercentCorrect[m_?MatrixQ,opts___?OptionQ]:=Module[{h,f},
	{h,f}=HitsAndFalseAlarms[m,opts];
	PercentCorrect[h,f,opts]]


(* Sensitivity *)
Options[Sensitivity]={Design->YesNo};

Sensitivity[h_Real,f_Real,opts___?OptionQ]:=Module[{design},
	design=Design/.{opts}/.Options[Sensitivity];
	(h-f)/(1-f)]

Sensitivity[m_?MatrixQ,opts___?OptionQ]:=Module[{design,h,f},
	design=Design/.{opts}/.Options[Sensitivity];
	{h,f}=HitsAndFalseAlarms[m,opts];
	Sensitivity[h,f,opts]]


(* d' *)
Options[DPrime]={Design->YesNo};

DPrime[h_Real,f_Real,opts___?OptionQ]:=Message[SDT::unity]/;h>1||f>1;
DPrime[pc_Real,opts___?OptionQ]:=Message[SDT::unity2]/;pc>1;
DPrime[pc_Real,m_Integer,Design->mAFC]:=Message[SDT::notimpl]/;m>12;
DPrime[pc_Real,m_Integer,Design->mAFC]:=Message[SDT::unity2]/;pc>1;


DPrime[pc_Real,opts___?OptionQ]:=Module[{design},
	design=Design/.{opts}/.Options[DPrime];
	Switch[design,
		YesNo,2z[pc],
		SameDifferent, 2z[(1+Sqrt[(2pc-1)])/2],
		_,Message[SDT::baddsgn,design]
	]]


DPrime[pc_Real,m_Integer,Design->mAFC]:=Module[{km},
	km=0.86-0.085 Log[m-1];
	km Log[((m-1)pc)/(1-pc)]
]/;m<12


DPrime[h_Real,f_Real,opts___?OptionQ]:=Module[{design},
	design=Design/.{opts}/.Options[DPrime];
	Switch[design,
		YesNo,(z[h]-z[f]),
		TwoAFC,1/Sqrt[2] (z[h]-z[f]),
		SameDifferent,DPrime[PercentCorrect[h,f,Design->SameDifferent],Design->SameDifferent],
		_,Message[SDT::baddsgn,design]
]]


DPrime[m_?MatrixQ,opts___?OptionQ]:=Module[{h,f},
	{h,f}=HitsAndFalseAlarms[m,opts];
	DPrime[h,f,opts]]


(* log(a) *)
Options[LogAlpha]={Design->YesNo};

LogAlpha[h_,f_,opts___?OptionQ]:=Message[SDT::unity]/;h>1||f>1;
LogAlpha[pc_,opts___?OptionQ]:=Message[SDT::unity2]/;pc>1;

LogAlpha[pc_,opts___?OptionQ]:=Module[{design},
	design=Design/.{opts}/.Options[LogAlpha];
	Switch[design,
		YesNo,Log[pc/(1-pc)],
		_,Message[SDT::baddsgn,design]
]]

LogAlpha[h_,f_,opts___?OptionQ]:=Module[{design},
	design=Design/.{opts}/.Options[LogAlpha];
	Switch[design,
		YesNo,Log[h/(1-h)]/2-Log[f/(1-f)]/2,
		TwoAFC,1/Sqrt[2] (z[h]-z[f]),
		_,Message[SDT::baddsgn,design]
]]


(* criterion *)
Options[Criterion]={Design->YesNo};

Criterion[h_Real,f_Real,opts___?OptionQ]:=Message[SDT::unity]/;h>1||f>1;
Criterion[pc_Real,m_Integer,Design->mAFC]:=Message[SDT::notimpl]/;m>12;
Criterion[pc_Real,m_Integer,Design->mAFC]:=Message[SDT::unity2]/;pc>1;

Criterion[h_Real,f_Real,opts___?OptionQ]:=Module[{design},
	design=Design/.{opts}/.Options[Criterion];
	Switch[design,
		YesNo,-0.5(z[h]+z[f]),
		TwoAFC,-0.5(z[h]+z[f]),
		_,Message[SDT::baddsgn,design]
]]

Criterion[pc_Real,m_Integer,Design->mAFC]:=(m-1) pc/(1-pc)

Criterion[m_?MatrixQ,opts___?OptionQ]:=Module[{design,dp,h,fa},
	design=Design/.{opts}/.Options[Criterion];
	If[design!=SameDifferent,Return[Criterion[N[m[[1,1]]/(Plus@@m[[1]])],N[m[[2,1]]/(Plus@@m[[2]])],opts]]];

	(* otherwise, we need these: *)
	dp=DPrime[m,Design->SameDifferent];
	{h,fa}=HitsAndFalseAlarms[m,Design->SameDifferent];

	Return[-Sqrt[2]z[fa/2]-(dp/2)]
]


(* ROC and Z plotting stuff *)
Options[ROC]={Design->YesNo};

ROC[dPrime_,f_,opts___?OptionQ]:=Message[SDT::unity]/;f>1;


ROC[dPrime_,f_,opts___?OptionQ]:=Module[{design},
	design=Design/.{opts}/.Options[ROC];
	Switch[design,
		YesNo,(1+Erf[(Sqrt[2]*dPrime+2*InverseErf[0,-1+2*f])/2])/2,
		TwoAFC,(1+Erf[dPrime+InverseErf[0,-1+2*f]])/2,
		_,Message[SDT::baddsgn,design]
]]

(*
Options[PlotROC]=Options[ROC];

PlotROC[dPrime_,opts___?OptionQ]:=Plot[ROC[dPrime,f,opts],{f,0,1},FrameLabel->{"p(False Alarm)","p(Hit)"},AspectRatio->1,Epilog->{Hue[1],Thickness[0.01`],Line[{{0,0},{1,1}}],Hue[0.5`],Thickness[0.0075`],Line[{{0,1},{0.5`,0.5`}}]},Frame->True,PlotRange->All]


ZZListPlot[data_,opts___]:=
	ScaledListPlot[data,ScaleFunction->{Z[#]&,Z[#]&},
	AxesOrigin->Automatic,opts,Options[ZZPlot]]


ZListPlot[data_,opts___]:=
ScaledListPlot[data,ScaleFunction->{#&,Z[#]&},
AxesOrigin->Automatic,opts,Options[ZZPlot]]
*)


Options[Make2x2]={Normalize->False};

Make2x2[srmat_,resps_:{"yes","no"},opts___?OptionQ]:=Module[{res,poss,norm},
	res=Table[0,{2},{2}];
	poss=Union[Flatten[srmat]];
	If[Sort[poss]!=Sort[resps],Message[Make2x2::badresp,resps];Return[Null]];

	norm=Normalize2x2/.{opts}/.Options[Make2x2];
	Map[If[#[[1]]==#[[2]],
		(* correct *)
		If[#[[1]]==resps[[1]],
			res[[1,1]]++,
			res[[2,2]]++],
		(* incorrect *)
		If[#[[1]]==resps[[1]],
			res[[1,2]]++,
			res[[2,1]]++]
			]&,srmat];
		(* If[norm,res=N[res/{Plus@@res[[1]],Plus@@res[[2]]}]]; *)
	
	If[norm,res=toprop[res]];
	Return[res]
];


(*
Z[p_]:=Message[SDT::unity3]/;p>1||p<0;
Z[p_]:=z[p]
*)


End[]


EndPackage[]
