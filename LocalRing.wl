(* ::Package:: *)

(* Spare matrix epsilon pole free *)


epCondition[expr_]:=(Expand[expr/.ep->0]!=0);
nonzeroCount[expr_]:=Length[ArrayRules[expr]]-1;

RowClear[row_]:=Module[{numlist,denlist},
	numlist=Numerator/@row;
	denlist=Denominator/@row;
	Return[row/(PolynomialGCD@@numlist)*(PolynomialGCD@@denlist)//Cancel];
];

Options[GaussianElimination]={"PivotStrategy"->"OriginalIfPossible"};
GaussianElimination[M1_,OptionsPattern[]]:=Module[{M=M1,m,n,rowpointer,colpointer,colIndex,subM,arrayrules,rowIndice,columnIndice,
rowCounting,columnCounting,MarkowitzRep,pivotPos,tempCol,temp,tempRow,pivot,head,relevantRow,i},
	{m,n}=Dimensions[M];
	If[m<=1,Return[M]];

	rowpointer=colpointer=1;
	colIndex=Range[n];
	
	While[rowpointer<=m&&colpointer<=n,
		(* Search for a pivot *)
		subM=M[[rowpointer;;,colpointer;;]];  (* This is a sub matrix *)
		arrayrules=Select[ArrayRules[subM],#[[1]]=!={_,_}&];
		arrayrules=Select[arrayrules,epCondition[#[[2]]]&];
		If[arrayrules=={},Break[]];
		
		(* Markowitz *)
		(* Markowitz Number computation *)
		rowIndice=#[[1]]&/@Keys[arrayrules]//Union;
		columnIndice=#[[2]]&/@Keys[arrayrules]//Union;
		rowCounting=#->nonzeroCount[subM[[#]]]&/@rowIndice;
		columnCounting=#->nonzeroCount[subM[[All,#]]]&/@columnIndice;
		MarkowitzRep=Dispatch[#->((#[[1]]/.rowCounting)-1)*((#[[2]]/.columnCounting)-1)&/@Keys[arrayrules]];
		
		Switch[OptionValue["PivotStrategy"],
			"OriginalIfPossible",
			pivotPos=SortBy[Keys[arrayrules],{colIndex[[#[[2]]+colpointer-1]],#/.MarkowitzRep}&][[1]], (* First compare the Original column number *)
			"Markowitz",
			pivotPos=SortBy[Keys[arrayrules],{#/.MarkowitzRep,colIndex[[#[[2]]+colpointer-1]]}&][[1]] (* First compare Markowitz number than compare the column number *)
		];
		
		
	
		pivotPos+={rowpointer-1,colpointer-1};    (* shift back to the original matrix's position *)
		
		(* Test *)
		(* ArrayRules[M[[rowpointer;;,colpointer]]]//Print;
		ArrayRules[M[[rowpointer;;,pivotPos[[2]]]]]//Print;
		Print[M[[pivotPos[[1]],pivotPos[[2]]]]]; *)
		
		
		(* Swapping *)
		tempCol=M[[All,colpointer]];
		M[[All,colpointer]]=M[[All,pivotPos[[2]]]];
		M[[All,pivotPos[[2]]]]=tempCol;
		
		temp=colIndex[[colpointer]];
		colIndex[[colpointer]]=colIndex[[pivotPos[[2]]]];
		colIndex[[pivotPos[[2]]]]=temp;
		
		tempRow=M[[rowpointer]];
		M[[rowpointer]]=M[[pivotPos[[1]]]];
		M[[pivotPos[[1]]]]=tempRow;
		
		(* Row Reduction *)
		pivot=M[[rowpointer,colpointer]];
		M[[rowpointer]]=M[[rowpointer]]/pivot//Cancel; (* Normalization Maybe slow *)
		
		relevantRow=Complement[Flatten[#[[1]]&/@Keys[ArrayRules[M[[rowpointer+1;;,colpointer]]]]],{_}]+rowpointer;
		For[i=1,i<=Length[relevantRow],i++,
			M[[relevantRow[[i]] ]]=M[[relevantRow[[i]] ]]-M[[relevantRow[[i]],colpointer]]*M[[rowpointer]]//Together;
			M[[relevantRow[[i]] ]]=M[[relevantRow[[i]] ]]//RowClear;
		];
		Print["Row ",rowpointer," with the column swap ",colpointer," <-> ",pivotPos[[2]],"   ",Ints[[colIndex[[pivotPos[[2]]]]]]," replaced by ", Ints[[colIndex[[colpointer]]]]];
		
		
		
		(* Next step *)
		rowpointer++;
		colpointer++;
		
	];	
	Return[{M,colIndex}];
];
TailReduction[M1__,OptionsPattern[]]:=Module[{M=M1,m,n,i,j,ar,pivotColIndex,pivot,localtimer,timer},
    timer=AbsoluteTime[];
	{m,n}=Dimensions[M];
	pivotColIndex=ArrayRules[#][[1,1,1]]&/@M;
	For[i=m,i>=1,i--,
	localtimer=AbsoluteTime[];
			pivot=M[[i,pivotColIndex[[i]]]];
			If[!(pivot===1),M[[i]]=Together[M[[i]]/pivot]];    (* Normalize the row *)
			For[j=m,j>i,j--,
					M[[i]]-=M[[i,pivotColIndex[[j]]]]*M[[j]]; (* tail reduction *)				
			];
			M[[i]]=Together[M[[i]]];
		Print["Row ",i, " finished ... ",AbsoluteTime[]-localtimer];
	];
	Print["total time ... ",AbsoluteTime[]-timer];
	Return[M];
];


