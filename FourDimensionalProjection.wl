(* ::Package:: *)

permutationMatrix[p_List]:=
	SparseArray[
		IdentityMatrix[Length[p]][[p]]
	] (*in earlier versions of Mathematica, PermutationMatrix was not implemented yet. This is an in-house implementation.*)


finiteCoefficients[ibps_List,G_(*integrals*),eps_(*dimreg parameter*)]:=
	Block[
		{
			localibps=ibps,
			nmis,
			kin,(*the kinematic variables*)
			allGs,(*all integrals appearing the IBPs*)
			reduced,
			reducible,
			mis={},
			shuffle
		},
		
		kin=DeleteCases[Variables[ibps],G[__]|eps];(*the Mandelstam invariants... *)
		kin=Thread[kin -> Prepend[RandomPrime[{10,20},Length[kin]-1],1]];(* ... are replaced by numerical values*)
			(*one of the variables can be set to 1 and restore by dimensional analysis*)
		
		localibps=localibps/.kin;
		localibps=If[!FreeQ[#,eps]&&Series[#,{eps,0,0}][[-4]]=={},Simplify[#/eps],#]&/@localibps;
			(*some of the previously derived ibp identity may have an overall eps, then we should remove this extra eps*)
		
		allGs=DeleteDuplicates@Cases[ibps,G[__],All];
		allGs=
			Reverse@
				SortBy[(*the integrals are sorted by their sectors and the total tensor rank*)
					allGs,
					{
						Length[DeleteCases[(List@@#),_?(#<=0&)]]&,
						Total[Abs[(List@@#)]]&
					}
				];
		nmis=Length[allGs]-Length[ibps];(*this is the number of MIs, provided that we use NeatIBP to generate the IBPs*)
		
		While[
			nmis!=Length[allGs],(*a while loop until the number of MIs with eps=0 is the same as the number of MIs with eps!=0*)
			
			(*The while loop is needed because, when we put eps=0, we may miss linear combinations of IBPs which are actually proportional to eps*)
		
			localibps=Transpose[Coefficient[localibps,#]&/@allGs];(*the IBP matrix*)
			localibps=Map[Simplify,localibps,{2}];
		
			reduced=RowReduce[localibps/.eps->0];(*IBP reduction with eps=0*)
		
			reducible=
				Part[#,1]&/@(*the position of the first 1 after a series of 0's in a row is a reducible integral*)
					DeleteCases[
						(Flatten/@
							(Position[#,1]&/@reduced)
						),
						{}
					];(*this integrals can be decomposed into MIs whose coefficients are finite as eps->0*)
			
			mis=Complement[Range[Length@allGs],reducible];(*the position of the MIs with respect to the IBPs with eps=0*)
		
			shuffle=permutationMatrix[Join[reducible,mis]];(*a matrix which shuffles the order of the integrals to prioritise the MIs with respect to the the identities at eps=0*)
			
			allGs=shuffle . allGs;
			localibps=localibps . Inverse[shuffle];
			
			(*there may be additional relations between the MIs at eps=0, which are linear combinations of the original ones but they come with an overall power of esp. We need to isolate such IBPs.*)
			
			If[
				(* !MatchQ[reduced,{}] *)
				Length[mis]-nmis>0,(*checking whether we missed linear relations*)

				reduced=RowReduce[localibps];(*IBP reduction with respect the newly introduced order*)
				reduced=reduced[[-Length[mis]+nmis;;]];(*isolating the subset of IBP which relates MIs at eps=0*)
				
				localibps=reduced . allGs;
				localibps=
					If[
						(!FreeQ[#,eps])&& (*it checks whether the IBP depend on eps...*)
							Series[#,{eps,0,-1}][[-4]]!={},(*... and if it is divergent as eps->0*)
						Power[
							eps,
							(-Series[#,{eps,0,-1}][[-3]])(*this is the leading power in eps in the IBP*)(*Series is performed twice, so this is not an optimal implementation, but this is not the bottleneck anyway.*)
						]*#,
						#(*if it is convergent as eps->0, then it does not touch the relation*)
					]&/@localibps;
			];
			
			allGs=allGs[[-Length[mis];;]](*we restricts to the subset of integrals which are MIs at eps=0 and run the algorithm again*)
		];
		
		Return[allGs]
	];
