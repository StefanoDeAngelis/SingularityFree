(* ::Package:: *)

permutationMatrix[p_List]:=
	SparseArray[
		IdentityMatrix[Length[p]][[p]]
	] (*in earlier versions of Mathematica, PermutationMatrix was not implemented yet. This is an in-house implementation.*)


gaussReductionMatrix[matrix_]:=
	Block[
		{
			n,m,
			reduced,
			gauss
		},
		
		{n,m}=Dimensions[matrix];(*storing the dimensions of the matrix*)
		
		reduced=ArrayFlatten[{{matrix,IdentityMatrix[n]}}];(*forming a matrix [matrix|identity]*)
		reduced=RowReduce[reduced];(*the identity keeps track of the steps in the Gaussian elimination*)
		
		gauss=reduced[[All, m+1;;]];
		reduced=reduced[[All, ;;m]];

		Return[{reduced,gauss}];
	]


finiteCoefficients[ibps_List,G_(*integrals*),eps_(*dimreg parameter*)]:=
	Block[
		{
			nmis,
			rank,
			kin,(*the kinematic variables*)
			allGs,(*all integrals appearing the IBPs*)
			
			localibps=ibps,
			localibps4D,
			reduced,
			gauss,(*the matrix performing the row reduction on `localibps4D`, giving `reduced` as a result*)
			pmatrix,(*the matrix to eliminate the overall powers of eps after Gauss elimination with gauss*)
			
			reducible,
			mis={},
			shuffle,
			
			count=0
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
		
		localibps=Transpose[Coefficient[localibps,#]&/@allGs];(*the IBP matrix*)
		localibps=Map[Simplify,localibps,{2}];
		localibps4D=localibps/.eps->0;
		
		rank=MatrixRank[localibps];
		
		count=1;
		
		{reduced,gauss}=gaussReductionMatrix[localibps4D];(*IBP reduction with eps=0*)
		
		While[
			MatrixRank[reduced] != rank,(*a while loop until the number of rank of the IBP matrix with eps=0 is the same as the rank of the IBP matrix with eps!=0*)
			
			(*The while loop is needed because, when we put eps=0, we may miss linear combinations of IBPs which are actually proportional to eps*)
			
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
			
			allGs = shuffle . allGs;
			localibps = gauss . localibps . Inverse[shuffle];
			
			(*There are additional relations between the MIs at eps=0, which are linear combinations of the original ones but they come with an overall power of esp.
			 The construction gauss . localibps . Inverse[shuffle] isolates such relations in the last MatrixRank[localibps]-MatrixRank[reduced]-nmis rows.
			 We need multiply by the appropiate power of 1/\[Epsilon] such that, when we set \[Epsilon]=0 in the next iteration, these relations are not discarded.*)
				
			localibps=Map[Simplify,localibps,{2}];(*the matrix product gives sums which should be simplified to identify the overall power of eps in the next step*)
			
			pmatrix = Min@Exponent[DeleteCases[#,0],eps,Min]&/@localibps;(*extracting the least common power of eps*)
			pmatrix = eps^(-pmatrix);
			pmatrix = DiagonalMatrix[pmatrix];
			
			localibps = pmatrix . localibps;
			
			localibps4D = localibps /. eps->0;
			
			{reduced,gauss} = gaussReductionMatrix[localibps4D];(*IBP reduction with eps=0 and Guass-elimination matrix*)
			
			count=count+1;
		];
		
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
			
		allGs = shuffle . allGs;
		
		Return[allGs[[-nmis;;]]]
	];
