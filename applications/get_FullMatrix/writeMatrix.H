
	Info<< "UEqn.symmetric(): " << UEqn.symmetric() << endl;

	const direction cmpt = 0; // x-direction
	const label nCells = UEqn.diag().size();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	// * * * * * * * * * *    Some Preparations    * * * * * * * * * * * * * //
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	// Copy some field for further necessary parameters
	// The following codes appear in fvMatrixSolve.C line 155
	// in the definition of function solveSegregated(...)
	
	Field<vector> source = UEqn.source();
	addBoundarySource(source, UEqn);

	scalarField psiCmpt(UEqn.psi().primitiveField().component(cmpt));
	scalarField sourceCmpt(source.component(cmpt));

	// ** addBoundary(UEqn.diag(), cmpt) -- Line 109 in fvMatrix.C
	// ** This function will change diag()

	addBoundaryDiag(UEqn, cmpt);

	FieldField<Field, scalar> bouCoeffsCmpt
	(
		UEqn.boundaryCoeffs().component(cmpt)
	);

	FieldField<Field, scalar> intCoeffsCmpt
	(
		UEqn.internalCoeffs().component(cmpt)
	);

	lduInterfaceFieldPtrsList interfaces = 
		UEqn.psi().boundaryField().scalarInterfaces();

	// Correct bouCoeffsCmpt for explicit part of 
	// the coupled boundary conditions
	UEqn.initMatrixInterfaces
	(
		true,
		bouCoeffsCmpt,
		interfaces,
		psiCmpt,
		sourceCmpt,
		cmpt
	);

	UEqn.updateMatrixInterfaces
	(
		true,
		bouCoeffsCmpt,
		interfaces,
		psiCmpt,
		sourceCmpt,
		cmpt
	);

	// writeFoamData1D(sourceCmpt, "sourceCmpt");

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	// *  Calculate fvMatrix Aij in full matrix form (could be large in size)  * //
	// *  BUT to write fvMatrix Aij in sparse COO matrix form (memory-saving)  * //
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	data1D temp1D(nCells, Zero);
	data2D Aij(nCells, temp1D);
	scalarField Aj(nCells, Zero); // The j-th column of the Matrix A_ij

	OFstream fsp(std::string("./A_Sparse_fromFull"));

	forAll(UEqn.diag(), cellI)
	{
		// A vector which all elements are zero except that the j-th element is 1
		scalarField Ij(nCells, Zero); 
		Ij[cellI] = 1.0; 
		UEqn.Amul(Aj, Ij, bouCoeffsCmpt, interfaces, cmpt);

		forAll(Aj, j)
		{
			Aij[j][cellI] = Aj[j];
			if (Aj[j] != 0)
			{
				fsp << j << ' ' << cellI << ' ' << Aj[j] << endl;
			}
		}
	}

	Info<< "Write fvMatrix to the sparse COO form of A_Sparse_fromFull completed..." << endl;
	// Info<< "Transform fvMatrix to the form of A_Full completed..." << endl;
	// writeData2D(Aij, "./A_Full");
