/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 AUTHOR,AFFILIATION
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
	get_FullMatrix

Description

	Explicit:  nu,             RS_DNS
	Implicit1: nu + nutL_DNS,  RS_DNS
	Implicit2: nu + nutL_RANS, RS_DNS
	Implicit3: nu + nut,       RS_DNS

	DynamicImplicitFoams:
		Imp1:	Baseline RANS correction
		Imp2:	Updated optimal Nut
		Imp3:	Linear version of Imp2
		Imp4:	Fixed optimal Nut
		Imp5:	Rescaled RANS Nut with a fixed factor C
		Imp6:	Rescaled RANS Nut with a updated factor C
		Imp7:	Rescaled RANS Nut with a fixed given factor C
				(provided as "FactorC" in constant/transportProperties)

\*---------------------------------------------------------------------------*/

#include <vector>
#include <iomanip>
#include <fstream>
#include "fvCFD.H"
#include "fvMatrix_ProtectedToGlobal.H"

typedef std::vector<Foam::scalar> data1D;
typedef std::vector<data1D> data2D;

template<typename Data> 
void writeFoamData1D(const Data& data, std::string fileName);

void writeData2D(const data2D &, std::string);

// create a copy of a transposed fvMatrix object
template<typename Type>
fvMatrix<Type> getTranspose(const fvMatrix<Type> & Mij);

template<typename Type>
void writeSparseFvMatrix
(
	const fvMatrix<Type> &, 
	const std::string fileName="A_Sparse",
	const direction cmpt=0
);

template <typename Type>
void correctBoundaryDiag(fvMatrix<Type> &, const direction cmpt);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
	#include "createMesh.H"
	// #include "setToLatestTime.H"
	#include "readFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	word propType(transportProperties.lookup("Prop_Type"));
	Info<< "Propagation type is: " << propType << endl << endl;

	if (propType == "Explicit_DNS")
	{
		#include "readForDNSProp.H"
		tmp<fvVectorMatrix> tUEqn
		(
			fvm::div(phi_DNS, U_DNS)
		  - fvm::laplacian(nu, U_DNS) // explicit treatment
		  + fvc::div(RS_DNS)
		);

		fvVectorMatrix & UEqn = tUEqn.ref();
		#include "writeMatrix.H"
	}

	else if (propType == "Implicit1_DNS")
	{
		#include "readForDNSProp.H"
		Info<< "Read field nutL_DNS ...\n" << endl;
		volScalarField nutL_DNS
		(
			IOobject
			(  
				"nutL_DNS",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);
		tmp<fvVectorMatrix> tUEqn
		(
			fvm::div(phi_DNS, U_DNS)
		  - fvm::laplacian(nu+nutL_DNS, U_DNS) // implicit treatment 1
		  + fvc::div(RS_DNS)
		);

		fvVectorMatrix & UEqn = tUEqn.ref();
		#include "writeMatrix.H"
	}

	else if (propType == "Implicit2_DNS")
	{
		#include "readForDNSProp.H"
		Info<< "Read field nutL_RANS ...\n" << endl;
		volScalarField nutL_RANS
		(
			IOobject
			(  
				"nutL_RANS",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);
		tmp<fvVectorMatrix> tUEqn
		(
			fvm::div(phi_DNS, U_DNS)
		  - fvm::laplacian(nu+nutL_RANS, U_DNS) // implicit treatment 2
		  + fvc::div(RS_DNS)
		);

		fvVectorMatrix & UEqn = tUEqn.ref();
		#include "writeMatrix.H"
	}

	else if (propType == "Implicit3_DNS")
	{
		#include "readForDNSProp.H"
		tmp<fvVectorMatrix> tUEqn
		(
			fvm::div(phi_DNS, U_DNS)
		  - fvm::laplacian(nu+nut, U_DNS) // implicit treatment 3
		  + fvc::div(RS_DNS)
		);

		fvVectorMatrix & UEqn = tUEqn.ref();
		#include "writeMatrix.H"
	}

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	else if (propType == "Implicit1")
	{
		Info<< "Read field nutL_DNS ...\n" << endl;
		volScalarField nutL_DNS
		(
			IOobject
			(  
				"nutL_DNS",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);
		tmp<fvVectorMatrix> tUEqn
		(
			fvm::div(phi, U)
		  - fvm::laplacian(nu+nutL_DNS, U) // implicit treatment 1
		  + fvc::div(RS_DNS)
		);

		fvVectorMatrix & UEqn = tUEqn.ref();
		#include "writeMatrix.H"
	}

	else if (propType == "Implicit2")
	{
		Info<< "Read field nutL_RANS ...\n" << endl;
		volScalarField nutL_RANS
		(
			IOobject
			(  
				"nutL_RANS",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);
		tmp<fvVectorMatrix> tUEqn
		(
			fvm::div(phi, U)
		  - fvm::laplacian(nu+nutL_RANS, U) // implicit treatment 2
		  + fvc::div(RS_DNS)
		);

		fvVectorMatrix & UEqn = tUEqn.ref();
		#include "writeMatrix.H"
	}

	else if (propType == "Implicit3")
	{
		tmp<fvVectorMatrix> tUEqn
		(
			fvm::div(phi, U)
		  - fvm::laplacian(nu+nut, U) // implicit treatment 3
		  + fvc::div(RS_DNS)
		);

		fvVectorMatrix & UEqn = tUEqn.ref();
		#include "writeMatrix.H"
	}

	else if (propType == "Imp7")
	{
		word FactorCWord(transportProperties.lookup("FactorC"));
		scalar C = readScalar(FactorCWord);
		volScalarField nut_New("nut_New", C*nut);

		tmp<fvVectorMatrix> tUEqn
		(
			fvm::div(phi, U)
		  - fvm::laplacian(nu+nut_New, U) // dynamicImplicitRSFoam7
		  + fvc::div(RS_DNS)
		);

		fvVectorMatrix & UEqn = tUEqn.ref();
		#include "writeMatrix.H"
	}

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
// ********************    Function Definitions   ************************** //
// ************************************************************************* //

template<typename Data> 
void writeFoamData1D(const Data& data, std::string fileName)
{
	// std::ofstream fout;
	// fout.open(fileName, std::ios_base::out);
	OFstream fout(fileName);
	// Foam::label N = data.size();
	forAll(data, ni)
	{
		fout << data[ni] << nl; // std::endl;
	}
	// fout.close();
}

template<typename Type>
fvMatrix<Type> getTranspose(const fvMatrix<Type> & Mij)
{
	fvMatrix<Type> MTij(Mij); // create a copy
	scalarField lower_temp = MTij.lower();
	MTij.lower() = MTij.upper();
	MTij.upper() = lower_temp;
	return MTij;
}

void writeData2D(const data2D & data, std::string fileName)
{
	// std::ofstream fout;
	// fout.open(fileName, std::ios_base::out);
	OFstream fout(fileName);
	for (auto A2D=data.cbegin(); A2D!=data.cend(); A2D++)
	{
		data1D Atemp = *A2D;
		for (auto A1D=Atemp.cbegin(); A1D!=Atemp.cend(); A1D++)
			fout << *A1D << ' ';
		fout << nl;
	}
	// fout.close();
}

// Functions below are global version of protected
// member functions in fvMatrix, which appear in 
// an important function -- solveSegregated(...)

template <typename Type>
void addToInternalField
(
	const labelUList & addr,
	const Field<Type> & pf,
	Field<Type>& intf
)
{
	if (addr.size() != pf.size())
	{
		FatalErrorInFunction
			<< "sizes of addressing and field are different"
			<< abort(FatalError);
	}

	forAll(addr, facei)
	{
		intf[addr[facei]] += pf[facei];
	}
}

template <typename Type>
void addToInternalField
(
	const labelUList & addr,
	const tmp<Field<Type>> & tpf,
	Field<Type>& intf
)
{
	addToInternalField(addr, tpf(), intf);
	tpf.clear();
}

template <typename Type>
void addBoundaryDiag
(
	fvMatrix<Type> & Eqn,
	const direction cmpt
)
{
	forAll(Eqn.internalCoeffs(), patchi)
	{
		addToInternalField
		(
			Eqn.lduAddr().patchAddr(patchi),
			Eqn.internalCoeffs()[patchi].component(cmpt),
			Eqn.diag()
		);
	}
}

template <typename Type>
void addBoundarySource
(
	Field<Type> & source,
	const fvMatrix<Type> & Eqn,
	const bool couples
)
{
	forAll(Eqn.psi().boundaryField(), patchi)
	{
		const fvPatchField<Type>& ptf = Eqn.psi().boundaryField()[patchi];
		const Field<Type>& pbc = Eqn.boundaryCoeffs()[patchi];

		if (!ptf.coupled())
		{
			addToInternalField(Eqn.lduAddr().patchAddr(patchi), pbc, source);
		}
		else if (couples)
		{
			const tmp<Field<Type>> tpnf = ptf.patchNeighbourField();
			const Field<Type>& pnf = tpnf();
			
			const labelUList& addr = Eqn.lduAddr().patchAddr(patchi);

			forAll(addr, facei)
			{
				source[addr[facei]] += cmptMultiply(pbc[facei], pnf[facei]);
			}
		}
	}
}

// ************************************************************************* //

template<typename Type>
void writeSparseFvMatrix
(
	const fvMatrix<Type> & UEqn, 
	const std::string fileName,
	const direction cmpt
)
{
	const label nCells = UEqn.diag().size();
	const label nFaces = UEqn.upper().size();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	// * * * * * * * * * Write fvMatrix UEqn in Sparse form  * * * * * * * * //
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	// std::ofstream fout;
	// fout.open(fileName, std::ios_base::out);
	OFstream fout(fileName);

	for (label cell=0; cell<nCells; cell++)
	{
		fout << cell << ' ' << cell << ' ' 
			 << UEqn.DD().ref()[cell].component(cmpt) << nl;
			 // << UEqn.diag()[cell] << nl;
	}
	
	for (label face=0; face<nFaces; face++)
	{
		const label indexI = UEqn.lduAddr().lowerAddr()[face];
		const label indexJ = UEqn.lduAddr().upperAddr()[face];

		// i<j corresponds to the upper part
		fout << indexI << ' ' << indexJ << ' '
			 << UEqn.upper()[face] << nl;

		// i>j corresponds to the lower part
		fout << indexJ << ' ' << indexI << ' '
			 << UEqn.lower()[face] << nl;
	}
	
	// fout.close();
	Info<< "Transform fvMatrix to the sparse form completed..." << endl;

}

template <typename Type>
void correctBoundaryDiag(fvMatrix<Type> & Eqn, const direction cmpt)
{
	Foam::tmp<Foam::Field<Type>> tdiag = Eqn.DD();
	Foam::Field<Type> & DD = tdiag.ref();
	Eqn.diag() = DD.component(cmpt);
}


// ************************************************************************* //
