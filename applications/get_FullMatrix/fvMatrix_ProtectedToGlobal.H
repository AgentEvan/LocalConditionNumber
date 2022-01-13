// Functions below are global version of protected
// member functions in fvMatrix, which appear in 
// an important function -- solveSegregated(...)

#ifndef fvMatrix_ProtectedToGlobal_H_
#define fvMatrix_ProtectedToGlobal_H_

#include "fvCFD.H"
#include "fvMatrix.H"
#include "Field.H"
#include "direction.H"

template <typename Type>
void addToInternalField
(
	const labelUList & addr,
	const Field<Type> & pf,
	Field<Type>& intf
);

template <typename Type>
void addToInternalField
(
	const labelUList & addr,
	const tmp<Field<Type>> & pf,
	Field<Type>& intf
);

template <typename Type>
void addBoundaryDiag
(
	fvMatrix<Type> & Eqn,
	const direction cmpt
);

template <typename Type>
void addBoundarySource
(
	Field<Type> & source,
	const fvMatrix<Type> & Eqn,
	const bool couples=true
);

#endif

