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
    get_InterpolatedDNS_Channel

Description

\*---------------------------------------------------------------------------*/

#include <string>
#include <fstream>
#include <vector>
#include "fvCFD.H"
#include "Functions.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
	#include "createMesh.H"
	#include "setToLatestTime.H"
	#include "readFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	

	// std::string caseName("0180");
	word caseName(transportProperties.lookup("Re_tau"));
	Info<< "Now the case is Re" << caseName << endl;

	// 6 columns in fileVel
	// y/delta		y^+		U		dUdy		W		P
	std::string fileVel("/data/guoxw/Data/DNS/Channel/data/Velocity_"+caseName);

	// 9 columns in fileRS
	// y/delta		y^+		uu    vv    ww    uv    uw    vw    k	
	std::string fileRS("/data/guoxw/Data/DNS/Channel/data/ReynoldsStress_"+caseName);
	
	int N_dns = get_ActualLines(fileVel);
	data1D temp1(6, Zero), temp2(9, Zero);
	data2D dataVel(N_dns, temp1), dataRS(N_dns, temp2);
	
	std::ifstream fin;
	fin.open(fileVel, std::ios_base::in);
	for(int i=0; i<N_dns; i++)
	{
		for(int j=0; j<6; j++)
			fin >> dataVel[i][j];
	}
	fin.close();

	fin.open(fileRS, std::ios_base::in);
	for(int i=0; i<N_dns; i++)
	{
		for(int j=0; j<9; j++)
			fin >> dataRS[i][j];
	}
	fin.close();

	scalarField y_DNS(N_dns, Zero);
	vectorField U_DNS_Field(N_dns, Zero);
	symmTensorField RS_DNS_Field(N_dns, Zero);

	forAll(y_DNS, ni)
	{
		y_DNS[ni] = dataVel[ni][0];

		U_DNS_Field[ni].x() = dataVel[ni][2];
		U_DNS_Field[ni].z() = dataVel[ni][4];

		RS_DNS_Field[ni].xx() = dataRS[ni][2];
		RS_DNS_Field[ni].yy() = dataRS[ni][3];
		RS_DNS_Field[ni].zz() = dataRS[ni][4];
		RS_DNS_Field[ni].xy() = dataRS[ni][5];
		RS_DNS_Field[ni].xz() = dataRS[ni][6];
		RS_DNS_Field[ni].yz() = dataRS[ni][7];
	}

	volVectorField U_DNS("U_DNS", U);
	volSymmTensorField RS_DNS("RS_DNS", RS_RANS);

	forAll(U_DNS, cellI)
	{
		const scalar y = mesh.C()[cellI].y();
	
		U_DNS.primitiveFieldRef()[cellI] = interpolateXY(y, y_DNS, U_DNS_Field);
		RS_DNS.primitiveFieldRef()[cellI] = interpolateXY(y, y_DNS, RS_DNS_Field);

		/*
		// For 2h case
		if (y<=0.0 && y>=-1.0)
		{
			U_DNS.primitiveFieldRef()[cellI] = interpolateXY(y+1, y_DNS, U_DNS_Field);
			RS_DNS.primitiveFieldRef()[cellI] = interpolateXY(y+1, y_DNS, RS_DNS_Field);
		}
		else
		{
			U_DNS.primitiveFieldRef()[cellI] = interpolateXY(1-y, y_DNS, U_DNS_Field);
			RS_DNS.primitiveFieldRef()[cellI] = interpolateXY(1-y, y_DNS, RS_DNS_Field);
		}
		*/
	}
	
	U_DNS.write();
	RS_DNS.write();

	volScalarField k_DNS("k_DNS", tr(RS_DNS));
	k_DNS.write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
