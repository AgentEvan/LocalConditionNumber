    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	// const volVectorField & PointsC = mesh.C();
	// const auto & Volume = mesh.V();

    Info<< "Read field RS_RANS ...\n" << endl;
    volSymmTensorField RS_RANS
    (
        IOobject
        (  
            "turbulenceProperties:R",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Read field U ...\n" << endl;
    volVectorField U
    (
        IOobject
        (  
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	Info<< "Read dictionary transportProperties ...\n" << endl;
	IOdictionary transportProperties
	(
		IOobject
		(
			"transportProperties",
			runTime.constant(),
			U.db(),
			IOobject::MUST_READ_IF_MODIFIED,
			IOobject::NO_WRITE
		)
	);

	/*
	Info<< "Create field nu from transportProperties ...\n" << endl;
	volScalarField nu
	(
		IOobject
		(
			"nu",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("nu", dimViscosity, transportProperties)
	);

	Info<< "Read/calculate face flux field phi\n" << endl;
	surfaceScalarField phi
	(
		IOobject
		(
			"phi",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		fvc::flux(U)
	);

    Info<< "Read field p ...\n" << endl;
    volScalarField p
    (
        IOobject
        (  
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Read field RS_DNS ...\n" << endl;
    volSymmTensorField RS_DNS
    (
        IOobject
        (  
            "RS_DNS",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
	
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

    Info<< "Read field RSN_DNS ...\n" << endl;
    volSymmTensorField RSN_DNS
    (
        IOobject
        (  
            "RSN_DNS",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
	
    Info<< "Read field U_DNS ...\n" << endl;
    volVectorField U_DNS
    (
        IOobject
        (  
            "U_DNS",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Read field p_DNS ...\n" << endl;
    volScalarField p_DNS
    (
        IOobject
        (  
            "p_DNS",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Read field nut ...\n" << endl;
    volScalarField nut
    (
        IOobject
        (  
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Read field k ...\n" << endl;
    volScalarField k
    (
        IOobject
        (  
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );	

    Info<< "Read field k_DNS ...\n" << endl;
    volScalarField k_DNS
    (
        IOobject
        (  
            "k_DNS",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Read field epsilon ...\n" << endl;
    volScalarField epsilon
    (
        IOobject
        (  
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Read field grad(p) ...\n" << endl;
    volVectorField gradp
    (
        IOobject
        (  
            "grad(p)",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Read field grad(k) ...\n" << endl;
    volVectorField gradk
    (
        IOobject
        (  
            "grad(k)",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Read field grad(U) ...\n" << endl;
    volTensorField gradU
    (
        IOobject
        (  
            "grad(U)",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
	*/

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
