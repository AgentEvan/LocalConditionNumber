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
	Info<< "Read/calculate face flux field phi_DNS\n" << endl;
	surfaceScalarField phi_DNS
	(
		IOobject
		(
			"phi_DNS",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		fvc::flux(U_DNS)
	);

