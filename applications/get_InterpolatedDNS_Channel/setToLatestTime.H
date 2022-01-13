	// Set runTime to the latestTime
	// According to $FOAM_SRC/OpenFOAM/include/checkTimeOptions.H
	
	Info<< "Set runTime to latestTime ...\n" << endl;

	instantList Times = runTime.times(); 
	
	Foam::label startTime = 0;

	startTime = Times.size() - 1;

	runTime.setTime(Times[startTime], startTime);

	Info<< "Now runTime.timeName() is " << runTime.timeName() << nl << nl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
