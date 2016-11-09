#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include "Work.h"
#include "ReadVariables.h"

/**
 Updates:
 25/08/2011
 - Incorporated MPI, if used.
 02/10/2011
 - Made all outputs to vector output_content critical so they are one at a time
 - A much better way to do this would be to open multiple file streams, and pass those around to be written too directly BUT
 I can't be bothered setting that out as it doesn't change the effectiveness of the algorithm. Whoever inherits this might
 want to sort that...
 06/10/2011
 - Also the run time is twice as long having changed to using vector-lists, so definitely not useful for big runs.
 Fortunately overall runtime is still short.
 24/10/2011
 - Moved MPI_INIT and MPI_FINALISE outside brackets
 **/

using namespace std;

void correctInputs()
{
	cout << "Inputs required:" << endl << endl;
	cout << "./zContrast input.file" << endl << endl;
	
}

int main(int argc, char *argv[])
// Here is our new program
// This class purely checks the input and decides what to do from there using the work class
// However it does return all the errors if inputs are incorrect
{
	// MPI  Initiate
	MPI_Init(&argc, &argv);
	
	// CHECK INPUTS //
	if (argc != 2) 
	{
		correctInputs();
		// return EXIT_FAILURE;
	}
	//////////////////
	
	ReadVariables red;
	bool check = red.openFile(argv[1]);
	
	if (check == EXIT_SUCCESS)
	{
		Work work;
		// Pass all variables to program
		work.setVariables(red);
		
		// MPI  Initiate
		// MPI_Init(&argc, &argv);
		
		// Run program as needed now settings are all configured
		work.run();
		
		// MPI Finialise (Stupid American spelling!)
		// MPI_Finalize();

	}
	else 
	{
		cout << "Error in input file. Exiting..." << endl;
		// return EXIT_FAILURE;
	}
	
	// MPI Finialise (Stupid American spelling!)
	MPI_Finalize();
	
	// return EXIT_SUCCESS;
}
