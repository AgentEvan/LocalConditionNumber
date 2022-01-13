//
// Useful Functions

#include "Functions.H"

// Counting total number of non-empty lines of one single text file

int get_ActualLines(const std::string &file)
{
	std::ifstream fin;
	fin.open(file, std::ios_base::in);
	if (!fin.is_open())
	{
		cerr << "Failed in opening file_dns!!\n";
		exit(EXIT_FAILURE);
	}
	
	int N_lines = 0, N_empty = 0;
	std::string temp;
	while (getline(fin, temp, '\n'))
	{
		N_lines++;
		if (temp[0]=='\0' || temp.size()==1) 
		{
			N_empty++;
			cout << "Line " << N_lines << " is empty!!\n";
		}

		// if (N_lines == 8712)
		//	  cout << "Length of Line " << N_lines << ": " << temp.size() << endl;
	}
	fin.close();

	std::cout << "There are total " << N_lines << " lines in " << file << std::endl;
	std::cout << "There are " << N_empty << " empty lines in " << file << std::endl;

	return (N_lines-N_empty);
}

void Print_data1D(const data1D & vec)
{
	data1D::const_iterator beg = vec.begin();
	for (; beg!=vec.end(); beg++)
		std::cout << *beg << std::endl;
	std::cout << std::endl;
}

//*******************************************************************************//

