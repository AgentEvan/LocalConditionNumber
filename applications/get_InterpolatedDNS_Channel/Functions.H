#ifndef Functions_H_
#define Functions_H_

#include <string>
#include <vector>
#include <fstream>
#include "fvCFD.H"

typedef std::vector<Foam::scalar> data1D;
typedef std::vector<data1D> data2D;

int get_ActualLines(const std::string &file);
void Print_data1D(const data1D &);

#endif
