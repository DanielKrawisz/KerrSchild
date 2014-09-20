#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

//These functions return true or false depending on whether an
//  input or output file was sucessfully opened to a stream.
//TODO: switch to exception handling. 
bool OutputFileOpen(ofstream& out, string location);

bool InputFileOpen(ifstream& in, string location);

//Turn an int into a string. 
string stringify(int x);
