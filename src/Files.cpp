#include "Files.h"

//using namespace std;

//These functions return true or false depending on whether an
//  input or output file was sucessfully opened to a stream.
bool OutputFileOpen(ofstream& out, string location)
{	out.open(location.c_str());
	if(!out.is_open())
	{
		return false;
	}
	return true;
}

bool InputFileOpen(ifstream& in, string location)
{	in.open(location.c_str());
	if(!in.is_open())
	{
		return false;
	}
	return true;
}

//Turn an int into a string. 
string stringify(int x)
{
	stringstream  jorf;
	jorf << x;
	return jorf.str();
}
