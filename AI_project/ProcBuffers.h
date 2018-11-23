#include "asiolib.h"

class ProcBuffers
{
private:

public:
	ProcBuffers();
	~ProcBuffers();
	static int Process(double**input, double**output);
};