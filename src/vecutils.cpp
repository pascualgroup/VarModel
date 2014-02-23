#include <cassert>
#include "vecutils.h"

using namespace std;

std::vector<size_t> makeRange(size_t size)
{
	return makeRange(0, size);
}

std::vector<size_t> makeRange(size_t from, size_t to)
{
	assert(to >= from);
	size_t size = to - from;
	vector<size_t> range(size);
	for(size_t i = 0; i < size; i++)
	{
		range[i] = from + i;
	}
	return range;
}
