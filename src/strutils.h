#ifndef __malariamodel__strutils__
#define __malariamodel__strutils__

#include <vector>
#include <string>
#include <sstream>

std::string strprintf(char const * format, ...);
std::string strjoin(std::vector<std::string> const & strings, std::string sep);
std::string escapeJsonString(std::string const & input);

#endif
