#include "strutils.h"

std::string strprintf(char const * format, ...)
{
	va_list args;
	
	std::unique_ptr<char[]> buf(new char[0]);
	va_start(args, format);
	int size = vsnprintf(buf.get(), 0, format, args) + 1;
	va_end(args);
	if(size <= 0) {
		throw std::runtime_error("Encoding error occurred in vsnprintf.");
	}
	
	buf = std::unique_ptr<char[]>(new char[size]);
	va_start(args, format);
	vsnprintf(buf.get(), size, format, args);
	va_end(args);
	
	return std::string(buf.get());
}

std::string strjoin(std::vector<std::string> const & strings, std::string sep)
{
	if(strings.size() == 0) {
		return "";
	}
	
	std::string str;
	
	for(auto itr = strings.begin(); itr != (strings.end() - 1); ++itr) {
		str += *itr;
		str += sep;
	}
	str += *(strings.end() - 1);
	return str;
}

// From Stack Overflow user simfoo:
// http://stackoverflow.com/questions/7724448/simple-json-string-escape-for-c
std::string escapeJsonString(std::string const & input)
{
    std::ostringstream ss;
    for (auto iter = input.cbegin(); iter != input.cend(); iter++) {
    //C++98/03:
    //for (std::string::const_iterator iter = input.begin(); iter != input.end(); iter++) {
        switch (*iter) {
            case '\\': ss << "\\\\"; break;
            case '"': ss << "\\\""; break;
            case '/': ss << "\\/"; break;
            case '\b': ss << "\\b"; break;
            case '\f': ss << "\\f"; break;
            case '\n': ss << "\\n"; break;
            case '\r': ss << "\\r"; break;
            case '\t': ss << "\\t"; break;
            default: ss << *iter; break;
        }
    }
    return ss.str();
}
