#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cctype>


inline std::string trim(
        const std::string& str,
        const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
    return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

inline std::string reduce(
        const std::string& str,
        const std::string& fill = " ",
        const std::string& whitespace = " \t")
{
    // trim first
    auto result = trim(str, whitespace);

    // replace sub ranges
    auto beginSpace = result.find_first_of(whitespace);
    while (beginSpace != std::string::npos)
    {
        const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
        const auto range = endSpace - beginSpace;

        result.replace(beginSpace, range, fill);

        const auto newStart = beginSpace + fill.length();
        beginSpace = result.find_first_of(whitespace, newStart);
    }
    return result;
}


//==================================== put into a new module, e.g., "fileUtils" (or "strUtils") ===============================================
// 150626: Hmm... We want, e.g., 'test\r\n' to be equal to 'test\r' (or was it \r and \n swapped?) but we don't want
//         'test' to be equal to '' or '\0'. The old version would do that.
//         Will this version work as we intend? Note that this is not order-independent wrt. arguments!!
//         Intended use: strMatchLeading( string_with_tagged_data, prefix_string )

inline bool strMatchLeading( const std::string &str, const std::string &prefix )
{
    if ( str.length() < prefix.length() ) {
        return false;
    }
    // assert( prefix.length() > 0 ); // Anything else doesn't really make sense.
    if ( prefix.length() == 0 ) {
        // But this should also be ok, and never cause a crash...
        return ( str.length() == 0 );
    }
    const size_t n = std::min( str.length(), prefix.length() );
    return ( prefix.substr(0, n) == str.substr(0, n) );
}

inline bool strMatchFull( const std::string &a, const std::string &b )
{
    const size_t n = std::min( a.length(), b.length() );
    if (n==0) {
        return ( a.length() == b.length() );
    }
    return (  ( a.length() == b.length() )  &&  ( a.substr(0, n) == b.substr(0, n) )  );
}



inline std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}





// If called without "init_line", it will always read at least one line from the file.
// If it should test on an already read line without necessarily reading more, call it with that line.
inline std::string scanForToken( std::istream &is, const std::string &str, const std::string init_line ="" )
{
    std::string line = init_line;
    while ( (!strMatchLeading(line, str)) && (is.good()) ) {
        getline( is, line );
        line = reduce(line);
    }
    if (is.good()) {
        return line;
    } else {
        return "";
    }
}

// Not sure about what this returns, but it seems customary to test on it for "goodness"...
inline std::istream& getline2( std::istream& is, std::string& str, std::string &prev_str, size_t &line_cntr )
{
    //        static int l = 0;
    //        l++;
    //        std::cout << "l=" << l << std::endl;
    line_cntr++;
    prev_str = str;
    return getline(is, str);
}


// Regarding occurences of tests like 's.substr(0, std::min(int(s.length()), 2)) == "v "' below:
// This is a way of getting Linux/Windows (i.e., LF/LF+CR) compatibility.


#if 0
inline bool strMatch( const std::string &a, const std::string &b )
{
    const size_t n = std::min( a.length(), b.length() );
    return ( a.substr(0, n) == b.substr(0, n) );
}
#endif




#if 0
inline bool strMatchFull( const std::string &a, const std::string &b )
{
    const size_t n = std::min( a.length(), b.length() );
    if (n==0) {
        return ( a.length() == b.length() );
    }
    return (  ( a.length() == b.length() )  &&  ( a.substr(0, n) == b.substr(0, n) )  );
}
#endif


// Definition of obj-file-comment
inline bool isComment(const std::string &s) {
    return ( (s.length()>0) && (s[0]=='#') );
}


inline std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos)
    {
        return filename;
    }
    return filename.substr(0, lastdot);
}

inline std::string get_path(const std::string& filename) {
    size_t lastslash = filename.find_last_of('/');
    if (lastslash == std::string::npos)
    {
        return filename;
    }
    return filename.substr(0, lastslash);
}

inline void trimLeft(std::istream &is) {
	int c = is.peek();
	while (std::isspace(c)) {
		is.ignore(1);
		c = is.peek();
	}
}


// Function reading a token from a potentially binary string. 
// A token contains printable, non-space characters
inline std::string readToken(std::istream &is) {
	int c = is.peek();
	std::string token;
	while (is.good() && std::isprint(c) && !std::isspace(c)) {
		char curr;
		is.read(&curr,1);
		token.append(1, curr);
		c = is.peek();
	}
	return token;
}

inline std::string readPrintableLine(std::istream &is) {
	int c = is.peek();
	std::string line;
	while (is.good() && std::isprint(c)) {
		char curr;
		is.read(&curr, 1);
		line.append(1, curr);
		c = is.peek();
	}
	return line;
}