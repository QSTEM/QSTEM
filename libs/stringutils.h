#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include <string>

// Next 3 functions from: http://www.cplusplus.com/forum/beginner/1962/
inline std::string ExtractDirectory( const std::string& path )
  {
  return path.substr( 0, path.find_last_of( "\\/") +1 );
  }

inline std::string ExtractFilename( const std::string& path )
  {
  return path.substr( path.find_last_of( "\\/") +1 );
  }

inline std::string ChangeExtension( const std::string& path, const std::string& ext )
  {
  std::string filename = ExtractFilename( path );
  return ExtractDirectory( path ) +filename.substr( 0, filename.find_last_of( '.' ) ) +ext;
  }

inline std::string StripSpaces( const std::string& path)
{
	std::string ext = path;
	ext.erase(std::remove_if(ext.begin(), ext.end(), isspace), ext.end());
	return ext;
}

// from http://stackoverflow.com/questions/11521914/splitting-a-c-string-with
inline std::string GetExtension( const std::string& path)
{
	std::string ext = path.substr(path.find_last_of(".")+1);
	//if(fn.substr(fn.find_last_of(".") + 1) == "conf")
	//return path.substr(path.rfind(".")+1);
	ext = StripSpaces(ext);
	return ext;
}

#endif