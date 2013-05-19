#include <string>

// Next 3 functions from: http://www.cplusplus.com/forum/beginner/1962/
std::string ExtractDirectory( const std::string& path )
  {
  return path.substr( 0, path.find_last_of( "\\/") +1 );
  }

std::string ExtractFilename( const std::string& path )
  {
  return path.substr( path.find_last_of( "\\/") +1 );
  }

std::string ChangeExtension( const std::string& path, const std::string& ext )
  {
  std::string filename = ExtractFilename( path );
  return ExtractDirectory( path ) +filename.substr( 0, filename.find_last_of( '.' ) ) +ext;
  }

// from http://stackoverflow.com/questions/11521914/splitting-a-c-string-with
std::string GetExtension( const std::string& path)
{
	return path.substr(path.rfind(".")+1);
}