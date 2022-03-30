#ifndef preprocessing_h    // To make sure you don't declare the function more than once by including the header multi>
#define preprocessing_h
//declare creating db and refs function
#include <vector>
#include <map>
//vector<vector <int>> createdb(int n , const char* filepath );

std::map <int, std::vector<std::vector<uint8_t>> >  createRefs(int n , string filepath );
std::vector<std::vector<uint8_t>> createRefs_spike(int n , string filepath );
std::vector<std::vector<uint8_t>> createdb(int n , std::vector<string> );
#endif