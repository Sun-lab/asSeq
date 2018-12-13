/*
   function returning read name for different platforms
*/

#include <stdint.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
using namespace std;
#include "alName.h"

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

int whichDelim(const std::string &s){
	size_t delim=0;
  size_t pos=-1;
	
  //cout << pos<<endl;
	pos=s.find(DELIM_LOOKUP[DELIM_SLASH]);
	//cout << pos<<endl;
	if(pos!=-1)delim=0;
	pos=-1;
	pos=s.find(DELIM_LOOKUP[DELIM_SPACE]);
	//cout << pos<<endl;
	
  if(pos!=-1)delim=1;
	return delim;
}
