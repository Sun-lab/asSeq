/*
   function returning read name for different platforms
*/
std::vector<std::string> split(const std::string &s, char delim);
int whichDelim(const std::string &s);

const char* const DELIM_LOOKUP = "/ ";
const uint8_t DELIM_SLASH    = 0;
const uint8_t DELIM_SPACE    = 1;
