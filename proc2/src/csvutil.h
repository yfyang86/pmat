#ifdef __HEAD__CSVUTIL__
#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

template <class T> std::vector<std::pair<std::string, std::vector<T>>> read_csv(std::string filename);

#define __HEAD_CSVUTIL__

#endif

