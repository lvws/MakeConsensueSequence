#include <string>
#include "reader.h"

using std::string;

string ReverseComplement(const string& str)
{
    string result;
    result.resize(str.size());
    size_t n = str.size() - 1;
    for (const auto& base : str)
    {
        switch (base)
        {
        case 'A':
            result[n] = 'T';
            break;
        case 'T':
            result[n] = 'A';
            break;
        case 'C':
            result[n] = 'G';
            break;
        case 'G':
            result[n] = 'C';
            break;
        default:
            result[n] = 'N';
            break;
        }
        n--;
    }
    return result;
}
