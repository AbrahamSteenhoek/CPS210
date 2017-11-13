#include <iostream>
#include <iostream>
#include <fstream>

using namespace std;

class CharToBits {
public:
	unsigned char *c;
	int shift;
	string bitString;
    CharToBits(string _bstr) {
        bitString = _bstr;
        shift = 0;
        c = (unsigned char*)calloc(1, sizeof(unsigned char));
    }
    
    void write(ofstream& out) {
        out << *c;
    }

	int insertBits(const char* filename) {
        ofstream out(filename);

        int totalBits = 0;
        int i = 0;
        while(i < bitString.size()) {
            if(bitString[i] == '1') {
                *c |= 1;
            }
            
            if(shift == 7) {
                shift = 0;
                write(out);
                free(c);
                c = (unsigned char*)calloc(1, sizeof(unsigned char));
            } else {
                *c <<= 1;
                ++shift;
                ++totalBits;
            }
            ++i;
        }

        return totalBits;
    }

    ~CharToBits() {
        if(c)
            free(c);
    }
	
};