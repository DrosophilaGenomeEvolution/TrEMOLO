#include <iostream>
#include <string>
/*
    Calcul number of chain
    For ID TrEMOLO
*/
using namespace std;
int main() {

    string chain="";
    for (string line; getline(cin, line);) {
        chain=chain+line;
    }

    int ind = 0;
    float s = 0;

    char ascii[4];
    for (int i = 0; i < chain.size(); ++i)
    {
        sprintf(ascii, "%d", chain[i]);
        s   += ( stoi(ascii) + ( (float)(7 + ind)/ (float)(5 - stoi(ascii)) ) );
        ind += 1;
    }
    s = ( (float)((float)s - (int)s) * 100000 ) + (int)s;
    printf("%.0f\n", s);
    return 0;
}
