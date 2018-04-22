#include "stdio.h"

inline void printer (size_t i) {
    printf("('mkdir' '%ld');\n",i);
    printf("('mv' 'BoostedZ-%ld.lhe' '%ld/');\n",i,i);
    printf("('cp' './LHENoISRMPI.in'  '%ld/');\n",i);
}

int main () {
    for(size_t i=1;i<=16;i++){
        printer(i);
    }
    return 0;
}
