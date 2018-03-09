#include "./all.cc"
const size_t MaxThreads = 4 ;
int main () {
    GenerateQCD slave ;
    #pragma omp parallel for
    for(size_t i=0;i<MaxThreads;i++) {slave.GenSim(i+1);}
    return 0;
}
