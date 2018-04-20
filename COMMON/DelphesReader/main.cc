#include "reader.cc"

int BenchmarkVectors () {
    for(double x=-1000.0;x<1000.0;x=x+1.0){
        for(double y=-1000;y<1000;y=y+1.0){
            for(double z=-1000;z<1000;z=z+1.0){
                double t = sqrt((x*x)+(y*y)+(z*z)+(1.0)) ;
                NewHEPHeaders::VECTORS::lorentz4vector <> myvector(x,y,z,t) ;
                NewHEPHeaders::VECTORS::lorentz4vector <> newvector ;
                newvector.SetPtEtaPhiM(myvector.pt(),myvector.eta(),myvector.phi(),myvector.m());
                double error =
                CPPFileIO::mymod(x-newvector[0])+
                CPPFileIO::mymod(y-newvector[1])+
                CPPFileIO::mymod(z-newvector[2])+
                CPPFileIO::mymod(t-newvector[3]);
                if(error>0.0001){printf("%e\n",error);}
            }
        }
    }
    return 0;

    for(double x=-5.0;x<5.0;x=x+0.15){
        for(double y=-5.0;y<5.0;y=y+0.15){
            for(double z=-5.0;z<5.0;z=z+0.15){
                NewHEPHeaders::VECTORS::euclid3vector <> myvector(x,y,z) ;
                NewHEPHeaders::VECTORS::euclid3vector <> newvector ;
                newvector.SetPtEtaPhi(myvector.pt(),myvector.eta(),myvector.phi());
                printf("(%e:%e) (%e:%e) (%e:%e)\n",x,x-newvector[0],y,y-newvector[1],z,z-newvector[2]);
            }
        }
    }
    return 0;
}

inline void ReadRootFilesAndAnalyzeWithMPI () {
    MainAnalyzer a("GraphsWithMPI");
    a("./BoostedZ/4/out.root");
    a("./BoostedZ/1/out.root");
    a("./BoostedZ/3/out.root");
    a("./BoostedZ/2/out.root");
    a("./BoostedZ/15/out.root");
    a("./BoostedZ/9/out.root");
    a("./BoostedZ/5/out.root");
    a("./BoostedZ/13/out.root");
    a("./BoostedZ/12/out.root");
    a("./BoostedZ/7/out.root");
    a("./BoostedZ/8/out.root");
    a("./BoostedZ/10/out.root");
    a("./BoostedZ/6/out.root");
    a("./BoostedZ/11/out.root");
    a("./BoostedZ/16/out.root");
    a("./BoostedZ/14/out.root");
}

inline void ReadRootFilesAndAnalyzeNoMPI () {
    MainAnalyzer a("GraphsNoMPI");
    a("./BoostedZ/4/NoISRout.root");
    a("./BoostedZ/1/NoISRout.root");
    a("./BoostedZ/3/NoISRout.root");
    a("./BoostedZ/2/NoISRout.root");
    a("./BoostedZ/15/NoISRout.root");
    a("./BoostedZ/9/NoISRout.root");
    a("./BoostedZ/5/NoISRout.root");
    a("./BoostedZ/13/NoISRout.root");
    a("./BoostedZ/12/NoISRout.root");
    a("./BoostedZ/7/NoISRout.root");
    a("./BoostedZ/8/NoISRout.root");
    a("./BoostedZ/10/NoISRout.root");
    a("./BoostedZ/6/NoISRout.root");
    a("./BoostedZ/11/NoISRout.root");
    a("./BoostedZ/16/NoISRout.root");
    a("./BoostedZ/14/NoISRout.root");
}

int main () {
    CPPFileIO::ForkMe forker ;
    if(forker.InKid()){ReadRootFilesAndAnalyzeWithMPI();}
    if(forker.InKid()){ReadRootFilesAndAnalyzeNoMPI();}
    return 0;
}
