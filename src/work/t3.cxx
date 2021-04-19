#include <iostream>
#include "TRandom.h"

using namespace std;

bool foo() {
   return 1;
}
int main(int argc, char **argv) 
{
   gRandom->SetSeed(0) ;
   for (int i=0; i < 10; i++)  
      cout << gRandom->Rndm()  << endl;

   foo() ;
}
