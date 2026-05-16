
//! o4_dispersion provides services to compute the modified dispersion curve
//! and modified distribution function.
class o4_dispersion
{
  public: // private
    double fT ;  //! Temperature
    double fLambda ; //!Cutoff in GeV
    double fMvacuum; //! Vacuum pion mass in GeV
    double fRatio ; //! Ratio of chiral condensate to vacuum one
  public:
    o4_dispersion(const double &zmin=, const double &Tpc=0.155) ;
    //! returns the modified distribution function
    double f(const double &p, const double &T) ;
    //! returns the vacuum distribution function
    double fvac(const double &p, const double &T) ;
    //! returns the modified dispersion curve E(p) in units of GeV
    double eofp(const double &p, const double &T) ;
    //! returns the modified pole mass squared in units of GeV**2j
    double m2pole(const double &p, const double &T) ;
    //! returns the chiral velocity  
    double v2(const double &p, const double &T) ;
    //!
    double zofT(const double &T)
    //! Returns the scaling function
    double fG(const double &z) ;

    //! Makes a plot of the dispersion curve
    void plot(const std::string &filename="pidispersion.out")  ;
} ;

o4_dispersion::o4_dispersion(const double &lambda_by_T, const double &ratio) :
  fT(0.155), fLambda(lambda_by_T*0.155), fMvacuum(0.1396), fRatio(ratio) 
{}

double o4_dispersion::eofp(const double &p)  { 
  return sqrt( v2(p) * p * p + m2pole(p) ) ;
}

double o4_dispersion::fvac(const double &p) {
  double x = sqrt(fMvacuum*fMvacuum + p*p)/fT ;
  return exp(-x)/(1. - exp(-x)) ;
}

double o4_dispersion::f(const double &p) {
  double x = eofp(p)/fT ;
  return exp(-x)/(1. - exp(-x)) ;
}

double o4_dispersion::m2pole(const double &p)  { 
  double mv2 = fMvacuum * fMvacuum ;
  double m02 = mv2 * fRatio ;
  double x = p*p/(fLambda*fLambda) ;
  return mv2 - (mv2 - m02)/(1 + x/2. + x *x) ;
}

double o4_dispersion::v2(const double &p)  { 
  double v02 = fRatio*fRatio ;
  double x = p*p/(fLambda*fLambda) ;
  return 1. - (1.-v02)/(1 + x/2. + x*x) ;
}

double o4_dispersion::fG(const double &z) 

void o4_dispersion::plot(const std::string &filename) 
{
   int    np   = 100;
   double pmin = 0.;
   double pmax = 1.5;
   double dp   = (pmax - pmin)/ (double) np ;
   int  ip ;
   FILE *fp = fopen(filename.c_str(), "w") ;
   for (ip = 0 ; ip < np ; ip++) {
     double p = pmin + ip*dp;
     fprintf(fp,"%15.5e ", p) ;
     fprintf(fp,"%15.5e ", v2(p)) ;
     fprintf(fp,"%15.5e ", m2pole(p)) ;
     fprintf(fp,"%15.5e ", eofp(p)) ;
     fprintf(fp,"%15.5e ", f(p)) ;
     fprintf(fp,"%15.5e ", fvac(p)) ;
     fprintf(fp,"\n") ;
   }
   fclose(fp) ;
}
