#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{
  H0 = Constants.H0_over_h * h;

  double OmegaR_numerator = 2 * pow(M_PI, 2) * pow(Constants.k_b * TCMB, 4) * 8 * M_PI * Constants.G;
  double OmegaR_denominator = 30 * pow(Constants.hbar, 3) * pow(Constants.c, 5) * 3 * pow(H0, 2);
  OmegaR = OmegaR_numerator / OmegaR_denominator;

  OmegaNu = Neff * 7/8 * pow(4.0/11.0, 4.0/3.0) * OmegaR;

  OmegaLambda = 1.0 - (OmegaK + OmegaB + OmegaCDM + OmegaR + OmegaNu);

  double x_val = 1e-3;
  std::cout << "OmegaB = " << get_OmegaB(x_val) << std::endl;
  std::cout << "OmegaR = " << get_OmegaR(x_val) << std::endl;
  std::cout << "OmegaNu = " << get_OmegaNu(x_val) << std::endl;
  std::cout << "OmegaCDM = " << get_OmegaCDM(x_val) << std::endl;
  std::cout << "OmegaLambda = " << get_OmegaLambda(x_val) << std::endl;
  std::cout << "OmegaK = " << get_OmegaK(x_val) << std::endl;
  std::cout << "TCMB = " << get_TCMB(x_val) << std::endl;
  std::cout << "r = " << get_comoving_distance_of_x(x_val) << std::endl;
  std::cout << "dL = " << get_luminosity_distance_of_x(x_val) << std::endl;
  //...
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  Vector x_array;

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...

    detadx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  // ...
  // ...
  // ...
  // ...

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

// H(x)
double BackgroundCosmology::H_of_x(double x) const{
  double term_b_CDM = (OmegaB + OmegaCDM) * exp(-3*x);
  double term_R_Nu = (OmegaR + OmegaNu) * exp(-4*x);
  double term_K = OmegaK * exp(-2*x);
  double H_sqrt = sqrt(term_b_CDM + term_R_Nu + term_K + OmegaLambda);

  return H0 * H_sqrt;
}

// Hp(x)
double BackgroundCosmology::Hp_of_x(double x) const{
  double H = H_of_x(x);

  return exp(x) * H;
}

// dHpdx(x)
double BackgroundCosmology::dHpdx_of_x(double x) const{
  double Hp = Hp_of_x(x);

  double term_b_CDM = -0.5 * (OmegaB + OmegaCDM) * exp(-x);
  double term_R_Nu = -(OmegaR + OmegaNu) * exp(-2*x);
  double term_Lambda = OmegaLambda * exp(2*x);
  double dHp_Bracket = term_b_CDM + term_R_Nu + term_Lambda;

  double dHpdx = pow(H0, 2) / Hp * dHp_Bracket;

  return dHpdx;
}

// ddHpddx(x)
double BackgroundCosmology::ddHpddx_of_x(double x) const{
  double Hp = Hp_of_x(x);
  double dHpdx = dHpdx_of_x(x);

  double term_b_CDM = 0.5 * (OmegaB + OmegaCDM) * exp(-x);
  double term_R_Nu = 2 * (OmegaR + OmegaNu) * exp(-2*x);
  double term_Lambda = 2 * OmegaLambda * exp(2*x);
  double ddHp_Bracket = term_b_CDM + term_R_Nu + term_Lambda;

  double ddHpddx = pow(H0, 2) / Hp * ddHp_Bracket - 1/Hp * pow(dHpdx, 2);

  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  double H = H_of_x(x);
  double denominator = exp(3*x) * pow(H, 2) / pow(H0, 2);

  return OmegaB / denominator;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  double H = H_of_x(x);
  double denominator = exp(4*x) * pow(H, 2) / pow(H0, 2);

  return OmegaR / denominator;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  double H = H_of_x(x);
  double denominator = exp(4*x) * pow(H, 2) / pow(H0, 2);

  return OmegaNu / denominator;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  double H = H_of_x(x);
  double denominator = exp(3*x) * pow(H, 2) / pow(H0, 2);

  return OmegaCDM / denominator;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  double H = H_of_x(x);
  double denominator = pow(H, 2) / pow(H0, 2);

  return OmegaLambda / denominator;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  double H = H_of_x(x);
  double denominator = exp(2*x) * pow(H, 2) / pow(H0, 2);

  return OmegaK / denominator;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  
  double r = get_comoving_distance_of_x(x);
  double dL = r * exp(-x);

  return dL;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  double eta = 0.0;
  double eta0 = 1.0;

  double chi = eta0 - eta;
  double u = sqrt(abs(OmegaK)) * H0 * chi / Constants.c; // sin/sinh - argument.
  double tol = 1e-10;
  double r;

  if (OmegaK < -tol) {
    r = chi * sin(u) / u;
  } else if (OmegaK > tol) {
    r = chi * sinh(u) / u;
  } else {
    r = chi;
  }

  return r;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

