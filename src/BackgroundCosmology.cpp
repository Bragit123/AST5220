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

  OmegaNu = Neff * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0) * OmegaR;

  OmegaLambda = 1.0 - (OmegaK + OmegaB + OmegaCDM + OmegaR + OmegaNu);
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  // Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  double npts = 100;
  std::cout << x_start << " : " << x_end << std::endl;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    double Hp = Hp_of_x(x);

    detadx[0] = Constants.c/Hp;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  Vector eta_ini{Constants.c/Hp_of_x(x_start)};
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ini);

  auto eta_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array, "eta");
  // Utils::EndTiming("Eta");


  // Cosmic time t
  // Utils::StartTiming("t");

  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    double H = H_of_x(x);

    dtdx[0] = 1/H;

    return GSL_SUCCESS;
  };

  Vector t_ini{1/(2*H_of_x(x_start))};
  ode.solve(dtdx, x_array, t_ini);

  auto t_array = ode.get_data_by_component(0);

  t_of_x_spline.create(x_array, t_array, "t");

  // Utils::EndTiming("t");
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

double BackgroundCosmology::get_curvature_scale_factor_of_chi(double chi) const{
  double r;
  double u = sqrt(abs(OmegaK)) * H0 * chi / Constants.c; // sin/sinh - argument.
  double tol = 1e-10;

  if (OmegaK < -tol) {
    r = chi * sin(u) / u;
  } else if (OmegaK > tol) {
    r = chi * sinh(u) / u;
  } else {
    r = chi;
  }

  return r;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  
  double chi = get_comoving_distance_of_x(x);
  
  double r = get_curvature_scale_factor_of_chi(chi);

  double dL = r * exp(-x);

  return dL;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  double eta0 = eta_of_x(0); // eta0 = eta today, which is when a=1, or x=log(a)=0
  double eta = eta_of_x(x);

  double chi = eta0 - eta;

  return chi;
}

double BackgroundCosmology::get_angular_distance_of_x(double x) const{
  double chi = get_comoving_distance_of_x(x);
  double r = get_curvature_scale_factor_of_chi(chi);
  return exp(x) * r;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
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

double BackgroundCosmology::get_z(double x) const{
  return exp(-x)-1;
}

double BackgroundCosmology::get_OmegaM(double x) const{
  return get_OmegaB(x) + get_OmegaCDM(x);
}

double BackgroundCosmology::get_OmegaRtot(double x) const{
  return get_OmegaR(x) + get_OmegaNu(x);
}


//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{
  double Gyr = 1e9 * 365 * 24 * 60 * 60;
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
  std::cout << "OmegaTotal   " << OmegaB+OmegaCDM+OmegaLambda+OmegaK+OmegaNu+OmegaR << "\n";


  const double x_min = x_start;
  const double x_max = x_end;
  const int    n_pts = 100;
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);
  
  Vector OmegaMR_array(n_pts, 0.0);
  Vector OmegaMDE_array(n_pts, 0.0);

  for (int i=0; i < n_pts; i++) {
    double xi = x_array.at(i);
    
    // Find Matter-radiation spline
    double OmegaMR_i = get_OmegaM(xi) - get_OmegaRtot(xi);
    if (OmegaMR_i < 0.0) {
      OmegaMR_i = -OmegaMR_i; // Want difference in absolute value
    }
    OmegaMR_array.at(i) = OmegaMR_i;

    // Find Matter-dark energy spline
    double OmegaMDE_i = get_OmegaM(xi) - get_OmegaLambda(xi);
    if (OmegaMDE_i < 0.0) {
      OmegaMDE_i = -OmegaMDE_i; // Want difference in absolute value
    }
    OmegaMDE_array.at(i) = OmegaMDE_i;
  }

  /*
  std::cout << "1" << std::endl;
  Spline OmegaMR(x_array, OmegaMR_array);
  std::cout << "2" << std::endl;
  std::cout << OmegaMR(exp(3.0*pow(10,-4))) << std::endl;
  std::cout << "3" << std::endl;
  double x_MR = Utils::binary_search_for_value(OmegaMR, 0.03);
  std::cout << "4" << std::endl;

  std::cout << "Matter-Radiation equality" << "\n";
  std::cout << "  x = " << x_MR << "\n";

  std::cout << "Universe begins to accelerate" << "\n";
  
  std::cout << "Matter-Dark energy equality" << "\n";
  */

  std::cout << "Today" << "\n";
  std::cout << "  a:       " << 1 << "\n";
  std::cout << "  z:       " << get_z(0) << "\n";
  std::cout << "  t (Gyr): " << t_of_x(0) / Gyr << "\n";
  std::cout << "  eta (Gyr):     " << eta_of_x(0) / Gyr << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = x_start;
  const double x_max = x_end;
  const int    n_pts = 100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);
  std::cout << x_min << std::endl;
  std::cout << x_max << std::endl;

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                                 << " ";
    fp << get_z(x)                          << " ";
    fp << eta_of_x(x)                       << " ";
    fp << t_of_x(x)                         << " ";
    fp << H_of_x(x)/get_H0()                << " ";
    fp << Hp_of_x(x)                        << " ";
    fp << dHpdx_of_x(x)                     << " ";
    fp << ddHpddx_of_x(x)                   << " ";
    fp << get_luminosity_distance_of_x(x)   << " ";
    fp << get_comoving_distance_of_x(x)     << " ";
    fp << get_angular_distance_of_x(x)      << " ";
    fp << get_OmegaB(x)                     << " ";
    fp << get_OmegaCDM(x)                   << " ";
    fp << get_OmegaLambda(x)                << " ";
    fp << get_OmegaR(x)                     << " ";
    fp << get_OmegaNu(x)                    << " ";
    fp << get_OmegaK(x); // Removed [<< " "] because it created an extra column of NaNs in pandas
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

