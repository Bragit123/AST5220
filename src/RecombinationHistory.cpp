#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{

}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");

  // Some values needed when computing n_e
  const double m_H = Constants.m_H;
  const double G = Constants.G;
  const double OmegaB0 = cosmo->get_OmegaB();
  const double H0 = cosmo->get_H0();
  const double rho_c0 = 3 * H0*H0 / (8 * M_PI * G);
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  const double x_min = x_start;
  const double x_max = x_end;
  
  Vector x_array = Utils::linspace(x_min, x_max, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);
  
  Vector Xe_saha_arr(npts_rec_arrays);
  Vector ne_saha_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){
    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================

    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    Xe_saha_arr[i] = Xe_current;
    ne_saha_arr[i] = ne_current;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      Vector Xe_ini{Xe_current};

      // Make an array containing only the remaining x-values
      double xi = x_array[i];
      Vector x_rest = Utils::linspace(xi, x_max, npts_rec_arrays - i);

      peebles_Xe_ode.solve(dXedx, x_rest, Xe_ini);
      auto Xe_rest = peebles_Xe_ode.get_data_by_component(0);

      for (int j=0; j < npts_rec_arrays - i; j++) {

        Xe_arr[i + j] = Xe_rest[j];

        // Compute n_e from x
        double a = exp(x_array[i+j]);
        double nb = OmegaB0 * rho_c0 / (m_H * pow(a, 3));
        ne_arr[i + j] = Xe_rest[j] * nb; // Ignore heavier elements than hydrogen, so nH = nb

        //// Also update saha arrays for each iteration
        auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i+j]);
        const double Xe_current = Xe_ne_data.first;
        const double ne_current = Xe_ne_data.second;
        Xe_saha_arr[i+j] = Xe_current;
        ne_saha_arr[i+j] = ne_current;
      }

      break; // Break the loop since we create all results from Peebles ODE in one go

    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  
  Vector log_Xe_arr = log(Xe_arr);
  Vector log_ne_arr = log(ne_arr);

  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "Xe");
  log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");
  
  Xe_saha_of_x_spline.create(x_array, Xe_saha_arr, "Xe_saha");
  ne_saha_of_x_spline.create(x_array, ne_saha_arr, "ne_saha");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double H0 = cosmo->get_H0();
  const double OmegaB0 = cosmo->get_OmegaB();
  const double TCMB0 = cosmo->get_TCMB();
  
  // Compute parameters
  const double rho_c0 = 3 * H0*H0 / (8 * M_PI * G);
  double nb = OmegaB0 * rho_c0 / (m_H * pow(a, 3));
  double Tb = TCMB0 / a;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================

  double coeff = 1.0/nb * pow(k_b*m_e*Tb / (2*M_PI*hbar*hbar), 3.0/2.0) * exp(-epsilon_0 / (k_b * Tb));
  
  double Xe = 0.0;
  double tol = 1e-4; // tolerance for taylor-approximation
  if (1.0/coeff >= tol) {
    // Quadratic formula with a=1, b=coeff and c=-coeff. Ignore negative solutions
    Xe = (-coeff + sqrt(coeff*coeff + 4*coeff)) / 2.0;
  } else {
    Xe = 1.0;
  }

  double ne = Xe * nb; // Ignore heavier elements than hydrogen, so nH = nb

  // Return electron fraction and number density
  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double a           = exp(x);
  const double X_e         = Xe[0];

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double H0 = cosmo->get_H0();
  double H = cosmo->H_of_x(x);
  const double TCMB0 = cosmo->get_TCMB();
  const double OmegaB0 = cosmo->get_OmegaB();

  double Tb = TCMB0 / a;

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  double nH = 3.0 * H0*H0 * OmegaB0 / (8.0 * M_PI * G * m_H * pow(a, 3.0));
  // double nH = 3.0 * H0*H0 * OmegaB0 / 8.0 / M_PI / G / m_H / pow(a, 3.0);
  double n_1s = (1.0 - X_e) * nH;

  double epsilon0_kb_Tb = epsilon_0 / (k_b * Tb);
  double phi2 = 0.448 * log(epsilon0_kb_Tb);
  double alpha2 = 8.0/sqrt(3.0*M_PI) * c * sigma_T * sqrt(epsilon0_kb_Tb) * phi2;
  
  double beta_coeff = alpha2 * pow(k_b*m_e*Tb / (2.0 * M_PI * hbar*hbar), 3.0/2.0);
  double beta = beta_coeff * exp(-epsilon0_kb_Tb);
  
  // Write beta2 explicity instead of in terms of beta to avoid overflow in exponent.
  double beta2 = beta_coeff * exp(-epsilon0_kb_Tb / 4.0);
  
  double lambda_alpha = H * pow(3.0*epsilon_0, 3.0) / (pow(8.0*M_PI, 2.0) * pow(c, 3) * pow(hbar, 3) * n_1s);
  double lambda_2s1s_alpha = lambda_2s1s + lambda_alpha;
  double Cr = (lambda_2s1s_alpha) / (lambda_2s1s_alpha + beta2);
  
  double fac = Cr/H;
  double first_term = beta*(1.0-X_e);
  double second_term = nH*alpha2*X_e*X_e;
  
  double rhs = Cr / H * (beta * (1.0 - X_e) - nH * alpha2 * X_e*X_e);
  
  dXedx[0] = rhs;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array_rev = Utils::linspace(x_end, x_start, npts); // Reversed to account for initial condition

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    const double c = Constants.c;
    const double sigma_T = Constants.sigma_T;
    double H = cosmo->H_of_x(x);
    double n_e = ne_of_x(x);

    // Set the derivative for photon optical depth
    double rhs = - c * n_e * sigma_T / H;
    
    dtaudx[0] = rhs;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================

  Vector tau_ini{0.0}; // Use some initialized value, then subtract value at x=0 to get correct.
  ODESolver tau_ode;
  tau_ode.solve(dtaudx, x_array_rev, tau_ini);

  auto tau_array_rev = tau_ode.get_data_by_component(0); // temporary tau array until reversing tau(x=0)
  Vector tau_array(npts, 0.0);
  for (int i=0; i < npts; i++) {
    tau_array[i] = tau_array_rev[npts-i];
  }
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  tau_of_x_spline.create(x_array, tau_array, "tau"); // Spline result


  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  Vector g_tilde_array(npts, 0.0);
  for (int i = 0; i < npts; i++) {
    double xi = x_array[i];
    double tau = tau_of_x(xi);
    double dtaudx = dtaudx_of_x(xi);

    g_tilde_array[i] = -dtaudx * exp(-tau);
  }

  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g");

  Utils::EndTiming("opticaldepth");

  //// Now solve for the sound-horizon
  const double c = Constants.c;
  const double OmegaR0 = cosmo->get_OmegaR();
  const double OmegaB0 = cosmo->get_OmegaB();
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    double Hp = cosmo->Hp_of_x(x);
    double a = exp(x);

    double R = 4 * OmegaR0 / (3 * OmegaB0 * a);
    double cs = c * sqrt(R / (3 * (1+R)));

    // Set the derivative for photon optical depth
    double rhs = cs / Hp;
    
    dsdx[0] = rhs;

    return GSL_SUCCESS;
  };

  double Hp_ini = cosmo->Hp_of_x(x_start);
  double a_ini = exp(x_start);
  double R_ini = 4 * OmegaR0 / (3 * OmegaB0 * a_ini);
  double cs_ini = c * sqrt(R_ini / (3.0 * (1+R_ini)));

  ODESolver s_ode;
  Vector s_ini{cs_ini/Hp_ini};
  s_ode.solve(dsdx, x_array, s_ini);

  auto s_array = s_ode.get_data_by_component(0);

  s_of_x_spline.create(x_array, s_array, "s");

}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::s_of_x(double x) const{
  return s_of_x_spline(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  double log_Xe = log_Xe_of_x_spline(x);
  double Xe = exp(log_Xe);

  return Xe;
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  double log_ne = log_ne_of_x_spline(x);
  double ne = exp(log_ne);

  return ne;
}

double RecombinationHistory::Xe_saha_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  double Xe_saha = Xe_saha_of_x_spline(x);

  return Xe_saha;
}

double RecombinationHistory::ne_saha_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  double ne_saha = ne_saha_of_x_spline(x);

  return ne_saha;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  // Compute decoupling times
  double x_decoupling = Utils::binary_search_for_value(tau_of_x_spline, 1);
  double z_decoupling = cosmo->get_z(x_decoupling);
  double t_decoupling = cosmo->t_of_x(x_decoupling);
  double rs = s_of_x(x_decoupling);

  double x_recombination = Utils::binary_search_for_value(log_Xe_of_x_spline, log(0.1));
  double z_recombination = cosmo->get_z(x_recombination);
  double t_recombination = cosmo->t_of_x(x_recombination);

  double Gyr = 1e9 * 365 * 24 * 60 * 60;
  double Mpc = Constants.Mpc;

  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << "rs (Mpc):    " << rs / Mpc    << "\n";
  std::cout << "Xe(x=0)      " << Xe_of_x(0)  << "\n";
  std::cout << "Xe_saha(x=0)      " << Xe_saha_of_x(0)  << "\n";
  std::cout << "Times at decoupling:\n";
  std::cout << "  a:         " << exp(x_decoupling) << "\n";
  std::cout << "  z:         " << z_decoupling << "\n";
  std::cout << "  t (Gyr):   " << t_decoupling / Gyr << "\n";
  std::cout << "Times at recombination:\n";
  std::cout << "  a:         " << exp(x_recombination) << "\n";
  std::cout << "  z:         " << z_recombination << "\n";
  std::cout << "  t (Gyr):   " << t_recombination / Gyr << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << Xe_saha_of_x(x)      << " ";
    fp << ne_saha_of_x(x)      << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x); // Removed [<< " "] because it created an extra column of NaNs in pandas
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

