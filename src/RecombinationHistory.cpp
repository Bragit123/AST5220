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
        auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
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
  
  Vector log_Xe_saha_arr = log(Xe_saha_arr);
  Vector log_ne_saha_arr = log(ne_saha_arr);
  
  log_Xe_saha_of_x_spline.create(x_array, log_Xe_saha_arr, "Xe");
  log_ne_saha_of_x_spline.create(x_array, log_ne_saha_arr, "ne");

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
  const double H0_over_h   = Constants.H0_over_h;

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

  double coeff = 1/nb * pow(k_b*m_e*Tb / (2*M_PI*hbar*hbar), 3/2) * exp(-epsilon_0 / (k_b * Tb));
  
  double Xe = 0.0;
  double tol = 1e-4; // tolerance for taylor-approximation
  // double tol = 0.3; // tolerance for taylor-approximation
  if (coeff >= tol) {
    // Quadratic formula with a=1, b=coeff and c=-coeff. Ignore negative solutions
    Xe = (-coeff + sqrt(coeff*coeff + 4*coeff)) / 2;
  } else {
    // Approximation for sqrt() if coeff << 1
    Xe = (-coeff + 2 * sqrt(coeff) * (1 + coeff/8)) / 2;
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
  const double X_e         = Xe[0];
  const double a           = exp(x);

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
  double nH = 3 * H0*H0 * OmegaB0 / (8 * M_PI * G * m_H * pow(a, 3));
  double n_1s = (1.0 - X_e) * nH;
  
  double alpha = 1 / 137.0359992; // Fine structure constant
  double phi2 = 0.448 * log(epsilon_0 / (k_b * Tb));
  // double alpha2 = hbar*hbar / (c * sqrt(k_b)) * 64*M_PI / sqrt(27*M_PI) * alpha*alpha / (m_e*m_e) * sqrt(epsilon_0/Tb) * phi2;
  double alpha2 = 8/sqrt(3*M_PI) * c * sigma_T * sqrt(epsilon_0 / (k_b * Tb)) * phi2;

  // double beta_coeff = pow(k_b, 3/2) / pow(hbar, 3) * alpha2 * pow(m_e*Tb / (2 * M_PI), 3/2);
  double beta_coeff = alpha2 * pow(k_b*m_e*Tb / (2 * M_PI * hbar*hbar), 3/2);
  double beta = beta_coeff * exp(-epsilon_0/(k_b*Tb));
  // Write beta2 explicity instead of in terms of beta to avoid overflow in exponent.
  double beta2 = beta_coeff * exp(-epsilon_0/(4*k_b*Tb));

  double lambda_alpha = 1/pow(hbar*c, 3) * H * pow(3*epsilon_0, 3) / (pow(8*M_PI, 2) * n_1s);
  double lambda_2s1s_alpha = lambda_2s1s + lambda_alpha;
  double Cr = (lambda_2s1s_alpha) / (lambda_2s1s_alpha + beta2);
  
  double rhs = Cr / H * (beta * (1.0 - X_e) - nH * alpha2 * X_e*X_e);
  dXedx[0] = rhs;
  // dXedx[0] = 0.0;

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
  Vector x_array = Utils::linspace(x_start, x_end, npts);

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

  Vector tau_ini{0}; // Use some initialized value, then subtract value at x=0 to get correct.
  ODESolver tau_ode;
  tau_ode.solve(dtaudx, x_array, tau_ini);

  auto tau_temp_array = tau_ode.get_data_by_component(0); // temporary tau array until subtracting tau(x=0)

  // Find tau_temp(x=0), and subtract that value from tau_temp to ensure tau(x=0)=0
  Spline tau_temp_spline(x_array, tau_temp_array);
  double tau_zero = tau_temp_spline(0);
  Vector tau_array(npts, 0.0); // Initialize tau array

  for (int i=0; i < npts; i++) {
    // Subtract tau_zero from each element in tau_temp_array, and store in tau_array
    tau_array[i] = tau_temp_array[i] - tau_zero;
  }

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
  double log_Xe_saha = log_Xe_saha_of_x_spline(x);
  double Xe_saha = exp(log_Xe_saha);

  return Xe_saha;
}

double RecombinationHistory::ne_saha_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  double log_ne_saha = log_ne_saha_of_x_spline(x);
  double ne_saha = exp(log_ne_saha);

  return ne_saha;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
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

