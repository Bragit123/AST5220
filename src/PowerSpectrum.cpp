#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  // Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  Vector k_array = Utils::linspace(k_min, k_max, n_k);
  Vector log_k_array = log(k_array);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  double z_min = 0.0;
  double z_max = k_max * cosmo->eta_of_x(0.0);
  double dz = 2.0*M_PI / 16.0;
  // std::cout << "z_min = " << z_min << std::endl;
  // std::cout << "z_max = " << z_max << std::endl;
  // std::cout << "dz = " << dz << std::endl;
  int n_z = (z_max - z_min) / dz;
  // std::cout << "n_z = " << n_z << std::endl;
  // std::cout << "----- 1 -----" << std::endl;
  Vector z_array = Utils::linspace(z_min, z_max, n_z);
  // std::cout << "----- 2 -----" << std::endl;

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];

    Vector j_ell_i(n_z, 0.0);

    for (int iz=0; iz < n_z; iz++) {
      double z = z_array[iz];
      j_ell_i[iz] = Utils::j_ell(ell, z);
    }

    Spline j_ell_spline_i(z_array, j_ell_i);
    j_ell_splines[i] = j_ell_spline_i;
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  int n_x = 1000;
  Vector x_array = Utils::linspace(Constants.x_start, 0.0, n_x);

  double eta0 = cosmo->eta_of_x(0.0);

  for (int i=0; i < ells.size(); i++) {
    for (size_t ik=0; ik < k_array.size(); ik++) {
      double k = k_array[ik];

      ODEFunction dThetaldx = [&](double x, const double *Thetal, double *dThetaldx){
        double S = source_function(x, k);
        double eta = cosmo->eta_of_x(x);
        double z = k * (eta0 - eta);
        double j_ell = j_ell_splines[i](z);
        
        dThetaldx[0] = S * j_ell;

        return GSL_SUCCESS;
      };
      
      Vector Thetal_ini{0.0};
      ODESolver ode;
      ode.solve(dThetaldx, x_array, Thetal_ini);

      auto Thetal_array = ode.get_data_by_component(0);
      
      result[i][ik] = Thetal_array[n_x - 1];
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

// Vector2D PowerSpectrum::line_of_sight_integration_single(
//     Vector & k_array, 
//     std::function<double(double,double)> &source_function){
//   Utils::StartTiming("lineofsight");

//   // Make storage for the results
//   Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

//   double eta0 = cosmo->eta_of_x(0.0);

//   int n_x = 1000;
//   Vector x_array = Utils::linspace(Constants.x_start, 0.0, n_x);

//   for(size_t ik = 0; ik < k_array.size(); ik++) {

//     for (int ell=0; ell < ells.size(); ell++) {
//       //=============================================================================
//       // TODO: Implement to solve for the general line of sight integral 
//       // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
//       // given value of k
//       //=============================================================================
//       double k = k_array[ik];

//       double x = x_array[0];
//       double S_xk = source_function(x, k);
//       double eta = cosmo->eta_of_x(x);
//       double z = k * (eta0 - eta);
//       double j_ell = j_ell_splines[ell](z);

//       double integrand = S_xk * j_ell;
//       double integrand_prev;

//       double integral = 0.0;
//       double integral_add;

//       for (int ix=1; ix < n_x; ix++) {
//         integrand_prev = integrand;
//         x = x_array[ix];
//         S_xk = source_function(x, k);

//         eta = cosmo->eta_of_x(x);
//         z = k * (eta0 - eta);
//         j_ell = j_ell_splines[ell](z);

//         integrand = S_xk * j_ell;

//         integral_add = (integrand + integrand_prev) / 2.0; // Trapezoidal method
//         integral = integral + integral_add;
//       }
      
//       // Store the result for Source_ell(k) in results[ell][ik]
//       result[ell][ik] = integral;
//     }
//   }

//   Utils::EndTiming("lineofsight");
//   return result;
// }

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);



  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  for (int i=0; i < nells; i++) {
    thetaT_ell_of_k_spline[i].create(k_array, thetaT_ell_of_k[i]);
  }

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  // if(Constants.polarization){

  // }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  
  Vector result = Vector(nells);
  
  double dlogk0 = log_k_array[1] - log_k_array[0];
  double dlogk1 = log_k_array[log_k_array.size()-1] - log_k_array[log_k_array.size()-2];
  Vector k_array = exp(log_k_array);
  double dk0 = k_array[1] - k_array[0];
  double dk1 = k_array[k_array.size()-1] - k_array[k_array.size()-2];
  double eta0 = cosmo->eta_of_x(0.0);
  double dk_wanted = 2*M_PI / (32 * eta0);
  double n0 = 2*M_PI / (dk0 * eta0);
  double n1 = 2*M_PI / (dk1 * eta0);
  double n_tot_need = (k_array[k_array.size()-1] - k_array[0]) / dk_wanted;

  for (int i=0; i < nells; i++) {
    ODEFunction dCldlogk = [&](double logk, const double *Cl, double *dCldlogk){
      double k = exp(logk);
      double P = primordial_power_spectrum(k);
      // double k_mpc = k * Constants.Mpc;
      // std::cout << "k = " << k << std::endl;
      // std::cout << "k_mpc = " << k_mpc << std::endl;
      // double P = get_matter_power_spectrum(0.0, k_mpc);
      // double P_mpc3 = get_matter_power_spectrum(0.0, k_mpc);
      // double P = P_mpc3 / (pow(Constants.Mpc, 3));
      double Theta_squared = f_ell_spline[i](k) * g_ell_spline[i](k);

      dCldlogk[0] = 4*M_PI * P * Theta_squared;

      return GSL_SUCCESS;
    };
    
    double k_ini = exp(log_k_array[0]);
    // Vector Cl_ini{primordial_power_spectrum(k_ini)};
    // Vector Cl_ini{1000.0};
    Vector Cl_ini{0.0};
    ODESolver ode;
    ode.solve(dCldlogk, log_k_array, Cl_ini);

    auto Cl_array = ode.get_data_by_component(0);
    
    int last_k_index = log_k_array.size() - 1;
    result[i] = Cl_array[last_k_index];
  }

  return result;
}

// Vector PowerSpectrum::solve_for_cell(
//     Vector & log_k_array,
//     std::vector<Spline> & f_ell_spline,
//     std::vector<Spline> & g_ell_spline){
//   const int nells      = ells.size();

//   //============================================================================
//   // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
//   // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
//   //============================================================================
  
//   Vector result = Vector(nells);

//   Vector k_array = exp(log_k_array);

//   for (int ell=0; ell < nells; ell++) {

//     double k = k_array[0];
//     double k_mpc = k * Constants.Mpc;
//     // double P = primordial_power_spectrum(k);
//     double P_mpc3 = get_matter_power_spectrum(0.0, k);
//     double P = P_mpc3 / pow(Constants.Mpc, 3);
//     double integrand = P * f_ell_spline[ell](k) * g_ell_spline[ell](k) / k;
//     double integrand_prev;
//     double integral_add;
//     double integral = 0.0;

//     for (int ik=0; ik < log_k_array.size(); ik++) {
//       integrand_prev = integrand;
//       k = k_array[ik];
//       k_mpc = k * Constants.Mpc;
//       // P = primordial_power_spectrum(k);
//       P_mpc3 = get_matter_power_spectrum(0.0, k_mpc);
//       P = P_mpc3 / pow(Constants.Mpc, 3);
//       integrand = P * f_ell_spline[ell](k) * g_ell_spline[ell](k) / k;

//       integral_add = (integrand_prev + integrand) / 2.0; // Trapezoidal method
//       integral = integral + integral_add;
//     }

//     result[ell] = 4*M_PI * integral;
//   }

//   return result;
// }

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  double k = k_mpc / Constants.Mpc;
  double c = Constants.c;
  double Phi = pert->get_Phi(x, k);
  double OmegaM = cosmo->get_OmegaM();
  double H0 = cosmo->get_H0();

  double num = c*c * k*k * Phi;
  double den = 3.0/2.0 * OmegaM * exp(-x) * H0*H0;
  double delta_M = num / den;

  // double pofk = abs(delta_M*delta_M) * primordial_power_spectrum(k) * 2*M_PI*M_PI / pow(k, 3.0);
  double pofk = abs(delta_M)*abs(delta_M) * primordial_power_spectrum(k) * 2*M_PI*M_PI / pow(k, 3.0);

  double pofk_mpc3 = pofk / pow(Constants.Mpc, 3.0);

  return pofk_mpc3;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    fp << cell_TT_spline( ell ) * normfactorN  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

void PowerSpectrum::output_theta(std::string filename) const{
  const double H0 = cosmo->get_H0();
  const double c = Constants.c;

  std::ofstream fp(filename.c_str());
  auto kvalues = Utils::linspace(k_min, k_max, n_k);
  auto print_data = [&] (const double k) {
    fp << c*k / H0                                 << " ";
    fp << thetaT_ell_of_k_spline[4](k)   << " ";
    fp << thetaT_ell_of_k_spline[19](k)  << " ";
    fp << thetaT_ell_of_k_spline[24](k)  << " ";
    fp << thetaT_ell_of_k_spline[32](k)  << " ";
    fp << thetaT_ell_of_k_spline[42](k)  << " ";
    fp << "\n";
  };
  std::for_each(kvalues.begin(), kvalues.end(), print_data);
}

void PowerSpectrum::output_matter(std::string filename) const{
  const double h = cosmo->get_h();
  const double Mpc = Constants.Mpc;

  std::ofstream fp(filename.c_str());
  auto kvalues = Utils::linspace(k_min, k_max, n_k);
  auto print_data = [&] (const double k) {
    double k_mpc = k * Mpc;
    fp << k_mpc/h                                 << " ";
    fp << get_matter_power_spectrum(0.0, k_mpc) * pow(h, 3.0)     << " ";
    fp << "\n";
  };
  std::for_each(kvalues.begin(), kvalues.end(), print_data);
}