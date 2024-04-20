#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // // Compute source functions and spline the result
  std::cout << "#### HUSK UNCOMMENT SOURCE FUNCTION ####" << std::endl;
  // compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================

  // Initialize vectors needed for splines
  Vector2D delta_cdm(n_x, Vector(n_k));
  Vector2D delta_b(n_x, Vector(n_k));
  Vector2D v_cdm(n_x, Vector(n_k));
  Vector2D v_b(n_x, Vector(n_k));
  Vector2D Phi(n_x, Vector(n_k));
  Vector2D Psi(n_x, Vector(n_k));
  std::vector<Vector2D> Theta(Constants.n_ell_theta, Vector2D(n_x, Vector(n_k)));

  // Find x_array
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make k_array
  double k_log_start = log(k_min);
  double k_log_end = log(k_max);
  Vector k_log_array = Utils::linspace(k_log_start, k_log_end, n_k);
  Vector k_array = exp(k_log_array);

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);

    // Find number of x-values for tight coupling and the full system
    double tight_fraction = (x_end_tight - x_start) / (x_end - x_start);
    double n_x_tc = floor(tight_fraction * n_x);
    double n_x_full = n_x - n_x_tc + 1;

    Vector x_array_tc(n_x_tc);
    Vector x_array_full(n_x_full);
    for (int i=0; i < n_x_tc; i++) {
      x_array_tc[i] = x_array[i];
    }
    for (int i=0; i < n_x_full; i++) {
      x_array_full[i] = x_array[n_x_tc + i - 1];
    }

    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    ODESolver ode_tc;
    ode_tc.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini);

    auto y_tight_coupling = ode_tc.get_data();

    //===================================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of
    // tight coupling)
    const int n_ell_tot_tc = Constants.n_ell_tot_tc;
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling.back(), x_end_tight, k);

    if (k*Constants.Mpc > 0.0093 && k*Constants.Mpc < 0.017) {
      std::cout << "k = " << k*Constants.Mpc << "  :  " << y_full_ini[Constants.ind_start_theta] << std::endl;
    }

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    ODESolver ode_full;
    ode_full.solve(dydx_full, x_array_full, y_full_ini);

    auto y_full = ode_full.get_data();

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================
    
    // Some values needed to compute Psi
    double H0 = cosmo->get_H0();
    double OmegaR0 = cosmo->get_OmegaR();
    double c = Constants.c;

    for (int ix = 0; ix < n_x_tc; ix++) {
      delta_cdm[ix][ik] = y_tight_coupling[ix][Constants.ind_deltacdm];
      delta_b[ix][ik] = y_tight_coupling[ix][Constants.ind_deltab];
      v_cdm[ix][ik] = y_tight_coupling[ix][Constants.ind_vcdm];
      v_b[ix][ik] = y_tight_coupling[ix][Constants.ind_vb];
      Phi[ix][ik] = y_tight_coupling[ix][Constants.ind_Phi];

      for (int l=0; l < Constants.n_ell_theta_tc; l++) {
        Theta[l][ix][ik] = y_tight_coupling[ix][Constants.ind_start_theta + l];
      }

      double Hp = cosmo->Hp_of_x(x_array[ix]);
      double dtaudx = rec->dtaudx_of_x(x_array[ix]);
      for (int l=Constants.n_ell_theta_tc; l < Constants.n_ell_theta; l++) {
        Theta[l][ix][ik] = -l / (2.0*l + 1.0) * c*k / (Hp*dtaudx) * Theta[l-1][ix][ik];
      }

      double a = exp(x_array[ix]);
      Psi[ix][ik] = -Phi[ix][ik] - 12*H0*H0 / (c*c * k*k * a*a) * OmegaR0 * Theta[2][ix][ik];
    }
    
    for (int ix_full = 0; ix_full < n_x_full; ix_full++) {
      double ix = n_x_tc - 1 + ix_full;
      delta_cdm[ix][ik] = y_full[ix_full][Constants.ind_deltacdm];
      delta_b[ix][ik] = y_full[ix_full][Constants.ind_deltab];
      v_cdm[ix][ik] = y_full[ix_full][Constants.ind_vcdm];
      v_b[ix][ik] = y_full[ix_full][Constants.ind_vb];
      Phi[ix][ik] = y_full[ix_full][Constants.ind_Phi];

      for (int l=0; l < Constants.n_ell_theta; l++) {
        Theta[l][ix][ik] = y_full[ix_full][Constants.ind_start_theta + l];
      }

      double a = exp(x_array[ix]);
      Psi[ix][ik] = -Phi[ix][ik] - 12.0*H0*H0 / (c*c * k*k * a*a) * OmegaR0 * Theta[2][ix][ik];
    }
  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  delta_cdm_spline.create(x_array, k_array, delta_cdm, "delta_cdm_spline");
  delta_b_spline.create(x_array, k_array, delta_b, "delta_b_spline");
  v_cdm_spline.create(x_array, k_array, v_cdm, "v_cdm_spline");
  v_b_spline.create(x_array, k_array, v_b, "v_b_spline");
  Phi_spline.create(x_array, k_array, Phi, "Phi_spline");
  Psi_spline.create(x_array, k_array, Psi, "Psi_spline");

  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  for (int l=0; l < Constants.n_ell_theta; l++) {
    Vector2D Theta_l = Theta[l];
    Theta_spline[l].create(x_array, k_array, Theta_l);
  }
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================

  double c = Constants.c;
  double Hp = cosmo->Hp_of_x(x);
  double ck_Hp = c*k/Hp;
  
  // Neglect neutrinos -> f_nu = 0
  double Psi = -2.0/3.0;
  Phi = -Psi;

  delta_cdm = -(3.0/2.0) * Psi;
  delta_b = -(3.0/2.0) * Psi;

  v_cdm = -0.5 * ck_Hp * Psi;
  v_b = -0.5 * ck_Hp * Psi;

  Theta[0] = -0.5 * Psi;
  Theta[1] = 1.0/6.0 * ck_Hp * Psi;

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;
  Phi = Phi_tc;
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];

  double c = Constants.c;
  double Hp = cosmo->Hp_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);
  double ck_Hpdtau = c*k / (Hp * dtaudx);

  Theta[2] = -20.0/45.0 * ck_Hpdtau * Theta[1];

  for (int l=3; l < n_ell_theta; l++) {
    Theta[l] = -l / (2.0*l + 1) * ck_Hpdtau * Theta[l-1];
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  double x_recomb = rec->get_x_recomb();
  double c = Constants.c;

  // Declare variables
  double Hp = 0.0;
  double dtau = 0.0;
  bool cond1 = true;
  bool cond2 = true;
  bool cond3 = true;
  double xi = 0.0;

  for (int i=0; i < n_x; i++) {
    xi = x_array[i];
    Hp = cosmo->Hp_of_x(xi);
    dtau = rec->dtaudx_of_x(xi);
    
    // Conditions for leaving tight coupling. (These are the reversed conditions
    // from Callin, since he gives the conditions that must be met to STAY in
    // the tight coupling regime, thus < or > are switched to >= and <= respectively.)
    cond1 = abs(c*k/(Hp*dtau)) >= 0.1;
    cond2 = abs(dtau) <= 10.0;
    cond3 = xi >= x_recomb;

    if ((cond1 || cond2) || cond3) {
      // If any of the three conditions are met, leave tight coupling regime
      return xi;
    }
  }

  std::cout << "WARNING: Never left tight coupling!" << std::endl;
  return x_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================

  // Find x_array
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make k_array
  double k_log_start = log(k_min);
  double k_log_end = log(k_max);
  Vector k_log_array = Utils::linspace(k_log_start, k_log_end, n_k);
  Vector k_array = exp(k_log_array);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      
      // From BackgroundCosmology
      const double Hp = cosmo->Hp_of_x(x);
      const double dHpdx = cosmo->dHpdx_of_x(x);

      // From Recombination
      const double tau = rec->tau_of_x(x);
      const double g_tilde = rec->g_tile_of_x(x);

      // From Perturbations
      const double Theta0 = get_Theta(x, k, 0);
      const double Psi = get_Psi(x, k);
      const double dPsidx = get_dPsidx(x, k);
      const double dPhidx = get_dPhidx(x, k);
      const double Pi = get_Theta(x, k, 2); // Ignore polarization terms
      const double dPidx = get_dThetadx(x, k, 2); // Ignore polarization terms
      const double dv_bdx = get_dv_bdx(x, k);

      /////////////// FORTSETT HER!!!!!!!! ////////////
      /////////////// FORTSETT HER!!!!!!!! ////////////
      /////////////// FORTSETT HER!!!!!!!! ////////////
      /////////////// FORTSETT HER!!!!!!!! ////////////
      /////////////// FORTSETT HER!!!!!!!! ////////////
      /////////////// FORTSETT HER!!!!!!!! ////////////
      /////////////// FORTSETT HER!!!!!!!! ////////////
      /////////////// FORTSETT HER!!!!!!!! ////////////

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  double a = exp(x);
  double c = Constants.c;

  // Values from BackgroundCosmology and Recombination
  double H0 = cosmo->get_H0();
  double OmegaR0 = cosmo->get_OmegaR();
  double OmegaB0 = cosmo->get_OmegaB();
  double OmegaCDM0 = cosmo->get_OmegaCDM();

  double Hp = cosmo->Hp_of_x(x);
  double dHpdx = cosmo->dHpdx_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);
  double ddtauddx = rec->ddtauddx_of_x(x);

  double ck_Hp = c*k/Hp; // Appears often, so this is just to simplify computations

  double R = 4.0 * OmegaR0 / (3.0 * OmegaB0 * a);
  double Theta2 = -20.0 * ck_Hp / (45.0 * dtaudx) * Theta[1];
  double Psi = -Phi - 12.0*H0*H0 / (c*c * k*k * a*a) * OmegaR0 * Theta2; // Ignore neutrinos

  // dPhidx
  double dPhidx_term1 = Psi - c*c * k*k / (3.0*Hp*Hp) * Phi;
  double dPhidx_term2_1 = OmegaCDM0 / a * delta_cdm + OmegaB0 / a * delta_b;
  double dPhidx_term2_2 = 4.0*OmegaR0 /(a*a) * Theta[0];
  dPhidx = dPhidx_term1 + H0*H0/(2.0*Hp*Hp) * (dPhidx_term2_1 + dPhidx_term2_2);

  // dThetadx0
  dThetadx[0] = -ck_Hp * Theta[1] - dPhidx;

  // q
  double q_term1 = -((1.0-R) * dtaudx + (1.0+R) * ddtauddx) * (3.0*Theta[1] + v_b);
  double q_term2 = -ck_Hp * Psi;
  double q_term3 = (1.0 - dHpdx/Hp) * ck_Hp * (-Theta[0] + 2.0*Theta2);
  double q_term4 = -ck_Hp * dThetadx[0];
  double q_numerator = q_term1 + q_term2 + q_term3 + q_term4;
  double q_denominator = (1.0 + R) * dtaudx + dHpdx/Hp - 1.0;
  double q = q_numerator / q_denominator;

  // delta_CDM, v_cdm, delta_b
  ddelta_cdmdx = ck_Hp * v_cdm - 3.0*dPhidx;
  ddelta_bdx = ck_Hp * v_b - 3.0*dPhidx;
  dv_cdmdx = -v_cdm - ck_Hp * Psi;

  // dv_bdx
  double dv_bdx_term1 = -v_b - ck_Hp * Psi;
  double dv_bdx_term2 = R * (q + ck_Hp * (-Theta[0] + 2.0*Theta2) - ck_Hp * Psi);
  dv_bdx = 1.0 / (1.0 + R) * (dv_bdx_term1 + dv_bdx_term2);

  // Theta1
  dThetadx[1] = 1.0/3.0 * (q - dv_bdx);

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];
  
  double a = exp(x);
  double c = Constants.c;

  // Cosmological parameters and variables
  double H0 = cosmo->get_H0();
  double Hp = cosmo->Hp_of_x(x);
  double OmegaB0 = cosmo->get_OmegaB();
  double OmegaR0 = cosmo->get_OmegaR();
  double OmegaCDM0 = cosmo->get_OmegaCDM();
  double eta = cosmo->eta_of_x(x);

  double R = 4.0*OmegaR0 / (3.0*OmegaB0 * a);

  // Recombination variables
  double dtaudx = rec->dtaudx_of_x(x);

  double ck_Hp = c*k/Hp;

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // Psi
  double Psi = -Phi - 12.0*H0*H0/(c*c * k*k * a*a) * OmegaR0 * Theta[2]; // Ignore neutrinos

  // Phi
  double Phi_term1 = Psi - 1.0/3.0 * ck_Hp * ck_Hp * Phi;
  double Phi_term2_1 = OmegaCDM0 / a * delta_cdm + OmegaB0 / a * delta_b;
  double Phi_term2_2 = 4.0 * OmegaR0 / (a*a) * Theta[0];
  dPhidx = Phi_term1 + H0*H0/(2.0*Hp*Hp) * (Phi_term2_1 + Phi_term2_2);

  // delta_cdm
  ddelta_cdmdx = ck_Hp * v_cdm - 3.0*dPhidx;
  dv_cdmdx = -v_cdm - ck_Hp * Psi;
  ddelta_bdx = ck_Hp * v_b - 3.0*dPhidx;
  dv_bdx = -v_b - ck_Hp * Psi + dtaudx * R * (3.0*Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -ck_Hp * Theta[1] - dPhidx;
  
  double Theta1_term1 = ck_Hp / 3.0 * Theta[0] - 2.0/3.0 * ck_Hp * Theta[2];
  double Theta1_term2 = ck_Hp / 3.0 * Psi + dtaudx * (Theta[1] + v_b / 3.0);
  dThetadx[1] = Theta1_term1 + Theta1_term2;

  double Theta2_term1 = 2.0/(2.0*2.0 + 1.0) * ck_Hp * Theta[1];
  double Theta2_term2 = -(2.0 + 1.0)/(2.0*2.0 + 1.0) * ck_Hp * Theta[3];
  double Theta2_term3 = dtaudx * (Theta[2] - 0.1 * Theta[2]);
  dThetadx[2] = Theta2_term1 + Theta2_term2 + Theta2_term3;

  double Thetal_term1;
  double Thetal_term2;
  double Thetal_term3;
  for (int l=3; l < n_ell_theta-1; l++) {
    Thetal_term1 = l/(2.0*l + 1.0) * ck_Hp * Theta[l-1];
    Thetal_term2 = -(l + 1.0)/(2.0*l + 1.0) * ck_Hp * Theta[l+1];
    Thetal_term3 = dtaudx * Theta[l];
    dThetadx[l] = Thetal_term1 + Thetal_term2 + Thetal_term3;
  }

  // l = l_max
  int l = n_ell_theta - 1;
  dThetadx[l] = ck_Hp * Theta[l-1] - c * (l+1.0)/(Hp*eta) * Theta[l] + dtaudx*Theta[l];

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_dv_bdx(const double x, const double k) const{
  return v_b_spline.deriv_x(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_dPhidx(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_dPsidx(const double x, const double k) const{
  return Psi_spline.deriv_x(x,y k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_dThetadx(const double x, const double k, const int ell) const{
  return Theta_spline[ell].deriv_x(x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_v_b(x,k)       << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    // fp << get_Pi(x,k)        << " ";
    // fp << get_Source_T(x,k)  << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}