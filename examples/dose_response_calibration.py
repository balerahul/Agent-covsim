'''
Dose-response model calibration example for COVASIM.

This example demonstrates how to calibrate the dose-response model parameters
to match observed epidemiological data and compare with traditional beta-based models.
'''

import numpy as np
import covasim as cv
import matplotlib.pyplot as plt
import sciris as sc
from scipy.optimize import minimize


def load_example_data():
    '''Load or generate example epidemiological data for calibration'''
    
    # For this example, we'll generate synthetic "observed" data
    # In practice, you would load real epidemic data
    
    days = np.arange(60)
    
    # Synthetic epidemic curve (daily new infections)
    peak_day = 25
    peak_infections = 150
    width = 8
    
    # Generate realistic epidemic curve
    new_infections = peak_infections * np.exp(-0.5 * ((days - peak_day) / width)**2)
    new_infections += np.random.poisson(5, len(days))  # Add noise
    new_infections = np.maximum(new_infections, 0)
    
    # Calculate cumulative infections
    cum_infections = np.cumsum(new_infections)
    
    return {
        'days': days,
        'new_infections': new_infections,
        'cum_infections': cum_infections
    }


def run_dose_response_sim(N0=100, infectivity_scale=1.0, viral_load_params=None):
    '''Run simulation with specified dose-response parameters'''
    
    if viral_load_params is None:
        viral_load_params = {'type': 'default', 'peak_load': 1e8}
    
    pars = {
        'pop_size': 5000,
        'n_days': 60,
        'use_dose_response': True,
        'N0': N0,
        'infectivity_scale': infectivity_scale,
        'viral_load_function': viral_load_params,
        'rand_seed': 42  # For reproducibility
    }
    
    sim = cv.Sim(pars)
    sim.run()
    
    return sim


def calculate_fit_metric(sim_data, obs_data, metric='rmse'):
    '''Calculate goodness of fit between simulation and observed data'''
    
    sim_new = np.diff(sim_data['cum_infections'])
    obs_new = obs_data['new_infections'][1:]  # Skip first day
    
    # Ensure same length
    min_len = min(len(sim_new), len(obs_new))
    sim_new = sim_new[:min_len]
    obs_new = obs_new[:min_len]
    
    if metric == 'rmse':
        return np.sqrt(np.mean((sim_new - obs_new)**2))
    elif metric == 'mae':
        return np.mean(np.abs(sim_new - obs_new))
    elif metric == 'mape':
        return np.mean(np.abs((sim_new - obs_new) / (obs_new + 1))) * 100
    else:
        raise ValueError(f"Unknown metric: {metric}")


def calibrate_dose_response_basic():
    '''Basic calibration of dose-response parameters'''
    
    print("=== Basic Dose-Response Calibration ===")
    
    # Load target data
    obs_data = load_example_data()
    
    # Define objective function
    def objective(params):
        N0, infectivity_scale = params
        
        # Run simulation
        sim = run_dose_response_sim(N0=N0, infectivity_scale=infectivity_scale)
        
        # Calculate fit
        fit = calculate_fit_metric(sim.results, obs_data)
        
        return fit
    
    # Initial guess and bounds
    initial_params = [100, 1.0]
    bounds = [(10, 1000), (0.1, 5.0)]
    
    print("Starting calibration...")
    result = minimize(objective, initial_params, bounds=bounds, method='L-BFGS-B')
    
    optimal_N0, optimal_infectivity = result.x
    print(f"Optimal parameters:")
    print(f"  N0: {optimal_N0:.1f}")
    print(f"  Infectivity scale: {optimal_infectivity:.3f}")
    print(f"  Final RMSE: {result.fun:.2f}")
    
    # Run final simulation with optimal parameters
    best_sim = run_dose_response_sim(N0=optimal_N0, infectivity_scale=optimal_infectivity)
    
    # Plot results
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 2, 1)
    sim_new = np.diff(best_sim.results['cum_infections'])
    plt.plot(obs_data['days'][1:], obs_data['new_infections'][1:], 'ko-', label='Observed', alpha=0.7)
    plt.plot(obs_data['days'][1:len(sim_new)+1], sim_new, 'r-', label='Dose-response model', linewidth=2)
    plt.xlabel('Day')
    plt.ylabel('New infections')
    plt.title('Daily New Infections')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.plot(obs_data['days'], obs_data['cum_infections'], 'ko-', label='Observed', alpha=0.7)
    plt.plot(obs_data['days'][:len(best_sim.results['cum_infections'])], 
             best_sim.results['cum_infections'], 'r-', label='Dose-response model', linewidth=2)
    plt.xlabel('Day')
    plt.ylabel('Cumulative infections')
    plt.title('Cumulative Infections')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return best_sim, result


def compare_models():
    '''Compare dose-response model with traditional beta-based model'''
    
    print("\n=== Model Comparison ===")
    
    obs_data = load_example_data()
    
    # Run dose-response model (using calibrated parameters from previous example)
    dose_sim = run_dose_response_sim(N0=120, infectivity_scale=1.2)
    
    # Run traditional beta-based model
    beta_pars = {
        'pop_size': 5000,
        'n_days': 60,
        'use_dose_response': False,
        'beta': 0.016,  # Calibrated beta value
        'rand_seed': 42
    }
    beta_sim = cv.Sim(beta_pars)
    beta_sim.run()
    
    # Calculate fit metrics
    dose_rmse = calculate_fit_metric(dose_sim.results, obs_data)
    beta_rmse = calculate_fit_metric(beta_sim.results, obs_data)
    
    print(f"Model comparison (RMSE):")
    print(f"  Dose-response model: {dose_rmse:.2f}")
    print(f"  Beta-based model: {beta_rmse:.2f}")
    
    # Plot comparison
    plt.figure(figsize=(15, 5))
    
    # Daily new infections
    plt.subplot(1, 3, 1)
    dose_new = np.diff(dose_sim.results['cum_infections'])
    beta_new = np.diff(beta_sim.results['cum_infections'])
    
    plt.plot(obs_data['days'][1:], obs_data['new_infections'][1:], 'ko-', 
             label='Observed', alpha=0.7, markersize=4)
    plt.plot(obs_data['days'][1:len(dose_new)+1], dose_new, 'r-', 
             label=f'Dose-response (RMSE={dose_rmse:.1f})', linewidth=2)
    plt.plot(obs_data['days'][1:len(beta_new)+1], beta_new, 'b--', 
             label=f'Beta-based (RMSE={beta_rmse:.1f})', linewidth=2)
    plt.xlabel('Day')
    plt.ylabel('New infections')
    plt.title('Daily New Infections')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Cumulative infections
    plt.subplot(1, 3, 2)
    plt.plot(obs_data['days'], obs_data['cum_infections'], 'ko-', 
             label='Observed', alpha=0.7, markersize=4)
    plt.plot(dose_sim.results['cum_infections'], 'r-', 
             label='Dose-response', linewidth=2)
    plt.plot(beta_sim.results['cum_infections'], 'b--', 
             label='Beta-based', linewidth=2)
    plt.xlabel('Day')
    plt.ylabel('Cumulative infections')
    plt.title('Cumulative Infections')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # R_effective comparison
    plt.subplot(1, 3, 3)
    dose_r_eff = dose_sim.compute_r_eff()
    beta_r_eff = beta_sim.compute_r_eff()
    
    plt.plot(dose_r_eff, 'r-', label='Dose-response', linewidth=2)
    plt.plot(beta_r_eff, 'b--', label='Beta-based', linewidth=2)
    plt.axhline(y=1, color='k', linestyle=':', alpha=0.5)
    plt.xlabel('Day')
    plt.ylabel('R_effective')
    plt.title('Effective Reproduction Number')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return dose_sim, beta_sim


def viral_load_calibration():
    '''Calibrate viral load parameters alongside dose-response parameters'''
    
    print("\n=== Viral Load Calibration ===")
    
    obs_data = load_example_data()
    
    def objective_viral_load(params):
        N0, infectivity_scale, peak_load, peak_day = params
        
        viral_load_params = {
            'type': 'default',
            'peak_load': peak_load,
            'peak_day': peak_day
        }
        
        sim = run_dose_response_sim(
            N0=N0, 
            infectivity_scale=infectivity_scale,
            viral_load_params=viral_load_params
        )
        
        return calculate_fit_metric(sim.results, obs_data)
    
    # Initial guess and bounds
    initial_params = [100, 1.0, 1e8, 5.0]
    bounds = [(10, 1000), (0.1, 5.0), (1e6, 1e9), (2.0, 10.0)]
    
    print("Calibrating viral load parameters...")
    result = minimize(objective_viral_load, initial_params, bounds=bounds, method='L-BFGS-B')
    
    optimal_N0, optimal_infectivity, optimal_peak_load, optimal_peak_day = result.x
    print(f"Optimal parameters:")
    print(f"  N0: {optimal_N0:.1f}")
    print(f"  Infectivity scale: {optimal_infectivity:.3f}")
    print(f"  Peak viral load: {optimal_peak_load:.2e}")
    print(f"  Peak day: {optimal_peak_day:.1f}")
    print(f"  Final RMSE: {result.fun:.2f}")
    
    # Show optimal viral load curve
    from covasim.viral_dynamics import DefaultViralLoad
    
    optimal_vl = DefaultViralLoad(
        peak_load=optimal_peak_load,
        peak_day=optimal_peak_day
    )
    
    times = np.linspace(0, 20, 200)
    loads = [optimal_vl.get_viral_load(t) for t in times]
    
    plt.figure(figsize=(8, 5))
    plt.plot(times, loads, 'g-', linewidth=2)
    plt.xlabel('Days since infection')
    plt.ylabel('Viral load')
    plt.title(f'Optimal Viral Load Curve\n(Peak: {optimal_peak_load:.1e} at day {optimal_peak_day:.1f})')
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    plt.show()
    
    return result


def uncertainty_analysis():
    '''Perform uncertainty analysis on calibrated parameters'''
    
    print("\n=== Uncertainty Analysis ===")
    
    # Use previously calibrated parameters with uncertainty
    base_params = {
        'N0': (120, 20),  # (mean, std)
        'infectivity_scale': (1.2, 0.2),
        'peak_load': (5e7, 1e7),
        'peak_day': (4.5, 1.0)
    }
    
    n_runs = 50
    results = []
    
    print(f"Running {n_runs} uncertainty samples...")
    
    for i in range(n_runs):
        # Sample parameters
        N0 = np.random.normal(base_params['N0'][0], base_params['N0'][1])
        N0 = max(10, N0)  # Ensure positive
        
        infectivity_scale = np.random.normal(
            base_params['infectivity_scale'][0], 
            base_params['infectivity_scale'][1]
        )
        infectivity_scale = max(0.1, infectivity_scale)
        
        peak_load = np.random.normal(
            base_params['peak_load'][0], 
            base_params['peak_load'][1]
        )
        peak_load = max(1e6, peak_load)
        
        peak_day = np.random.normal(
            base_params['peak_day'][0], 
            base_params['peak_day'][1]
        )
        peak_day = max(2.0, peak_day)
        
        # Run simulation
        viral_load_params = {
            'type': 'default',
            'peak_load': peak_load,
            'peak_day': peak_day
        }
        
        sim = run_dose_response_sim(
            N0=N0,
            infectivity_scale=infectivity_scale,
            viral_load_params=viral_load_params
        )
        
        results.append({
            'N0': N0,
            'infectivity_scale': infectivity_scale,
            'peak_load': peak_load,
            'peak_day': peak_day,
            'final_infections': sim.results['cum_infections'][-1],
            'peak_day_sim': np.argmax(np.diff(sim.results['cum_infections'])),
            'cum_infections': sim.results['cum_infections']
        })
    
    # Analyze results
    final_infections = [r['final_infections'] for r in results]
    peak_days_sim = [r['peak_day_sim'] for r in results]
    
    print(f"Uncertainty analysis results:")
    print(f"  Final infections: {np.mean(final_infections):.0f} ± {np.std(final_infections):.0f}")
    print(f"  Peak day: {np.mean(peak_days_sim):.1f} ± {np.std(peak_days_sim):.1f}")
    
    # Plot uncertainty bands
    cum_infections_array = np.array([r['cum_infections'] for r in results])
    
    percentiles = [5, 25, 50, 75, 95]
    cum_percentiles = np.percentile(cum_infections_array, percentiles, axis=0)
    
    plt.figure(figsize=(10, 6))
    days = np.arange(len(cum_percentiles[0]))
    
    # Plot uncertainty bands
    plt.fill_between(days, cum_percentiles[0], cum_percentiles[4], 
                     alpha=0.2, color='blue', label='90% CI')
    plt.fill_between(days, cum_percentiles[1], cum_percentiles[3], 
                     alpha=0.4, color='blue', label='50% CI')
    plt.plot(days, cum_percentiles[2], 'b-', linewidth=2, label='Median')
    
    plt.xlabel('Day')
    plt.ylabel('Cumulative infections')
    plt.title('Uncertainty Analysis - Cumulative Infections')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
    
    return results


if __name__ == '__main__':
    print("=== Dose-Response Calibration Examples ===\n")
    
    # Basic calibration
    best_sim, calib_result = calibrate_dose_response_basic()
    
    # Model comparison
    dose_sim, beta_sim = compare_models()
    
    # Viral load calibration
    vl_result = viral_load_calibration()
    
    # Uncertainty analysis
    uncertainty_results = uncertainty_analysis()
    
    print("\n=== Calibration examples completed ===")