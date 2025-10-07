'''
Basic dose-response model example for COVASIM.

This example demonstrates how to configure and run a simulation using the 
dose-response transmission model instead of the default beta-based approach.
'''

import numpy as np
import covasim as cv
import matplotlib.pyplot as plt

# Example 1: Basic dose-response simulation
def basic_dose_response_example():
    '''Simple example with default dose-response parameters'''
    
    print("Running basic dose-response simulation...")
    
    pars = {
        'pop_size': 10000,
        'pop_type': 'hybrid',  # Use hybrid to get h, s, w, c layers
        'n_days': 180,
        'use_dose_response': True,
        'N0': 50,  # Infectious dose scaling parameter
        'infectivity_scale': 2.0,  # Additional infectivity scaling
      #   'viral_load_function': {
      #       'type': 'default',
      #       'peak_load': 1e8,
      #       'peak_day': 4.0
      #   }
      #   'viral_load_function': {
      #       'type': 'hcd',
      #   }
        'viral_load_function': {
            'type': 'hcd_database',
            'n_profiles': 20,  # Number of profiles to sample
            # 'seed': 32,  # Random seed for reproducibility
        }                
    }
    
    # Create and run simulation
    sim = cv.Sim(pars)
    sim.run()
    
    print(f"Final infections: {sim.results['cum_infections'][-1]}")
    
    # Plot results
    sim.plot(['cum_infections', 'new_infections'])
    plt.show()
    
    return sim


# Comparison: Traditional beta-based transmission model
def basic_traditional_example():
    '''Equivalent example using traditional beta-based transmission for comparison'''
    
    print("Running basic traditional (beta-based) simulation...")
    
    pars = {
        'pop_size': 10000,
        'pop_type': 'hybrid',  # Use hybrid to get h, s, w, c layers
        'n_days': 180,
        'use_dose_response': False,  # Use traditional beta-based model
        # Traditional parameters - these are COVASIM defaults
        'beta': 0.016,  # Base transmission probability per contact
        # Layer-specific beta multipliers (can be customized)
      #   'beta_layer': {'h': 1.0, 's': 0.5, 'w': 0.5, 'c': 0.1}
    }
    
    # Create and run simulation
    sim = cv.Sim(pars)
    sim.run()
    
    print(f"Final infections: {sim.results['cum_infections'][-1]}")
    
    # Plot results
    sim.plot(['cum_infections', 'new_infections'])
    plt.show()
    
    return sim


# Example 2: Custom deposition rates for different environments
def custom_deposition_rates_example():
    '''Example with custom deposition rates for different contact layers'''
    
    print("Running simulation with custom deposition rates...")
    
    # Define custom deposition rates based on environment characteristics
    # Units: viral particles deposited per minute of contact
    custom_direct_rates = {
            'a': 6.0e-9,   # All contacts (random population) - average value
            'h': 10.e-9,   # Household - close proximity, no masks typically
            's': 5.0e-9,   # School - moderate proximity
            'w': 5.0e-9,   # Work - moderate proximity  
            'c': 3.0e-9,   # Community - variable proximity
            'l': 8.0e-9,   # Long-term care facilities - close proximity care
    }
    
    custom_indirect_rates = {
            'a': 1.5e-9,   # All contacts (random population) - average value
            'h': 2.0e-9,   # Household - smaller space, less ventilation
            's': 1.5e-9,   # School - larger space, moderate ventilation
            'w': 1.5e-9,   # Work - office spaces
            'c': 0.5e-9,   # Community - highly variable, often well-ventilated
            'l': 2.5e-9,   # Long-term care facilities - confined spaces
    }
    
    pars = {
        'pop_size': 1000,
        'pop_type': 'hybrid',  # Use hybrid to get h, s, w, c layers
        'n_days': 180,
        'use_dose_response': True,
        'N0': 150,
        'infectivity_scale': 1.2,
        'deposition_rates_direct': custom_direct_rates,
        'deposition_rates_indirect': custom_indirect_rates,
        'viral_load_function': {
            'type': 'default',
            'peak_load': 5e7,
            'peak_day': 4.5
        }
    }
    
    sim = cv.Sim(pars)
    sim.run()
    
    print(f"Final infections with custom rates: {sim.results['cum_infections'][-1]}")
    
    return sim


# Example 3: Different viral load models
def viral_load_comparison_example():
    '''Compare different viral load models'''
    
    print("Comparing different viral load models...")
    
    # Define different viral load configurations
    viral_load_configs = {
        'Default': {
            'type': 'default',
            'peak_load': 1e8,
            'peak_day': 5.0
        },
        'Gaussian': {
            'type': 'gaussian',
            'peak_load': 1e8,
            'peak_day': 5.0,
            'sigma': 2.5
        },
        'Piecewise': {
            'type': 'piecewise',
            'times': [0, 2, 4, 7, 14, 21],
            'values': [0, 1e6, 1e8, 5e7, 1e5, 0]
        },
        'HCD Model': {
            'type': 'hcd',
            'sample_params': False  # Use mean parameters
        }
    }
    
    base_pars = {
        'pop_size': 1000,
        'pop_type': 'hybrid',  # Use hybrid to get h, s, w, c layers
        'n_days': 180,
        'use_dose_response': True,
        'N0': 100,
        'infectivity_scale': 1.0
    }
    
    results = {}
    
    for name, vl_config in viral_load_configs.items():
        print(f"  Running with {name} viral load...")
        
        pars = base_pars.copy()
        pars['viral_load_function'] = vl_config
        
        try:
            sim = cv.Sim(pars)
            sim.run()
            results[name] = sim.results['cum_infections'][-1]
        except Exception as e:
            print(f"    Error with {name}: {e}")
            results[name] = None
    
    print("\nFinal infection counts by viral load model:")
    for name, count in results.items():
        if count is not None:
            print(f"  {name}: {count}")
        else:
            print(f"  {name}: Failed")
    
    return results


# Example 4: Custom viral load function
def custom_viral_load_example():
    '''Example using a user-defined viral load function'''
    
    print("Running simulation with custom viral load function...")
    
    # Define custom viral load function
    def bimodal_viral_load(t):
        '''Bimodal viral load with two peaks'''
        # First peak around day 3
        peak1 = 5e7 * np.exp(-0.5 * ((t - 3) / 1.5)**2)
        # Second peak around day 8
        peak2 = 3e7 * np.exp(-0.5 * ((t - 8) / 2.0)**2)
        return np.maximum(peak1 + peak2, 0)
    
    pars = {
        'pop_size': 800,
        'pop_type': 'hybrid',  # Use hybrid to get h, s, w, c layers
        'n_days': 50,
        'use_dose_response': True,
        'N0': 120,
        'infectivity_scale': 0.8,
        'viral_load_function': bimodal_viral_load  # Pass function directly
    }
    
    sim = cv.Sim(pars)
    sim.run()
    
    print(f"Final infections with bimodal viral load: {sim.results['cum_infections'][-1]}")
    
    # Plot viral load function
    times = np.linspace(0, 20, 200)
    loads = [bimodal_viral_load(t) for t in times]
    
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(times, loads)
    plt.xlabel('Days since infection')
    plt.ylabel('Viral load')
    plt.title('Custom Bimodal Viral Load')
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    plt.plot(sim.results['new_infections'])
    plt.xlabel('Day')
    plt.ylabel('New infections')
    plt.title('Simulation Results')
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()
    
    return sim


# Example 5: Sensitivity analysis
def sensitivity_analysis_example():
    '''Perform sensitivity analysis on dose-response parameters'''
    
    print("Running dose-response sensitivity analysis...")
    
    # Parameter ranges to test
    N0_values = [50, 100, 200, 500]
    infectivity_scales = [0.5, 1.0, 1.5, 2.0]
    
    base_pars = {
        'pop_size': 500,
        'pop_type': 'hybrid',  # Use hybrid to get h, s, w, c layers
        'n_days': 40,
        'use_dose_response': True,
        'viral_load_function': {
            'type': 'default',
            'peak_load': 1e8
        }
    }
    
    results_grid = np.zeros((len(N0_values), len(infectivity_scales)))
    
    print("Testing parameter combinations:")
    for i, N0 in enumerate(N0_values):
        for j, inf_scale in enumerate(infectivity_scales):
            print(f"  N0={N0}, infectivity_scale={inf_scale}")
            
            pars = base_pars.copy()
            pars['N0'] = N0
            pars['infectivity_scale'] = inf_scale
            
            sim = cv.Sim(pars)
            sim.run()
            results_grid[i, j] = sim.results['cum_infections'][-1]
    
    # Plot heatmap
    plt.figure(figsize=(8, 6))
    plt.imshow(results_grid, aspect='auto', cmap='viridis')
    plt.colorbar(label='Final infections')
    plt.xlabel('Infectivity scale')
    plt.ylabel('N0 (infectious dose)')
    plt.xticks(range(len(infectivity_scales)), infectivity_scales)
    plt.yticks(range(len(N0_values)), N0_values)
    plt.title('Dose-Response Parameter Sensitivity')
    
    # Add text annotations
    for i in range(len(N0_values)):
        for j in range(len(infectivity_scales)):
            plt.text(j, i, f'{int(results_grid[i, j])}', 
                    ha='center', va='center', color='white' if results_grid[i, j] < results_grid.max()/2 else 'black')
    
    plt.tight_layout()
    plt.show()
    
    return results_grid


# Example 6: Mask usage and deposition rate modifications
def mask_and_deposition_example():
    '''
    Example demonstrating how to simulate mask usage and play with deposition rates.

    Note: Mask usage is simulated by modifying deposition rates based on mask efficacy.
    This is a workaround since explicit mask interventions aren't implemented yet.
    '''

    print("Running mask usage and deposition rate example...")

    # Helper function to apply mask effects to deposition rates
    def apply_mask_effect(base_rates, mask_efficacy, adoption_rate=1.0):
        '''
        Apply mask effect to deposition rates.

        Args:
            base_rates (dict): Base deposition rates by layer
            mask_efficacy (float): Mask efficacy (0-1), e.g., 0.5 = 50% reduction
            adoption_rate (float): Fraction of population wearing masks (0-1)

        Returns:
            dict: Modified deposition rates
        '''
        # Effective reduction = efficacy × adoption rate
        reduction_factor = 1.0 - (mask_efficacy * adoption_rate)

        masked_rates = {}
        for layer, rate in base_rates.items():
            masked_rates[layer] = rate * reduction_factor

        return masked_rates

    # Define baseline deposition rates (no masks)
    baseline_direct = {
        'h': 10.0e-9,   # Household
        's': 8.0e-9,    # School
        'w': 6.0e-9,    # Work
        'c': 4.0e-9,    # Community
    }

    baseline_indirect = {
        'h': 2.0e-9,    # Household
        's': 1.5e-9,    # School
        'w': 1.2e-9,    # Work
        'c': 0.8e-9,    # Community
    }

    # Define mask efficacies (based on literature)
    mask_types = {
        'No masks': {'efficacy': 0.0, 'adoption': 1.0},
        'Cloth masks (50% adoption)': {'efficacy': 0.30, 'adoption': 0.5},
        'Surgical masks (70% adoption)': {'efficacy': 0.50, 'adoption': 0.7},
        'N95 masks (90% adoption)': {'efficacy': 0.80, 'adoption': 0.9},
    }

    # Base simulation parameters
    base_pars = {
        'pop_size': 5000,
        'pop_type': 'hybrid',
        'n_days': 120,
        'use_dose_response': True,
        'N0': 100,
        'infectivity_scale': 1.0,
        'viral_load_function': {
            'type': 'default',
            'peak_load': 1e8,
            'peak_day': 5.0
        }
    }

    # Run scenarios
    results = {}
    sims = {}

    for scenario_name, mask_params in mask_types.items():
        print(f"\n  Running: {scenario_name}")

        # Apply mask effect to deposition rates
        masked_direct = apply_mask_effect(
            baseline_direct,
            mask_params['efficacy'],
            mask_params['adoption']
        )
        masked_indirect = apply_mask_effect(
            baseline_indirect,
            mask_params['efficacy'],
            mask_params['adoption']
        )

        # Create simulation with modified deposition rates
        pars = base_pars.copy()
        pars['deposition_rates_direct'] = masked_direct
        pars['deposition_rates_indirect'] = masked_indirect

        sim = cv.Sim(pars)
        sim.run()

        results[scenario_name] = sim.results['cum_infections'][-1]
        sims[scenario_name] = sim

        print(f"    Final infections: {results[scenario_name]:,.0f}")

    # Plot comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Time series comparison
    ax1 = axes[0]
    for scenario_name, sim in sims.items():
        ax1.plot(sim.results['cum_infections'], label=scenario_name, linewidth=2)

    ax1.set_xlabel('Day')
    ax1.set_ylabel('Cumulative Infections')
    ax1.set_title('Impact of Mask Usage on Infection Spread')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Bar chart of final infections
    ax2 = axes[1]
    scenario_names = list(results.keys())
    final_infections = list(results.values())
    colors = ['red', 'orange', 'yellow', 'green']

    bars = ax2.bar(range(len(scenario_names)), final_infections, color=colors, alpha=0.7)
    ax2.set_xticks(range(len(scenario_names)))
    ax2.set_xticklabels(scenario_names, rotation=15, ha='right')
    ax2.set_ylabel('Final Infections')
    ax2.set_title('Final Infection Count by Scenario')
    ax2.grid(True, alpha=0.3, axis='y')

    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, final_infections)):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(val):,}',
                ha='center', va='bottom')

    plt.tight_layout()
    plt.show()

    # Print summary
    print("\n" + "="*60)
    print("MASK USAGE IMPACT SUMMARY")
    print("="*60)
    baseline_infections = results['No masks']
    for scenario_name, infections in results.items():
        reduction = (baseline_infections - infections) / baseline_infections * 100
        print(f"{scenario_name:40s}: {infections:7,.0f} infections ({reduction:+5.1f}%)")
    print("="*60)

    return sims, results


# Example 7: Enhanced intervention scenario with time-dependent changes
def intervention_scenario_example():
    '''
    Demonstrate time-dependent interventions: mask mandate, improved ventilation.

    This shows how to simulate interventions that start on specific days
    by running multiple simulations with different parameter sets.
    '''

    print("Running intervention scenario example...")

    # Scenario: Mask mandate introduced on day 30
    baseline_direct = {
            'h': 5.0e-9,   # Household - close proximity, no masks typically
            's': 5.0e-9,   # School - moderate proximity
            'w': 5.0e-9,   # Work - moderate proximity  
            'c': 5.0e-9,   # Community - variable proximity
    }

    baseline_indirect = {
            'h': 0.5e-9,   # Household - smaller space, less ventilation
            's': 0.1e-9,   # School - larger space, moderate ventilation
            'w': 0.1e-9,   # Work - office spaces
            'c': 0.1e-9,   # Community - highly variable, often well-ventilated
    }

    # After intervention: masks + improved ventilation
    # Masks: 60% efficacy, 80% adoption in schools/work/community
    # (households remain unchanged)
    intervention_direct = {
        'h': 5.0e-9,   # No change at home
        's': 5.0e-9* (1 - 0.8 * 0.8),   # 60% efficacy × 80% adoption
        'w': 5.0e-9* (1 - 0.8 * 0.8),
        'c': 5.0e-9* (1 - 0.8 * 0.8),
    }

    # Improved ventilation reduces indirect transmission by 40%
    intervention_indirect = {
        'h': 0.5e-9,     # No change at home
        's': 0.1e-9 * 0.6,   # 40% reduction
        'w': 0.1e-9 * 0.6,
        'c': 0.1e-9 * 0.6,
    }

    base_pars = {
        'pop_size': 10000,
        'pop_type': 'hybrid',
        'n_days': 180,
        'use_dose_response': True,
        'N0': 50,
        'infectivity_scale': 2.0,
      #   'viral_load_function': {
      #       'type': 'default',
      #       'peak_load': 1e8,
      #       'peak_day': 5.0
      #   }
        'viral_load_function': {
            'type': 'hcd_database',
            'n_profiles': 20,  # Number of profiles to sample
            # 'seed': 32,  # Random seed for reproducibility
        }      
    }

    # Scenario 1: No intervention
    print("\n  Scenario 1: No intervention (baseline)")
    pars1 = base_pars.copy()
    pars1['deposition_rates_direct'] = baseline_direct
    pars1['deposition_rates_indirect'] = baseline_indirect
    sim1 = cv.Sim(pars1)
    sim1.run()

    # Scenario 2: Intervention from day 30
    # Note: True time-dependent parameter changes would require
    # modifying COVASIM's intervention system. Here we approximate
    # by using average rates weighted by time period.
    print("  Scenario 2: Masks + ventilation from day 30")

    # Weighted average: 30 days baseline + 60 days intervention
    weighted_direct = {}
    weighted_indirect = {}

    for layer in baseline_direct.keys():
        weighted_direct[layer] = (
            baseline_direct[layer] * 30 +
            intervention_direct[layer] * 60
        ) / 90

    for layer in baseline_indirect.keys():
        weighted_indirect[layer] = (
            baseline_indirect[layer] * 30 +
            intervention_indirect[layer] * 60
        ) / 90

    pars2 = base_pars.copy()
    pars2['deposition_rates_direct'] = weighted_direct
    pars2['deposition_rates_indirect'] = weighted_indirect
    sim2 = cv.Sim(pars2)
    sim2.run()

    # Scenario 3: Intervention from day 0
    print("  Scenario 3: Masks + ventilation from day 0")
    pars3 = base_pars.copy()
    pars3['deposition_rates_direct'] = intervention_direct
    pars3['deposition_rates_indirect'] = intervention_indirect
    sim3 = cv.Sim(pars3)
    sim3.run()

    # Plot comparison
    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    plt.plot(sim1.results['cum_infections'], label='No intervention', linewidth=2, color='red')
    plt.plot(sim2.results['cum_infections'], label='Intervention from day 30', linewidth=2, color='orange')
    plt.plot(sim3.results['cum_infections'], label='Intervention from day 0', linewidth=2, color='green')
    plt.axvline(x=30, color='gray', linestyle='--', alpha=0.5, label='Intervention start')
    plt.xlabel('Day')
    plt.ylabel('Cumulative Infections')
    plt.title('Impact of Intervention Timing')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(1, 2, 2)
    plt.plot(sim1.results['new_infections'], label='No intervention', linewidth=2, color='red')
    plt.plot(sim2.results['new_infections'], label='Intervention from day 30', linewidth=2, color='orange')
    plt.plot(sim3.results['new_infections'], label='Intervention from day 0', linewidth=2, color='green')
    plt.axvline(x=30, color='gray', linestyle='--', alpha=0.5)
    plt.xlabel('Day')
    plt.ylabel('Daily New Infections')
    plt.title('Daily Infection Rate')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    # Print summary
    print("\n" + "="*60)
    print("INTERVENTION TIMING ANALYSIS")
    print("="*60)
    print(f"No intervention:          {sim1.results['cum_infections'][-1]:7,.0f} infections")
    print(f"Intervention from day 30: {sim2.results['cum_infections'][-1]:7,.0f} infections")
    print(f"Intervention from day 0:  {sim3.results['cum_infections'][-1]:7,.0f} infections")

    early_reduction = (sim1.results['cum_infections'][-1] - sim3.results['cum_infections'][-1]) / sim1.results['cum_infections'][-1] * 100
    delayed_reduction = (sim1.results['cum_infections'][-1] - sim2.results['cum_infections'][-1]) / sim1.results['cum_infections'][-1] * 100

    print(f"\nEarly intervention impact:   {early_reduction:.1f}% reduction")
    print(f"Delayed intervention impact: {delayed_reduction:.1f}% reduction")
    print("="*60)

    return {'no_intervention': sim1, 'day_30': sim2, 'day_0': sim3}


# Example 8: Simple ventilation comparison
def ventilation_comparison_example():
    '''
    Compare two ventilation scenarios: poor vs improved ventilation.

    This demonstrates how changing both direct and indirect transmission rates
    affects infection spread. Poor ventilation increases both face-to-face
    exposure (due to stagnant air) and indirect airborne transmission.
    '''

    print("Running ventilation comparison example...")

    # Scenario 1: Poor ventilation (baseline)
    # Higher rates for both direct and indirect transmission
    poor_ventilation_direct = {
        'h': 5.0e-9,   # Household - stagnant air increases direct exposure
        's': 5.0e-9,   # School
        'w': 5.0e-9,   # Work
        'c': 5.0e-9,   # Community
    }

    poor_ventilation_indirect = {
        'h': 0.5e-9 ,   # Household - poor ventilation, high airborne concentration
        's': 0.1e-9 ,   # School
        'w': 0.1e-9 ,  # Work
        'c': 0.1e-9 ,   # Community
    }

    # Scenario 2: Improved ventilation
    # 40% reduction in direct (clearer air reduces droplet exposure)
    # 60% reduction in indirect (better air circulation removes airborne particles)
    improved_ventilation_direct = {
        'h': 5.0e-9 * 0.5,   # 40% reduction
        's': 5.0e-9 * 0.5,
        'w': 5.5e-9 * 0.5,
        'c': 5.5e-9 * 0.5,
    }

    improved_ventilation_indirect = {
        'h': 0.5e-9 * 0.5,   # 60% reduction
        's': 0.1e-9 * 0.5,
        'w': 0.1e-9 * 0.5,
        'c': 0.1e-9 * 0.5,
    }

    base_pars = {
        'pop_size': 10000,
        'pop_type': 'hybrid',
        'n_days': 180,
        'use_dose_response': True,
        'N0': 50,
        'infectivity_scale': 2.0,
        'viral_load_function': {
            'type': 'hcd_database',
            'n_profiles': 20,
        }
    }

    # Scenario 1: Poor ventilation
    print("\n  Scenario 1: Poor ventilation (high direct and indirect transmission)")
    pars1 = base_pars.copy()
    pars1['deposition_rates_direct'] = poor_ventilation_direct
    pars1['deposition_rates_indirect'] = poor_ventilation_indirect
    sim1 = cv.Sim(pars1)
    sim1.run()

    # Scenario 2: Improved ventilation
    print("  Scenario 2: Improved ventilation (40% reduction direct, 60% reduction indirect)")
    pars2 = base_pars.copy()
    pars2['deposition_rates_direct'] = improved_ventilation_direct
    pars2['deposition_rates_indirect'] = improved_ventilation_indirect
    sim2 = cv.Sim(pars2)
    sim2.run()

    # Plot comparison
    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    plt.plot(sim1.results['cum_infections'], label='Poor ventilation', linewidth=2, color='red')
    plt.plot(sim2.results['cum_infections'], label='Improved ventilation', linewidth=2, color='green')
    plt.xlabel('Day')
    plt.ylabel('Cumulative Infections')
    plt.title('Impact of Ventilation Improvement')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(1, 2, 2)
    plt.plot(sim1.results['new_infections'], label='Poor ventilation', linewidth=2, color='red')
    plt.plot(sim2.results['new_infections'], label='Improved ventilation', linewidth=2, color='green')
    plt.xlabel('Day')
    plt.ylabel('Daily New Infections')
    plt.title('Daily Infection Rate')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    # Print summary
    poor_vent_infections = sim1.results['cum_infections'][-1]
    improved_vent_infections = sim2.results['cum_infections'][-1]
    reduction = (poor_vent_infections - improved_vent_infections) / poor_vent_infections * 100

    print("\n" + "="*60)
    print("VENTILATION COMPARISON RESULTS")
    print("="*60)
    print(f"Poor ventilation:      {poor_vent_infections:7,.0f} infections")
    print(f"Improved ventilation:  {improved_vent_infections:7,.0f} infections")
    print(f"Reduction:             {reduction:7.1f}%")
    print("="*60)

    return {'poor_ventilation': sim1, 'improved_ventilation': sim2}


# Comparison example: Run both models side by side
def compare_transmission_models():
    '''Compare dose-response vs traditional beta-based transmission models'''
    
    print("=== Model Comparison ===")
    print("Running identical simulations with both transmission models...\n")
    
    # Run dose-response model
    print("1. Dose-Response Model:")
    sim_dose = basic_dose_response_example()
    dose_infections = sim_dose.results['cum_infections'][-1]
    
    print("\n" + "="*50 + "\n")
    
    # Run traditional model
    print("2. Traditional Beta Model:")
    sim_traditional = basic_traditional_example()
    traditional_infections = sim_traditional.results['cum_infections'][-1]
    
    # Compare results
    print("\n" + "="*50)
    print("COMPARISON RESULTS:")
    print(f"Dose-Response Model:  {dose_infections:,} infections")
    print(f"Traditional Model:    {traditional_infections:,} infections")
    print(f"Difference:           {dose_infections - traditional_infections:,} infections")
    print(f"Ratio (Dose/Traditional): {dose_infections/traditional_infections:.2f}")
    
    return {'dose_response': sim_dose, 'traditional': sim_traditional}


if __name__ == '__main__':
    # Run examples
    print("=== Dose-Response Model Examples ===\n")

   #  # Example 1: Dose-Response Model
   #  sim1 = basic_dose_response_example()
   #  print()

   #  # Traditional Model Comparison
   #  sim_traditional = basic_traditional_example()
   #  print()

    # Side-by-side comparison (uncomment to run)
    # comparison_results = compare_transmission_models()
    # print()

   #  # Example 2: Custom deposition rates
   #  sim2 = custom_deposition_rates_example()
   #  print()

   #  # Example 3: Viral load comparison
   #  results3 = viral_load_comparison_example()
   #  print()

   #  # Example 4: Custom viral load
   #  sim4 = custom_viral_load_example()
   #  print()

   #  # Example 5: Sensitivity analysis
   #  results5 = sensitivity_analysis_example()
   #  print()

   #  # Example 6: Mask usage and deposition rates (NEW!)
   #  sims6, results6 = mask_and_deposition_example()
   #  print()

   #  # Example 7: Intervention scenarios (NEW!)
   #  sims7 = intervention_scenario_example()
   #  print()

    # Example 7: Intervention scenarios (NEW!)
    sims7 = ventilation_comparison_example()
    print()

    print("\n=== All examples completed ===")