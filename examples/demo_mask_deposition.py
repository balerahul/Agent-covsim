'''
Demo: Playing with Mask Usage and Deposition Rates in COVASIM Dose-Response Model

This script demonstrates how users can experiment with:
1. Deposition rates for direct and indirect transmission
2. Mask usage effects (simulated via deposition rate modifications)
3. Intervention scenarios with different timings

Run this script to see interactive examples of how mask usage and ventilation
interventions affect infection spread in the dose-response model.
'''

import os
os.environ['COVASIM_INTERACTIVE'] = '0'  # Disable interactive plots

import numpy as np
import covasim as cv
import matplotlib.pyplot as plt


def simple_mask_example():
    '''
    Simplest example: Compare no masks vs surgical masks
    '''
    print("\n" + "="*70)
    print("SIMPLE EXAMPLE: No Masks vs Surgical Masks")
    print("="*70)

    # Base deposition rates (no masks)
    no_mask_direct = {
        'h': 10.0e-9,  # Household
        's': 8.0e-9,   # School
        'w': 6.0e-9,   # Work
        'c': 4.0e-9,   # Community
    }

    no_mask_indirect = {
        'h': 2.0e-9,
        's': 1.5e-9,
        'w': 1.2e-9,
        'c': 0.8e-9,
    }

    # With surgical masks (50% efficacy, 70% adoption)
    # Effective reduction = 0.5 × 0.7 = 35%
    mask_reduction = 1.0 - (0.5 * 0.7)

    with_mask_direct = {k: v * mask_reduction for k, v in no_mask_direct.items()}
    with_mask_indirect = {k: v * mask_reduction for k, v in no_mask_indirect.items()}

    # Common parameters
    base_pars = {
        'pop_size': 5000,
        'pop_type': 'hybrid',
        'n_days': 90,
        'use_dose_response': True,
        'N0': 100,
        'infectivity_scale': 1.0,
        'viral_load_function': {'type': 'default', 'peak_load': 1e8},
    }

    # Scenario 1: No masks
    print("\nRunning scenario 1: No masks...")
    pars1 = base_pars.copy()
    pars1['deposition_rates_direct'] = no_mask_direct
    pars1['deposition_rates_indirect'] = no_mask_indirect
    sim1 = cv.Sim(pars1)
    sim1.run()

    # Scenario 2: With masks
    print("Running scenario 2: Surgical masks (50% eff, 70% adoption)...")
    pars2 = base_pars.copy()
    pars2['deposition_rates_direct'] = with_mask_direct
    pars2['deposition_rates_indirect'] = with_mask_indirect
    sim2 = cv.Sim(pars2)
    sim2.run()

    # Results
    no_mask_infections = sim1.results['cum_infections'][-1]
    with_mask_infections = sim2.results['cum_infections'][-1]
    reduction = (no_mask_infections - with_mask_infections) / no_mask_infections * 100

    print(f"\nRESULTS:")
    print(f"  No masks:              {no_mask_infections:6,.0f} infections")
    print(f"  With surgical masks:   {with_mask_infections:6,.0f} infections")
    print(f"  Reduction:             {reduction:6.1f}%")
    print("="*70)

    return {'no_masks': sim1, 'with_masks': sim2}


def ventilation_example():
    '''
    Example: Compare different ventilation scenarios
    '''
    print("\n" + "="*70)
    print("VENTILATION EXAMPLE: Impact of Improved Air Circulation")
    print("="*70)

    # Baseline: Poor ventilation (high indirect transmission)
    poor_vent_indirect = {
        'h': 3.0e-9,
        's': 2.5e-9,
        'w': 2.0e-9,
        'c': 1.0e-9,
    }

    # Good ventilation: 50% reduction in indirect transmission
    good_vent_indirect = {k: v * 0.5 for k, v in poor_vent_indirect.items()}

    # Excellent ventilation: 70% reduction in indirect transmission
    excellent_vent_indirect = {k: v * 0.3 for k, v in poor_vent_indirect.items()}

    # Direct transmission stays the same
    direct_rates = {
        'h': 10.0e-9,
        's': 8.0e-9,
        'w': 6.0e-9,
        'c': 4.0e-9,
    }

    base_pars = {
        'pop_size': 5000,
        'pop_type': 'hybrid',
        'n_days': 90,
        'use_dose_response': True,
        'N0': 100,
        'infectivity_scale': 1.0,
        'viral_load_function': {'type': 'default', 'peak_load': 1e8},
        'deposition_rates_direct': direct_rates,
    }

    scenarios = {
        'Poor ventilation': poor_vent_indirect,
        'Good ventilation': good_vent_indirect,
        'Excellent ventilation': excellent_vent_indirect,
    }

    results = {}
    for name, indirect_rates in scenarios.items():
        print(f"\nRunning: {name}...")
        pars = base_pars.copy()
        pars['deposition_rates_indirect'] = indirect_rates
        sim = cv.Sim(pars)
        sim.run()
        results[name] = sim.results['cum_infections'][-1]

    print(f"\nRESULTS:")
    baseline = results['Poor ventilation']
    for name, infections in results.items():
        reduction = (baseline - infections) / baseline * 100 if name != 'Poor ventilation' else 0
        print(f"  {name:25s}: {infections:6,.0f} infections ({reduction:+5.1f}%)")
    print("="*70)

    return results


def combined_interventions_example():
    '''
    Example: Combine masks + ventilation
    '''
    print("\n" + "="*70)
    print("COMBINED INTERVENTIONS: Masks + Ventilation")
    print("="*70)

    # Baseline deposition rates
    baseline_direct = {'h': 10.0e-9, 's': 8.0e-9, 'w': 6.0e-9, 'c': 4.0e-9}
    baseline_indirect = {'h': 2.0e-9, 's': 1.5e-9, 'w': 1.2e-9, 'c': 0.8e-9}

    scenarios = {
        'No intervention': {
            'direct': baseline_direct,
            'indirect': baseline_indirect,
        },
        'Masks only (60% eff, 70% adopt)': {
            'direct': {k: v * (1 - 0.6*0.7) for k, v in baseline_direct.items()},
            'indirect': {k: v * (1 - 0.6*0.7) for k, v in baseline_indirect.items()},
        },
        'Ventilation only (50% improvement)': {
            'direct': baseline_direct,
            'indirect': {k: v * 0.5 for k, v in baseline_indirect.items()},
        },
        'Masks + Ventilation': {
            'direct': {k: v * (1 - 0.6*0.7) for k, v in baseline_direct.items()},
            'indirect': {k: v * 0.5 * (1 - 0.6*0.7) for k, v in baseline_indirect.items()},
        },
    }

    base_pars = {
        'pop_size': 5000,
        'pop_type': 'hybrid',
        'n_days': 90,
        'use_dose_response': True,
        'N0': 100,
        'infectivity_scale': 1.0,
        'viral_load_function': {'type': 'default', 'peak_load': 1e8},
    }

    results = {}
    for name, rates in scenarios.items():
        print(f"\nRunning: {name}...")
        pars = base_pars.copy()
        pars['deposition_rates_direct'] = rates['direct']
        pars['deposition_rates_indirect'] = rates['indirect']
        sim = cv.Sim(pars)
        sim.run()
        results[name] = sim.results['cum_infections'][-1]

    print(f"\nRESULTS:")
    baseline = results['No intervention']
    for name, infections in results.items():
        reduction = (baseline - infections) / baseline * 100 if name != 'No intervention' else 0
        print(f"  {name:40s}: {infections:6,.0f} ({reduction:+5.1f}%)")
    print("="*70)

    return results


def custom_parameters_template():
    '''
    Template showing how users can customize all parameters
    '''
    print("\n" + "="*70)
    print("TEMPLATE: Customize Your Own Scenario")
    print("="*70)
    print("""
This template shows all the parameters you can adjust:

```python
# 1. DEPOSITION RATES (viral particles deposited per minute)
custom_direct_rates = {
    'h': 10.0e-9,   # Household: face-to-face contact
    's': 8.0e-9,    # School
    'w': 6.0e-9,    # Work
    'c': 4.0e-9,    # Community
}

custom_indirect_rates = {
    'h': 2.0e-9,    # Household: shared air
    's': 1.5e-9,    # School
    'w': 1.2e-9,    # Work
    'c': 0.8e-9,    # Community
}

# 2. MASK EFFECTS (modify deposition rates)
mask_efficacy = 0.5      # 50% reduction (cloth=0.3, surgical=0.5, N95=0.8)
mask_adoption = 0.7      # 70% of people wear masks
mask_factor = 1.0 - (mask_efficacy * mask_adoption)

# Apply masks to specific layers (e.g., not at home)
masked_direct = custom_direct_rates.copy()
masked_direct['s'] *= mask_factor  # Masks at school
masked_direct['w'] *= mask_factor  # Masks at work
masked_direct['c'] *= mask_factor  # Masks in community

# 3. VENTILATION (affects indirect transmission only)
ventilation_improvement = 0.6  # 40% reduction in airborne transmission

ventilated_indirect = {
    k: v * ventilation_improvement for k, v in custom_indirect_rates.items()
}

# 4. RUN SIMULATION
pars = {
    'pop_size': 5000,
    'pop_type': 'hybrid',
    'n_days': 90,
    'use_dose_response': True,
    'N0': 100,                           # Infectious dose parameter
    'infectivity_scale': 1.0,            # Overall infectivity
    'deposition_rates_direct': masked_direct,
    'deposition_rates_indirect': ventilated_indirect,
    'viral_load_function': {
        'type': 'default',
        'peak_load': 1e8,
        'peak_day': 5.0
    }
}

sim = cv.Sim(pars)
sim.run()
```

Try modifying these values to explore different scenarios!
    """)
    print("="*70)


if __name__ == '__main__':
    print("\n" + "="*70)
    print("DOSE-RESPONSE MODEL: MASK AND DEPOSITION RATE EXAMPLES")
    print("="*70)
    print("\nThis demo shows how to use COVASIM's dose-response model to")
    print("experiment with mask usage and ventilation interventions.\n")

    # Run examples
    simple_mask_example()
    ventilation_example()
    combined_interventions_example()
    custom_parameters_template()

    print("\n" + "="*70)
    print("DEMO COMPLETE!")
    print("="*70)
    print("\nTo run the full interactive examples with plots, see:")
    print("  - examples/dose_response_basic.py")
    print("    Functions: mask_and_deposition_example()")
    print("               intervention_scenario_example()")
    print("="*70 + "\n")