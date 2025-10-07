'''
Configuration templates for dose-response models in COVASIM.

This file provides ready-to-use configuration templates for different
scenarios and environments when using the dose-response transmission model.
'''

import numpy as np
import covasim as cv


# =============================================================================
# BASIC CONFIGURATION TEMPLATES
# =============================================================================

def get_household_config():
    '''Configuration optimized for household transmission studies'''
    return {
        'N0': 80,  # Lower threshold for close household contact
        'infectivity_scale': 1.3,  # Higher infectivity in close quarters
        'deposition_rates_direct': {
            'h': 18.0,  # High face-to-face deposition in households
            's': 8.0,   # Standard school rates
            'w': 6.0,   # Standard workplace rates
            'c': 4.0,   # Standard community rates
        },
        'deposition_rates_indirect': {
            'h': 4.0,   # Poor ventilation in households
            's': 2.0,   # Moderate ventilation in schools
            'w': 1.5,   # Good ventilation in workplaces
            'c': 0.8,   # Variable community ventilation
        },
        'viral_load_function': {
            'type': 'default',
            'peak_load': 8e7,  # Slightly lower peak for household setting
            'peak_day': 4.5,   # Earlier peak
            'decay_rate': 0.4  # Slower decay
        }
    }


def get_school_config():
    '''Configuration for school environment transmission'''
    return {
        'N0': 120,  # Moderate threshold
        'infectivity_scale': 1.1,  # Moderate infectivity
        'deposition_rates_direct': {
            'h': 15.0,  # Household rates
            's': 10.0,  # Enhanced school rates (close classroom contact)
            'w': 6.0,   # Workplace rates
            'c': 4.0,   # Community rates
        },
        'deposition_rates_indirect': {
            'h': 3.0,   # Household rates
            's': 2.5,   # Moderate school ventilation
            'w': 1.5,   # Workplace rates
            'c': 0.8,   # Community rates
        },
        'viral_load_function': {
            'type': 'default',
            'peak_load': 1e8,
            'peak_day': 5.0,
            'decay_rate': 0.5
        }
    }


def get_workplace_config():
    '''Configuration for workplace transmission (offices, etc.)'''
    return {
        'N0': 150,  # Higher threshold due to better ventilation/distancing
        'infectivity_scale': 0.9,  # Lower infectivity
        'deposition_rates_direct': {
            'h': 15.0,  # Household rates
            's': 8.0,   # School rates
            'w': 7.0,   # Moderate workplace rates
            'c': 4.0,   # Community rates
        },
        'deposition_rates_indirect': {
            'h': 3.0,   # Household rates
            's': 2.0,   # School rates
            'w': 1.2,   # Good workplace ventilation
            'c': 0.8,   # Community rates
        },
        'viral_load_function': {
            'type': 'default',
            'peak_load': 1e8,
            'peak_day': 5.0,
            'decay_rate': 0.5
        }
    }


def get_healthcare_config():
    '''Configuration for healthcare settings with high-risk contacts'''
    return {
        'N0': 60,   # Low threshold due to high viral exposure risk
        'infectivity_scale': 1.5,  # High infectivity in healthcare
        'deposition_rates_direct': {
            'h': 15.0,  # Household rates
            's': 8.0,   # School rates
            'w': 12.0,  # High workplace rates for healthcare
            'c': 4.0,   # Community rates
            'l': 20.0,  # Very high long-term care rates
        },
        'deposition_rates_indirect': {
            'h': 3.0,   # Household rates
            's': 2.0,   # School rates
            'w': 2.0,   # Higher workplace rates
            'c': 0.8,   # Community rates
            'l': 3.5,   # High long-term care rates
        },
        'viral_load_function': {
            'type': 'default',
            'peak_load': 1.2e8,  # Higher peak viral load
            'peak_day': 4.0,     # Earlier peak
            'decay_rate': 0.4    # Slower decay
        }
    }


# =============================================================================
# VARIANT-SPECIFIC CONFIGURATIONS
# =============================================================================

def get_alpha_variant_config():
    '''Configuration for Alpha variant (higher transmissibility)'''
    return {
        'N0': 100,
        'infectivity_scale': 1.4,  # ~40% more transmissible
        'deposition_rates_direct': {
            'h': 15.0, 's': 8.0, 'w': 6.0, 'c': 4.0
        },
        'deposition_rates_indirect': {
            'h': 3.0, 's': 2.0, 'w': 1.5, 'c': 0.8
        },
        'viral_load_function': {
            'type': 'default',
            'peak_load': 1.3e8,  # Higher peak viral load
            'peak_day': 4.5,
            'decay_rate': 0.45
        }
    }


def get_delta_variant_config():
    '''Configuration for Delta variant (much higher transmissibility)'''
    return {
        'N0': 80,   # Lower threshold
        'infectivity_scale': 2.0,  # ~100% more transmissible
        'deposition_rates_direct': {
            'h': 18.0, 's': 10.0, 'w': 8.0, 'c': 5.0
        },
        'deposition_rates_indirect': {
            'h': 4.0, 's': 2.5, 'w': 2.0, 'c': 1.0
        },
        'viral_load_function': {
            'type': 'default',
            'peak_load': 2e8,    # Much higher peak
            'peak_day': 4.0,     # Earlier peak
            'decay_rate': 0.4    # Slower decay
        }
    }


def get_omicron_variant_config():
    '''Configuration for Omicron variant (high transmissibility, different kinetics)'''
    return {
        'N0': 70,   # Low threshold due to high transmissibility
        'infectivity_scale': 2.2,  # Very high transmissibility
        'deposition_rates_direct': {
            'h': 20.0, 's': 12.0, 'w': 9.0, 'c': 6.0
        },
        'deposition_rates_indirect': {
            'h': 4.5, 's': 3.0, 'w': 2.2, 'c': 1.2
        },
        'viral_load_function': {
            'type': 'gaussian',  # Different kinetics
            'peak_load': 1.5e8,
            'peak_day': 3.5,     # Much earlier peak
            'sigma': 2.0         # Narrower distribution
        }
    }


# =============================================================================
# INTERVENTION-SPECIFIC CONFIGURATIONS
# =============================================================================

def get_masked_environment_config():
    '''Configuration for environments with widespread mask use'''
    base_config = get_workplace_config()
    
    # Reduce deposition rates due to masks
    mask_reduction = 0.3  # 70% reduction
    
    for layer in base_config['deposition_rates_direct']:
        base_config['deposition_rates_direct'][layer] *= mask_reduction
    
    for layer in base_config['deposition_rates_indirect']:
        base_config['deposition_rates_indirect'][layer] *= mask_reduction
    
    return base_config


def get_improved_ventilation_config():
    '''Configuration for environments with improved ventilation'''
    base_config = get_workplace_config()
    
    # Significantly reduce indirect transmission
    ventilation_reduction = 0.4  # 60% reduction in indirect transmission
    
    for layer in base_config['deposition_rates_indirect']:
        base_config['deposition_rates_indirect'][layer] *= ventilation_reduction
    
    return base_config


def get_social_distancing_config():
    '''Configuration for environments with social distancing measures'''
    base_config = get_workplace_config()
    
    # Reduce direct transmission due to increased distance
    distancing_reduction = 0.5  # 50% reduction in direct transmission
    
    for layer in base_config['deposition_rates_direct']:
        base_config['deposition_rates_direct'][layer] *= distancing_reduction
    
    return base_config


# =============================================================================
# ADVANCED VIRAL LOAD CONFIGURATIONS
# =============================================================================

def get_hcd_model_config():
    '''Configuration using HCD (Human Challenge Data) viral load model'''
    return {
        'N0': 100,
        'infectivity_scale': 1.0,
        'deposition_rates_direct': {
            'h': 15.0, 's': 8.0, 'w': 6.0, 'c': 4.0
        },
        'deposition_rates_indirect': {
            'h': 3.0, 's': 2.0, 'w': 1.5, 'c': 0.8
        },
        'viral_load_function': {
            'type': 'hcd',
            'sample_params': True,  # Sample from parameter distribution
            'seed': None  # Random seed for parameter sampling
        }
    }


def get_custom_viral_load_config():
    '''Configuration with custom viral load function for special scenarios'''
    
    def prolonged_shedding_viral_load(t):
        '''Viral load with prolonged shedding (e.g., immunocompromised)'''
        # Initial peak similar to normal
        main_peak = 1e8 * np.exp(-0.5 * ((t - 5) / 2.5)**2)
        
        # Prolonged low-level shedding
        prolonged = 1e6 * np.exp(-0.1 * np.maximum(t - 10, 0))
        
        return np.maximum(main_peak + prolonged, 0)
    
    return {
        'N0': 100,
        'infectivity_scale': 1.0,
        'deposition_rates_direct': {
            'h': 15.0, 's': 8.0, 'w': 6.0, 'c': 4.0
        },
        'deposition_rates_indirect': {
            'h': 3.0, 's': 2.0, 'w': 1.5, 'c': 0.8
        },
        'viral_load_function': prolonged_shedding_viral_load
    }


# =============================================================================
# COMPLETE SIMULATION CONFIGURATIONS
# =============================================================================

def create_household_study_sim():
    '''Create a complete simulation for household transmission study'''
    config = get_household_config()
    
    pars = {
        'pop_size': 2000,
        'n_days': 90,
        'use_dose_response': True,
        **config
    }
    
    return cv.Sim(pars)


def create_school_outbreak_sim():
    '''Create a simulation for school outbreak analysis'''
    config = get_school_config()
    
    pars = {
        'pop_size': 5000,
        'n_days': 60,
        'use_dose_response': True,
        **config
    }
    
    return cv.Sim(pars)


def create_variant_comparison_sims():
    '''Create simulations for comparing different variants'''
    variants = {
        'Original': get_workplace_config(),
        'Alpha': get_alpha_variant_config(),
        'Delta': get_delta_variant_config(),
        'Omicron': get_omicron_variant_config()
    }
    
    sims = {}
    for variant_name, config in variants.items():
        pars = {
            'pop_size': 3000,
            'n_days': 80,
            'use_dose_response': True,
            'label': variant_name,
            **config
        }
        sims[variant_name] = cv.Sim(pars)
    
    return sims


def create_intervention_comparison_sims():
    '''Create simulations for comparing different interventions'''
    interventions = {
        'Baseline': get_workplace_config(),
        'Masks': get_masked_environment_config(),
        'Ventilation': get_improved_ventilation_config(),
        'Distancing': get_social_distancing_config()
    }
    
    sims = {}
    for intervention_name, config in interventions.items():
        pars = {
            'pop_size': 3000,
            'n_days': 80,
            'use_dose_response': True,
            'label': intervention_name,
            **config
        }
        sims[intervention_name] = cv.Sim(pars)
    
    return sims


# =============================================================================
# CONFIGURATION UTILITIES
# =============================================================================

def validate_dose_response_config(config):
    '''Validate a dose-response configuration'''
    required_keys = ['N0', 'infectivity_scale']
    
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required key: {key}")
    
    # Validate ranges
    if config['N0'] <= 0:
        raise ValueError("N0 must be positive")
    
    if config['infectivity_scale'] <= 0:
        raise ValueError("infectivity_scale must be positive")
    
    print("Configuration validation passed")
    return True


def print_config_summary(config):
    '''Print a summary of a dose-response configuration'''
    print("=== Dose-Response Configuration Summary ===")
    
    print(f"N0 (infectious dose): {config.get('N0', 'Not set')}")
    print(f"Infectivity scale: {config.get('infectivity_scale', 'Not set')}")
    
    if 'deposition_rates_direct' in config:
        print("\nDirect deposition rates:")
        for layer, rate in config['deposition_rates_direct'].items():
            print(f"  {layer}: {rate:.1f}")
    
    if 'deposition_rates_indirect' in config:
        print("\nIndirect deposition rates:")
        for layer, rate in config['deposition_rates_indirect'].items():
            print(f"  {layer}: {rate:.1f}")
    
    vl_func = config.get('viral_load_function')
    if vl_func:
        if isinstance(vl_func, dict):
            print(f"\nViral load model: {vl_func.get('type', 'unknown')}")
            if 'peak_load' in vl_func:
                print(f"Peak load: {vl_func['peak_load']:.1e}")
            if 'peak_day' in vl_func:
                print(f"Peak day: {vl_func['peak_day']}")
        else:
            print(f"\nViral load model: custom function")


if __name__ == '__main__':
    # Example usage
    print("=== Dose-Response Configuration Templates ===\n")
    
    # Test basic configurations
    configs = {
        'Household': get_household_config(),
        'School': get_school_config(),
        'Workplace': get_workplace_config(),
        'Healthcare': get_healthcare_config()
    }
    
    for name, config in configs.items():
        print(f"\n{name} Configuration:")
        validate_dose_response_config(config)
        print_config_summary(config)
        print("-" * 50)
    
    # Test variant configurations
    print("\n=== Variant Configurations ===")
    variant_configs = {
        'Alpha': get_alpha_variant_config(),
        'Delta': get_delta_variant_config(),
        'Omicron': get_omicron_variant_config()
    }
    
    for name, config in variant_configs.items():
        print(f"\n{name} Variant:")
        print_config_summary(config)
    
    print("\n=== Configuration templates ready for use ===")