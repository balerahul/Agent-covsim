'''
Contact duration sampling module for COVASIM dose-response model.
Implements duration sampling from distributions with age-based modifiers.
Based on the interaction-time.py implementation.
'''

import numpy as np
import sciris as sc
from . import defaults as cvd


class ContactDurationSampler(sc.prettyobj):
    '''
    Sample contact durations based on layer-specific distributions and age interactions.
    
    This class implements the exact duration sampling logic from interaction-time.py,
    including mixture lognormal distributions for different layers and age-based modifiers.
    
    Args:
        duration_config (dict): Configuration for duration distributions per layer
        seed (int): Random seed for reproducibility
    '''
    
    def __init__(self, duration_config=None, seed=None):
        self.duration_config = sc.mergedicts(self._default_duration_config(), duration_config)
        self.seed = seed
        if seed is not None:
            self.rng = np.random.default_rng(seed)
        else:
            self.rng = np.random.default_rng()
        self.initialized = True
        return
    
    def _default_duration_config(self):
        '''
        Default duration distribution parameters from interaction-time.py.
        '''
        return {
            'h': {  # Home
                'type': 'mixture_lognormal',
                'weights': [0.70, 0.30],
                'medians': [60, 180],  # minutes
                'sigmas': [0.9, 0.6],
                'cap': 360
            },
            's': {  # School
                'type': 'mixture_lognormal',
                'weights': [0.55, 0.45],
                'medians': [50, 2],  # minutes
                'sigmas': [0.25, 0.6],
                'cap': 60
            },
            'w': {  # Work
                'type': 'mixture_lognormal',
                'weights': [0.65, 0.25, 0.10],
                'medians': [1.5, 7, 45],  # minutes
                'sigmas': [0.7, 0.6, 0.4],
                'cap': 120
            },
            'c': {  # Community
                'type': 'weibull',
                'k': 0.6,
                'lambda': 5.0,
                'cap': 60
            },
            'l': {  # Long-term care facilities (default fallback)
                'type': 'mixture_lognormal',
                'weights': [0.8, 0.2],
                'medians': [2, 30],
                'sigmas': [0.6, 0.7],
                'cap': 120
            }
        }
    
    def lognormal_from_median_sigma(self, median, sigma, size=None):
        '''
        Generate lognormal samples from median and sigma parameters.
        
        Args:
            median (float): Median of the lognormal distribution
            sigma (float): Sigma parameter of the underlying normal distribution
            size (int or tuple): Output shape
        
        Returns:
            samples (float or array): Lognormal samples
        '''
        mu = np.log(median)
        return self.rng.lognormal(mean=mu, sigma=sigma, size=size)
    
    def mixture_lognormal(self, weights, medians, sigmas, cap=None):
        '''
        Sample from a mixture of lognormal distributions.
        
        Args:
            weights (list): Mixture weights (must sum to 1)
            medians (list): Medians for each component
            sigmas (list): Sigmas for each component
            cap (float): Maximum value cap
        
        Returns:
            sample (float): Duration sample
        '''
        # Ensure weights sum to 1
        weights = np.array(weights)
        weights = weights / weights.sum()
        
        # Choose component
        idx = self.rng.choice(len(weights), p=weights)
        
        # Sample from chosen component
        t = self.lognormal_from_median_sigma(medians[idx], sigmas[idx])
        
        # Apply cap if specified
        if cap is not None:
            t = min(t, cap)
        
        return t
    
    def weibull_sample(self, k, lam, cap=None):
        '''
        Sample from Weibull distribution using inverse transform.
        
        Args:
            k (float): Shape parameter
            lam (float): Scale parameter
            cap (float): Maximum value cap
        
        Returns:
            sample (float): Duration sample
        '''
        u = self.rng.uniform()
        t = lam * (-np.log(1 - u)) ** (1.0 / k)
        
        if cap is not None:
            t = min(t, cap)
        
        return t
    
    def age_interaction_modifier(self, age_i, age_j, flat_band=5.0, sigma_decay=4.0, hard_cut=20.0):
        '''
        Calculate age-based interaction modifier for duration.
        Returns a multiplier in [0,1] based on age difference.
        
        Args:
            age_i (float): Age of person i
            age_j (float): Age of person j
            flat_band (float): Age difference within which modifier = 1.0
            sigma_decay (float): Gaussian decay parameter
            hard_cut (float): Age difference beyond which modifier = 0.0
        
        Returns:
            modifier (float): Duration modifier in [0, 1]
        '''
        d = abs(float(age_i) - float(age_j))
        
        if d <= flat_band:
            return 1.0
        if d >= hard_cut:
            return 0.0
        
        x = d - flat_band
        return float(np.exp(-0.5 * (x / sigma_decay) ** 2))
    
    def home_age_modifier(self, age_i, age_j,
                         flat_band=5.0, sigma_decay=4.0, hard_cut=20.0,
                         cross_peak_adult=(25, 50), cross_peak_child=(0, 20),
                         cross_sigma=5.0, cross_boost=1.0):
        '''
        Special age modifier for HOME contacts with cross-generation boost.
        
        Args:
            age_i, age_j (float): Ages of individuals
            flat_band (float): Same-age flat band
            sigma_decay (float): Same-age Gaussian decay
            hard_cut (float): Same-age hard cutoff
            cross_peak_adult (tuple): Adult age range for cross-generation boost
            cross_peak_child (tuple): Child age range for cross-generation boost
            cross_sigma (float): Cross-generation Gaussian width
            cross_boost (float): Cross-generation boost factor
        
        Returns:
            modifier (float): Duration modifier
        '''
        # Same-age kernel
        d = abs(float(age_i) - float(age_j))
        if d <= flat_band:
            base_w = 1.0
        elif d >= hard_cut:
            base_w = 0.0
        else:
            x = d - flat_band
            base_w = float(np.exp(-0.5 * (x / sigma_decay) ** 2))
        
        # Cross-generation boost
        def in_band(age, band):
            return band[0] <= age <= band[1]
        
        cross_w = 0.0
        if (in_band(age_i, cross_peak_adult) and in_band(age_j, cross_peak_child)) or \
           (in_band(age_j, cross_peak_adult) and in_band(age_i, cross_peak_child)):
            # Taper boost with distance from center of bands
            mid_adult = 0.5 * (cross_peak_adult[0] + cross_peak_adult[1])
            mid_child = 0.5 * (cross_peak_child[0] + cross_peak_child[1])
            dist_adult = min(abs(age_i - mid_adult), abs(age_j - mid_adult))
            dist_child = min(abs(age_i - mid_child), abs(age_j - mid_child))
            taper = np.exp(-0.5 * ((dist_adult/cross_sigma)**2 + (dist_child/cross_sigma)**2))
            cross_w = cross_boost * taper
        
        # Combine effects: keep the stronger of the two
        return max(base_w, cross_w)
    
    def sample_duration(self, layer, age_i=None, age_j=None, overlap_minutes=np.inf, theta_ab=1.0):
        '''
        Sample contact duration for a given layer and age pair.
        
        Args:
            layer (str): Contact layer (h, s, w, c, l)
            age_i (float): Age of person i (optional)
            age_j (float): Age of person j (optional)
            overlap_minutes (float): Maximum possible overlap time
            theta_ab (float): Dyad-specific scaling factor
        
        Returns:
            duration (float): Contact duration in minutes
        '''
        # Get configuration for this layer
        if layer not in self.duration_config:
            # Use default (l) configuration if layer not found
            config = self.duration_config['l']
        else:
            config = self.duration_config[layer]
        
        # Sample base duration based on distribution type
        if config['type'] == 'mixture_lognormal':
            t = self.mixture_lognormal(
                config['weights'],
                config['medians'],
                config['sigmas'],
                config.get('cap')
            )
        elif config['type'] == 'weibull':
            t = self.weibull_sample(
                config['k'],
                config['lambda'],
                config.get('cap')
            )
        else:
            # Default fallback
            t = 10.0  # Default 10 minutes
        
        # Apply dyad factor
        t *= float(theta_ab)
        
        # Apply age-based modifier if ages provided
        if age_i is not None and age_j is not None:
            if layer == 'h':  # Home layer uses special modifier
                w = self.home_age_modifier(age_i, age_j)
            else:
                w = self.age_interaction_modifier(age_i, age_j)
            t *= w
        
        # Truncate by available overlap
        t = max(0.1, min(t, overlap_minutes))
        
        return t
    
    def sample_durations_vectorized(self, layer, n_contacts, ages_i=None, ages_j=None, 
                                   overlap_minutes=None, theta_ab=None):
        '''
        Sample multiple contact durations efficiently.
        
        Args:
            layer (str): Contact layer
            n_contacts (int): Number of contacts to sample
            ages_i (array): Ages of source individuals (optional)
            ages_j (array): Ages of target individuals (optional)
            overlap_minutes (float or array): Maximum overlap times
            theta_ab (float or array): Dyad-specific scaling factors
        
        Returns:
            durations (array): Contact durations in minutes
        '''
        durations = np.zeros(n_contacts)
        
        # Handle scalar inputs
        if overlap_minutes is None:
            overlap_minutes = np.full(n_contacts, np.inf)
        elif np.isscalar(overlap_minutes):
            overlap_minutes = np.full(n_contacts, overlap_minutes)
        
        if theta_ab is None:
            theta_ab = np.ones(n_contacts)
        elif np.isscalar(theta_ab):
            theta_ab = np.full(n_contacts, theta_ab)
        
        # Sample durations for each contact
        for i in range(n_contacts):
            age_i = ages_i[i] if ages_i is not None else None
            age_j = ages_j[i] if ages_j is not None else None
            overlap = overlap_minutes[i]
            theta = theta_ab[i]
            
            durations[i] = self.sample_duration(layer, age_i, age_j, overlap, theta)
        
        return durations