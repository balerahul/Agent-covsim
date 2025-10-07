'''
Dose-response model for COVASIM infection risk calculation.
Implements P_dose = 1 - exp(-N/N₀ · I) formulation based on inhaled viral particles.
'''

import numpy as np
import sciris as sc
from . import utils as cvu
from . import defaults as cvd


class DoseResponseModel(sc.prettyobj):
    '''
    Dose-response model for computing infection probability based on inhaled viral particles.
    
    This model replaces the beta-based transmission calculation with a dose-response function:
    P_dose = 1 - exp(-N/N₀ · I)
    
    where:
        N = N_direct + N_indirect (total inhaled viral particles)
        N_direct = Ṅᵛ_direct × τ_direct × λ(t)
        N_indirect = Ṅᵛ_indirect × τ_indirect × λ(t)
        N₀ = infectious dose scaling parameter
        I = infectivity scaling factor
    
    Args:
        N0 (float): Infectious dose scaling parameter (viral particles)
        infectivity_scale (float): Additional scaling factor (I) for infection probability
        deposition_rates_direct (dict): Ṅᵛ_direct - Rate of viral particle deposition for direct contact per layer/interaction type
        deposition_rates_indirect (dict): Ṅᵛ_indirect - Rate of viral particle deposition for indirect contact per layer
    '''
    
    def __init__(self, N0=100, infectivity_scale=1.0, 
                 deposition_rates_direct=None, deposition_rates_indirect=None):
        self.N0 = N0
        self.infectivity_scale = infectivity_scale
        
        # Initialize deposition rates - these should be provided by the user
        # Default values are placeholders and should be replaced with actual data
        self.deposition_rates_direct = sc.mergedicts(
            self._default_deposition_rates_direct(), 
            deposition_rates_direct
        )
        self.deposition_rates_indirect = sc.mergedicts(
            self._default_deposition_rates_indirect(), 
            deposition_rates_indirect
        )

        self.initialized = False
        return
    
    def _default_deposition_rates_direct(self):
        '''
        Default viral particle deposition rates for direct (face-to-face) contact.
        These are placeholder values and should be replaced with user-provided data.
        Units: particles/minute (rate of deposition on upper respiratory tract)
        '''
        return {
            'a': 6.0e-9,   # All contacts (random population) - average value
            'h': 5.0e-9,   # Household - close proximity, no masks typically
            's': 5.0e-9,   # School - moderate proximity
            'w': 5.0e-9,   # Work - moderate proximity  
            'c': 5.0e-9,   # Community - variable proximity
            'l': 8.0e-9,   # Long-term care facilities - close proximity care
        }
    
    def _default_deposition_rates_indirect(self):
        '''
        Default viral particle deposition rates for indirect (shared air) contact.
        These are placeholder values and should be replaced with CFD-derived data.
        Units: particles/minute (rate of deposition on upper respiratory tract)
        '''
        return {
            'a': 0.5e-9,   # All contacts (random population) - average value
            'h': 0.5e-9,   # Household - smaller space, less ventilation
            's': 0.1e-9,   # School - larger space, moderate ventilation
            'w': 0.1e-9,   # Work - office spaces
            'c': 0.1e-9,   # Community - highly variable, often well-ventilated
            'l': 1.e-9,   # Long-term care facilities - confined spaces
        }
    
    def compute_direct_dose(self, deposition_rate, duration, viral_load):
        '''
        Calculate N_direct from face-to-face contact.
        N_direct = Ṅᵛ_direct × τ_direct × λ(t)
        
        Args:
            deposition_rate (float or array): Ṅᵛ_direct - Rate of viral particle deposition
            duration (float or array): τ_direct - Contact duration in minutes
            viral_load (float or array): λ(t) - Viral load at time t
        
        Returns:
            N_direct (float or array): Direct dose of viral particles
        '''
        N_direct = deposition_rate * duration * viral_load
        return N_direct
    
    def compute_indirect_dose(self, deposition_rate, duration, viral_load):
        '''
        Calculate N_indirect from shared air exposure.
        N_indirect = Ṅᵛ_indirect × τ_indirect × λ(t)
        
        Args:
            deposition_rate (float or array): Ṅᵛ_indirect - Rate of viral particle deposition
            duration (float or array): τ_indirect - Exposure duration in minutes
            viral_load (float or array): λ(t) - Viral load at time t
        
        Returns:
            N_indirect (float or array): Indirect dose from shared air
        '''
        N_indirect = deposition_rate * duration * viral_load
        return N_indirect
    
    def compute_infection_probability(self, N_total):
        '''
        Calculate infection probability using dose-response function.
        P_dose = 1 - exp(-N/N₀ · I)
        
        Args:
            N_total (float or array): Total inhaled viral particles
        
        Returns:
            P_dose (float or array): Infection probability
        '''
        # Avoid numerical issues with very small or negative N
        N_total = np.maximum(N_total, 0)
        
        # P_dose = 1 - exp(-N/N₀ · I)
        exponent = -(N_total / self.N0) * self.infectivity_scale
        P_dose = 1.0 - np.exp(exponent)
        
        # Ensure probability is in [0, 1]
        P_dose = np.clip(P_dose, 0.0, 1.0)
        
        return P_dose
    
    def compute_dose_response_transmission(self, sources, targets, viral_loads, 
                                          direct_durations, indirect_durations,
                                          layer='c', rel_trans=None, rel_sus=None):
        '''
        Main method to compute transmission using dose-response model.
        
        Args:
            sources (array): Source individual indices
            targets (array): Target individual indices  
            viral_loads (array): Viral loads λ(t) for all individuals
            direct_durations (array): τ_direct - Direct contact durations for each source-target pair
            indirect_durations (array): τ_indirect - Indirect exposure durations
            layer (str): Contact layer type (h, s, w, c, l)
            rel_trans (array): Relative transmissibility modifiers
            rel_sus (array): Relative susceptibility modifiers
        
        Returns:
            trans_probs (array): Transmission probabilities for each contact
        '''
        n_contacts = len(sources)
        
        # Get deposition rates for this layer
        deposition_rate_direct = self.deposition_rates_direct.get(layer, 5.0)
        deposition_rate_indirect = self.deposition_rates_indirect.get(layer, 1.0)
        
        # Get viral loads for source individuals
        source_viral_loads = viral_loads[sources]
        
        # Apply relative transmissibility if provided
        if rel_trans is not None:
            source_viral_loads = source_viral_loads * rel_trans[sources]
        # Compute direct dose: N_direct = Ṅᵛ_direct × τ_direct × λ(t)
        N_direct = self.compute_direct_dose(
            deposition_rate_direct,
            direct_durations, 
            source_viral_loads
        )
        
        # Compute indirect dose: N_indirect = Ṅᵛ_indirect × τ_indirect × λ(t)
        N_indirect = self.compute_indirect_dose(
            deposition_rate_indirect,
            indirect_durations,
            source_viral_loads
        )
        
        # Total dose: N = N_direct + N_indirect
        N_total = N_direct + N_indirect
        
        # Compute base infection probability: P_dose = 1 - exp(-N/N₀ · I)
        P_dose = self.compute_infection_probability(N_total)
        
        # Apply relative susceptibility if provided
        if rel_sus is not None:
            P_dose = P_dose * rel_sus[targets]
        
        return P_dose
    
    def integrate_viral_load(self, viral_load_func, t_start, t_end, n_points=10):
        '''
        Integrate viral load over exposure duration for more accurate dose calculation.
        This is used when viral load changes significantly during the contact period.
        
        Args:
            viral_load_func (callable): Function that returns viral load given time since infection
            t_start (float): Start time of exposure (days since infection)
            t_end (float): End time of exposure (days since infection)
            n_points (int): Number of integration points
        
        Returns:
            avg_viral_load (float): Average viral load over the exposure period
        '''
        if t_start >= t_end:
            return 0.0
        
        # Simple trapezoidal integration
        t_points = np.linspace(t_start, t_end, n_points)
        viral_loads = np.array([viral_load_func(t) for t in t_points])
        avg_viral_load = np.trapz(viral_loads, t_points) / (t_end - t_start)
        
        return avg_viral_load