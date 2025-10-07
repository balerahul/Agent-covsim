'''
Viral dynamics module for COVASIM dose-response model.
Provides pluggable interface for user-supplied viral load functions.
'''

import numpy as np
import sciris as sc
from abc import ABC, abstractmethod
from scipy.integrate import odeint
import json
import os


class ViralLoadFunction(ABC):
    '''
    Abstract base class for viral load functions.
    
    Users should inherit from this class and implement the get_viral_load method
    to provide custom viral load dynamics.
    '''
    
    @abstractmethod
    def get_viral_load(self, t_since_infection):
        '''
        Get viral load at time t since infection.
        
        Args:
            t_since_infection (float or array): Time since infection in days
        
        Returns:
            viral_load (float or array): Viral load (particles/mL or user-defined units)
        '''
        pass
    
    def get_integrated_load(self, t_start, t_end, n_points=10):
        '''
        Get average viral load over a time period.
        
        Args:
            t_start (float): Start time (days since infection)
            t_end (float): End time (days since infection)
            n_points (int): Number of integration points
        
        Returns:
            avg_load (float): Average viral load over the period
        '''
        if t_start >= t_end:
            return 0.0
        
        t_points = np.linspace(t_start, t_end, n_points)
        loads = np.array([self.get_viral_load(t) for t in t_points])
        avg_load = np.trapz(loads, t_points) / (t_end - t_start)
        
        return avg_load


class DefaultViralLoad(ViralLoadFunction):
    '''
    Default viral load implementation for testing.
    
    Uses a simple lognormal-like curve with:
    - Incubation period with low/no viral load
    - Rapid rise to peak
    - Gradual decline
    
    Args:
        incubation_period (float): Days before viral load becomes detectable
        peak_day (float): Day of peak viral load
        peak_load (float): Maximum viral load
        decay_rate (float): Rate of viral load decay after peak
    '''
    
    def __init__(self, incubation_period=2.0, peak_day=5.0, peak_load=1e8, decay_rate=0.5):
        self.incubation_period = incubation_period
        self.peak_day = peak_day
        self.peak_load = peak_load
        self.decay_rate = decay_rate
    
    def get_viral_load(self, t_since_infection):
        '''
        Get viral load using a lognormal-like curve.
        
        Args:
            t_since_infection (float or array): Time since infection in days
        
        Returns:
            viral_load (float or array): Viral load
        '''
        t = np.asarray(t_since_infection)
        viral_load = np.zeros_like(t, dtype=float)
        
        # No viral load during early incubation
        mask_active = t > self.incubation_period
        
        if np.any(mask_active):
            t_active = t[mask_active] - self.incubation_period
            
            # Lognormal-like curve
            sigma = 2.0  # Width of the curve
            mu = np.log(self.peak_day - self.incubation_period)
            
            # Calculate viral load
            log_t = np.log(np.maximum(t_active, 0.1))
            exponent = -0.5 * ((log_t - mu) / sigma) ** 2
            vl = self.peak_load * np.exp(exponent)
            
            # Apply decay after peak
            mask_decay = t[mask_active] > self.peak_day
            if np.any(mask_decay):
                days_after_peak = t[mask_active][mask_decay] - self.peak_day
                vl[mask_decay] *= np.exp(-self.decay_rate * days_after_peak)
            
            viral_load[mask_active] = vl
        
        # Return scalar if input was scalar
        if np.isscalar(t_since_infection):
            return float(viral_load)
        
        return viral_load


class GaussianViralLoad(ViralLoadFunction):
    '''
    Gaussian-shaped viral load curve.
    
    Args:
        peak_day (float): Day of peak viral load
        peak_load (float): Maximum viral load
        sigma (float): Standard deviation of the Gaussian
    '''
    
    def __init__(self, peak_day=5.0, peak_load=1e8, sigma=2.0):
        self.peak_day = peak_day
        self.peak_load = peak_load
        self.sigma = sigma
    
    def get_viral_load(self, t_since_infection):
        '''
        Get viral load using a Gaussian curve.
        
        Args:
            t_since_infection (float or array): Time since infection in days
        
        Returns:
            viral_load (float or array): Viral load
        '''
        t = np.asarray(t_since_infection)
        scalar_input = np.isscalar(t_since_infection)
        
        exponent = -0.5 * ((t - self.peak_day) / self.sigma) ** 2
        viral_load = self.peak_load * np.exp(exponent)
        
        # Zero out negative times
        viral_load = np.where(t >= 0, viral_load, 0)
        
        # Return scalar if input was scalar
        if scalar_input:
            return float(viral_load)
        
        return viral_load


class PiecewiseLinearViralLoad(ViralLoadFunction):
    '''
    Piecewise linear viral load curve defined by time-value pairs.
    
    Args:
        times (list): Time points (days since infection)
        values (list): Viral load values at each time point
    '''
    
    def __init__(self, times=None, values=None):
        if times is None:
            # Default profile
            times = [0, 2, 5, 7, 14, 21]
            values = [0, 1e4, 1e8, 5e7, 1e5, 0]
        
        self.times = np.array(times)
        self.values = np.array(values)
    
    def get_viral_load(self, t_since_infection):
        '''
        Get viral load using linear interpolation.
        
        Args:
            t_since_infection (float or array): Time since infection in days
        
        Returns:
            viral_load (float or array): Viral load
        '''
        viral_load = np.interp(t_since_infection, self.times, self.values)
        return viral_load


class UserSuppliedViralLoad(ViralLoadFunction):
    '''
    Wrapper for user-supplied viral load function.
    
    Args:
        viral_load_func (callable): Function that takes time since infection and returns viral load
    '''
    
    def __init__(self, viral_load_func):
        if not callable(viral_load_func):
            raise ValueError("viral_load_func must be callable")
        self.viral_load_func = viral_load_func
    
    def get_viral_load(self, t_since_infection):
        '''
        Get viral load from user-supplied function.
        
        Args:
            t_since_infection (float or array): Time since infection in days
        
        Returns:
            viral_load (float or array): Viral load
        '''
        return self.viral_load_func(t_since_infection)


class HCDViralLoad(ViralLoadFunction):
    '''
    HCD (Human Challenge Data) viral load model based on within-host ODE dynamics.
    
    Samples parameters from population distribution for each agent to create
    individual viral load trajectories. For pre-computed representative profiles,
    use HCDViralLoadDatabase instead.
    
    Implements the ODE system:
        dU/dt = -beta * U * V
        dI/dt = beta * U * V - gamma * I
        dV/dt = (p / V_m) * I - c * V
    
    Where:
        U: Uninfected target cells
        I: Infected cells
        V: Free virions (viral load)
        beta: Infection rate constant
        gamma: Infected cell clearance rate
        p: Virion production scale
        c: Virion clearance rate
        V_m: Scaling constant
    
    Args:
        params_file (str): Path to JSON file with HCD population parameters
        seed (int): Random seed for parameter sampling
    '''
    
    def __init__(self, params_file=None, seed=None):
        # Default to package data if no file specified
        if params_file is None:
            # Get the path to the data file in the covasim package
            package_dir = os.path.dirname(os.path.abspath(__file__))
            params_file = os.path.join(package_dir, 'data', 'hcd_population_summary.json')
        
        # Load population parameters
        with open(params_file, 'r') as f:
            self.population_data = json.load(f)
        
        # Extract model constants
        constants = self.population_data['model_constants']
        self.U0 = constants['U0']
        self.I0 = constants['I0']
        self.V0 = constants['V0']
        self.c = constants['c']
        self.V_m = constants['V_m']
        
        # Always sample parameters for individual variation
        self._use_mean_parameters()
        # Cache for computed viral loads
        self._cache = {}
        self._max_time = 30.0  # Maximum time to simulate (days)
        self._time_resolution = 0.1  # Time step for ODE integration
        
        # Pre-compute viral load trajectory
        self._precompute_trajectory()
    

    
    def _use_mean_parameters(self):
        '''Use population mean parameters'''
        mu_log = self.population_data['mu_log']
        self.beta = np.exp(mu_log['beta'])
        self.gamma = np.exp(mu_log['gamma'])
        self.p = np.exp(mu_log['p'])
    
    def _hcd_ode(self, y, t):
        '''HCD ODE system'''
        U, I, V = y
        
        dU_dt = -self.beta * U * V
        dI_dt = self.beta * U * V - self.gamma * I
        dV_dt = (self.p / self.V_m) * I - self.c * V
        
        return [dU_dt, dI_dt, dV_dt]
    
    def _precompute_trajectory(self):
        '''Pre-compute viral load trajectory by solving ODE'''
        # Time points to evaluate
        t_points = np.arange(0, self._max_time, self._time_resolution)
        
        # Initial conditions
        y0 = [self.U0, self.I0, self.V0]
        
        # Solve ODE
        solution = odeint(self._hcd_ode, y0, t_points)
        
        # Extract viral load (V)
        self._t_points = t_points
        self._viral_loads = solution[:, 2]  # V is the third component
        
        # Ensure non-negative values
        self._viral_loads = np.maximum(self._viral_loads, 0)
    
    def get_viral_load(self, t_since_infection):
        '''
        Get viral load at time t since infection.
        
        Args:
            t_since_infection (float or array): Time since infection in days
        
        Returns:
            viral_load (float or array): Viral load (virions)
        '''
        # Handle scalar and array inputs
        t = np.asarray(t_since_infection)
        scalar_input = np.isscalar(t_since_infection)
        
        # Interpolate from pre-computed trajectory
        viral_load = np.interp(t, self._t_points, self._viral_loads)
        
        # Zero out negative times
        viral_load = np.where(t >= 0, viral_load, 0)
        
        # Return scalar if input was scalar
        if scalar_input:
            return float(viral_load)
        
        return viral_load
    
    def get_parameters(self):
        '''Return the current parameter values'''
        return {
            'beta': self.beta,
            'gamma': self.gamma,
            'p': self.p,
            'c': self.c,
            'V_m': self.V_m,
            'U0': self.U0,
            'I0': self.I0,
            'V0': self.V0
        }


class HCDViralLoadDatabase(ViralLoadFunction):
    '''
    Database of pre-computed HCD viral load profiles.
    
    This class pre-computes a representative set of viral load trajectories
    and randomly assigns them to infected agents, ensuring:
    - Controlled variation (avoiding extreme outliers)
    - Computational efficiency (pre-computed trajectories)
    - Consistent viral load per agent throughout infection
    
    Args:
        params_file (str): Path to JSON file with HCD population parameters
        n_profiles (int): Number of viral load profiles to pre-compute (default 12)
        seed (int): Random seed for reproducibility
    '''
    
    def __init__(self, params_file=None, n_profiles=12, seed=None):
        # Default to package data if no file specified
        if params_file is None:
            package_dir = os.path.dirname(os.path.abspath(__file__))
            params_file = os.path.join(package_dir, 'data', 'hcd_population_summary.json')
        
        self.n_profiles = n_profiles
        self.params_file = params_file
        self.seed = seed
        
        # Load population parameters
        with open(params_file, 'r') as f:
            self.population_data = json.load(f)
        
        # Pre-compute viral load profiles with consistent generation
        self._create_profile_database()
        
        # Track agent-profile assignments (to be managed externally)
        self.agent_profiles = {}
    
    def _create_profile_database(self):
        '''Create database of pre-computed viral load profiles'''
        
        # Extract model constants
        constants = self.population_data['model_constants']
        
        # Get population parameters
        mu_log = self.population_data['mu_log']
        mu_log_arr = np.array([mu_log['beta'], mu_log['gamma'], mu_log['p']])
        
        # Generate parameter sets
        # Use a stratified approach: include mean + representative samples
        parameter_sets = []
        
        # 1. Always include the population mean
        mean_params = np.exp(mu_log_arr)
        parameter_sets.append(mean_params)
        
        # 2. Generate representative samples
        # Use quantiles to ensure good coverage of the distribution
        n_remaining = self.n_profiles - 1
        
        if n_remaining > 0:
            # Use covariance matrix if available
            Sigma_log = self.population_data.get('Sigma_log_matrix')
            if Sigma_log is not None:
                Sigma_log = np.array(Sigma_log)
                jitter = self.population_data.get('Sigma_log_regularization_jitter', 1e-8)
                Sigma_log += np.eye(3) * jitter
                
                # Generate samples using stratified approach
                # Create a local random generator for consistent profile generation
                rng = np.random.RandomState(42)  # Fixed seed for profile generation
                
                # Generate samples
                for i in range(n_remaining):
                    # Use stratified sampling to ensure good coverage
                    # Scale factor controls distance from mean (0.8 limits extremes)
                    scale_factor = rng.randn(3) * 0.8
                    log_params = mu_log_arr + np.dot(np.linalg.cholesky(Sigma_log), scale_factor)
                    params = np.exp(log_params)
                    parameter_sets.append(params)
            else:
                # Sample independently using marginal standard deviations
                sigma_log = self.population_data['sigma_log_marginal']
                sigma_log_arr = np.array([sigma_log['beta'], sigma_log['gamma'], sigma_log['p']])
                
                # Create a local random generator for consistent profile generation
                rng = np.random.RandomState(42)  # Fixed seed for profile generation
                
                for i in range(n_remaining):
                    scale_factor = rng.randn(3) * 0.8  # Limit to avoid extremes
                    log_params = mu_log_arr + sigma_log_arr * scale_factor
                    params = np.exp(log_params)
                    parameter_sets.append(params)
        
        # Pre-compute viral load trajectories for each parameter set
        self.profiles = []
        self._max_time = 30.0
        self._time_resolution = 0.1
        self._t_points = np.arange(0, self._max_time, self._time_resolution)
        
        for i, (beta, gamma, p) in enumerate(parameter_sets):
            # Solve ODE for this parameter set
            trajectory = self._solve_ode(beta, gamma, p, constants)
            
            profile = {
                'id': i,
                'beta': beta,
                'gamma': gamma,
                'p': p,
                'viral_loads': trajectory,
                'peak_load': np.max(trajectory),
                'peak_day': self._t_points[np.argmax(trajectory)]
            }
            self.profiles.append(profile)
        
        print(f"Created viral load database with {self.n_profiles} profiles")
        print(f"Peak loads range: {min(p['peak_load'] for p in self.profiles):.2e} - {max(p['peak_load'] for p in self.profiles):.2e}")
        print(f"Peak days range: {min(p['peak_day'] for p in self.profiles):.1f} - {max(p['peak_day'] for p in self.profiles):.1f}")
    
    def _solve_ode(self, beta, gamma, p, constants):
        '''Solve HCD ODE for given parameters'''
        
        def hcd_ode(y, t):
            U, I, V = y
            dU_dt = -beta * U * V
            dI_dt = beta * U * V - gamma * I
            dV_dt = (p / constants['V_m']) * I - constants['c'] * V
            return [dU_dt, dI_dt, dV_dt]
        
        # Initial conditions
        y0 = [constants['U0'], constants['I0'], constants['V0']]
        
        # Solve ODE
        solution = odeint(hcd_ode, y0, self._t_points)
        
        # Extract viral load (V) and ensure non-negative
        viral_loads = np.maximum(solution[:, 2], 0)
        
        return viral_loads
    
    def assign_profile(self, agent_id, profile_id=None):
        '''
        Assign a viral load profile to an agent.
        
        Args:
            agent_id: Unique identifier for the agent
            profile_id: Specific profile to assign (if None, randomly selected)
        
        Returns:
            profile_id: The assigned profile ID
        '''
        if profile_id is None:
            # Use the instance seed if provided for consistent assignment
            if self.seed is not None:
                # Create a deterministic profile ID based on agent_id and seed
                rng = np.random.RandomState(self.seed + agent_id)
                profile_id = rng.randint(0, self.n_profiles)
            else:
                # Use global random state
                profile_id = np.random.randint(0, self.n_profiles)
        
        self.agent_profiles[agent_id] = profile_id
        return profile_id
    
    def get_agent_profile(self, agent_id):
        '''
        Get the profile ID assigned to an agent.
        
        Args:
            agent_id: Unique identifier for the agent
        
        Returns:
            profile_id: The profile ID, or None if not assigned
        '''
        return self.agent_profiles.get(agent_id)
    
    def get_viral_load(self, t_since_infection, agent_id=None, profile_id=None):
        '''
        Get viral load at time t since infection.
        
        Args:
            t_since_infection (float or array): Time since infection in days
            agent_id: Agent identifier (to look up assigned profile)
            profile_id: Specific profile to use (overrides agent_id)
        
        Returns:
            viral_load (float or array): Viral load (virions)
        '''
        # Determine which profile to use
        if profile_id is None:
            if agent_id is not None:
                profile_id = self.agent_profiles.get(agent_id)
            
            if profile_id is None:
                # No specific profile requested, use random
                profile_id = np.random.randint(0, self.n_profiles)
        
        # Get the profile
        profile = self.profiles[profile_id]
        
        # Handle scalar and array inputs
        t = np.asarray(t_since_infection)
        scalar_input = np.isscalar(t_since_infection)
        
        # Interpolate from pre-computed trajectory
        viral_load = np.interp(t, self._t_points, profile['viral_loads'])
        
        # Zero out negative times
        viral_load = np.where(t >= 0, viral_load, 0)
        
        # Return scalar if input was scalar
        if scalar_input:
            return float(viral_load)
        
        return viral_load
    
    def get_profile_stats(self):
        '''Get statistics about the viral load profiles'''
        stats = {
            'n_profiles': self.n_profiles,
            'peak_loads': [p['peak_load'] for p in self.profiles],
            'peak_days': [p['peak_day'] for p in self.profiles],
            'parameters': [(p['beta'], p['gamma'], p['p']) for p in self.profiles]
        }
        return stats


def create_viral_load_function(config):
    '''
    Factory function to create viral load function from configuration.
    
    Args:
        config (dict or callable or ViralLoadFunction): Configuration for viral load function
    
    Returns:
        viral_load_func (ViralLoadFunction): Viral load function instance
    '''
    # If already a ViralLoadFunction instance, return as-is
    if isinstance(config, ViralLoadFunction):
        return config
    
    # If callable, wrap in UserSuppliedViralLoad
    if callable(config):
        return UserSuppliedViralLoad(config)
    
    # If dict, create based on type
    if isinstance(config, dict):
        vl_type = config.get('type', 'default')
        
        if vl_type == 'default':
            return DefaultViralLoad(**{k: v for k, v in config.items() if k != 'type'})
        elif vl_type == 'gaussian':
            return GaussianViralLoad(**{k: v for k, v in config.items() if k != 'type'})
        elif vl_type == 'piecewise':
            return PiecewiseLinearViralLoad(**{k: v for k, v in config.items() if k != 'type'})
        elif vl_type == 'hcd':
            return HCDViralLoad(**{k: v for k, v in config.items() if k != 'type'})
        elif vl_type == 'hcd_database':
            return HCDViralLoadDatabase(**{k: v for k, v in config.items() if k != 'type'})
        else:
            raise ValueError(f"Unknown viral load type: {vl_type}")
    
    # Default
    return DefaultViralLoad()