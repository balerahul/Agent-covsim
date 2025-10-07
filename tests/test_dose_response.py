'''
Unit tests for dose-response model functionality.

Tests both the DoseResponseModel class and ViralLoadFunction implementations
to ensure correct behavior and integration with COVASIM.
'''

import pytest
import numpy as np
import sciris as sc
import covasim as cv
from covasim.dose_response import DoseResponseModel
from covasim.viral_dynamics import (
    DefaultViralLoad, GaussianViralLoad, PiecewiseLinearViralLoad, 
    HCDViralLoad, UserSuppliedViralLoad, create_viral_load_function
)


class TestDoseResponseModel:
    '''Test the DoseResponseModel class'''
    
    def test_initialization(self):
        '''Test model initialization with default and custom parameters'''
        # Default initialization
        model = DoseResponseModel()
        assert model.N0 == 100
        assert model.infectivity_scale == 1.0
        assert len(model.deposition_rates_direct) == 5
        assert len(model.deposition_rates_indirect) == 5
        
        # Custom initialization
        custom_direct = {'h': 15.0, 'w': 8.0}
        custom_indirect = {'h': 3.0, 'w': 2.0}
        model = DoseResponseModel(
            N0=200, 
            infectivity_scale=1.5,
            deposition_rates_direct=custom_direct,
            deposition_rates_indirect=custom_indirect
        )
        assert model.N0 == 200
        assert model.infectivity_scale == 1.5
        assert model.deposition_rates_direct['h'] == 15.0
        assert model.deposition_rates_indirect['h'] == 3.0
    
    def test_compute_direct_dose(self):
        '''Test direct dose calculation'''
        model = DoseResponseModel()
        
        # Test scalar inputs
        N_direct = model.compute_direct_dose(10.0, 5.0, 1e6)
        expected = 10.0 * 5.0 * 1e6
        assert N_direct == expected
        
        # Test array inputs
        deposition_rates = np.array([10.0, 15.0])
        durations = np.array([5.0, 3.0])
        viral_loads = np.array([1e6, 2e6])
        N_direct = model.compute_direct_dose(deposition_rates, durations, viral_loads)
        expected = np.array([10.0 * 5.0 * 1e6, 15.0 * 3.0 * 2e6])
        np.testing.assert_array_equal(N_direct, expected)
    
    def test_compute_indirect_dose(self):
        '''Test indirect dose calculation'''
        model = DoseResponseModel()
        
        # Test scalar inputs
        N_indirect = model.compute_indirect_dose(2.0, 60.0, 1e5)
        expected = 2.0 * 60.0 * 1e5
        assert N_indirect == expected
        
        # Test array inputs
        deposition_rates = np.array([2.0, 3.0])
        durations = np.array([60.0, 45.0])
        viral_loads = np.array([1e5, 5e5])
        N_indirect = model.compute_indirect_dose(deposition_rates, durations, viral_loads)
        expected = np.array([2.0 * 60.0 * 1e5, 3.0 * 45.0 * 5e5])
        np.testing.assert_array_equal(N_indirect, expected)
    
    def test_compute_infection_probability(self):
        '''Test infection probability calculation'''
        model = DoseResponseModel(N0=100, infectivity_scale=1.0)
        
        # Test zero dose
        P_dose = model.compute_infection_probability(0.0)
        assert P_dose == 0.0
        
        # Test known values
        N_total = 100.0  # Equal to N0
        P_dose = model.compute_infection_probability(N_total)
        expected = 1.0 - np.exp(-1.0)  # 1 - exp(-1)
        assert abs(P_dose - expected) < 1e-10
        
        # Test array inputs
        N_total = np.array([0.0, 50.0, 100.0, 200.0])
        P_dose = model.compute_infection_probability(N_total)
        expected = 1.0 - np.exp(-N_total / 100.0)
        np.testing.assert_array_almost_equal(P_dose, expected)
        
        # Test probability bounds
        assert np.all(P_dose >= 0.0)
        assert np.all(P_dose <= 1.0)
        
        # Test negative input handling
        P_dose_neg = model.compute_infection_probability(-10.0)
        assert P_dose_neg == 0.0
    
    def test_compute_dose_response_transmission(self):
        '''Test full transmission calculation'''
        model = DoseResponseModel(N0=100, infectivity_scale=1.0)
        
        # Simple test case
        sources = np.array([0, 1])
        targets = np.array([2, 3])
        viral_loads = np.array([0, 1e6, 0, 0])  # Only source 1 has viral load
        direct_durations = np.array([5.0, 3.0])
        indirect_durations = np.array([60.0, 45.0])
        
        trans_probs = model.compute_dose_response_transmission(
            sources, targets, viral_loads, direct_durations, indirect_durations,
            layer='h'
        )
        
        # Check that results are valid probabilities
        assert len(trans_probs) == 2
        assert np.all(trans_probs >= 0.0)
        assert np.all(trans_probs <= 1.0)
        
        # First contact (source 0) should have zero transmission (no viral load)
        assert trans_probs[0] == 0.0
        
        # Second contact should have positive transmission
        assert trans_probs[1] > 0.0
    
    def test_with_relative_factors(self):
        '''Test transmission with relative transmissibility and susceptibility'''
        model = DoseResponseModel(N0=1000)  # Higher N0 to avoid saturation
        
        sources = np.array([0])
        targets = np.array([1])
        viral_loads = np.array([1e5, 0])  # Lower viral load to avoid saturation
        direct_durations = np.array([2.0])  # Shorter duration
        indirect_durations = np.array([20.0])
        
        # Test with relative transmissibility
        rel_trans = np.array([2.0, 1.0])  # Double transmissibility for source
        trans_probs_high = model.compute_dose_response_transmission(
            sources, targets, viral_loads, direct_durations, indirect_durations,
            rel_trans=rel_trans
        )
        
        # Test without relative factors
        trans_probs_base = model.compute_dose_response_transmission(
            sources, targets, viral_loads, direct_durations, indirect_durations
        )
        
        # High transmissibility should give higher probability
        assert trans_probs_high[0] > trans_probs_base[0]
        
        # Test with relative susceptibility
        rel_sus = np.array([1.0, 0.5])  # Reduced susceptibility for target
        trans_probs_low = model.compute_dose_response_transmission(
            sources, targets, viral_loads, direct_durations, indirect_durations,
            rel_sus=rel_sus
        )
        
        # Lower susceptibility should give lower probability
        assert trans_probs_low[0] < trans_probs_base[0]
    
    def test_integrate_viral_load(self):
        '''Test viral load integration over time'''
        model = DoseResponseModel()
        
        # Simple constant viral load function
        def constant_viral_load(t):
            return 1e6
        
        avg_load = model.integrate_viral_load(constant_viral_load, 0.0, 1.0)
        assert abs(avg_load - 1e6) < 1e-10
        
        # Linear viral load function
        def linear_viral_load(t):
            return 1e6 * t
        
        avg_load = model.integrate_viral_load(linear_viral_load, 0.0, 2.0)
        expected = 1e6 * 1.0  # Average of linear function from 0 to 2
        assert abs(avg_load - expected) < 1e3  # Small tolerance for numerical integration
        
        # Test edge case: t_start >= t_end
        avg_load = model.integrate_viral_load(constant_viral_load, 1.0, 1.0)
        assert avg_load == 0.0


class TestViralLoadFunctions:
    '''Test viral load function implementations'''
    
    def test_default_viral_load(self):
        '''Test DefaultViralLoad implementation'''
        vl_func = DefaultViralLoad(incubation_period=2.0, peak_day=5.0, peak_load=1e8)
        
        # Test before incubation
        assert vl_func.get_viral_load(1.0) == 0.0
        
        # Test at peak
        peak_load = vl_func.get_viral_load(5.0)
        assert peak_load > 0
        
        # Test array input
        times = np.array([0, 2, 5, 10])
        loads = vl_func.get_viral_load(times)
        assert len(loads) == 4
        assert loads[0] == 0.0  # Before incubation
        assert loads[2] > loads[1]  # Increasing to peak
    
    def test_gaussian_viral_load(self):
        '''Test GaussianViralLoad implementation'''
        vl_func = GaussianViralLoad(peak_day=5.0, peak_load=1e8, sigma=2.0)
        
        # Test at peak
        peak_load = vl_func.get_viral_load(5.0)
        assert abs(peak_load - 1e8) < 1e6
        
        # Test symmetry around peak
        load_before = vl_func.get_viral_load(3.0)
        load_after = vl_func.get_viral_load(7.0)
        assert abs(load_before - load_after) < 1e6
        
        # Test negative time handling
        assert vl_func.get_viral_load(-1.0) == 0.0
    
    def test_piecewise_linear_viral_load(self):
        '''Test PiecewiseLinearViralLoad implementation'''
        times = [0, 2, 5, 10]
        values = [0, 1e6, 1e8, 1e4]
        vl_func = PiecewiseLinearViralLoad(times=times, values=values)
        
        # Test exact points
        assert vl_func.get_viral_load(0) == 0
        assert vl_func.get_viral_load(5) == 1e8
        
        # Test interpolation
        mid_load = vl_func.get_viral_load(3.5)  # Midpoint between 2 and 5
        expected = (1e6 + 1e8) / 2
        assert abs(mid_load - expected) < 1e6
    
    def test_user_supplied_viral_load(self):
        '''Test UserSuppliedViralLoad wrapper'''
        def custom_func(t):
            return 1e6 * np.exp(-0.1 * t)
        
        vl_func = UserSuppliedViralLoad(custom_func)
        
        # Test that it calls the user function
        load = vl_func.get_viral_load(5.0)
        expected = 1e6 * np.exp(-0.5)
        assert abs(load - expected) < 1e3
        
        # Test error for non-callable
        with pytest.raises(ValueError):
            UserSuppliedViralLoad("not_callable")
    
    def test_hcd_viral_load(self):
        '''Test HCDViralLoad implementation'''
        # Create test parameters file
        test_params = {
            "model_constants": {
                "U0": 4e8,
                "I0": 0,
                "V0": 1,
                "c": 5.2,
                "V_m": 1e-3
            },
            "mu_log": {
                "beta": -9.0,
                "gamma": 0.5,
                "p": 6.0
            },
            "sigma_log_marginal": {
                "beta": 0.5,
                "gamma": 0.3,
                "p": 0.4
            }
        }
        
        # Write temporary parameter file
        import tempfile
        import json
        import os
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(test_params, f)
            temp_file = f.name
        
        try:
            # Test with mean parameters
            vl_func = HCDViralLoad(params_file=temp_file, sample_params=False)
            
            # Test basic functionality
            load_early = vl_func.get_viral_load(1.0)
            load_peak = vl_func.get_viral_load(5.0)
            load_late = vl_func.get_viral_load(15.0)
            
            # Should show typical viral load trajectory
            assert load_early >= 0
            assert load_peak > load_early
            assert load_late < load_peak
            
            # Test parameter retrieval
            params = vl_func.get_parameters()
            assert 'beta' in params
            assert 'gamma' in params
            assert 'p' in params
            
        finally:
            os.unlink(temp_file)
    
    def test_create_viral_load_function(self):
        '''Test viral load function factory'''
        # Test with existing instance
        existing = DefaultViralLoad()
        result = create_viral_load_function(existing)
        assert result is existing
        
        # Test with callable
        def custom_func(t):
            return 1e6
        
        result = create_viral_load_function(custom_func)
        assert isinstance(result, UserSuppliedViralLoad)
        
        # Test with dict configs
        config_default = {'type': 'default', 'peak_load': 2e8}
        result = create_viral_load_function(config_default)
        assert isinstance(result, DefaultViralLoad)
        assert result.peak_load == 2e8
        
        config_gaussian = {'type': 'gaussian', 'sigma': 3.0}
        result = create_viral_load_function(config_gaussian)
        assert isinstance(result, GaussianViralLoad)
        assert result.sigma == 3.0
        
        # Test unknown type
        with pytest.raises(ValueError):
            create_viral_load_function({'type': 'unknown'})


class TestIntegrationWithCovasim:
    '''Test integration of dose-response model with COVASIM simulation'''
    
    def test_dose_response_simulation(self):
        '''Test running a simulation with dose-response model enabled'''
        # Create simple simulation with dose-response
        pars = {
            'pop_size': 100,
            'n_days': 10,
            'use_dose_response': True,
            'N0': 100,
            'infectivity_scale': 1.0,
            'viral_load_function': {
                'type': 'default',
                'peak_load': 1e7
            }
        }
        
        sim = cv.Sim(pars)
        sim.run()
        
        # Check that simulation completed
        assert hasattr(sim, 'results')
        assert len(sim.results['cum_infections']) == pars['n_days'] + 1
        
        # Should have some infections
        final_infections = sim.results['cum_infections'][-1]
        assert final_infections > 0
    
    def test_backward_compatibility(self):
        '''Test that dose-response doesn't break standard simulations'''
        # Standard simulation
        sim_standard = cv.Sim(pop_size=100, n_days=10, use_dose_response=False)
        sim_standard.run()
        
        # Should complete without errors
        assert hasattr(sim_standard, 'results')
        assert len(sim_standard.results['cum_infections']) == 11


if __name__ == '__main__':
    # Run tests
    pytest.main([__file__, '-v'])