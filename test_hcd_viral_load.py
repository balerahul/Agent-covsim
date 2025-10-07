#!/usr/bin/env python
'''
Test script for HCD viral load model implementation
'''

import numpy as np
import matplotlib.pyplot as plt
from covasim.viral_dynamics import HCDViralLoad, create_viral_load_function

def test_hcd_viral_load():
    '''Test the HCD viral load model'''
    
    print("Testing HCD Viral Load Model")
    print("=" * 50)
    
    # Test 1: Create HCD viral load with mean parameters
    print("\n1. Testing with mean population parameters...")
    vl_mean = HCDViralLoad(sample_params=False)
    
    # Get parameters
    params = vl_mean.get_parameters()
    print(f"   Beta: {params['beta']:.2e}")
    print(f"   Gamma: {params['gamma']:.4f}")
    print(f"   p: {params['p']:.4f}")
    
    # Test viral load at specific time points
    test_times = [0, 1, 3, 5, 7, 10, 14, 21]
    print("\n   Viral load at different time points:")
    for t in test_times:
        vl = vl_mean.get_viral_load(t)
        print(f"   Day {t:2d}: {vl:.2e} virions")
    
    # Test 2: Create HCD viral load with sampled parameters
    print("\n2. Testing with sampled parameters (3 agents)...")
    for i in range(3):
        vl_sampled = HCDViralLoad(sample_params=True, seed=i)
        params = vl_sampled.get_parameters()
        print(f"\n   Agent {i+1}:")
        print(f"   Beta: {params['beta']:.2e}, Gamma: {params['gamma']:.4f}, p: {params['p']:.4f}")
        
        # Get peak viral load
        times = np.linspace(0, 14, 100)
        loads = vl_sampled.get_viral_load(times)
        peak_idx = np.argmax(loads)
        print(f"   Peak viral load: {loads[peak_idx]:.2e} on day {times[peak_idx]:.1f}")
    
    # Test 3: Test with factory function
    print("\n3. Testing factory function...")
    config = {'type': 'hcd', 'sample_params': False}
    vl_factory = create_viral_load_function(config)
    vl_test = vl_factory.get_viral_load(5)
    print(f"   Viral load at day 5: {vl_test:.2e}")
    
    # Test 4: Plot viral load trajectories
    print("\n4. Plotting viral load trajectories...")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot mean trajectory
    ax1 = axes[0]
    times = np.linspace(0, 21, 200)
    vl_mean = HCDViralLoad(sample_params=False)
    loads_mean = vl_mean.get_viral_load(times)
    ax1.semilogy(times, loads_mean, 'b-', linewidth=2, label='Population mean')
    ax1.set_xlabel('Days since infection')
    ax1.set_ylabel('Viral load (virions)')
    ax1.set_title('HCD Viral Load - Population Mean')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim([1e0, 1e10])
    
    # Plot multiple sampled trajectories
    ax2 = axes[1]
    np.random.seed(42)
    for i in range(10):
        vl_sampled = HCDViralLoad(sample_params=True)
        loads_sampled = vl_sampled.get_viral_load(times)
        ax2.semilogy(times, loads_sampled, alpha=0.5, linewidth=1)
    
    ax2.set_xlabel('Days since infection')
    ax2.set_ylabel('Viral load (virions)')
    ax2.set_title('HCD Viral Load - 10 Sampled Agents')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([1e0, 1e10])
    
    plt.tight_layout()
    plt.savefig('hcd_viral_load_test.png', dpi=150)
    print(f"   Plot saved as 'hcd_viral_load_test.png'")
    
    # Test 5: Test array input
    print("\n5. Testing array input...")
    times_array = np.array([0, 5, 10, 15])
    loads_array = vl_mean.get_viral_load(times_array)
    print(f"   Input shape: {times_array.shape}, Output shape: {loads_array.shape}")
    print(f"   Viral loads: {loads_array}")
    
    print("\n" + "=" * 50)
    print("All tests completed successfully!")

if __name__ == '__main__':
    test_hcd_viral_load()