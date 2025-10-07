'''
Test script for HCDViralLoadDatabase implementation
'''

import numpy as np
import matplotlib.pyplot as plt
from covasim.viral_dynamics import HCDViralLoadDatabase, HCDViralLoad
import covasim as cv

def test_database_creation():
    '''Test the HCDViralLoadDatabase creation and profile generation'''
    
    print("Testing HCD Viral Load Database")
    print("=" * 50)
    
    # Create database with 12 profiles
    print("\n1. Creating viral load database with 12 profiles...")
    db = HCDViralLoadDatabase(n_profiles=12, seed=42)
    
    # Get profile statistics
    stats = db.get_profile_stats()
    print(f"   Number of profiles: {stats['n_profiles']}")
    print(f"   Peak loads range: {min(stats['peak_loads']):.2e} - {max(stats['peak_loads']):.2e}")
    print(f"   Peak days range: {min(stats['peak_days']):.1f} - {max(stats['peak_days']):.1f}")
    
    # Test profile assignment
    print("\n2. Testing profile assignment to agents...")
    for agent_id in range(5):
        profile_id = db.assign_profile(agent_id)
        print(f"   Agent {agent_id} assigned profile {profile_id}")
    
    # Test consistent viral load retrieval
    print("\n3. Testing consistent viral load retrieval...")
    agent_id = 0
    times = [1, 5, 10, 15]
    loads = []
    for t in times:
        load = db.get_viral_load(t, agent_id=agent_id)
        loads.append(load)
        print(f"   Agent {agent_id}, Day {t}: {load:.2e} virions")
    
    # Verify consistency - should get same values for same agent
    print("\n4. Verifying consistency for same agent...")
    load_check = db.get_viral_load(5, agent_id=0)
    assert abs(load_check - loads[1]) < 1e-10, "Inconsistent viral load for same agent!"
    print(f"   ✓ Consistent viral load confirmed: {load_check:.2e}")
    
    return db


def plot_database_profiles(db):
    '''Plot all viral load profiles in the database'''
    
    print("\n5. Plotting viral load profiles...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    times = np.linspace(0, 20, 200)
    
    # Plot all profiles
    for profile in db.profiles:
        loads = [db.get_viral_load(t, profile_id=profile['id']) for t in times]
        if profile['id'] == 0:
            # Highlight the mean profile
            ax1.semilogy(times, loads, 'r-', linewidth=2, alpha=1.0, label='Mean profile')
        else:
            ax1.semilogy(times, loads, alpha=0.6, linewidth=1)
    
    ax1.set_xlabel('Days since infection')
    ax1.set_ylabel('Viral load (virions)')
    ax1.set_title('Database: 12 Pre-computed Profiles')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_ylim([1e0, 1e10])
    
    # Compare with individual sampling approach
    for i in range(12):
        vl = HCDViralLoad(sample_params=True, seed=i)
        loads = [vl.get_viral_load(t) for t in times]
        ax2.semilogy(times, loads, alpha=0.6, linewidth=1)
    
    ax2.set_xlabel('Days since infection')
    ax2.set_ylabel('Viral load (virions)')
    ax2.set_title('Old Approach: Individual Sampling')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([1e0, 1e10])
    
    plt.tight_layout()
    plt.savefig('hcd_database_comparison.png', dpi=150)
    print(f"   Plot saved as 'hcd_database_comparison.png'")
    plt.show()


def test_in_simulation():
    '''Test the database in an actual COVASIM simulation'''
    
    print("\n6. Testing in COVASIM simulation...")
    
    # Create simulation with dose-response and HCD database
    sim = cv.Sim(
        pop_size=100,
        n_days=30,
        use_dose_response=True,
        viral_load_function={'type': 'hcd_database', 'n_profiles': 12, 'seed': 42},
        verbose=0
    )
    
    # Run simulation
    sim.run()
    
    # Check that viral load profiles were assigned
    if hasattr(sim.viral_load_func, 'agent_profiles'):
        n_assigned = len(sim.viral_load_func.agent_profiles)
        print(f"   ✓ {n_assigned} agents were assigned viral load profiles")
        
        # Show distribution of profile assignments
        profile_counts = {}
        for agent_id, profile_id in sim.viral_load_func.agent_profiles.items():
            profile_counts[profile_id] = profile_counts.get(profile_id, 0) + 1
        
        print("   Profile distribution:")
        for profile_id in sorted(profile_counts.keys()):
            print(f"     Profile {profile_id}: {profile_counts[profile_id]} agents")
    
    print("   ✓ Simulation completed successfully")
    
    return sim


if __name__ == '__main__':
    # Run tests
    db = test_database_creation()
    plot_database_profiles(db)
    sim = test_in_simulation()
    
    print("\n" + "=" * 50)
    print("All tests passed successfully!")