# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## ⚠️ IMPORTANT: Start Every Session Here

**1. CHECK TODO LIST FIRST**: Always review `TODO.md` for current implementation status and pending tasks.
**2. PROJECT OBJECTIVES**: Refer to `Project_Decription.md` for detailed requirements.
**3. IMPLEMENTATION PLAN**: See `Plan.md` for technical design and architecture.

## Current Project Focus

**TASK TRACKING**: The `TODO.md` file contains the complete list of completed and pending tasks. Always check this file at the start of each session to understand the current implementation status.

**IMPORTANT**: Please refer to `Project_Decription.md` for the current project objectives. We are modifying COVASIM's infection risk calculation to use a dose-response formulation based on inhaled viral particles, replacing the baseline-beta approach. This involves:
- Implementing a dose-response function: P_dose = 1 - exp(-N/N₀ · I)
- Separating direct (face-to-face) and indirect (shared air) transmission
- Integrating user-supplied viral load functions and CFD-based inhalation profiles
- Maintaining backward compatibility with existing functionality

See `Project_Decription.md` for full technical specifications and requirements.

## Current Status Summary (from TODO.md)
- **✅ 11/15 tasks completed (73%)**
- **Core dose-response functionality**: Fully implemented and tested
- **Examples and tests**: Complete and working
- **Pending**: CFD profiles, orientation factors, logging, documentation

**IMPLEMENTATION PLAN**: Refer to `Plan.md` for the detailed implementation plan, including:
- Current implementation status
- CFD-based inhalation profiles design
- Face-to-face orientation factors
- Enhanced environment mapping
- Testing and validation strategies

## About Covasim

Covasim is a stochastic agent-based COVID-19 simulator for analyzing infection dynamics and evaluating interventions like social distancing, testing, contact tracing, and vaccination. It was developed by the Institute for Disease Modeling with contributions from multiple institutions.

## Key Commands

### Installation
```bash
# Standard install (core dependencies only)
pip install -e .

# Full install (includes optional dependencies like Plotly, Optuna, SynthPops)
pip install -e .[full]
```

### Running Tests
```bash
# Run all tests with parallelization
cd tests
./run_tests

# Run tests with coverage report
./check_coverage

# Run a single test file
pytest test_baselines.py

# Run tests with specific options
pytest test_*.py -n auto --durations=0
```

### Code Quality
```bash
# Run style checks with pylint
cd tests
./check_style
# or
pylint ../covasim
```

### Documentation
```bash
# Build documentation
cd docs
make html
```

### Quick Start
```python
import covasim as cv
sim = cv.Sim()
sim.run()
sim.plot()
```

## Architecture Overview

### Core Classes

1. **Sim** (`sim.py`): Main simulation class that orchestrates the epidemic dynamics
   - Handles initialization, running, and results
   - Manages time steps and interventions
   - Integrates all other components

2. **People** (`people.py`): Manages the population of agents
   - Tracks health states (susceptible, exposed, infected, recovered, etc.)
   - Handles state transitions and disease progression
   - Manages immunity and vaccination status

3. **Interventions** (`interventions.py`): Implements various control measures
   - Testing strategies
   - Contact tracing
   - Social distancing
   - Vaccination campaigns
   - Custom interventions can be created by inheriting base classes

### Module Dependencies

The modules build on each other in this order:
1. `version.py`, `requirements.py` - Version info and dependency checking
2. `utils.py`, `misc.py` - Helper functions and utilities
3. `settings.py`, `defaults.py` - Configuration and default values
4. `parameters.py` - Parameter creation and validation
5. `plotting.py` - Visualization functions
6. `base.py` - Base classes (ParsObj, BaseSim, BasePeople)
7. `people.py` - Population management
8. `population.py` - Population generation (age, contacts, networks)
9. `interventions.py` - Intervention implementations
10. `immunity.py` - Immunity, variants, and neutralizing antibodies
11. `sim.py` - Main simulation class
12. `run.py` - Multi-sim runs, scenarios, and parallel execution
13. `analysis.py` - Analysis tools, fitting, and transmission trees

### Key Data Structures

- **Parameters**: Dict-like object containing all simulation parameters
- **People arrays**: NumPy arrays storing agent properties (one value per person)
- **Contact layers**: Different types of contacts (household, school, work, community)
- **Results**: Time series data stored as arrays in sim.results

### State Management

Agent states are tracked as boolean arrays:
- `susceptible`, `exposed`, `infectious`, `symptomatic`, `severe`, `critical`
- `recovered`, `dead`, `diagnosed`, `quarantined`, `vaccinated`

Dates are tracked as float arrays (day numbers):
- `date_exposed`, `date_infectious`, `date_symptomatic`, `date_diagnosed`
- `date_recovered`, `date_dead`, `date_vaccinated`

## Development Tips

### Running Single Tests
```python
# To run a specific test function
pytest tests/test_baselines.py::test_baseline -v
```

### Updating Baselines
When intentionally changing model behavior:
```bash
cd tests
./update_baseline
```

### Common Simulation Patterns
```python
# Custom parameters
sim = cv.Sim(pop_size=10000, n_days=180, rand_seed=1)

# With interventions
tp = cv.test_prob(symp_prob=0.1, asymp_prob=0.01)
ct = cv.contact_tracing(trace_probs=0.8)
sim = cv.Sim(interventions=[tp, ct])

# Loading data for calibration
sim = cv.Sim(datafile='data.csv')

# Running scenarios
scens = cv.Scenarios(base_sim=sim, scenarios={
    'baseline': {},
    'masks': {'interventions': cv.change_beta(days=30, changes=0.5)}
})
scens.run()
```

### Environment Variables
- `COVASIM_INTERACTIVE=0` - Disable plot display (useful for tests)
- `COVASIM_WARNINGS='error'` - Convert warnings to errors (for debugging)

## Python Version Requirements

Covasim requires Python 3.7-3.9 (Python 3.10+ not yet supported due to Numba compatibility).

## Python Environment 
local-covasim-env

## Working directory 
/Users/rahul/data/rahul/Work/drafts/eccomas-2024/data/agent_model/covasim-gituhub/Agent-covsim


