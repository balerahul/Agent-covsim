# Project: Modify COVASIM Infection Risk Calculation to Use Dose–Response Formulation

## Objective
Refactor the infection risk computation in COVASIM to use a **dose–response model** based on the inhaled number of viral particles, replacing the existing baseline-beta-based approach for infection probability.

## Background
In the current COVASIM model (see original formulation), infection probability is computed as:

$$ P_T = \beta \cdot \gamma \cdot f_{quar} \cdot f_{asymp} \cdot f_{iso} \cdot \beta_{layer} \cdot \lambda $$

Here:
- `β`, `β_layer`, and `λ` collectively account for base infection probability, location effects, and viral load variations.

In the new approach, **these factors are consolidated** into a **dose–response function**:

$$
P_{\text{dose}} = 1 - \exp\!\left(-\frac{N}{N_0} \cdot I\right)
$$

where:
- $ N = N_{\text{direct}} + N_{\text{indirect}} $
- $ N_{\text{direct}} = \dot{N}^{v}_{direct} \tau_{direct} \lambda(t) $ : Viral particles inhaled during face-to-face contact (depends on exposure duration, viral load, proximity, etc.)
   - $\dot{N}^{v}_{direct}$ : The rate of deposition of viral particles on the upper respiratory tract. This is be setup as a list for different  types of interactions to be filled in by the user. 
   - $\tau_{direct}$: Contact duration during face-to-face direct contact. 
   - $\lambda(t)$ : Viral load as function of time $t$ after infection of the infected host. 
- $ N_{\text{indirect}} = \dot{N}^{v}_{indirect} \tau_{indirect} \lambda(t) $ : Viral particles inhaled via shared-air exposure in enclosed spaces (depends on CFD-derived inhalation rates for each environment type)
   - $\dot{N}^{v}_{indirect}$ : The rate of deposition of viral particles on the upper respiratory tract due to passive indirect interaction (being present in the same room).  This is be setup as a list for different  types of interactions to be filled in by the user. 
   - $\tau_{direct}$: Contact duration during indirect contact. 
   - $\lambda(t)$ : Viral load as function of time $t$ after infection of the infected host. 
- $ N_0 $ : Infectious dose scaling parameter  
- $ I $ : Infectivity scaling factor (can be set to 1 if not required)

The final transmission probability becomes:

$$
P_T = P_{\text{dose}} \cdot \gamma \cdot f_{quar} \cdot f_{asymp} \cdot f_{iso}
$$

## Key Requirements

### 1. Replace Infection Probability Core
- Remove `β`, `β_layer`, and `λ` from the infection probability formula when using the new dose–response option.
- Implement `P_dose` using the dose–response function above.
- Keep `γ`, `f_quar`, `f_asymp`, and `f_iso` as multiplicative modifiers.

### 2. Inhaled Dose Components
- **Direct Dose (`N_direct`)**  
  - Computed from face-to-face contact duration, viral load of the infector at time of contact, and other modifying factors (mask usage, distance, breathing rate, etc.).
  - Contact duration should be sampled from a distribution per contact type (Gaussian or user-defined). A sample code is provided in tmp-files/interaction-time.py. This can be used for the time being. 
- **Indirect Dose (`N_indirect`)**  
  - Computed from shared-environment exposure (e.g., same room) using precomputed CFD-based inhalation rates per environment type.
  - Scales with the number of infectors present and their viral loads.
- For both of these a tabulated or list of N = theta * exposure_duration * viral_load * face2face_factor will be provided by the user. 
  - Face to face factor is used to scale up or down the total droplets attaching to upper resporatory tract (URT) depending on the interaction is head on or sideways, etc. 

### 3. Viral Load Function
- Remove the existing viral load model in COVASIM.
- Add a new, user-supplied `viral_load_function(t_since_infection)` to be used for all dose calculations.
- This function will be provided externally; implementation should allow plug-and-play.
- The details of the viral load model is provided in temp-files/README_HCD.md file and the related data for the viral load model is contained in temp-files/hcd_population_summary.json. 

### 4. Integration
- Introduce a simulation parameter (e.g., `use_dose_response`) to enable the new calculation path.
- Ensure backward compatibility by keeping the old method available when the flag is `False`.
- Integrate new dose computations into the per-contact infection loop in COVASIM.

### 5. Inputs & Configuration
- `N0` (infectious dose scale)
- `viral_load_function` (user-defined)
- Contact duration distributions per contact type
- CFD-derived inhalation profiles per environment type

### 6. Output & Debug
- Optionally log per-contact `N_direct`, `N_indirect`, total `N`, and resulting `P_dose` for debugging or analysis.
- Ensure reproducibility with fixed random seeds for duration sampling.

## Deliverables
1. Modified COVASIM codebase with:
   - Dose–response infection risk calculation
   - New inputs and configuration options
   - Backward compatibility support
2. Example configuration file demonstrating:
   - Sample `viral_load_function`
   - Contact duration distributions
   - CFD inhalation profiles
3. Documentation outlining:
   - How the new infection probability is computed
   - How to provide viral load and CFD data
   - Calibration guidance for `N0` and other parameters
4. Unit tests covering:
   - Direct and indirect dose calculations
   - Dose–response function outputs
   - Integration with COVASIM’s transmission pipeline

## Notes
- The AI model implementing this should determine:
  - Data structures for storing per-contact durations and inhalation profiles
  - Whether dose calculations are vectorized or computed per-contact
  - How viral load is integrated over the exposure duration
  - Exact handling of units and scaling
- CFD profiles and viral load function will be provided externally by the user.

---
