from helpers import *
from ice_albedo import *

def run_model(
                # model timing
                dt = 1e3/1e6,     # time step [Myr]
                t_max = 80,       # model run time [Myr]
    
                # ice-albedo feedback
                Ti = 260,         # ice-covered threshold temperature [K]
                To = 295,         # ice-free threshold temperature [K]
                ai = 0.6,         # ice-covered albedo
                ao = 0.28,        # ice-free albedo
    
                # temperatures
                T0 = 288,         # initial temperature [K]
                Teq = 288,        # equilibrium temperature [K]
                
                # volcanic degassing and weathering
                V_C = 8.5,        # volcanic degassing [examol/Myr]
                V_red = 1.7,     # volcanic reduced gases [examol/Myr]
                W_sea = 1.6,      # seafloor weathering [examol/Myr]
                n = 0.2,          # silicate weathering feedback strength
    
                # productivity
                forg = 0.2,       # fraction of buried carbon that is organic
                CPsed = 250,      # sedimentary organic C:P ratio
                nb=1,             # exponent for organic burial dependence on P
    
                # oxygen
                O20 = 0.1,        # O2 as % of PAL
    
                # phosphorus
                P_conc = 2.2,     # seawater phosphate concentration [uM]
                W_pho0 = 4e-2,    # phosphorus weathering [examol/Myr]
                Pinorg_mode = "scale", # inorganic burial mode
    
                # perturbation
                C_imb = 0,        # imposed generalized carbon cycle imbalance [examol/Myr]
                tau = np.inf,     # e-folding timescale of forcing decay [Myr]
                beta = 0,         # fraction of carbon cycle imbalance happening organically
                W_LIP = 0,        # imposed LIP C sequestration [examol/Myr] (replaced C_imb)
                n_LIP = False,      # LIP silicate weathering feedback strength
                PC_LIP = 0.009,   # molar ratio of weatherable P to Ca+Mg in LIP
                suppress_Borg = False, # allows burial enhancement due to P input to be suppressed to isolate effect of generalized C imbalance
        
                # output
                verbose = False,
            ):
    
    # --- setup ---
    t = 0
    snowball = False
    
    # --- initial conditions ---
    if not n_LIP:
        n_LIP = n
    
    # climate
    T = T0
    pCO2, a = get_conditions(T0, Ti=Ti, To=To, ai=ai, ao=ao)
    
    # equilibrium climate
    pCO20, a0 = get_conditions(Teq, Ti=Ti, To=To, ai=ai, ao=ao)
        # getting pCO20 ensures that weathering is still scaled appropriately if the system starts in a snowball
    
    # carbon inventory    
    Nat = pCO2/1e6 * Na/1e18  # convert to examol
    N = 2.83 * np.sqrt(pCO2/280)  # surficial carbon [examol]
    
    # oxygen inventory
    O2 = O20 * 0.2 * 1.8e20/1e18  # examol
        
    # phosphorus inventory
    P = P0 = P_conc * Voc / 1e6 / 1e18  # convert uM to examol
    
    if verbose:
        print("Initial conditions")
        print(f"    T0 = {T0:0.0f} K")
        print(f"    albedo = {a:0.2f}")
        print(f"    pCO2 = {pCO2:0.0f} ppm")
        print(f"    pO2 = {O20:0.1f} PAL")
        print(f"    {Nat:0.1f} Emol C in atmos")
        print(f"    {N:0.1f} Emol C in surficial system")
        
    if verbose and Teq != T0:
        print("Equilibrium conditions")
        print(f"    Teq = {Teq:0.0f} K")
        print(f"    albedo = {a0:0.2f}")
        print(f"    pCO2 = {pCO20:0.1e} ppm")
        
    # --- initial rates ---
    
    # carbon cycle
    WC_total = V_C-V_red
    W_sil = W_sil0 = WC_total - W_sea
    
    # balancing C and O cycles to find W_org and B_org
    W_org = W_org0 = ((V_C*forg)-V_red)/(1-forg)
    BorgC = BorgC0 = (V_C+W_org)*forg
    
    # phosphorus cycle
    W_pho = W_pho0
    Pbur_org = Pbur_org0 = BorgC/CPsed
    Pbur_inorg = Pbur_inorg0 = W_pho - Pbur_org
    
    if Pinorg_mode == "scale":
        f_Porg = Pbur_org / W_pho0
    elif Pinorg_mode == "constant":
        pass
    else:
        raise ValueError("Invalid input for Pinorg_mode: must be 'constant' or 'scale'")
        
    if verbose:
        print(f"Carbon cycle")
        print(f"    CO2 degassing = {V_C:0.1f} Emol/Myr")
        print(f"    Sil. weathering = {W_sil:0.1f} Emol/Myr")
        print(f"    Sea. weathering = {W_sea:0.1f} Emol/Myr")
        print(f"    Org. weathering = {W_org:0.1f} Emol/Myr")
        print(f"    Org. burial = {BorgC:0.1f} Emol/Myr")
        
        print(f"Oxygen cycle")
        print(f"    Reduced gases = {V_red:0.1f} Emol/Myr")
        print(f"    Org. weathering = {W_org:0.1f} Emol/Myr")
        print(f"    Org. C burial = {BorgC:0.1f} Emol/Myr")
        
        print(f"Phosphorus cycle")
        print(f"    Weathering = {W_pho:0.2f} Emol/Myr")
        print(f"    Org. burial = {Pbur_org:0.2f} Emol/Myr")
        print(f"    Inorg. burial = {Pbur_inorg:0.2f} Emol/Myr")
    
    # --- checks and balances (pun intended) ---
    
    # check that nothing is negative
    rates_to_check = ['W_org','W_sil','BorgC','W_pho','Pbur_org','Pbur_inorg']
    for rate in rates_to_check:
        val = locals()[rate]
        if val < 0:
            raise ValueError(f"initial {rate} is negative.")
            
    # double check that everything balances
    floating_point_err = 1e-15
    dNdt = V_C + W_org - W_sil - W_sea - BorgC
    if np.abs(dNdt) > floating_point_err:
        warnings.warn(f"Carbon cycle out of balance by {np.abs(dNdt):0.2e}")
        
    dOdt = BorgC - W_org - V_red
    if np.abs(dOdt) > floating_point_err:
        warnings.warn(f"Oxygen cycle out of balance by {np.abs(dOdt):0.2e}")
  
    dPdt = W_pho - Pbur_org - Pbur_inorg
    if np.abs(dPdt) > floating_point_err:
        warnings.warn(f"Phosphorus cycle out of balance by {np.abs(dPdt):0.2e}")
        
    # --- setup for perturbation ---
    
    C_imb0 = C_imb
    W_LIP0 = W_LIP
    if C_imb and W_LIP:
        warnings.warn(f"Model running with both a generalized carbon cycle imbalance and LIP weathering.")
    
    # --- set up arrays to store values for tracking ---
    
    all_params = [k for k, v in locals().items() if isinstance(v, (int, float, bool, np.number))]
    all_params += ['dNatdt', 'pH', 'f']
        # add parameters that don't exist yet because they're created inside loop
    arrays = {param: [] for param in all_params}
    for param in all_params:
        if param in locals():
            arrays[param].append(locals()[param])
        else:
            arrays[param].append(np.nan)
        
    # --- run model ---
    while t < t_max:
        
        # check if in a snowball
        snowball = a >= (ai - 0.01)
        
        # --- update sources and sinks ---
        
        # weathering
        y = pCO2 / pCO20
        W_sil = W_sil0 * y**n if not snowball else 0
        W_org = W_org0 * y**n if not snowball else 0
        W_pho = W_pho0 * y**n if not snowball else 0
        
        # perturbation 
        C_imb = C_imb0 * np.exp(-t/tau) if not snowball else 0
        O_imb = C_imb * beta
        W_LIP = W_LIP0 * np.exp(-t/tau) * y**n_LIP if not snowball else 0
        Wpho_LIP = W_LIP * PC_LIP
        
        # organic burial
        BorgC = BorgC0 * (P/P0)**nb
        if suppress_Borg:
            # special modification to isolate generalized carbon cycle imbalance
            BorgC = BorgC0 if not snowball else 0
        
        # phosphorus sinks
        if Pinorg_mode == "scale":
            Pbur_inorg = Pbur_inorg0 * P/P0
        elif Pinorg_mode == "constant":
            Pbur_inorg = Pbur_inorg0
        Pbur_org = BorgC/CPsed
        
        # --- update reservoirs ---
        
        # oxygen
        dOdt = BorgC - V_red - W_org                # main equation
        dOdt += O_imb                               # add extra burial if applicable
        O2 = max(0, O2 + dOdt * dt)                 # update reservoir
        
        # phosphorus
        dPdt = W_pho - Pbur_org - Pbur_inorg        # main equation
        dPdt += Wpho_LIP                            # add LIP weathering if appliable
        P = max(0, P + dPdt * dt)                   # update reservoir
        
        # carbon
        dNdt = V_C + W_org - W_sil - W_sea - BorgC  # main equation
        dNdt -= C_imb                               # substract generalized imbalance if applicable
        dNdt -= W_LIP                               # subtract LIP weathering if applicable

        # partition carbon between reservoirs
        F = Nat / N
        pH = pH_fn(F)
        f = f_fn(pH)
        dNatdt = dNdt * f
        
        # update carbon reservoirs
        N = max(0, N + dNdt * dt)
        Nat = max(0, Nat + dNatdt * dt)
        pCO2 = Nat * 1e18 / Na * 1e6
        
        # --- final updates ---
        
        # update albedo and temp
        a = albedo(T, Ti=Ti, To=To, ai=ai, ao=ao)
        T = get_T(pCO2, a)
        
        # increment time step
        t += dt
        
        # record parameters for tracking
        for param in all_params:
            arrays[param].append(locals().get(param))
    
    # --- outputs ----
    
    for param in all_params:
        arrays[param] = np.array(arrays[param])
        
    # convert O2 reservoir from Emol to % of PAL
    arrays['pO2'] = arrays['O2']*1e18/1.8e20/0.2
    
    return arrays
