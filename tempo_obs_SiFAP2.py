#!/usr/bin/env python3
import numpy as np

def main():
    print("-" * 50)
    print(" SiFAP2 Observation Time Estimation Tool ")
    print("-" * 50)

    # 1. Magnitude Input
    try:
        mag = float(input("Enter magnitude (V band): "))
    except ValueError:
        print("Error: Please enter a valid number for magnitude.")
        return

    # 2. Background-Subtracted Amplitude Input
    try:
        a_sub_percent = float(input("Enter target bkg-subtracted amplitude (in %): "))
        A_sub = a_sub_percent / 100.0  # Convert percentage to fraction
    except ValueError:
        print("Error: Please enter a valid number for amplitude.")
        return

    # 3. Background Rate Confirmation
    bkg_default = 1.5e4  # 15 kcps
    bkg_choice = input(f"Assume default bkg rate of 15 kcps ({bkg_default:.0f} c/s)? [Y/n]: ").strip().lower()
    
    if bkg_choice in ['', 'y', 'yes', 's', 'si']:
        R_bkg = bkg_default
    else:
        try:
            R_bkg_kcps = float(input("Enter desired bkg rate in kcps: "))
            R_bkg = R_bkg_kcps * 1000.0
        except ValueError:
            print("Error: Invalid input. Falling back to default value.")
            R_bkg = bkg_default

    # 4. Trials Selection (Single vs Multiple)
    trials_choice = input("Is this a Single Trial search? [Y/n]: ").strip().lower()
    if trials_choice in ['', 'y', 'yes', 's', 'si']:
        N_trials = 1
    else:
        try:
            N_trials = int(input("Enter the number of independent trials (N_trials): "))
            if N_trials < 1:
                print("Error: N_trials must be at least 1. Setting to 1.")
                N_trials = 1
        except ValueError:
            print("Error: Invalid integer. Falling back to Single Trial.")
            N_trials = 1

    # --- PHYSICAL CALCULATIONS ---
    
    # SiFAP2 calibration parameters
    m = 0.942  
    q = -18.67 
    
    # R_net (Net source count rate) in c/s
    y_kcps = 10**(-((m * mag + q) / 2.512))
    R_net = y_kcps * 1000.0
    
    # R_tot (Total count rate) in c/s
    R_tot = R_net + R_bkg
    
    # Threshold Power P (Vaughan 1994)
    P_single = 2.24 * np.log(1 / 2.7e-3)
    if N_trials == 1:
        P = P_single
    else:
        P = P_single + 2.24 * np.log(N_trials)
    
    # Calculate amplitude A (not subtracted) from A_sub
    A = A_sub * (R_net / R_tot)
    
    # Calculate observation time T_obs
    T_obs = (1.61**2 * P) / (R_tot * A**2)
    
    # Final integrated counts
    N_net = R_net * T_obs
    N_bkg = R_bkg * T_obs
    N_tot = N_net + N_bkg
    
    # Double check consistency of the achieved bkg-subtracted amplitude
    A_achieved_sub = 1.61*np.sqrt(P/N_tot)*N_tot/N_net

    # --- PRINT RESULTS ---
    print("\n" + "="*50)
    print("                    RESULTS                    ")
    print("="*50)
    print(f"Net Source Rate (R_net)     : {R_net:.2f} c/s")
    print(f"Background Rate (R_bkg)     : {R_bkg:.2f} c/s")
    print(f"Total Count Rate (R_tot)    : {R_tot:.2f} c/s")
    print("-" * 50)
    print(f"Number of Trials            : {N_trials}")
    print(f"Threshold Power (P)         : {P:.2f}")
    print(f"Source counts (N_net)       : {N_net:.2e}")
    print(f"Background counts (N_bkg)   : {N_bkg:.2e}")
    print(f"Total counts (N_tot)        : {N_tot:.2e}")
    print("-" * 50)
    
    # Time formatting
    if T_obs < 60:
        print(f"REQUIRED OBS. TIME (Tobs)   : {T_obs:.2f} seconds")
    elif T_obs < 3600:
        print(f"REQUIRED OBS. TIME (Tobs)   : {T_obs/60:.2f} minutes ({T_obs:.1f} s)")
    else:
        print(f"REQUIRED OBS. TIME (Tobs)   : {T_obs/3600:.2f} hours ({T_obs/60:.1f} min)")
    
    print("-" * 50)
    print(f"Bkg-subtracted amplitude: {A_sub * 100:.4f} %")
    print(f"Amplitude A (including bkg)   : {A * 100:.4f} %")
    print(f"(Double check) Bkg-subtracted amplitude: {A_achieved_sub * 100:.4f} %")
    print("="*50 + "\n")

if __name__ == "__main__":
    main()
