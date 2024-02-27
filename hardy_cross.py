def calculate_moments(FEM_DC, FEM_DA, FEM_AD, FEM_AB, L_DC, L_DA, L_AD, L_AB):
    # Calculate stiffness using moment of inertia and length
    def K(L1, L2):
        return (1 / L1) / ((1 / L1) + (1 / L2))

    # Calculate the Distribution Factors (DF) for each joint
    DF_DC = K(L_DC, L_DA)  # Distribution factor at joint D
    DF_DA = K(L_DA, L_DC)
    DF_AD = K(L_AD, L_AB) # Distribution factor at joint A
    DF_AB = K(L_AB, L_AD)

    # Initialize lists for each member
    DC_list = [FEM_DC]
    DA_list = [FEM_DA]
    AD_list = [FEM_AD]
    AB_list = [FEM_AB]


    # Set a threshold for the unbalanced moment
    threshold = 0.01

    # Iterate until the unbalanced moment is less than the threshold
    while True:
        # Calculate unbalanced moments at each joint
        unbalanced_moment_D = DC_list[-1] + DA_list[-1]   # Unbalanced moment at joint D 
        unbalanced_moment_A = AD_list[-1] + AB_list[-1]   # Unbalanced moment at joint A 

        # Apply distribution factors to calculate balanced moments 
        Mdc_balanced = abs(unbalanced_moment_D * DF_DC)   # Absolute value of balanced moment for member DC 
        Mda_balanced = abs(unbalanced_moment_D * DF_DA)   # Absolute value of balanced moment for member DA 

        Mab_balanced = abs(unbalanced_moment_A * DF_AB)   # Absolute value of balanced moment for member AB
        Mad_balanced = abs(unbalanced_moment_A * DF_AD)   # Absolute value of balanced moment for member AD 

        # Assign the sign of the FEM with the smaller absolute value to the balanced moments
        if abs(FEM_DC) < abs(FEM_DA):
            Mdc_balanced *= -1 if FEM_DC < 0 else 1
            Mda_balanced *= -1 if FEM_DC < 0 else 1
        else:
            Mdc_balanced *= -1 if FEM_DA < 0 else 1
            Mda_balanced *= -1 if FEM_DA < 0 else 1

        if abs(FEM_AD) < abs(FEM_AB):
            Mad_balanced *= -1 if FEM_AD < 0 else 1
            Mab_balanced *= -1 if FEM_AD < 0 else 1
        else:
            Mad_balanced *= -1 if FEM_AB < 0 else 1
            Mab_balanced *= -1 if FEM_AB < 0 else 1

        # Calculate carry over moments
        CM_DA = 0.5 * Mad_balanced  # Carry over moment for member DA
        CM_AD = 0.5 * Mda_balanced  # Carry over moment for member AD
        CM_DC = 0  # Carry over moment for member DC
        CM_AB = 0  # Carry over moment for member AB

        # Add the balanced moments and carry over moments to the lists
        DC_list.extend([Mdc_balanced, CM_DC])
        DA_list.extend([Mda_balanced, CM_DA])
        AD_list.extend([Mad_balanced, CM_AD])
        AB_list.extend([Mab_balanced, CM_AB])

        # Check if the unbalanced moment is less than the threshold
        if max(abs(unbalanced_moment_D), abs(unbalanced_moment_A)) < threshold:
            break

    # Return the lists
    return DC_list, DA_list, AD_list, AB_list
