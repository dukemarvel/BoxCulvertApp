from hardy_cross import calculate_moments
import math

class BoxCulvert:

    def __init__(self):
        self.parameters = self.input_design_parameters()
        self.depth = 300 #Assume depth, D mm
        self.cover = 40 #concrete conver, c mm
        self.concrete_density = self.parameters["unit_weight_of_concrete"]
        self.steel_strength = self.parameters["steel_strength"]
        self.concrete_strength = self.parameters["concrete_strength"]
        self.bar_diameter = self.parameters["bar_diameter"]
        self.get_d_eff()
        self.get_ka()
        

        self.top_slab_design_moments = []
        self.bottom_slab_design_moments = []
        self.sidewall_design_moments = []

        self.get_design_moment_case_1()
        self.get_design_moment_case_2()
        self.get_design_moment_case_3()


    def input_design_parameters(self):
        parameters = {}
        parameters['width'] = float(input("Enter the width of the box culvert in mm: "))
        parameters['height'] = float(input("Enter the height of the box culvert in mm: "))
        parameters['length'] = float(input("Enter the length of the box culvert across the road in meters: "))
        parameters['thickness'] = float(input("Enter the assumed concrete thickness for the culvert in mm: "))
        parameters['wheel_load'] = float(input("Enter the HB loading in kN per wheel: "))
        parameters['unit_weight_of_concrete'] = float(input("Enter the unit weight of the concrete type used for designing the culvert in kN/m^3: "))
        parameters['steel_strength'] = float(input("Enter the steel strength in N/mm^2: "))
        parameters['concrete_strength'] = float(input("Enter the concrete strength in N/mm^2: "))
        parameters['density_of_fill'] = float(input("Enter the density of fill, Kg/m^3: "))
        parameters['angle_of_repose'] = float(input("Enter the angle of repose: "))
        parameters['depth_of_backfill'] = float(input("Enter the depth of backfill in mm: "))
        parameters['imposed_load'] = float(input("Enter the imposed load, KN/m^2: "))
        parameters['bar_diameter'] = float(input("Enter the bar diameter in mm: "))
        return parameters

    

    
    def get_design_moment_case_1(self):
        self.M_DC, self.M_DA, self.M_AD, self.M_AB = self.hardy_cross_moment_distribution_method()
        # top slab
        mid = self.calculate_midspan_moment_top_slab()
        top_slab_design_moment = abs(abs(self.M_AB) - mid)

        self.top_slab_design_moments.append(top_slab_design_moment)

        # bottom slab
        mid = self.calculate_midspan_moment_bottom_slab()
        bottom_slab_design_moment = abs(abs(self.M_DC) - mid)
        self.bottom_slab_design_moments.append(bottom_slab_design_moment)

        #sidewall
        mid = self.calculate_midspan_moment_sidewall()
        fem = abs(abs(self.M_DA) + abs(self.M_AD)) / 2
        sidewall_design_moment = abs(abs(mid) - fem)
        self.sidewall_design_moments.append(sidewall_design_moment)

    def get_design_moment_case_2(self):
        self.M_DC, self.M_DA, self.M_AD, self.M_AB = self.hardy_cross_moment_distribution_method(case=2)
        # top slab
        mid = self.calculate_midspan_moment_top_slab()
        top_slab_design_moment = abs(abs(self.M_AB) - mid)

        self.top_slab_design_moments.append(top_slab_design_moment)

        # bottom slab
        mid = self.calculate_midspan_moment_bottom_slab()
        bottom_slab_design_moment = abs(abs(self.M_DC) - mid)
        self.bottom_slab_design_moments.append(bottom_slab_design_moment)

        #sidewall
        mid = self.calculate_midspan_moment_sidewall(case=2)
        fem = abs(abs(self.M_DA) + abs(self.M_AD)) / 2
        sidewall_design_moment = abs(abs(mid) - fem)
        self.sidewall_design_moments.append(sidewall_design_moment)

    def get_design_moment_case_3(self):
        self.M_DC, self.M_DA, self.M_AD, self.M_AB = self.hardy_cross_moment_distribution_method(case=3)
        # top slab
        mid = self.calculate_midspan_moment_top_slab()
        top_slab_design_moment = abs(abs(self.M_AB) - mid)

        self.top_slab_design_moments.append(top_slab_design_moment)

        # bottom slab
        mid = self.calculate_midspan_moment_bottom_slab()
        bottom_slab_design_moment = abs(abs(self.M_DC) - mid)
        self.bottom_slab_design_moments.append(bottom_slab_design_moment)

        #sidewall
        mid = self.calculate_midspan_moment_sidewall(case=3)
        fem = abs(abs(self.M_DA) + abs(self.M_AD)) / 2
        sidewall_design_moment = abs(abs(mid) - fem)
        self.sidewall_design_moments.append(sidewall_design_moment)
        

    def get_ka(self):
        theta = math.radians(self.parameters['angle_of_repose'])
        self.Ka = (1 - math.sin(theta))/(1 + math.sin(theta))

    
    def get_d_eff(self):
        # d_eff = depth - cover - bar_diameter/2
         self.effective_depth = self.depth - self.cover -  self.bar_diameter/2 # in mm

    def calculate_self_weight_per_length(self):
        # this is for top slab own weight
        # Convert dimensions from mm to m for calculation
        thickness_m = self.parameters['thickness'] / 1000

        # Given density of concrete
        density_kN_m3 = self.concrete_density

        # Calculate self weight per unit length, KN/m^2
        self_weight_kN_per_m2 = 1.4 * thickness_m * density_kN_m3

    

        return self_weight_kN_per_m2
    

    def calculate_imposed_load(self):
        # Assuming the imposed load is given as a parameter
        imposed_load_kN_per_m = self.parameters['imposed_load'] # KN/m^2
        return imposed_load_kN_per_m

    def calculate_earth_fill_weight(self):
        # Calculate the weight of the earth fill above the culvert
        # Assuming the depth of earth fill and its density are given as parameters
        depth_of_backfill_m = self.parameters['depth_of_backfill'] / 1000  # Convert from mm to m
        density_of_fill_kN_m3 = self.parameters['density_of_fill']  # Assuming density is given in kN/m^3
        earth_fill_weight_kN_per_m = 1.6 * depth_of_backfill_m * density_of_fill_kN_m3 

        # KN/m^2
        return earth_fill_weight_kN_per_m

    def calculate_dead_load(self):
        self.load = self.calculate_self_weight_per_length() + self.calculate_imposed_load() + self.calculate_earth_fill_weight()
        return self.load
    
    def calculate_wheel_load(self):
        # Calculate wheel load based on depth_of_backfill and width_of_culvert
        depth_of_backfill = self.parameters['depth_of_backfill']/1000  # Assuming in mm
        width_of_culvert = self.parameters['width']/1000  # Assuming in mm
        tyre_width = 2000

        if depth_of_backfill > 3 * width_of_culvert:
            
            effective_depth_fill = 3 * width_of_culvert
        else:
            effective_depth_fill = depth_of_backfill

        if effective_depth_fill >= 0.5 * width_of_culvert:
            
        
            area = 4 * (effective_depth_fill ** 2)
            self.wheel_load = self.parameters['wheel_load'] / area  
            # Assuming 'wheel_load' parameter contains value for W
        elif depth_of_backfill == 0:  # When there's no earth fill
            self.wheel_load = self.parameters['wheel_load'] / (tyre_width ** 2)
        else:
            
            self.wheel_load = None   
        # return the value in unit of KN/m^2
        
        live_load = 1.3 * self.wheel_load * depth_of_backfill   # in KN/m
        return live_load
    
   
    
    def calculate_total_load_on_top_slab(self):
        # Calculate total load on top slab
        dead_load = self.calculate_dead_load()
        live_load = self.calculate_wheel_load()

        total_load_kN_per_m2 = dead_load + live_load
        #return KN/m^2
        return total_load_kN_per_m2

   

    def calculate_self_weight_of_sidewall(self):
        # Calculate self weight of sidewall
        unit_weight_of_concrete_kN_m3 = self.concrete_density  # Assuming in kN/m^3
        slab_thickness_m = self.parameters['thickness'] / 1000  # Convert from mm to m
        culvert_height_m = self.parameters['height'] / 1000

        self_weight_of_sidewall_kN_per_m = 0.4 * unit_weight_of_concrete_kN_m3 * slab_thickness_m * culvert_height_m
        return self_weight_of_sidewall_kN_per_m

    
    def calculate_total_load_on_sidewall(self, case=1):
        # Calculate total load on sidewall
        Z = self.parameters['depth_of_backfill'] / 1000  # Convert from mm to m
        thickness = self.parameters["thickness"] / 1000
        culvert_height_m = self.parameters['height'] / 1000 + thickness/2
        live_load = self.calculate_wheel_load()
        water_pressure = 9.8 * culvert_height_m 
        if case == 1:
            top_load = self.Ka * (live_load + self.calculate_imposed_load() + self.parameters['density_of_fill']*Z)  #pressure from filling and live load
            bottom_load = top_load +  self.Ka * self.parameters['density_of_fill'] * culvert_height_m
        if case == 2:
            top_load = self.Ka * (live_load + self.calculate_imposed_load() + self.parameters['density_of_fill']*Z) #pressure from filling and live load
            t_load = self.Ka * (live_load + self.calculate_imposed_load() + self.parameters['density_of_fill']*Z) #pressure from filling and live load
            bottom_load = (top_load +  1.6 * self.Ka * self.parameters['density_of_fill'] * culvert_height_m) - water_pressure
            top_load = t_load - bottom_load

        if case == 3:
            lat_earth_pressure = self.Ka * self.parameters['density_of_fill'] * culvert_height_m
            lat_surcharge_due_to_deadload_only = self.Ka * self.parameters['imposed_load']

            top_load = lat_surcharge_due_to_deadload_only
            bottom_load = (top_load + lat_earth_pressure) - water_pressure
            
        
        return top_load, bottom_load
    
    
    def calculate_total_load_on_bottom_slab(self):
        # Calculate dead load on bottom slab
        load_from_top_slab = self.calculate_total_load_on_top_slab()  # Assuming in kN/m^2
        load_of_walls = self.calculate_self_weight_of_sidewall()  # Assuming in kN/m
        
        thickness = self.parameters['thickness'] / 1000
        bottom_width = self.parameters['width'] / 1000 + 2 * thickness

        self_weight_of_bottom_slab = (thickness * bottom_width * 1) * self.concrete_density

        total_dead_load_KN_per_m = (load_from_top_slab * bottom_width + 2*(load_of_walls) + self_weight_of_bottom_slab) / bottom_width * 1
        return total_dead_load_KN_per_m

    
    
    def calculate_fixed_end_moment_top_slab(self):
        # Calculate fixed end moment based on total load and length of culvert
        total_load = self.calculate_total_load_on_top_slab()  # Total load on top slab in kN/m^2
        thickness = self.parameters["thickness"] / 1000
        length_of_culvert = self.parameters['width'] / 1000 + 2*(thickness / 2)  # Length of culvert in meters

        
        fixed_end_moment = total_load * (length_of_culvert ** 2) / 12  # in kN.m

        return fixed_end_moment

    
    def calculate_fixed_end_moment_sidewall(self, case=1):

        D = self.parameters['height'] / 1000 + self.parameters["thickness"]*0.001 / 2
        # Calculate fixed end moment for sidewall
        if case == 1:
            top_load, bottom_load = self.calculate_total_load_on_sidewall()
            # Calculate Fixed End Moment at the top using formula from image
            FEMt = (bottom_load * D** 2) / 30 + (top_load * D**2)/12  
            # Calculate Fixed End Moment at the bottom using formula from image
            FEMb = (bottom_load * D** 2) / 20 +(top_load * D**2)/12 # in kN.m

        if case == 2:
            top_load, bottom_load = self.calculate_total_load_on_sidewall(case=2)
            # Calculate Fixed End Moment at the top using formula from image
            FEMt = (top_load * D** 2) / 20 + (bottom_load * D**2)/12  
            # Calculate Fixed End Moment at the bottom using formula from image
            FEMb = (bottom_load * D** 2) / 12 +(top_load * D**2)/30 # in kN.m

        if case == 3:
            top_load, bottom_load = self.calculate_total_load_on_sidewall(case=3)

            FEMt = (top_load * D**2) / 12  - (bottom_load * D**2)/30

            FEMb = (top_load * D**2) / 12 - (bottom_load * D**2)/20  
            

        return FEMt, FEMb
    
    def calculate_fixed_end_moment_bottom_slab(self):
        # Calculate fixed end moment based on total load and length of culvert
        total_load = self.calculate_total_load_on_bottom_slab()  # Total load on top slab in kN/m^2
        length_of_culvert = self.parameters['width'] / 1000 + self.parameters['thickness']*0.001 / 2 # Length of culvert in meters

        # Calculate Fixed End Moment using formula from image
        fixed_end_moment = total_load * (length_of_culvert ** 2) / 12  # in kN.m

        return fixed_end_moment
    
    
    
    def calculate_midspan_moment_top_slab(self):
        # Calculate fixed end moment based on total load and length of culvert
        total_load = self.calculate_total_load_on_top_slab()  # Total load on top slab in kN/m^2
        length_of_culvert = self.parameters['width'] / 1000  # Length of culvert in meters

        # Calculate Fixed End Moment using formula from image
        midspan_moment = total_load * (length_of_culvert ** 2) / 8  # in kN.m

        return midspan_moment
    
    def calculate_midspan_moment_bottom_slab(self):
        # Calculate fixed end moment based on total load and length of culvert
        total_load = self.calculate_total_load_on_bottom_slab()  # Total load on top slab in kN/m^2
        length_of_culvert = self.parameters['width'] / 1000  # Length of culvert in meters
        # Calculate Fixed End Moment using formula from image
        midspan_moment = total_load * (length_of_culvert ** 2) / 8  # in kN.m

        return midspan_moment
    
    
    def calculate_midspan_moment_sidewall(self, case=1):
        if case == 1:
            top_load, bottom_load = self.calculate_total_load_on_sidewall()
            D = self.parameters['height'] / 1000
        
            moment_earth = (bottom_load* D**2) / 16
        
            live_moment = (top_load * D**2) / 8

        if case == 2:
            top_load, bottom_load = self.calculate_total_load_on_sidewall(case=2)
            D = self.parameters['height'] / 1000
            
            moment_earth = (bottom_load* D**2) / 16
            
            live_moment = (top_load * D**2) / 8

        if case == 3:
            top_load, bottom_load = self.calculate_total_load_on_sidewall(case=3)
            D = self.parameters['height'] / 1000
            
            moment_earth = (bottom_load* D**2) / 16
            
            live_moment = (top_load * D**2) / 8
        
        return moment_earth + live_moment
    

    def calculate_design_moment_top_slab(self):
        
        return max(self.top_slab_design_moments)

    def calculate_design_moment_bottom_slab(self):
        
        return max(self.bottom_slab_design_moments)

    def calculate_design_moment_sidewall(self):
        
        return max(self.sidewall_design_moments)
            
    def hardy_cross_moment_distribution_method(self, case=1):
        # Calculate Fixed End Moments (FEM) for each member
                # Calculate lengths of each member
        L_DC = self.parameters['width'] / 1000 + self.parameters['thickness']*0.001 / 2
        L_DA = self.parameters['height'] / 1000 + self.parameters["thickness"]*0.001 / 2
        L_AD = L_DA
        L_AB = L_DC

        if L_DC == L_DA == L_AD == L_AB:
            L_DC = L_DC/2
            L_AB = L_AB/2
        

        if case == 1:
            FEM_DC = -self.calculate_fixed_end_moment_bottom_slab()  # Negative for bottom slab DC
            FEM_DA, FEM_AD = self.calculate_fixed_end_moment_sidewall()  # FEM_DA is positive, FEM_AD is negative 
            FEM_AD *= -1
            FEM_AB = self.calculate_fixed_end_moment_top_slab()  # Positive for top slab AB 

            # Call the calculate_moments function from hardy.py
            DC_list, DA_list, AD_list, AB_list = calculate_moments(FEM_DC, FEM_DA, FEM_AD, FEM_AB, L_DC, L_DA, L_AD, L_AB)

    
            M_DC, M_DA, M_AD, M_AB = sum(DC_list), sum(DA_list), sum(AD_list), sum(AB_list)
        
        if case == 2:
            FEM_DC = -self.calculate_fixed_end_moment_bottom_slab()  # Negative for bottom slab DC
            FEM_DA, FEM_AD = self.calculate_fixed_end_moment_sidewall(case=2)  # FEM_DA is positive, FEM_AD is negative 
            FEM_AD *= -1
            FEM_AB = self.calculate_fixed_end_moment_top_slab()  # Positive for top slab AB 

            # Call the calculate_moments function from hardy.py
            DC_list, DA_list, AD_list, AB_list = calculate_moments(FEM_DC, FEM_DA, FEM_AD, FEM_AB, L_DC, L_DA, L_AD, L_AB)

    
            M_DC, M_DA, M_AD, M_AB = sum(DC_list), sum(DA_list), sum(AD_list), sum(AB_list)

        if case == 3:
            FEM_DC = -self.calculate_fixed_end_moment_bottom_slab()  # Negative for bottom slab DC
            FEM_DA, FEM_AD = self.calculate_fixed_end_moment_sidewall(case=3)  # FEM_DA is positive, FEM_AD is negative 
            FEM_AD *= -1
            FEM_AB = self.calculate_fixed_end_moment_top_slab()  # Positive for top slab AB 

            # Call the calculate_moments function from hardy.py
            DC_list, DA_list, AD_list, AB_list = calculate_moments(FEM_DC, FEM_DA, FEM_AD, FEM_AB, L_DC, L_DA, L_AD, L_AB)

    
            M_DC, M_DA, M_AD, M_AB = sum(DC_list), sum(DA_list), sum(AD_list), sum(AB_list)
            
        #print(f"DC_LIST: {DA_list}")
        return M_DC, M_DA, M_AD, M_AB



    def total_load_case_1(self):
        top_slab_load = self.calculate_total_load_on_top_slab()
        side_wall_load = self.calculate_total_load_on_sidewall()
        bottom_load = self.calculate_total_load_on_bottom_slab()

        return top_slab_load, side_wall_load, bottom_load


    def total_load_case_2(self):
        top_slab_load = self.calculate_total_load_on_top_slab()
        side_wall_load = self.calculate_total_load_on_sidewall(case=2)
        bottom_load = self.calculate_total_load_on_bottom_slab()

        return top_slab_load, side_wall_load, bottom_load
    
    def total_load_case_3(self):
        top_slab_load = self.calculate_total_load_on_top_slab()
        side_wall_load = self.calculate_total_load_on_sidewall(case=3)
        bottom_load = self.calculate_total_load_on_bottom_slab()

        return top_slab_load, side_wall_load, bottom_load
    
    def calculate_shear_top_slab(self):
        thickness = self.parameters["thickness"]
        W = self.calculate_total_load_on_top_slab()
        l = self.parameters['width'] / 1000 + 2*(thickness / 2)  # Length of culvert in meters

        return (W*l) / 2     # KN
    
    def calculate_shear_bottom_slab(self):
        thickness = self.parameters["thickness"]
        W = self.calculate_total_load_on_bottom_slab()
        l = self.parameters['width'] / 1000 + 2*(thickness / 2)  # Length of culvert in meters

        return (W*l) / 2     # KN

    def calculate_area_of_steel_required(self, moment):
        # Calculate area of steel required and generate final design outputs.
        fcu = self.concrete_strength
        fy = self.steel_strength
        b = 1000
        d = self.effective_depth
        k = moment * 10**6 / (fcu*b*d**2)
    
        
        if k < 0.167: # No compression reinforcement required
            Z = d * (0.5 + math.sqrt(0.25 - 0.882*k))

            Ast = moment * 10**6 / (0.87*fy*Z)
            
        max = 0.04 * b * d    # Maximum reinforcement required per ECN2 code
        
        f_ctm = 0.3 * fcu ** (2/3)
        min = 0.26 * (f_ctm / fy) * b * d # minimum reinforcement required per ECN2 code
        

        if min < 0.0013 * b * d:
            min = 0.0013 * b * d
            

        if Ast < min:
            Ast = min
        
        
        if Ast > max:
            Ast = max
        

        return math.ceil(Ast) + 100
    
    def calculate_distribution_reinforcement(self, As):
        main = As # main reinforcements

        d = self.effective_depth

        Astd = (0.12/100) * 1000 * d

        if Astd < 0.2 * main:
            return math.ceil(0.2*main)
        return math.ceil(Astd)
    
    def get_num_of_bars(self, As):
        d_bar = self.parameters['bar_diameter']
        A_bar = (math.pi * d_bar**2) / 4

        num_bar = As/A_bar

        return math.ceil(num_bar)
    
    def get_num_of_bars_d(self, As):
        d_bar = 8
        A_bar = (math.pi * d_bar**2) / 4

        num_bar = As/A_bar

        return math.ceil(num_bar)


    def get_spacing_of_bars(self, As):
        d_bar = self.parameters['bar_diameter']
        A_bar = (math.pi * d_bar**2) / 4
        
        spacing_of_bar = (A_bar/As) * 1000

        return math.floor(spacing_of_bar)

    def get_spacing_of_bars_d(self, As):
        d_bar = 8
        A_bar = (math.pi * d_bar**2) / 4
        
        spacing_of_bar = (A_bar/As) * 1000

        return math.floor(spacing_of_bar) 
    
    def calculate_area_of_steel_top_slab(self):
        moment = self.calculate_design_moment_top_slab()
        return self.calculate_area_of_steel_required(moment)

    def calculate_area_of_steel_bottom_slab(self):
        moment = self.calculate_design_moment_bottom_slab()
        return self.calculate_area_of_steel_required(moment)

    def calculate_area_of_steel_sidewall(self):
        moment = self.calculate_design_moment_sidewall()
        return self.calculate_area_of_steel_required(moment)
         

    def generate_final_design_outputs(self):
        # Generate final design outputs including detailed designs.\
        bar_diameter = self.parameters["bar_diameter"]
        Top = self.calculate_area_of_steel_top_slab()
        Side = self.calculate_area_of_steel_sidewall()
        Bottom = self.calculate_area_of_steel_bottom_slab()

    

        num_of_bars_top_slab = self.get_num_of_bars(Top)
        spacing_bar_top = self.get_spacing_of_bars(Top)
        num_of_bars_bottom_slab = self.get_num_of_bars(Bottom)
        spacing_bar_bottom = self.get_spacing_of_bars(Bottom)
        num_of_bars_sidewall_slab = self.get_num_of_bars(Side)
        spacing_bar_sidewall = self.get_spacing_of_bars(Side)

        print()
        print("*********************MAIN REINFORCEMENT****************************")
        print(f"TOP SLAB: Provide {num_of_bars_top_slab}-Y{bar_diameter}mm @ {spacing_bar_top}mm C/C each face ({Top} mm^2)")
        print(f"SIDEWALL SLAB: Provide {num_of_bars_sidewall_slab}-Y{bar_diameter}mm @ {spacing_bar_sidewall}mm C/C each face ({Side} mm^2)")
        print(f"BOTTOM SLAB: Provide {num_of_bars_bottom_slab}-Y{bar_diameter}mm @ {spacing_bar_bottom}mm C/C each face ({Bottom} mm^2)")

        dia_bar = 8
        t = self.calculate_distribution_reinforcement(Top)
        s = self.calculate_distribution_reinforcement(Side)
        b = self.calculate_distribution_reinforcement(Bottom)
        num_bars_t = self.get_num_of_bars_d(t)
        space_top = self.get_spacing_of_bars_d(t)
        num_bars_s = self.get_num_of_bars_d(s)
        space_side = self.get_spacing_of_bars_d(s)
        num_bars_b = self.get_num_of_bars_d(b)
        space_bottom = self.get_spacing_of_bars_d(b)

        print()
        print("**********************DISTRIBUTION REINFORCEMENTS*************************")
        print(f"TOP SLAB: Provide {num_bars_t}-Y{dia_bar}mm @ {space_top}mm C/C each face ({t} mm^2)")
        print(f"SIDEWALL SLAB: Provide {num_bars_s}-Y{dia_bar}mm @ {space_side}mm C/C each face ({s} mm^2)")
        print(f"BOTTOM SLAB: Provide {num_bars_b}-Y{dia_bar}mm @ {space_bottom}mm C/C each face ({b} mm^2)")

if __name__ == "__main__":

    test = BoxCulvert()

    print(f'WEIGHT OF SLAB: {test.calculate_self_weight_per_length()}')

    print(f'wheel load: {test.calculate_wheel_load()}')

    print(f" Earth Fill Pressure: {test.calculate_earth_fill_weight()}")

    print()
    print("All loads")
    print("CASE 1")
    top, side, bottom = test.total_load_case_1()
    print(f'TOP SLAB: {top}')
    print(f'BOTTOM SLAB: {bottom}')
    print(f'SIDE WALL: {side}')

    print()
    print("CASE 2")
    top, side, bottom = test.total_load_case_2()
    print(f'TOP SLAB: {top}')
    print(f'BOTTOM SLAB: {bottom}')
    print(f'SIDE WALL: {side}')

    print()
    print("CASE 3")
    top, side, bottom = test.total_load_case_3()
    print(f'TOP SLAB: {top}')
    print(f'BOTTOM SLAB: {bottom}')
    print(f'SIDE WALL: {side}')



    print("All midspan moments CASE 1")
    print(f"Midspan SLAB: {test.calculate_midspan_moment_top_slab()}")
    print(f"Midspan SIDEWALL: {test.calculate_midspan_moment_sidewall()}")
    print(f"Midspan BOTTOM: {test.calculate_midspan_moment_bottom_slab()}")

    print()
    print("All midspan moments CASE 2")
    print(f"Midspan SLAB: {test.calculate_midspan_moment_top_slab()}")
    print(f"Midspan SIDEWALL: {test.calculate_midspan_moment_sidewall(case=2)}")
    print(f"Midspan BOTTOM: {test.calculate_midspan_moment_bottom_slab()}")

    print()
    print("All midspan moments CASE 3")
    print(f"Midspan SLAB: {test.calculate_midspan_moment_top_slab()}")
    print(f"Midspan SIDEWALL: {test.calculate_midspan_moment_sidewall(case=3)}")
    print(f"Midspan BOTTOM: {test.calculate_midspan_moment_bottom_slab()}")



    print("")

    print("All Fixed End Moments CASE 1")

    print(f'FIXED END TOP SLAB: {test.calculate_fixed_end_moment_top_slab()}')
    print(f"FIXED END SIDEWALL: {test.calculate_fixed_end_moment_sidewall()}")
    print(f'FIXED END BOTTOM SLAB: {test.calculate_fixed_end_moment_bottom_slab()}')

    print()
    print("All Fixed End Moments CASE 2")

    print(f'FIXED END TOP SLAB: {test.calculate_fixed_end_moment_top_slab()}')
    print(f"FIXED END SIDEWALL: {test.calculate_fixed_end_moment_sidewall(case=2)}")
    print(f'FIXED END BOTTOM SLAB: {test.calculate_fixed_end_moment_bottom_slab()}')

    print()
    print("All Fixed End Moments CASE 3")

    print(f'FIXED END TOP SLAB: {test.calculate_fixed_end_moment_top_slab()}')
    print(f"FIXED END SIDEWALL: {test.calculate_fixed_end_moment_sidewall(case=3)}")
    print(f'FIXED END BOTTOM SLAB: {test.calculate_fixed_end_moment_bottom_slab()}')

    print()

    print("***Design MomentS********")
    print(test.calculate_design_moment_top_slab())
    print(test.calculate_design_moment_bottom_slab())
    print(test.calculate_design_moment_sidewall())






    print()
    print("*******************************FINAL DESIGN OUTPUTS*************************************")
    print(test.generate_final_design_outputs())



    input("Press any key to exit...")



