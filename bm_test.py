import cantera as ct
import math 
import yaml

# def bm_rate(bm_string, deltaH, T):
#     bm_list = bm_string.split()
#     A = float(bm_list[0][15:-1])
#     b = float(bm_list[1][2:-1])
#     E0 = float(bm_list[2][12:-1])
#     E = float(bm_list[3][2:-1])
#     w = float(bm_list[4][2:-1]) 
#     R = 8314.46261815324 # J/ (kmol.K)
    
#     if deltaH < -4 * E0:
#         E = 0 
#     elif deltaH > 4 * E0:
#         E = deltaH
#     else:
#         vp = 2 * w * (w+E0) / (w - E0)
#         E = (w + deltaH / 2) * (vp - 2 * w + deltaH) ** 2 / (vp ** 2 - 4 * w ** 2 + deltaH ** 2)
    
#     rate = A * T ** b * math.exp(-E / (R * T))
#     # print(f'A: {A}, b:{b}, E:{E}, w:{w}')
#     return rate

def bm_rate(A, b, E0, w, deltaH, T):

    R = 8314.46261815324 # J/ (kmol.K)
#     print(A, b, E, w)
    if deltaH < -4 * E0:
        E = 0 
    elif deltaH > 4 * E0:
        E = deltaH
    else:
        vp = 2 * w * (w+E0) / (w - E0)
        E = (w + deltaH / 2) * (vp - 2 * w + deltaH) ** 2 / (vp ** 2 - 4 * w ** 2 + deltaH ** 2)
    
    rate = A * T ** b * math.exp(-E / (R * T))
#     print(f'A: {A}, b:{b}, E:{E}, w:{w}') 
    return rate


def change_species_enthalpy(species_name, dH):
    """
    Find the species by name and change it's enthlapy by dH (in J/kmol)
    """
    index = gas.species_index(species_name)

    species = gas.species(index)
    print(f"Initial H(298) = {species.thermo.h(298)/1e6:.1f} kJ/mol")
    dx = dH / ct.gas_constant  # 'dx' is in fact (delta H / R). Note that R in cantera is 8314.462 J/kmol
    assert isinstance(species.thermo, ct.NasaPoly2)
    # print(species.thermo.coeffs)
    perturbed_coeffs = species.thermo.coeffs.copy()
    perturbed_coeffs[6] += dx
    perturbed_coeffs[13] += dx
    
    species.thermo = ct.NasaPoly2(species.thermo.min_temp, species.thermo.max_temp, 
                            species.thermo.reference_pressure, perturbed_coeffs)
    #print(species.thermo.coeffs)
    gas.modify_species(index, species)
    print(f"Modified H(298) = {species.thermo.h(298)/1e6:.1f} kJ/mol")

yaml_file = '/home/chao/bm_test/bm_example.yml'
gas = ct.Solution(yaml_file)
temp = gas.T

# print(gas.forward_rate_constants, gas.reaction(0).rate)
# bm_str1 = str(gas.reaction(0).rate)
with open(yaml_file, 'r') as f:
    content = yaml.load(f, Loader=yaml.FullLoader)
rate_constant = content['reactions'][0]['rate-constant']
A = rate_constant['A'] / 1000
b = rate_constant['b']
E0 = rate_constant['Ea'] * 4.184e3
w = float(rate_constant['w']) * 4.184e3

r1 = bm_rate(A, b, E0, w, gas.delta_enthalpy[0], temp)


print("Cantera output rate constant:",gas.forward_rate_constants[0], "test rate cosntant:", r1)
print("The rate parameters before change the enthalpy", gas.reaction(0).rate)

change_species_enthalpy('OH', +100e6)

r2 = bm_rate(A, b, E0, w, gas.delta_enthalpy[0], temp)

# print(gas.forward_rate_constants, gas.reaction(0).rate)

print("Cantera output rate constant:", gas.forward_rate_constants[0], "test rate cosntant:", r2)
print("The rate parameters after change the enthalpy", gas.reaction(0).rate)