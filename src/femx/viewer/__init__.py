from enum import Enum


class Output(Enum):
    displacement_x = 1
    displacement_y = 2
    stress_x = 3
    stress_y = 4
    stress_xy = 5
    strain_x = 6
    strain_y = 7
    strain_xy = 8
    p_stress_1 = 9
    p_stress_2 = 10
    averaged_stress_x = 11
    averaged_stress_y = 12
    averaged_stress_xy = 13
    averaged_strain_x = 14
    averaged_strain_y = 15
    averaged_strain_xy = 16
    von_mises_stress = 17
