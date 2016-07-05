import math


def write_in_script():
    filename = "in.kick"
    f = open(filename, "w")
    f.write("units           lj\natom_style      full\n#read_restart    old_config_1.dat\nread_data       Kick_in.dat\n")
    f.write("neighbor        0.3 bin\n")
    f.write("neigh_modify    every 1 delay 1\n")

    f.write("bond_style      fene\nbond_coeff      1 30.0 1.5 1.0 1.0\nspecial_bonds   lj 0.0 1.0 1.0\n")
    f.write("angle_style     harmonic\nangle_coeff      1 10.0 180\n")
    f.write("pair_style      table linear 1000\n")
    f.write("pair_coeff      2 2 tabulated_potential TAIL_TAIL\n")
    f.write("pair_coeff      1 2 tabulated_potential HEAD_TAIL\n")
    f.write("pair_coeff      1 1 tabulated_potential HEAD_TAIL\n")
    f.write("pair_coeff      1 3 tabulated_potential MONOMER_HEAD\n")
    f.write("pair_coeff      1 4 tabulated_potential MONOMER_HEAD\n")
    f.write("pair_coeff      1 5 tabulated_potential MONOMER_HEAD\n")
    f.write("pair_coeff      2 3 tabulated_potential MONOMER_TAIL\n")
    f.write("pair_coeff      2 4 tabulated_potential MONOMER_TAIL\n")
    f.write("pair_coeff      2 5 tabulated_potential MONOMER_TAIL\n")
    f.write("pair_coeff      3 3 tabulated_potential MONOMER_MONOMER\n")
    f.write("pair_coeff      3 4 tabulated_potential MONOMER_MONOMER\n")
    f.write("pair_coeff      3 5 tabulated_potential MONOMER_MONOMER\n")
    f.write("pair_coeff      4 4 tabulated_potential MONOMER_MONOMER\n")
    f.write("pair_coeff      4 5 tabulated_potential MONOMER_MONOMER\n")
    f.write("pair_coeff      5 5 tabulated_potential MONOMER_MONOMER\n")

    f.write("pair_modify        shift yes\n")
    f.write("neigh_modify exclude molecule all\n\n")
    f.write("# Additional variables to check\n")
    f.write("group           head type 1\n")
    f.write("group           tail type 2\n")
    f.write("group           bilayer union head tail\n")
    f.write("variable        top equal bound(bilayer,zmax)\n")
    f.write("variable        bottom equal bound(bilayer,zmin)\n")
    f.write("variable        Lz      equal zhi-zlo\n\n")
    f.write("# Fixes\n")

    f.write("velocity        all create 1.0 1 \n")

    f.write("fix             MCSWP all atom/swap 10 1 1234 1 types 4 5 \n")

    f.write("dump            id all xyz 100 kick_out.xyz\n")
    f.write("restart         10000 old_config_1.dat  old_config_2.dat\n")
    f.write("fix             0 all rigid/nve molecule\n")
    f.write("fix             2 all langevin 1 1 1 12345\n")
    f.write("fix             1 all press/berendsen x 0.0 0.0 1000.0 y 0.0 0.0 1000.0\n")
    f.write("dump            dDUMPALL all custom 500 \"data.lammpstrj\" id type x y z \n")
    f.write("dump_modify     dDUMPALL sort id \n")

    f.write("thermo          3000\n")
    f.write("thermo_style    custom step temp press etotal epair vol v_top v_bottom v_Lz\n")
    f.write("thermo_modify   flush yes\n")
    f.write("timestep        0.008\n")
    f.write("run             15000000\n")
    f.close()

# def make_monomers(xindex, yindex, Vx, Vy, Vz, )

if __name__ == "__main__":

    side_length = 1.1
    Lx = Ly = 15
    N_ref = 6

    '''Radius of spheres for proteins.'''
    radius = 2
    N_skip = (radius + 1)*(radius + 1)

    '''Box dimensions.'''
    Lz = 3*side_length*Lx
    Lx *= side_length
    Ly *= side_length

    positions = ["\nAtoms \n \n"]
    bonds = ["\nBonds \n \n"]
    angles = ["\nAngles \n \n"]

    k = 0
    m = 0
    Nbonds = 0
    bond_id = 0
    Nangles = 0
    angle_id = 0
    for j in range(1, int(math.floor(Ly)-1)):
        for s in range(1, int(math.floor(Lx)-1)):
            k += 1
            ID = k
            Vx = s * side_length
            Vy = j * side_length
            Vz = .5 * Lz

            if (10 > j > 7) and (10 > s > 7):
                if j == 8 and s == 8:
                    for q in range(0, 2):
                        m += 1
                        if q == 0:
                            positions.append("\t " + str(m) + " " + str(k) + " 3 0 " + str(Vx) + " " + str(Vy) +
                                     " " + str(Vz + (q - 2) * side_length) + " 0 0 0 \n")
                        if q == 1:
                            positions.append("\t " + str(m) + " " + str(k) + " 4 0 " + str(Vx) + " " + str(Vy) +
                                     " " + str(Vz + (q - 2) * side_length) + " 0 0 0 \n")
                    for q in range(1,3):
                        m += 1
                        if q == 1:
                            positions.append("\t " + str(m) + " " + str(k) + " 4 0 " + str(Vx) + " " + str(Vy) +
                                     " " + str(Vz + q * side_length) + " 0 0 0 \n")
                        if q == 2:
                            positions.append("\t " + str(m) + " " + str(k) + " 3 0 " + str(Vx) + " " + str(Vy) +
                                     " " + str(Vz + q * side_length) + " 0 0 0 \n")

                    # Building Reservoir for exchanging with membrane
                    for q in range(0, 3):
                        m += 1
                        if q == 0:
                            positions.append("\t " + str(m) + " " + str(k) + " 5 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Lz + (q - 2) * side_length) + " 0 0 0 \n")
                            # bond_id += 1
                            # bonds.append("\t " + str(bond_id) + " 1 " + str(m) + " " + str(m+1) + " \n")
                            # Nbonds += 1
                        if q == 1:
                            positions.append("\t " + str(m) + " " + str(k) + " 5 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Lz + (q - 2) * side_length) + " 0 0 0 \n")
                            # bond_id += 1
                            # bonds.append("\t " + str(bond_id) + " 1 " + str(m) + " " + str(m+1) + " \n")
                            # Nbonds += 1
                        if q == 2:
                            positions.append("\t " + str(m) + " " + str(k) + " 5 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Lz + (q - 2) * side_length) + " 0 0 0 \n")

            elif (5 > j > 2) and (5 > s > 2):
                if j == 3 and s == 3:
                    for q in range(0, 2):
                        m += 1
                        if q == 0:
                            positions.append("\t " + str(m) + " " + str(k) + " 3 0 " + str(Vx) + " " + str(Vy) +
                                     " " + str(Vz + (q - 2) * side_length) + " 0 0 0 \n")
                        if q == 1:
                            positions.append("\t " + str(m) + " " + str(k) + " 4 0 " + str(Vx) + " " + str(Vy) +
                                     " " + str(Vz + (q - 2) * side_length) + " 0 0 0 \n")
                    for q in range(1,3):
                        m += 1
                        if q == 1:
                            positions.append("\t " + str(m) + " " + str(k) + " 4 0 " + str(Vx) + " " + str(Vy) +
                                     " " + str(Vz + q * side_length) + " 0 0 0 \n")
                        if q == 2:
                            positions.append("\t " + str(m) + " " + str(k) + " 3 0 " + str(Vx) + " " + str(Vy) +
                                     " " + str(Vz + q * side_length) + " 0 0 0 \n")

                    # Building Reservoir for exchanging with membrane
                    for q in range(0, 3):
                        m += 1
                        if q == 0:
                            positions.append("\t " + str(m) + " " + str(k) + " 5 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Lz + (q - 2) * side_length) + " 0 0 0 \n")
                            # bond_id += 1
                            # bonds.append("\t " + str(bond_id) + " 1 " + str(m) + " " + str(m+1) + " \n")
                            # Nbonds += 1
                        if q == 1:
                            positions.append("\t " + str(m) + " " + str(k) + " 5 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Lz + (q - 2) * side_length) + " 0 0 0 \n")
                            # bond_id += 1
                            # bonds.append("\t " + str(bond_id) + " 1 " + str(m) + " " + str(m+1) + " \n")
                            # Nbonds += 1
                        if q == 2:
                            positions.append("\t " + str(m) + " " + str(k) + " 5 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Lz + (q - 2) * side_length) + " 0 0 0 \n")

            else:

                for q in range(0, 3):
                    m += 1
                    if q == 0:
                        positions.append("\t " + str(m) + " " + str(k) + " 1 0 " + str(Vx) + " " + str(Vy) +
                                         " " + str(Vz + (q - 2) * side_length) + " 0 0 0 \n")
                        bond_id += 1
                        bonds.append("\t " + str(bond_id) + " 1 " + str(m) + " " + str(m+1) + " \n")
                        Nbonds += 1
                    if q == 1:
                        positions.append("\t " + str(m) + " " + str(k) + " 2 0 " + str(Vx) + " " + str(Vy) +
                                         " " + str(Vz + (q - 2) * side_length) + " 0 0 0 \n")
                        bond_id += 1
                        bonds.append("\t " + str(bond_id) + " 1 " + str(m) + " " + str(m+1) + " \n")
                        Nbonds += 1
                    if q == 2:
                        positions.append("\t " + str(m) + " " + str(k) + " 2 0 " + str(Vx) + " " + str(Vy) +
                                         " " + str(Vz + (q - 2) * side_length) + " 0 0 0 \n")
                        angle_id += 1
                        angles.append("\t " + str(angle_id) + " 1 " + str(m-2) + " " + str(m-1) + " " + str(m) + "\n")
                        Nangles += 1

                for q in range(1, 4):
                    m += 1
                    if q == 1 or q == 2:
                        positions.append("\t " + str(m) + " " + str(k) + " 2 0 " + str(Vx) + " " + str(Vy) +
                                         " " + str(Vz + q * side_length) + " 0 0 0 \n")
                        bond_id += 1
                        bonds.append("\t " + str(bond_id) + " 1 " + str(m) + " " + str(m+1) + " \n")
                        Nbonds += 1
                    if q == 3:
                        positions.append("\t " + str(m) + " " + str(k) + " 1 0 " + str(Vx) + " " + str(Vy) +
                                         " " + str(Vz + q * side_length) + " 0 0 0 \n")
                        angle_id += 1
                        angles.append("\t " + str(angle_id) + " 1 " + str(m-2) + " " + str(m-1) + " " + str(m) + "\n")
                        Nangles += 1

    Natoms = m
    header = ["LAMMPS Description \n \n",
              "\t " + str(Natoms) + " atoms \n \t " + str(Nbonds) +
              " bonds \n \t " + str(Nangles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t 5 atom types \n \t 1 bond types \n \t 1 angle types \n \t 0 dihedral types \n \t 0 improper types \n",
              "\n \t " + str(0.0) + " " + str(Lx) + " xlo xhi\n \t", str(0.0) + " " + str(Ly) + " ylo yhi \n \t",
              str(0.0) + " " + str(Lz) + " zlo zhi\n", "\nMasses \n \n", "\t 1 1.0000 \n", "\t 2 1.0000 \n",
              "\t 3 1.0000 \n", "\t 4 1.0000 \n", "\t 5 1.0000 \n"]

    f = open("Kick_in.dat", "w")
    for item in header:
        f.write("%s " % item)

    for item in positions:
        f.write("%s " % item)

    for item in bonds:
        f.write("%s" % item)

    for item in angles:
        f.write("%s" % item)

    write_in_script()
