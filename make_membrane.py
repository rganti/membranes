import numpy as np


def write_in_script():
    filename = "in.kick"
    f = open(filename, "w")
    f.write("units           lj\natom_style      full\nread_data       Kick_in.dat\n")
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
    # f.write("pair_coeff      1 5 tabulated_potential MONOMER_HEAD\n")
    f.write("pair_coeff      2 3 tabulated_potential MONOMER_TAIL\n")
    f.write("pair_coeff      2 4 tabulated_potential MONOMER_TAIL\n")
    # f.write("pair_coeff      2 5 tabulated_potential MONOMER_TAIL\n")
    f.write("pair_coeff      3 3 tabulated_potential MONOMER_MONOMER\n")
    f.write("pair_coeff      3 4 tabulated_potential MONOMER_MONOMER\n")
    # f.write("pair_coeff      3 5 tabulated_potential MONOMER_MONOMER\n")
    f.write("pair_coeff      4 4 tabulated_potential MONOMER_MONOMER\n")
    # f.write("pair_coeff      4 5 tabulated_potential MONOMER_MONOMER\n")
    # f.write("pair_coeff      5 5 tabulated_potential MONOMER_MONOMER\n")
    f.write("pair_modify        shift yes\n\n")

    f.write("# Additional variables to check\n")
    f.write("group           monomer type 3 4\n")

    f.write("group           head type 1\n")
    f.write("group           tail type 2\n")
    f.write("group           bilayer union head tail\n")
    f.write("variable        top equal bound(monomer,zmax)\n")
    f.write("variable        bottom equal bound(monomer,zmin)\n")
    f.write("variable        Lz      equal zhi-zlo\n")
    f.write("neigh_modify    exclude molecule monomer\n\n")
    f.write("# Fixes\n")

    f.write("velocity        all create 1.0 1 \n")

    # f.write("fix             MCSWP all atom/swap 10 1 1234 1 types 4 5 \n")

    f.write("dump            id all xyz 100 kick_out.xyz\n")
    f.write("restart         10000 restart.dat\n")
    # f.write("fix             0 all rigid/nve molecule\n")

    f.write("fix             fNPT bilayer npt temp 1.0 1.0 1 x 0.0 0.0 10 y 0.0 0.0 10 couple xy\n")
    # f.write("fix             fMon monomer rigid/npt molecule temp 1.0 1.0 1 x 0.0 0.0 10 y 0.0 0.0 10 couple xy\n")
    # f.write("fix             fMon monomer rigid/nvt molecule temp 1.0 1.0 5\n\n")
    f.write("fix             fMon monomer rigid/nve molecule\n\n")
    f.write("thermo          3000\n")
    f.write("thermo_style    custom step temp press etotal epair vol v_top v_bottom v_Lz\n")
    f.write("thermo_modify   flush yes\n")
    f.write("timestep        0.008\n")
    f.write("run             15000000\n")
    f.close()


class MakeMembrane(object):

    def __init__(self):
        self.side_length = 1.1
        self.Lx_index = self.Ly_index = 20
        self.Lz = 3 * self.side_length * self.Lx_index
        self.Lx = self.Lx_index * self.side_length
        self.Ly = self.Ly_index * self.side_length

        self.positions = ["\nAtoms \n \n"]
        self.bonds = ["\nBonds \n \n"]
        self.angles = ["\nAngles \n \n"]

        # m is number of atoms and k is number of 6 atom-molecules
        self.m = 0
        self.k = 0

        self.Nbonds = 0
        self.bond_id = 0

        self.Nangles = 0
        self.angle_id = 0

    def make_monomer(self, Vx, Vy, Vz):
        for q in range(0, 2):
            self.m += 1
            if q == 0:
                self.positions.append("\t " + str(self.m) + " " + str(self.k) + " 3 0 " + str(Vx) + " "
                                      + str(Vy) + " " + str(Vz + (q - 2) * self.side_length) + " 0 0 0 \n")
            if q == 1:
                self.positions.append("\t " + str(self.m) + " " + str(self.k) + " 4 0 " + str(Vx) + " "
                                      + str(Vy) + " " + str(Vz + (q - 2) * self.side_length) + " 0 0 0 \n")
        for q in range(1, 3):
            self.m += 1
            if q == 1:
                self.positions.append("\t " + str(self.m) + " " + str(self.k) + " 4 0 " + str(Vx) + " "
                                      + str(Vy) + " " + str(Vz + q * self.side_length) + " 0 0 0 \n")
            if q == 2:
                self.positions.append("\t " + str(self.m) + " " + str(self.k) + " 3 0 " + str(Vx) + " "
                                      + str(Vy) + " " + str(Vz + q * self.side_length) + " 0 0 0 \n")

    def monomer_group(self, min, max, j, s, Vx, Vy, Vz):
        if (j == min + 1 and s == min + 1) or (j == max - 1 and s == max - 1) or \
                (j == min + 1 and s == max - 1) or (j == max - 1 and s == min + 1):
            self.make_monomer(Vx, Vy, Vz)

    # def make_monomer_group(self, j, s, Vx, Vy, Vz):
    #     min = 26
    #     max = 31
    #     if (max > j > min) and (max > s > min):
    #         self.monomer_group(min, max, j, s, Vx, Vy, Vz)
    #
    #     if (max-2 > j > min+1) and (min > s > min-2):
    #         j_center = np.mean([max-2, min+1])
    #         print("j_ctr = " + str(j_center))
    #         s_center = np.mean([min, min-2])
    #         print("s_ctr = " + str(s_center))
    #         if j == j_center and s == s_center:
    #             self.make_monomer(Vx, Vy, Vz)
    #
    #     if (max > j > max-2) and (max+2 > s > max):
    #         j_center = np.mean([max, max-2])
    #         s_center = np.mean([max+2, max])
    #         if j == j_center and s == s_center:
    #             self.make_monomer(Vx, Vy, Vz)

    def make_membrane(self):
        for j in range(self.Ly_index):
            for s in range(self.Lx_index):
                self.k += 1
                Vx = s * self.side_length
                Vy = j * self.side_length
                Vz = .5 * self.Lz
                min1 = 6
                max1 = 11
                # min2 = 26
                # max2 = 31
                # min3 = 49
                # max3 = 54

                # if (max1 >= j >= min1) and (max1 >= s >= min1):
                #     self.monomer_group(min1, max1, j, s, Vx, Vy, Vz)
                #
                if (max1 - 2 >= j >= min1 + 1) and (min1 >= s >= min1 - 2):
                    j_center = np.mean([max1 - 2, min1 + 1])
                    s_center = np.mean([min1, min1 - 2])
                    if j == j_center and s == s_center:
                        self.make_monomer(Vx, Vy, Vz)
                #
                # elif (max1-1 >= j >= min1 + 2) and (max1 + 2 >= s >= max1):
                #     j_center = np.mean([max1 - 1, min1 + 2])
                #     s_center = np.mean([max1 + 2, max1])
                #     if j == j_center and s == s_center:
                #         self.make_monomer(Vx, Vy, Vz)

                elif (max1 - 1 >= j >= min1 + 2) and (max1 + 2 >= s >= max1):
                    j_center = np.mean([max1 - 1, min1 + 2])
                    s_center = np.mean([max1 + 2, max1])
                    if j == j_center and s == s_center:
                        self.make_monomer(Vx, Vy, Vz)

                # elif (max2 >= j >= min2) and (max2 >= s >= min2):
                #     self.monomer_group(min2, max2, j, s, Vx, Vy, Vz)
                #
                # elif (max2 - 2 >= j >= min2 + 1) and (min2 >= s >= min2 - 2):
                #     j_center = np.mean([max2 - 2, min2 + 1])
                #     s_center = np.mean([min2, min2 - 2])
                #     if j == j_center and s == s_center:
                #         self.make_monomer(Vx, Vy, Vz)
                #
                # elif (max2-1 >= j >= min2 + 2) and (max2 + 2 >= s >= max2):
                #     j_center = np.mean([max2 - 1, min2 + 2])
                #     s_center = np.mean([max2 + 2, max2])
                #     if j == j_center and s == s_center:
                #         self.make_monomer(Vx, Vy, Vz)
                #
                # elif (max3 >= j >= min3) and (max3 >= s >= min3):
                #     self.monomer_group(min3, max3, j, s, Vx, Vy, Vz)
                #
                # elif (max3 - 2 >= j >= min3 + 1) and (min3 >= s >= min3 - 2):
                #     j_center = np.mean([max3 - 2, min3 + 1])
                #     s_center = np.mean([min3, min3 - 2])
                #     if j == j_center and s == s_center:
                #         self.make_monomer(Vx, Vy, Vz)
                #
                # elif (max3 - 1 >= j >= min3 + 2) and (max3 + 2 >= s >= max3):
                #     j_center = np.mean([max3 - 1, min3 + 2])
                #     s_center = np.mean([max3 + 2, max3])
                #     if j == j_center and s == s_center:
                #         self.make_monomer(Vx, Vy, Vz)

                else:
                    for q in range(0, 3):
                        self.m += 1
                        if q == 0:
                            self.positions.append("\t " + str(self.m) + " " + str(self.k) + " 1 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Vz + (q - 2) * self.side_length) + " 0 0 0 \n")
                            self.bond_id += 1
                            self.bonds.append("\t " + str(self.bond_id) + " 1 " + str(self.m) + " " + str(self.m + 1) + " \n")
                            self.Nbonds += 1
                        if q == 1:
                            self.positions.append("\t " + str(self.m) + " " + str(self.k) + " 2 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Vz + (q - 2) * self.side_length) + " 0 0 0 \n")
                            self.bond_id += 1
                            self.bonds.append("\t " + str(self.bond_id) + " 1 " + str(self.m) + " " + str(self.m + 1) + " \n")
                            self.Nbonds += 1
                        if q == 2:
                            self.positions.append("\t " + str(self.m) + " " + str(self.k) + " 2 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Vz + (q - 2) * self.side_length) + " 0 0 0 \n")
                            self.angle_id += 1
                            self.angles.append(
                                "\t " + str(self.angle_id) + " 1 " + str(self.m - 2) + " " + str(self.m - 1) + " " + str(self.m) + "\n")
                            self.Nangles += 1

                    for q in range(1, 4):
                        self.m += 1
                        if q == 1 or q == 2:
                            self.positions.append("\t " + str(self.m) + " " + str(self.k) + " 2 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Vz + q * self.side_length) + " 0 0 0 \n")
                            self.bond_id += 1
                            self.bonds.append("\t " + str(self.bond_id) + " 1 " + str(self.m) + " " + str(self.m + 1) + " \n")
                            self.Nbonds += 1
                        if q == 3:
                            self.positions.append("\t " + str(self.m) + " " + str(self.k) + " 1 0 " + str(Vx) + " " + str(Vy) +
                                             " " + str(Vz + q * self.side_length) + " 0 0 0 \n")
                            self.angle_id += 1
                            self.angles.append(
                                "\t " + str(self.angle_id) + " 1 " + str(self.m - 2) + " " + str(self.m - 1) + " " + str(self.m) + "\n")
                            self.Nangles += 1

if __name__ == "__main__":

    membrane = MakeMembrane()
    membrane.make_membrane()

    header = ["LAMMPS Description \n \n",
              "\t " + str(membrane.m) + " atoms \n \t " + str(membrane.Nbonds) +
              " bonds \n \t " + str(membrane.Nangles) + " angles \n \t 0 dihedrals \n \t 0 impropers \n",
              "\n \t 4 atom types \n \t 1 bond types \n \t 1 angle types \n \t 0 dihedral types \n \t 0 improper types \n",
              "\n \t " + str(0.0) + " " + str(membrane.Lx) + " xlo xhi\n \t", str(0.0) + " " + str(membrane.Ly) + " ylo yhi \n \t",
              str(0.0) + " " + str(membrane.Lz) + " zlo zhi\n", "\nMasses \n \n", "\t 1 1.0000 \n", "\t 2 1.0000 \n",
              "\t 3 1.0000 \n", "\t 4 1.0000 \n"]  # "\t 5 1.0000 \n"]

    f = open("Kick_in.dat", "w")
    for item in header:
        f.write("%s " % item)

    for item in membrane.positions:
        f.write("%s " % item)

    for item in membrane.bonds:
        f.write("%s" % item)

    for item in membrane.angles:
        f.write("%s" % item)

    write_in_script()
