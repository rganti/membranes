import numpy as np
import argparse


def write_file(file_object, r, force, v, label):
    file_object.write('\n')
    file_object.write(str(label) + ' \n')
    file_object.write('N {0} R {1} {2} \n\n'.format(len(r), r[0], r[len(r)-1]))
    for i in range(0, len(r)):
        f.write('{} {:06.6f} {:06.4f} {:06.4f} \n'.format(i+1, r[i], v[i], force[i]))


class Potential(object):

    def __init__(self, sigma, w_c):
        self.sigma = sigma
        self.w_c = w_c
        self.r_c = sigma*(2.0**(1.0/6.0))

    def attractive(self):

        # Potential energy and force r < self.r_c
        r_low = np.arange(0, self.r_c, 0.02)
        with np.errstate(all='ignore'):
            v_att_low = np.zeros_like(r_low) - 1
            v_rep_low = np.zeros_like(r_low) + 4.0*((self.sigma/r_low) ** 12 - (self.sigma/r_low) ** 6 + (1.0 / 4.0))
            force_low = np.zeros_like(r_low) + 4.0*(12*(self.sigma**12)/(r_low**13) - (6*(self.sigma**6)/(r_low**7)))

        # Potential energy and force r_c <= r <= r_c + w_c
        r_mid = np.arange(r_low[len(r_low)-1] + 0.02, self.r_c + self.w_c, 0.02)
        v_att_mid = np.zeros_like(r_mid) - (np.cos(np.pi*(r_mid - self.r_c)/(2.0*self.w_c)))**2
        v_rep_mid = np.zeros_like(r_mid)
        force_mid = np.zeros_like(r_mid) - np.cos(np.pi*(r_mid - self.r_c)/(2.0*self.w_c))*np.sin(np.pi*(r_mid - self.r_c)/(2.0*self.w_c))*(np.pi/self.w_c)

        # For r > r_c + w_c
        r_hi = np.arange(r_mid[len(r_mid)-1] + 0.02, 4.02, 0.02)
        v_att_hi = np.zeros_like(r_hi)
        v_rep_hi = np.zeros_like(r_hi)
        force_hi = np.zeros_like(r_hi)

        # Concatenate for full attractive forces
        r = np.append(np.append(r_low, r_mid), r_hi)
        v_attractive = np.append(np.append((v_att_low + v_rep_low), (v_att_mid + v_rep_mid)), v_att_hi+v_rep_hi)
        force_attractive = np.append(np.append(force_low, force_mid), force_hi)

        force_attractive = np.nan_to_num(force_attractive)
        v_attractive = np.nan_to_num(v_attractive)
        r[0] = 1.0e-6
        return r, force_attractive, v_attractive

    def repulsive(self):

        # For r < r_c
        r_low = np.arange(0, self.r_c, 0.02)
        with np.errstate(all='ignore'):
            v_rep_low = np.zeros_like(r_low) + 4.0*((self.sigma/r_low)**(12) - (self.sigma/r_low)**(6) + (1.0/4.0))
            force_low = np.zeros_like(r_low) + 4.0*(12*(self.sigma**12)/(r_low**13) - (6*(self.sigma**6)/(r_low**7)))

        r_hi = np.arange(r_low[len(r_low)-1] + 0.02, 4.02, 0.02)
        v_rep_hi = np.zeros_like(r_hi)
        force_hi = np.zeros_like(r_hi)

        # Concatenate for head-head and head-tail
        r = np.append(r_low, r_hi)
        v_repulsive = np.append(v_rep_low, v_rep_hi)
        force_repulsive = np.append(force_low, force_hi)

        force_repulsive = np.nan_to_num(force_repulsive)
        v_repulsive = np.nan_to_num(v_repulsive)
        r[0] = 1.0e-6
        return r, force_repulsive, v_repulsive

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script for generating tabulated potentials.")
    parser.add_argument('--wc', '-wc', dest='wc', action='store', type=float, default=1.4, help='range of attraction')
    args = parser.parse_args()

    f = open("tabulated_potential", 'w')
    f.write('# Potential Interactions \n')

    sigma = 1
    w_c = args.wc

    parameter = open("wc_parameter", 'w')
    parameter.write("w_c = {0}".format(w_c))
    parameter.close()

    # Potential Values Tail-Tail
    tail_tail = Potential(sigma, w_c)
    r_full, force_tail_tail, v_tail_tail = tail_tail.attractive()
    write_file(f, r_full, force_tail_tail, v_tail_tail, "TAIL_TAIL")

    # Monomer interactions with head and tail
    sigma = 2
    mon = Potential(sigma, w_c)
    r_full, force_mon_tail, v_mon_tail = mon.attractive()
    write_file(f, r_full, force_mon_tail, v_mon_tail, "MONOMER_TAIL")

    r_full, force_mon_head, v_mon_head = mon.repulsive()
    write_file(f, r_full, force_mon_head, v_mon_head, "MONOMER_HEAD")

    # Head - Head and Head - Tail interactions
    sigma = 0.95

    head_tail = Potential(sigma, w_c)
    r_full, force_head_tail, v_head_tail = head_tail.repulsive()
    write_file(f, r_full, force_head_tail, v_head_tail, "HEAD_TAIL")

    # Monomer-Monomer
    sigma = 2
    mon_mon = Potential(sigma, w_c)
    r_full, force_mon_mon, v_mon_mon = mon_mon.attractive()
    force_mon_mon.fill(0)
    v_mon_mon.fill(0)
    write_file(f, r_full, force_mon_mon, v_mon_mon, "MONOMER_MONOMER")

    f.close()
