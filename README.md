# membranes
Membrane model as described in Cooke et. al (2008) Efficient tunable generic model for fluid bilayer membranes

To create lammps crystal configuration Kick_in.dat and lmp file in.kick, execute the following commands

mkdir Simulation

cd Simulation/

python ../make_membrane.py

To create tabulated potential that is used by in.kick, run

python ../tabulated_energies.py --wc (Set wc anywhere from 0.8 - 1.8) {Default argument is 1.4, should be stable at T=1, see Figure 1 in Cooke et al.}

This will also create a file wc_parameter that stores the current value of wc that the tabulated potential uses. To run lammps file,

lammps_executable -in in.kick

To see plots of the potential interactions, a jupyter notebook file Membrane.ipynb runs the code written in make_membrane.py in each cell and then plots the potential and the forces. You are welcome to change the parameters and see how the potential changes. 
To run the notebook, install and run jupyter as explained in http://jupyter.readthedocs.io/en/latest/install.html and then launch by clicking on Membrane.ipynb

