Apply Redshift space distorsions to a catalogue in a box along the x-axis

usage: 
	./apply_rsd catalog Lbox z Omega_m output

- catalog: ascii file with columns: X, Y, Z, Vx ...
- LBox: size of the box in Mpc/h
- z: Redshift of the snapshot
- Omega_m: cosmological parameter Omega Matter 
- output: name of the desired output file, wil only write: X Y Z
