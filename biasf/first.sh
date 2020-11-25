#!/bin/bash


echo Calculates the power spectrum of 76 files. Each one contains X, Y, Z information of halos with a certain mass. Then it compares all power spectra with the theoretical power spectrum and then is extracted for each one of the files a particular bias.


cd positions_for_each_mass_bin/
for i in *.txt  
do
  echo "Looping ... i is set to $i"
  python3 ../bias_bernhard.py ../positions_for_each_mass_bin/$i > ../bias/bias_$i
  
done

