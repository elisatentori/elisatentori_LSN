#!/bin/bash

mkdir simulation

a=(0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2.0)

for ((i=0;i<31; i++))
do
    cp config.final config.0
	./executable ${a[i]}
    
	mv output_ene.dat ${i}_output_ene.dat
    mv output_heat.dat ${i}_output_heat.dat
    mv output_mag.dat ${i}_output_mag.dat
    mv output_chi.dat ${i}_output_chi.dat
    
    mv ${i}_output_ene.dat simulation
    mv ${i}_output_heat.dat simulation
    mv ${i}_output_mag.dat simulation
    mv ${i}_output_chi.dat simulation
    
    cp config.final ${i}_config.final
    
    mv ${i}_config.final simulation
    rm config.0

done

# READ ME: tieni il codice cpp di questa cartella, metti nstep=1000 e cambia T

