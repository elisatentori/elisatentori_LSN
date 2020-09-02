#!/bin/bash

mkdir thermalisation


for ((i=1;i<20; i++))
do
    cp config.final config.0
    cp old.final old.0
	./therm.exe

	mv output_etot.dat ${i}_etot.dat
    mv output_temp.dat ${i}_temp.dat
    mv output_ekin.dat ${i}_ekin.dat
    mv output_epot.dat ${i}_epot.dat
    
    mv ${i}_etot.dat thermalisation
    mv ${i}_epot.dat thermalisation
    mv ${i}_ekin.dat thermalisation
    mv ${i}_temp.dat thermalisation
    
    cp config.final ${i}_config.final
    cp old.final ${i}_old.final
    
    mv ${i}_config.final thermalisation
    mv ${i}_old.final thermalisation
    
    #mkdir thermalisation/${i}_frames
    #mv frames/* thermalisation/${i}_frames

done

# READ ME: tieni il codice cpp di questa cartella, metti nstep=1000 e cambia T

