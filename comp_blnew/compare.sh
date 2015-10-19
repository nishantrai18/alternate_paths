#!bin/bash

make all

./bdv_sst Romein romebdv < rome99.gr
./as_sst Romein romeas < rome99.gr
./bl_plt Romein romebl < rome99.gr
./compare romebl romebdv romeas
