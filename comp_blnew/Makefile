all : as_sstplot bdv_sstplot bl_plot compare

as_sstplot : as_sstplotter.cpp
	g++ -O2 as_sstplotter.cpp -o as_sst

bdv_sstplot : bdv_sstplotter.cpp
	g++ -O2 bdv_sstplotter.cpp -o bdv_sst
	
bl_plot : bl_plotter.cpp
	g++ -O2 bl_plotter.cpp -o bl_plt

compare : compare.cpp
	g++ compare.cpp -o compare
