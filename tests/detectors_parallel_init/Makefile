default: input

input: clean
	gmsh -2 -o src/square.msh src/square.geo
	cp src/square.msh .

clean: clean-mesh clean-run
clean-mesh:
	rm -f *.ele *.edge *.face *.node *.halo *.pvtu *.vtu *.detectors *.detectors.dat *.stat
	rm -f gmsh.log gmsh_err.log
clean-run:
	rm -f matrixdump matrixdump.info
	rm -f fluidity.err-* fluidity.log-*
	rm -rf detectors_[0-4]* *flredecomp*
