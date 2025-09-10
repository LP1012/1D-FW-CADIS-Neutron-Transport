#pragma once

#include <stdio.h>
#include "Neutron.h"
#include "Point.h"

class MCSlab
{
public:
	/// method to run simulation
	void k_eigenvalue();

	
protected:
	/// mesh size for Shannon entropy calcs
	const unsigned int _x_cells;
	
	/// number of neutrons per generation
	const unsigned int _n_particles;
	/// number of generations
	const unsigned int _n_generations;

	/// vector of slab data: size, Sigma_a, Sigma_s, nu_Sigma_f
	const std::vector<std::vector<double>> _slabs;
	
	/// bank of source sites
	std::vector<Neutron> _fission_bank;

	/// flux at each point in mesh
	std::vector<double> _scalar_flux;

	
};
