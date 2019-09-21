/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./embed.h"
#include "../modules/PhysiCell_settings.h"

// declare cell definitions here
Cell_Definition embed_cell;   // embedded
Cell_Definition env_cell;     // enveloping

void create_cell_types( void )
{
	// housekeeping 
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// turn the default cycle model to live, 
	// so it's easier to turn off proliferation
	
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;  
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// use default proliferation and death 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 
	
	// set default uptake and secretion 
	static int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0
	
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_ID] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_ID] = 10; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_ID] = 38; 

	// set the default cell type to no phenotype updates 
	// cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_oncoprotein; 
	cell_defaults.functions.update_phenotype = switch_celltype_check; 

	// not motile 
	cell_defaults.phenotype.motility.is_motile = false; 

	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
	cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = 5.0;

	// add custom data 
	// cell_defaults.custom_data.add_variable( "oncoprotein" , "dimensionless", 1.0 ); 


	
	//-------  cell type 0 ------------
	embed_cell = cell_defaults;
	embed_cell.name = "embed cell"; 
	embed_cell.type = 0; 
	embed_cell.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0; 
	embed_cell.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.01; 
	
	//-------  cell type 1 ------------
	env_cell = cell_defaults;
	env_cell.name = "envelop cell"; 
	env_cell.type = 1; 
	env_cell.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0; 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters  -- now done in config .xml

	// Ensure a 2D model
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	initialize_microenvironment(); 	
	return; 
}	

void setup_tissue( void )
{
	double x,y,z;
	int cell_type;

	Cell* pC;
	
	// std::string cell_file = parameters.strings( "cell_file" ); 
	std::string cell_file = "cells1552.dat"; 
	std::cout << "-------- reading " << cell_file << std::endl;

	std::ifstream infile(cell_file);
	// while ((infile >> x >> sep >> y >> sep >> z >> sep >> cell_type ) && (sep == ' '))
	while (infile >> x >> y >> z >> cell_type ) 
	{
//		std::cout << "cell type " << cell_type <<":  x,y= "<< x << "," << y ;

		if (cell_type == 0)
			pC = create_cell(embed_cell); 
		else
			pC = create_cell(env_cell); 

		pC->assign_position( x,y, 0.0 );
	}
}

void switch_celltype_check( Cell* pCell, Phenotype& phenotype, double dt )
{
	// std::cout << __FUNCTION__ << "------- " << std::endl;
	float t_thresh = 120.0;
	if ((PhysiCell_globals.current_time > t_thresh)&&(PhysiCell_globals.current_time < t_thresh+0.1))
	{
		if (pCell->type == 0)
		{
			std::cout << __FUNCTION__ << "------- switch cell type: time=" << PhysiCell_globals.current_time << std::endl;
			// void Cell::convert_to_cell_definition( Cell_Definition& cd )
			pCell->convert_to_cell_definition(env_cell);
			// pCell->type = 1;
		}
	}
	return; 
}
	
// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt )
{
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// if cell is dead, don't bother with future phenotype changes. 
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

	// multiply proliferation rate by the oncoprotein 
	
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 

	phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ; 
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector<std::string> output = { "green" , "black" , "cyan", "black" }; 
	
	if( pCell->type == 0 )  // embedded
	{ 
		output[0] = "cyan"; 
		output[2] = "gray"; 
	} 
	else if( pCell->type == 1 )  // enveloping
	{ 
		output[0] = "red"; 
		output[2] = "gray"; 
	} 
	return output; 
}
