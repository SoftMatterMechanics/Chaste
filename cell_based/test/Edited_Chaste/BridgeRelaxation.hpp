/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef BridgeRelaxation_HPP_
#define BridgeRelaxation_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"
#include "FaceValueAndStressStateModifier.hpp"
#include "PolarityModifier.hpp"

#include "FakePetscSetup.hpp"

#include "PlaneBoundaryCondition.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "MyXToroidalHoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "MyNagaiHondaForceWithStripesAdhesion.hpp"
#include "DiffusionForce.hpp"
#include "PlaneBasedCellKiller.hpp"

// #include <ctime>

class BridgeRelaxation : public AbstractCellBasedTestSuite
{
public:

    void TestStripSubstrateAdhesion()
    {
        // assert(false);

        /*------------------------------START: Basic Settings----------------------------*/
        double target_shape_index = 4.0;//p0
        double reference_area = M_PI;
        double initial_area = reference_area;

        double nagai_honda_membrane_surface_energy_parameter = 0.2/(M_PI/reference_area);//Gamma
        double Km_for_myosin_feedback = 1.0; // 1.0 for defaut
        double cell_cell_adhesion_parameter = -nagai_honda_membrane_surface_energy_parameter*(target_shape_index*sqrt(reference_area));
        double Ks_for_adhesion_feedback = 0.5; // 1.0 for defaut

        double time_for_changing = 150.0; // DOUBLE_UNSET;
        double change_factor = 0.8;
        double changed_nagai_honda_membrane_surface_energy_parameter = 0.2/(M_PI/reference_area)*change_factor;//Gamma
        double changed_Km_for_myosin_feedback = 1.0*change_factor;
        double changed_cell_cell_adhesion_parameter = -changed_nagai_honda_membrane_surface_energy_parameter*(target_shape_index*sqrt(reference_area))*change_factor;
        double changed_Ks_for_adhesion_feedback = 0.5*change_factor; // 1.0 for defaut

/* Strip Structure & Cell Mesh */
        bool   strip_width_doubled_for_multiple_leading_cells = false;
        double strip_width_mutiple = 8.0;
        int    move_mesh_right_for_N_periods = 1; // for display of multiple periods
        bool   one_strip_only_in_a_period = true;

        unsigned num_ele_cross = 6; // must be even number
        unsigned num_ele_up = 16;

        double center_of_width = 0.0;       // change made by Chao
        double width = num_ele_cross*sqrt(initial_area/(sqrt(3)/2));   //width of reservoir, change made by Chao

        double strip_width = 0.5*sqrt(initial_area/(sqrt(3)/2)); // default =0.9523 (1/2 cell width)
        if  (strip_width_doubled_for_multiple_leading_cells)
            strip_width = strip_width*strip_width_mutiple;

        double strip_distance = width;

        double strip_start_x_location = 0.0;
        if  (strip_width_doubled_for_multiple_leading_cells)
            strip_start_x_location += 1/2*sqrt(initial_area/(sqrt(3)/2));
        double strip_start_y_location = 1.5*num_ele_up*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2));

        /* Energy equation form: 1/2*KA*(A-A0)^2 + 1/2*Gamma*P^2 + 1.0*Lambda*P - alpha*Sa. ( Lambda=-Ga*P0, p0=P0/sqrt(A0) ) */
/* 1. Area expansion resistance */
        double nagai_honda_deformation_energy_parameter = 1; // KA
        bool   use_fixed_target_area = true; // fix A0

/* 2. Myosin activity */
        bool   if_consider_feedback_of_face_values = true;
        double feedback_rate_for_myosin_activity = 0.1/(M_PI/reference_area);//beta
        double hill_power_for_myosin_activity = 8.0; // 8.0 for default

        // change feedback after a time.
        double changed_feedback_rate = feedback_rate_for_myosin_activity;

/* 3. Cell-Cell adhesion */
        double cell_boundary_adhesion_parameter = 0.0; // cell-cell adhesion at boundary
        
        bool   if_consider_feedback_of_cell_cell_adhesion = true;
        double feedback_rate_for_adhesion = 0.1;
        double hill_power_for_adhesion = 8.0;
        double reference_stress_for_cc_adhesion = 2.0; // sigma0.

/* 4. Substrate Ahesion */
        bool   if_consider_strip_substrate_adhesion = true;
        bool   if_consider_reservoir_substrate_adhesion = true;// true for default
        double basic_SSA = -1.0/(M_PI/reference_area);
        double reservoir_substrate_adhesion_parameter = 2.0*basic_SSA;
        double homogeneous_substrate_adhesion_parameter = 2.0*basic_SSA;
        // Strip substrate adhesion form:
        bool   consider_consistency_for_SSA = true;
        bool   if_strip_substrate_adhesion_is_homogeneous = true;
        double small_change_for_area_calculation = 0.25/sqrt((M_PI/reference_area));
        // some basic mechanisms:
        bool   classify_elements_with_group_numbers = true;
        bool   mark_leading_cells = true;
        bool   kill_cells_group_number_from1 = true;
        bool   if_check_for_T4_swaps = false;
        
/* 5. Pulling Force */
        bool   multiple_leading_cells = false;
        unsigned leading_cell_number = 1;
        if (!multiple_leading_cells)
           leading_cell_number = 1;
        double pulling_force_on_leading_cell = 12.0;// Fy
        // homogeneous SSA case:
        bool   add_pulling_force_on_node_individually = false;
        bool   add_pulling_force_evenly_on_nodes_of_leading_cell = true;
        // non-homogeneous SSA case: several leading tops
        bool   use_my_detach_pattern_method = false;
        
        // Strip substrate adhesion (SSA):
        assert( !(if_strip_substrate_adhesion_is_homogeneous && (use_my_detach_pattern_method) ) );
        assert( !(use_my_detach_pattern_method) );

/* 6. Random force */
        bool   add_random_force = true;

        // Brownian random force
        bool   has_brownian_random_force = false;
        double translational_diffusion_constant = 0.0;
        double set_node_radius = 1.0/pow(1.0,2.0)*50*1e2*((M_PI/reference_area)*(M_PI/reference_area)); // effective radius of node in Einstein relation

        // Cell polarity
        bool   has_polarity = true;
        double polarity_magnitude = 0.05;
        double polarity_magnitude_equilibrium = 0.2;
        bool   seed_manually = true;
        unsigned seed_for_initial_random_polarity = 1u;
        double rotational_diffusion_constant = 0.4/(M_PI/reference_area);

        if (polarity_magnitude==0.0)
           has_polarity = false;        
        if ((!has_brownian_random_force)&&(!has_polarity))
           add_random_force = false;
        if (polarity_magnitude!=0.0)
           assert(add_random_force == true);

/* 7. Cell division */
        bool   run_with_birth = false;
        double growth_rate_for_target_area_after_division = 0.1/(M_PI/reference_area);

/* 8. Time */
        bool   if_equilibrate_for_a_while = false;
        double end_time_for_equilibrium = 0.0;
        if (end_time_for_equilibrium <= 0.0)
           if_equilibrate_for_a_while = false;
        
        double dt = 0.025*(M_PI/reference_area); // Previously 0.025
        double end_time = 400.0*(M_PI/reference_area);
        double max_movement_per_timestep = 0.025/sqrt((M_PI/reference_area)); // Previously 0.05

        bool   apply_my_change_to_make_timestep_adaptive = true;
        bool   restrict_vertex_movement = false; // 
        if (apply_my_change_to_make_timestep_adaptive)
           assert(restrict_vertex_movement == false);

        double sampling_time = 1.0*(M_PI/reference_area);

/* 9. Cell rearrangement */
        double cell_rearrangement_threshold = 0.01/sqrt((M_PI/reference_area)); // previously: 0.05. // the minimum threshold distance for element T1 rearrangement 

/* 10. Output & display */
        bool   output_concise_swap_information_when_remesh = false;
        bool   output_detailed_swap_information_when_remesh = false;
        bool   output_numerical_method_information = false;
        bool   output_information_for_nagai_honda_force = false;     // for output and display each component in Paraview

        /*------------------------------END: Basic Settings----------------------------*/


        /*------------------------------START: Mesh Structure------------------------------*/
        // Strips structure of substrate adhesion
        bool   if_update_face_elements_in_mesh = if_consider_feedback_of_face_values;
        
        MyXToroidalHoneycombVertexMeshGenerator generator(num_ele_cross, num_ele_up, initial_area, cell_rearrangement_threshold, 0.001/sqrt((M_PI/reference_area)));
        MyXToroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        p_mesh->SetDistanceForT3SwapChecking(5.0/sqrt((M_PI/reference_area)));
        p_mesh->SetMoveMeshRightForNPeriods(move_mesh_right_for_N_periods);        
        p_mesh->SetUpdateFaceElementsInMeshBoolean(if_update_face_elements_in_mesh);
        p_mesh->SetIfClassifyElementsWithGroupNumbers(classify_elements_with_group_numbers);
        p_mesh->SetMarkLeadingCells(mark_leading_cells);
        p_mesh->SetIfCheckForT4Swaps(if_check_for_T4_swaps);
        p_mesh->SetOutputConciseSwapInformationWhenRemesh(output_concise_swap_information_when_remesh);
        p_mesh->SetOutputDetailedSwapInformationWhenRemesh(output_detailed_swap_information_when_remesh);
        /*------------------------------END: Mesh Structure------------------------------*/


        /*------------------------------START: Cells and Cell population--------------------------*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<BernoulliTrialCellCycleModel, 2> cells_generator;

        bool use_bernoulli_trial_cell_cycle_model = true;
        cells_generator.SetUseBernoulliTrialCellCycleModel(use_bernoulli_trial_cell_cycle_model);

        double divide_probability_for_a_cell = 0.0;
        double minimum_division_age = -0.01;
        cells_generator.SetDivideProbability(divide_probability_for_a_cell);
        cells_generator.SetMinimumDivisionAge(minimum_division_age);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetRestrictVertexMovementBoolean(restrict_vertex_movement);
        /*------------------------------END: Cells and Cell population--------------------------*/


        OffLatticeSimulation<2> simulator(cell_population);
        // output cell velocity:
        bool my_output_cell_velocity = true;
        simulator.SetMyOutputCellVelocities(my_output_cell_velocity);
        if (my_output_cell_velocity && seed_manually)
          simulator.SetMySeed(seed_for_initial_random_polarity);
        
        bool output_cell_velocity = true;
        simulator.SetOutputCellVelocities(output_cell_velocity);

        simulator.SetNoBirth(!run_with_birth);


        /*--------------------------------START: TargetAreaModifier------------------------------*/
        MAKE_PTR_ARGS(TargetAreaLinearGrowthModifier<2>, p_growth_modifier, ());
        p_growth_modifier->SetUseUseMyOwnRuleForUpdateTargetAreaOfCell(true);
        p_growth_modifier->SetReferenceTargetArea(reference_area);
        p_growth_modifier->SetGrowthRate(growth_rate_for_target_area_after_division);
        p_growth_modifier->SetAgeToStartGrowing(0.0);
        simulator.AddSimulationModifier(p_growth_modifier);
        /*--------------------------------END: TargetAreaModifier------------------------------*/


        /*---------------------------------START: Add Numerical Method-----------------------------*/
        boost::shared_ptr<ForwardEulerNumericalMethod<2,2>> p_numerical_method(new ForwardEulerNumericalMethod<2,2>());
        p_numerical_method->SetCellPopulation(&cell_population);
        //Note: Here the address of std::vector of force collections are different from that in simulator!
        //But the corresponding shared_ptrs in the vectors should be the same.
        std::vector< boost::shared_ptr<AbstractForce<2, 2>> > force_collection = simulator.rGetForceCollection();
        p_numerical_method->SetForceCollection(&force_collection);
        p_numerical_method->SetMaxMovementPerTimestep(max_movement_per_timestep);
        p_numerical_method->SetOutputNumericalMethodInformation(output_numerical_method_information);
        simulator.SetNumericalMethod(p_numerical_method);
        /*---------------------------------END: Add Numerical Method-----------------------------*/


        /*-----------------------------START: MyNagaiHondaForceWithStripesAdhesion---------------------*/
        MAKE_PTR(MyNagaiHondaForceWithStripesAdhesion<2>, p_force);
        
        p_force->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        p_force->SetTimeForEquilibrium(end_time_for_equilibrium);

        // structure info
        p_force->SetCenterOfWidth(center_of_width);
        p_force->SetWidth(width);//check later!
        p_force->SetStripWidth(strip_width);
        p_force->SetStripDistance(strip_distance);
        p_force->SetStripStartXLocation(strip_start_x_location);
        p_force->SetStripStartYLocation(strip_start_y_location);

        // basic parameters
        p_force->SetNagaiHondaDeformationEnergyParameter(nagai_honda_deformation_energy_parameter); // KA
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter); // Gamma
        bool if_use_face_element_to_get_adhesion_parameter = if_consider_feedback_of_face_values;
        p_force->SetUseFaceElementToGetAdhesionParameterBoolean(if_use_face_element_to_get_adhesion_parameter);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(cell_cell_adhesion_parameter); // Lambda
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(cell_boundary_adhesion_parameter);
        p_force->SetIfConsiderSubstrateAdhesion(if_consider_strip_substrate_adhesion);
        p_force->SetIfSubstrateAdhesionIsHomogeneous(if_strip_substrate_adhesion_is_homogeneous);
        p_force->SetHomogeneousSubstrateAdhesionParameter(homogeneous_substrate_adhesion_parameter);
        p_force->SetIfConsiderReservoirSubstrateAdhesion(if_consider_reservoir_substrate_adhesion);
        p_force->SetReservoirSubstrateAdhesionParameter(reservoir_substrate_adhesion_parameter);
        p_force->SetUseFixedTargetArea(use_fixed_target_area); // used in the case where there is no target area modifier! (no division)
        p_force->SetFixedTargetArea(reference_area);
        p_force->SetTargetShapeIndex(target_shape_index);
        p_force->SetConsiderConsistencyForSSA(consider_consistency_for_SSA);
        p_force->SetSmallChangeForAreaCalculation(small_change_for_area_calculation);

        // pulling force
        p_force->SetAddPullingForceOnNodeIndividually(add_pulling_force_on_node_individually);
        p_force->SetAddPullingForceEvenlyOnNodesOfLeadingCell(add_pulling_force_evenly_on_nodes_of_leading_cell);
        p_force->SetPullingForceOnLeadingCell(pulling_force_on_leading_cell);
        p_force->SetLeadingCellNumber(leading_cell_number);

        // changed basic parameters
        p_force->SetTimeForChanging(time_for_changing);
        p_force->SetChangedNagaiHondaMembraneSurfaceEnergyParameter(changed_nagai_honda_membrane_surface_energy_parameter); // Gamma
        p_force->SetChangedNagaiHondaCellCellAdhesionEnergyParameter(changed_cell_cell_adhesion_parameter); // Lambda                

        p_force->SetOutputInformationForNagaiHondaForce(output_information_for_nagai_honda_force);
        simulator.AddForce(p_force);
        /*-----------------------------END: MyNagaiHondaForceWithStripesAdhesion---------------------*/


        /*-----------------------------START: RandomForce and PolarityModifier--------------------------*/
        // Brownian diffusion
        if (add_random_force)
        {
          MAKE_PTR(DiffusionForce<2>, p_force2);
          bool   use_the_same_node_radius = true;
          double node_radius = 50*1e2*((M_PI/reference_area)*(M_PI/reference_area));
          if  (translational_diffusion_constant >0.0)
              node_radius = p_force2->GetDiffusionScalingConstant()/translational_diffusion_constant;
          else
              node_radius = set_node_radius; // 50*1e2*((M_PI/reference_area)*(M_PI/reference_area));
          translational_diffusion_constant = p_force2->GetDiffusionScalingConstant()/node_radius;

          p_force2->SetHasBrownianRandomForce(has_brownian_random_force);
          p_force2->SetUseTheSameNodeRadius(use_the_same_node_radius);
          p_force2->SetTheSameNodeRadius(node_radius);
          p_force2->SetHasPolarity(has_polarity);
          p_force2->SetOnePeriodOnly(one_strip_only_in_a_period);

          p_force2->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
          p_force2->SetTimeForEquilibrium(end_time_for_equilibrium);

          simulator.AddForce(p_force2);
        }

        // Polarity modifier
        if (has_polarity)
        {
          MAKE_PTR_ARGS(PolarityModifier<2>, p_polarity_modifier, ());
          p_polarity_modifier->SetPolarityMagnitude(polarity_magnitude);
          p_polarity_modifier->SetPolarityMagnitudeEquilibrium(polarity_magnitude_equilibrium);
          p_polarity_modifier->SetSeedManually(seed_manually);
          p_polarity_modifier->SetSeedForInitialRandomPolarity(seed_for_initial_random_polarity);
          p_polarity_modifier->SetD(rotational_diffusion_constant);
          simulator.AddSimulationModifier(p_polarity_modifier);
        }
        /*-----------------------------END: RandomForce and PolarityModifier---------------------------*/


        /*------------------------START: !!!!!Feedback: FaceValueAndStressStateModifier: need modification---------------*/
        MAKE_PTR_ARGS(FaceValueAndStressStateModifier<2>, p_face_value_and_stress_state_modifier, ());
        
        if  (feedback_rate_for_myosin_activity == 0.0)
            if_consider_feedback_of_face_values = false;
        if  (feedback_rate_for_adhesion == 0.0)
            if_consider_feedback_of_cell_cell_adhesion = false;
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValues(if_consider_feedback_of_face_values);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfCellCellAdhesion(if_consider_feedback_of_cell_cell_adhesion);

        // equilibrium related
        p_face_value_and_stress_state_modifier->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        p_face_value_and_stress_state_modifier->SetEndTimeForEquilibrium(end_time_for_equilibrium);

        // structure info
        p_face_value_and_stress_state_modifier->SetStripWidth(strip_width);
        p_face_value_and_stress_state_modifier->SetStripDistance(strip_distance); // change made by Chao
        p_face_value_and_stress_state_modifier->SetCenterOfWidth(center_of_width); // change made by Chao
        p_face_value_and_stress_state_modifier->SetWidth(width); // change made by Chao
        p_face_value_and_stress_state_modifier->SetStripStartXLocation(strip_start_x_location);
        p_face_value_and_stress_state_modifier->SetStripStartYLocation(strip_start_y_location);

        // basic parameters
        p_face_value_and_stress_state_modifier->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        p_face_value_and_stress_state_modifier->SetNagaiHondaCellCellAdhesionEnergyParameter(cell_cell_adhesion_parameter);
        p_face_value_and_stress_state_modifier->SetStripSubstrateAdhesionParameter(homogeneous_substrate_adhesion_parameter);
        p_face_value_and_stress_state_modifier->SetReservoirSubstrateAdhesionParameter(reservoir_substrate_adhesion_parameter);
        p_face_value_and_stress_state_modifier->SetFixedTargetArea(reference_area);
        p_face_value_and_stress_state_modifier->SetFixedTargetPerimeter(target_shape_index*sqrt(reference_area));
        p_face_value_and_stress_state_modifier->SetConsiderConsistencyForSSA(consider_consistency_for_SSA);
        p_face_value_and_stress_state_modifier->SetSmallChangeForAreaCalculation(small_change_for_area_calculation); // change made by Chao
        
        // feedback parameters
        p_face_value_and_stress_state_modifier->SetKmForMyosinFeedback(Km_for_myosin_feedback);
        p_face_value_and_stress_state_modifier->SetFeedbackRateForMyosinActivity(feedback_rate_for_myosin_activity);
        p_face_value_and_stress_state_modifier->SetHillPowerForMyosinActivity(hill_power_for_myosin_activity);
        p_face_value_and_stress_state_modifier->SetKsForAdhesionFeedback(Ks_for_adhesion_feedback);
        p_face_value_and_stress_state_modifier->SetFeedbackRateForAdhesion(feedback_rate_for_adhesion);
        p_face_value_and_stress_state_modifier->SetHillPowerForAdhesion(hill_power_for_adhesion);
        p_face_value_and_stress_state_modifier->SetCalculateStressStateBoolean(true);  
        p_face_value_and_stress_state_modifier->SetReferenceStress(reference_stress_for_cc_adhesion);

        // changed feedback parameters
        p_face_value_and_stress_state_modifier->SetTimeForChanging(time_for_changing);
        p_face_value_and_stress_state_modifier->SetChangedNagaiHondaMembraneSurfaceEnergyParameter(changed_nagai_honda_membrane_surface_energy_parameter);
        p_face_value_and_stress_state_modifier->SetChangedNagaiHondaCellCellAdhesionEnergyParameter(changed_cell_cell_adhesion_parameter);
        p_face_value_and_stress_state_modifier->SetChangedKmForMyosinFeedback(changed_Km_for_myosin_feedback);
        p_face_value_and_stress_state_modifier->SetChangedKsForAdhesionFeedback(changed_Ks_for_adhesion_feedback);
        p_face_value_and_stress_state_modifier->SetChangedFeedbackRate(changed_feedback_rate);

        // my group number modifier
        p_face_value_and_stress_state_modifier->SetWriteGroupNumberToCell(use_my_detach_pattern_method);

        // leading cell marker
        p_face_value_and_stress_state_modifier->SetMarkLeadingCells(mark_leading_cells);
        p_face_value_and_stress_state_modifier->SetMultipleLeadingCells(multiple_leading_cells);
        p_face_value_and_stress_state_modifier->SetLeadingCellNumber(leading_cell_number);

        // outout information
        p_face_value_and_stress_state_modifier->SetOutputModifierInformationBoolean(false);
        simulator.AddSimulationModifier(p_face_value_and_stress_state_modifier);
        /*-----------------------END: !!!!!!Feedback: FaceValueAndStressStateModifier: need modification---------------*/


        /*------------------------------------START: CellKiller---------------------------------------*/
        c_vector<double, 2> point_vec = zero_vector<double>(2);
        double kill_height = 10000.0;
        point_vec[1] = kill_height;
        c_vector<double, 2> norm_vec = zero_vector<double>(2);
        norm_vec[1] = 1.0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point_vec, norm_vec) );
        p_killer->SetKillCellsGroupNumberFrom1(kill_cells_group_number_from1);
        simulator.AddCellKiller(p_killer);
        /*------------------------------------END: CellKiller---------------------------------------*/


        /*------------------------------------START: Timestep---------------------------------------*/
        unsigned sampling_timestep_multiple = (unsigned) round(sampling_time/dt);

        simulator.SetApplyMyChangesToMakeTimestepAdaptive(apply_my_change_to_make_timestep_adaptive);
        simulator.SetDt(dt);

        simulator.SetApplySamplingTimeInsteadOfSamplingTimestep(apply_my_change_to_make_timestep_adaptive);
        simulator.SetSamplingTimestepMultiple(sampling_timestep_multiple);
        simulator.SetSamplingTime(sampling_time);

        simulator.SetEndTime(end_time);
        /*------------------------------------END: Timestep---------------------------------------*/


        /*--------------------------------START: Boundary condition-----------------------------*/
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        double stop_time = DOUBLE_UNSET;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal, stop_time));

        c_vector<double,2> point2 = zero_vector<double>(2);
        c_vector<double,2> normal2 = zero_vector<double>(2);
        point2(1) = strip_start_y_location;
        normal2(1) = 1.0;
        double stop_time2 = end_time_for_equilibrium;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2, stop_time2));

        simulator.AddCellPopulationBoundaryCondition(p_bc);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        /*--------------------------------END: Boundary condition-----------------------------*/
        

        /*--------------------------START: Output Directory and Simulation Information File---------------------*/
        // Output directory:
        std::ostringstream oss;
        time_t raw_time = time(0);
        struct tm * now = localtime(& raw_time);

        std::string output_directory = 
            "BridgeRelaxation/Date: ";
        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << '/';
        output_directory += oss.str();

        oss.str("");
        oss << "NumUp=" << num_ele_up;
        oss << "_NumCr=" << num_ele_cross;

        if (if_strip_substrate_adhesion_is_homogeneous)
        {
          if (add_pulling_force_on_node_individually)
            oss << "_Fy_lead_node=" << std::fixed << setprecision(1) << pulling_force_on_leading_cell;
          if (add_pulling_force_evenly_on_nodes_of_leading_cell)
            oss << "_Fy=" << std::fixed << setprecision(1) << pulling_force_on_leading_cell;
        }

        oss << "_p0=" << std::fixed << setprecision(2) << target_shape_index;        
        oss << "_Km=" << std::fixed << setprecision(2) << Km_for_myosin_feedback;
        oss << "_Ks=" << std::fixed << setprecision(2) << Ks_for_adhesion_feedback;
        oss << "_Ga=" << ((nagai_honda_membrane_surface_energy_parameter>=0.01 || nagai_honda_membrane_surface_energy_parameter==0.0)? std::fixed : std::scientific) 
                << setprecision(2) << nagai_honda_membrane_surface_energy_parameter;
        if (time_for_changing<end_time)
        {
          oss << "_cp0=" << -changed_cell_cell_adhesion_parameter/changed_nagai_honda_membrane_surface_energy_parameter/sqrt(reference_area);
          oss << "_cKm=" << changed_Km_for_myosin_feedback;
          oss << "_cKs=" << changed_Ks_for_adhesion_feedback;
          oss << "_cGa=" << changed_nagai_honda_membrane_surface_energy_parameter;
        }  

        oss << "_Fp=" << ((polarity_magnitude>=0.01 || polarity_magnitude==0.0)? std::fixed : std::scientific) << setprecision(2) << polarity_magnitude;
        if (polarity_magnitude!=0.0)
        {
          oss << "_Dr=" << ((rotational_diffusion_constant>=0.01)? std::fixed : std::scientific) << setprecision(2) << rotational_diffusion_constant;
          if (seed_manually)
            oss << "_PSeed=" << seed_for_initial_random_polarity;
          else
            oss << "_PSeed=N";
        }

        output_directory += oss.str();

        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday 
            << ", " << now->tm_hour << ':' << now->tm_min << ':' << now->tm_sec;
        output_directory += "_Timestamp=" + oss.str();

        std::string concise_output_directory = output_directory;
        
        // Concise information written to directoory.
        oss.str("");
        oss << ((nagai_honda_membrane_surface_energy_parameter>=0.01 || nagai_honda_membrane_surface_energy_parameter==0.0)? std::fixed : std::scientific) 
                << setprecision(2) << nagai_honda_membrane_surface_energy_parameter;
        output_directory += "_|Ga=" + oss.str();
        oss.str("");
        oss << std::fixed << setprecision(2) << target_shape_index;
        output_directory += "_p0=" + oss.str();

        // Detailed information  written to directoory.
        output_directory += "/";

        if (add_random_force)
        {
          if (has_brownian_random_force)
          {
            oss.str("");
            oss << std::scientific << setprecision(2) << translational_diffusion_constant;
            output_directory += "_|D_trans=" + oss.str();
          }
          if (has_polarity)
          {
            oss.str("");
            oss << ((polarity_magnitude>=0.01 || polarity_magnitude==0.0)? std::fixed : std::scientific) << setprecision(2) << polarity_magnitude;
            output_directory += "_fp=" + oss.str();
            oss.str("");
            oss << ((rotational_diffusion_constant>=0.01)? std::fixed : std::scientific) << setprecision(2) << rotational_diffusion_constant;
            output_directory += "_Dr=" + oss.str();
          }
        }   

        if (if_equilibrate_for_a_while)
        {
          oss.str("");
          oss << "End_time_equilib=" << std::fixed << setprecision(1) << end_time_for_equilibrium;
          output_directory += oss.str();
        }

        bool omit_file_name_results_from_time_X = true;
        bool output_simulatin_information_to_file = true;
        simulator.SetOmitFileNameResultsFromTimeX(omit_file_name_results_from_time_X);
        simulator.SetOutputSimulationInformationToFile(output_simulatin_information_to_file);
        if (output_simulatin_information_to_file)
          simulator.InputSimulationInformation(output_directory);

        bool use_concise_output_directory = true;
        if (use_concise_output_directory)
        {
          simulator.SetOutputDirectory(concise_output_directory);
          std::cout << std::endl << "Concise output directory is set: " << concise_output_directory << std::endl;
        }
        else
        {
          simulator.SetOutputDirectory(output_directory);
          std::cout << std::endl << "Output directory is set: " << output_directory << std::endl;
        }
        /*--------------------------END: Output Directory and Simulation Information File---------------------*/

        simulator.Solve();
    }

};

#endif /* BRIDGERELAXATION_HPP_ */
