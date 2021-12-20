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

#include "MyNagaiHondaForceWithStripesAdhesion.hpp"
#include "MyXToroidal2dVertexMesh.hpp"

template<unsigned DIM>
MyNagaiHondaForceWithStripesAdhesion<DIM>::MyNagaiHondaForceWithStripesAdhesion()
   : AbstractForce<DIM>(),
     mNagaiHondaDeformationEnergyParameter(1.0),
     mNagaiHondaMembraneSurfaceEnergyParameter(0.1),
     mNagaiHondaCellCellAdhesionEnergyParameter(0.0),
     mNagaiHondaCellBoundaryAdhesionEnergyParameter(0.0),
     mTimeForChanging(150.0),
     mChangedNagaiHondaMembraneSurfaceEnergyParameter(0.1),
     mChangedNagaiHondaCellCellAdhesionEnergyParameter(0.0),
     mUseFixedTargetArea(true),
     mIfUseFaceElementToGetAdhesionParameter(false),
     mOutputInformationForNagaiHondaForce(false),
     mLeadingCellNumber(1),
     mPullingForceOnLeadingCell(0.0),
     mTimeForEquilibrium(150.0)
{
}

template<unsigned DIM>
MyNagaiHondaForceWithStripesAdhesion<DIM>::~MyNagaiHondaForceWithStripesAdhesion()
{
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{    
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("NagaiHondaForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
        try
        {
            // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
            // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
            // understandable message. There is a slight chance that the exception is thrown although the error is not about the
            // target areas.
            if (mUseFixedTargetArea)
            {
                target_areas[elem_index] = mFixedTargetArea;
            }
            else
                target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");

        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use NagaiHondaForce");
        }

    }

    unsigned node_index_for_leading_top_of_the_group0 = 0;
    if (mAddPullingForceOnNodeIndividually)
        node_index_for_leading_top_of_the_group0 = p_cell_population->rGetMesh().GetNodeIndexForLeadingTopOfTheGroup(0);

    unsigned num_nodes_leading_cell_top = 0;
    if (mAddPullingForceEvenlyOnNodesOfLeadingCell)
        num_nodes_leading_cell_top = p_cell_population->rGetMesh().GetNumNodesOfLeadingCellTopOfGroup(0);

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all CellPtrs in
         * the cell population. The free energy of each CellPtr is comprised of three
         * parts - a cell deformation energy, a membrane surface tension energy
         * and an adhesion energy.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> membrane_surface_tension_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> adhesion_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> strip_substrate_adhesion_contribution = zero_vector<double>(DIM);
        // bool node_has_SSA = false;
        // c_vector<double, DIM> averaged_inner_strip_substrate_adhesion_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> reservoir_substrate_adhesion_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's deformation energy (note the minus sign)
            c_vector<double, DIM> element_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
            deformation_contribution -= GetNagaiHondaDeformationEnergyParameter()*(element_areas[elem_index] - target_areas[elem_index])*element_area_gradient;
            
            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

            // Compute the adhesion parameter for each of these edges
            // changes to be made:
            double previous_edge_adhesion_parameter = GetAdhesionParameter(p_previous_node, p_this_node, *p_cell_population);
            double next_edge_adhesion_parameter = GetAdhesionParameter(p_this_node, p_next_node, *p_cell_population);

            // Compute the gradient of each these edges, computed at the present node
            c_vector<double, DIM> previous_edge_gradient = -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
            c_vector<double, DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary adhesion (note the minus sign)
            adhesion_contribution -= previous_edge_adhesion_parameter*previous_edge_gradient + next_edge_adhesion_parameter*next_edge_gradient;
            
            // Add the force contribution from this cell's membrane surface tension (note the minus sign)
            c_vector<double, DIM> element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
            double element_myosin_activity = p_element->GetElementMyosinActivity();                
            membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*element_myosin_activity*element_perimeters[elem_index]*element_perimeter_gradient;
            
            /*---------------------------------Start of substrate adhesion contribution---------------------------------*/
            //  --B
            //    |
            // C--A
            Node<DIM>* pNodeC = p_element->GetNode((local_index-1+p_element->GetNumNodes())%p_element->GetNumNodes());
            Node<DIM>* pNodeA = p_element->GetNode(local_index);
            Node<DIM>* pNodeB = p_element->GetNode((local_index+1)%p_element->GetNumNodes());
            c_vector<double,DIM> location_c = pNodeC->rGetLocation();
            c_vector<double,DIM> location_a = pNodeA->rGetLocation();
            c_vector<double,DIM> location_b = pNodeB->rGetLocation();
            c_vector<double,DIM> vec_ca =  location_a - location_c;
            c_vector<double,DIM> vec_ab =  location_b - location_a;

            // consider periodicity: modify the node location!
            if (dynamic_cast<MyXToroidal2dVertexMesh*>(&p_cell_population->rGetMesh()) != nullptr)
            {
                bool triangle_straddles_left_right_boundary = false;

                if (fabs(vec_ca[0]) > 0.5*mWidth)
                    triangle_straddles_left_right_boundary = true;
                if (fabs(vec_ab[0]) > 0.5*mWidth)
                    triangle_straddles_left_right_boundary = true;
                if (triangle_straddles_left_right_boundary)
                {
                    if (location_c[0] < mCenterOfWidth)
                        location_c[0] += mWidth;
                    if (location_a[0] < mCenterOfWidth)
                        location_a[0] += mWidth;
                    if (location_b[0] < mCenterOfWidth)
                        location_b[0] += mWidth;
                }
            }

            double triangle_box_bottom = std::min(std::min(location_c[1],location_a[1]),location_b[1]);
            double triangle_box_top = std::max(std::max(location_c[1],location_a[1]),location_b[1]);
            double triangle_box_left = std::min(std::min(location_c[0],location_a[0]),location_b[0]);
            double triangle_box_right = std::max(std::max(location_c[0],location_a[0]),location_b[0]);

            bool full_strip_substrate_adhesion = false;
            bool partial_strip_substrate_adhesion = false;
            bool full_reservoir_substrate_adhesion = false;
            bool partial_reservoir_substrate_adhesion = false;
            bool no_substrate_adhesion = false;

            double strip_distance = this->mStripDistance;
            double strip_start_y_location = this->mStripStartYLocation;
            double strip_right = this->mStripStartXLocation + this->mStripWidth/2.0;
            double strip_left = this->mStripStartXLocation - this->mStripWidth/2.0;
            // consider periodicity. Important!!!
            if (fabs(triangle_box_left-triangle_box_right) > strip_distance/2)
                partial_strip_substrate_adhesion = false; 

            if (triangle_box_bottom >= strip_start_y_location)
            {
                if ((triangle_box_right <= strip_left) || (triangle_box_left >= strip_right))
                    no_substrate_adhesion = true;   
                else if ((triangle_box_left >= strip_left) && (triangle_box_right <= strip_right) )
                    full_strip_substrate_adhesion = true;
                else
                    partial_strip_substrate_adhesion = true;   
            }
            else if ((triangle_box_top >= strip_start_y_location) && (triangle_box_bottom <= strip_start_y_location))
            {
                if ((triangle_box_right <= strip_left) || (triangle_box_left >= strip_right))
                    partial_reservoir_substrate_adhesion = true;  
                else
                {   
                    partial_strip_substrate_adhesion = true;
                    partial_reservoir_substrate_adhesion = true;
                }
            }
            else if ((triangle_box_top <= strip_start_y_location) && (triangle_box_bottom >= 0.0))
                full_reservoir_substrate_adhesion = true;
            else if ((triangle_box_top <= strip_start_y_location) && (triangle_box_top >= 0.0) && (triangle_box_bottom <= 0.0))
                partial_reservoir_substrate_adhesion = true;
            else
                no_substrate_adhesion = true;   

            c_vector<double, DIM> strip_substrate_adhesion_area_gradient = zero_vector<double>(DIM);
            c_vector<double, DIM> reservoir_substrate_adhesion_area_gradient = zero_vector<double>(DIM);
            if (!no_substrate_adhesion)
            {
                if (full_strip_substrate_adhesion)
                    strip_substrate_adhesion_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
                else if (full_reservoir_substrate_adhesion)
                    reservoir_substrate_adhesion_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);     
                else 
                {
                    if (partial_strip_substrate_adhesion)
                        strip_substrate_adhesion_area_gradient = this->GetStripSubstrateAdhesionAreaGradientOfElementAtNode(rCellPopulation, p_element, local_index);
                    if (partial_reservoir_substrate_adhesion)
                        reservoir_substrate_adhesion_area_gradient = this->GetReservoirSubstrateAdhesionAreaGradientOfElementAtNode(rCellPopulation, p_element, local_index);
                }
            }

            if (mIfConsiderSubstrateAdhesion)
                strip_substrate_adhesion_contribution -= mHomogeneousSubstrateAdhesionParameter*strip_substrate_adhesion_area_gradient;

            if (mIfConsiderReservoirSubstrateAdhesion)
                reservoir_substrate_adhesion_contribution -= mReservoirSubstrateAdhesionParameter*reservoir_substrate_adhesion_area_gradient;

            /*---------------------------------End of substrate adhesion contribution---------------------------------*/
        }// end of 'Iterate over these elements'

        c_vector<double, DIM> force_on_node = deformation_contribution + membrane_surface_tension_contribution + adhesion_contribution + strip_substrate_adhesion_contribution + reservoir_substrate_adhesion_contribution;
        
        if (mAddPullingForceOnNodeIndividually && p_this_node->GetIndex()==node_index_for_leading_top_of_the_group0)
            force_on_node[1] += mPullingForceOnLeadingCell;
        else if (mAddPullingForceEvenlyOnNodesOfLeadingCell)
        {
            bool node_belongs_to_leading_cell = false;
            for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
                iter != containing_elem_indices.end();
                ++iter)
            {
                if (p_cell_population->GetElement(*iter)->GetIsLeadingCellTop())
                    node_belongs_to_leading_cell = true;
            }

            if (node_belongs_to_leading_cell)
            {
                assert(num_nodes_leading_cell_top!=0);
                double t_now = SimulationTime::Instance()->GetTime();
                if ( !(mIfEquilibrateForAWhile && t_now<=mTimeForEquilibrium) )
                {
                    if (mLeadingCellNumber==1)
                        force_on_node[1] += mPullingForceOnLeadingCell/num_nodes_leading_cell_top;
                    else
                    {
                        assert(mLeadingCellNumber >= 2);
                        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
                            iter != containing_elem_indices.end();
                            ++iter)
                        {
                            if (p_cell_population->GetElement(*iter)->GetIsLeadingCellTop())
                            {
                                unsigned num_nodes_of_this_leading_cell = p_cell_population->GetElement(*iter)->GetNumNodes();
                                force_on_node[1] += (mPullingForceOnLeadingCell/mLeadingCellNumber)/num_nodes_of_this_leading_cell;
                            }
                        }

                    }
                }
            }
        }

        // tmp: possilbe output for testing new SSA:
        bool output_while_testing_new_SSA = false;
        if (output_while_testing_new_SSA)
        {
            VertexElement<DIM, DIM>* pElement = p_cell_population->GetElement(* p_this_node->rGetContainingElementIndices().begin());
            if (pElement->GetIsLeadingCell() && pElement->GetGroupNumber()>0)
            {
                double t_now = SimulationTime::Instance()->GetTime();
                if (p_this_node->rGetContainingElementIndices().size()==1)
                {
                    std::cout << "Time: " << t_now << std::endl;
                    std::cout << "Node of LeadingCell, node index=" << p_this_node->GetIndex() 
                            << ", SSA case number=0, node location: x=" << p_this_node->rGetLocation()[0] << ", y=" << p_this_node->rGetLocation()[1]
                            << std::endl << "Lamellipodium strength=" << pElement->GetLamellipodiumStrength()
                            << std::endl << "SSA: x=" << strip_substrate_adhesion_contribution[0] << ", y=" << strip_substrate_adhesion_contribution[1] << std::endl;
                }
                else if (p_this_node->rGetContainingElementIndices().size()==2)
                {
                    std::cout << "Time: " << t_now << std::endl;
                    std::cout << "Node of LeadingCell, node index=" << p_this_node->GetIndex() 
                            << ", SSA case number=1, node location: x=" << p_this_node->rGetLocation()[0] << ", y=" << p_this_node->rGetLocation()[1]
                            << std::endl << "Lamellipodium strength1=" << pElement->GetLamellipodiumStrength() 
                            << ", Lamellipodium strength2=" << p_cell_population->GetElement(* p_this_node->rGetContainingElementIndices().begin()++)->GetLamellipodiumStrength()
                            << std::endl << "SSA: x=" << strip_substrate_adhesion_contribution[0] << ", y=" << strip_substrate_adhesion_contribution[1] << std::endl;
                }
            }
        }

        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);

        // tmp
        if (mOutputInformationForNagaiHondaForce)
        {
            if (norm_2(force_on_node)>2)
            {
                std::cout << "Weird! Force is too large! Node Index: " << p_this_node->GetIndex() << std::endl;
                std::cout << "Node location: " << p_this_node->rGetLocation()[0] << ", " << p_this_node->rGetLocation()[1] << std::endl;
                std::cout << "X direction:" << std::endl;
                std::cout << "force on node=" << force_on_node[0]
                        << " deformation_contribution=" << deformation_contribution[0]
                        << " membrane_surface_tension_contribution=" << membrane_surface_tension_contribution[0]
                        << " adhesion_contribution=" << adhesion_contribution[0]
                        << " strip_substrate_adhesion_contribution=" << strip_substrate_adhesion_contribution[0]
                        << " reservoir_substrate_adhesion_contribution=" << reservoir_substrate_adhesion_contribution[0] << std::endl;
                std::cout << "Y direction:" << std::endl;
                std::cout << "force on node=" << force_on_node[1]
                        << " deformation_contribution=" << deformation_contribution[1]
                        << " membrane_surface_tension_contribution=" << membrane_surface_tension_contribution[1]
                        << " adhesion_contribution=" << adhesion_contribution[1]
                        << " strip_substrate_adhesion_contribution=" << strip_substrate_adhesion_contribution[1]
                        << " reservoir_substrate_adhesion_contribution=" << reservoir_substrate_adhesion_contribution[1] << std::endl;
            }
        }

    }// end of 'Iterate over nodes(vertices) in the cell population'

}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    double adhesion_parameter = 0.0;

    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        adhesion_parameter = GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
        return adhesion_parameter;
    }

    // my changes
    if (mIfUseFaceElementToGetAdhesionParameter)
    {
        unsigned element_global_index = *shared_elements.begin();
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rVertexCellPopulation);
        VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(element_global_index);
        if ((p_element->GetNodeLocalIndex(pNodeB->GetIndex()) - p_element->GetNodeLocalIndex(pNodeA->GetIndex())+p_element->GetNumNodes())%p_element->GetNumNodes()==1)
        {
            VertexElement<DIM-1, DIM>* pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeA->GetIndex(), pNodeB->GetIndex()));
            adhesion_parameter =GetNagaiHondaCellCellAdhesionEnergyParameter() * pFace->GetUnifiedCellCellAdhesionEnergyParameter();
        }
        else
        {
            if ( (p_element->GetNodeLocalIndex(pNodeA->GetIndex()) - p_element->GetNodeLocalIndex(pNodeB->GetIndex())+p_element->GetNumNodes())%p_element->GetNumNodes()!=1)
            {
                std::cout << std::endl << "ERROR: Method MyNagaiHondaForceWithStripesAdhesion::GetAdhesionParameter";
                std::cout << std::endl << "Not Reachable!" << std::endl;
            }
            assert((p_element->GetNodeLocalIndex(pNodeA->GetIndex()) - p_element->GetNodeLocalIndex(pNodeB->GetIndex())+p_element->GetNumNodes())%p_element->GetNumNodes()==1);
            VertexElement<DIM-1, DIM>* pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeB->GetIndex(), pNodeA->GetIndex()));
            adhesion_parameter = GetNagaiHondaCellCellAdhesionEnergyParameter() * pFace->GetUnifiedCellCellAdhesionEnergyParameter();    
        }
    }
    else
        adhesion_parameter = GetNagaiHondaCellCellAdhesionEnergyParameter();

    return adhesion_parameter;
}

template <unsigned DIM>
c_vector<double, DIM> MyNagaiHondaForceWithStripesAdhesion<DIM>::GetStripSubstrateAdhesionAreaGradientOfElementAtNode(AbstractCellPopulation<DIM>& rCellPopulation, VertexElement<DIM, DIM>* pElement, unsigned localIndex)
{
    assert(DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

//    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    c_vector<double, DIM> strip_substrate_adhesion_area_gradient = zero_vector<double>(DIM);

    // Parameters
    double strip_width = this->mStripWidth;
    double strip_distance = this->mStripDistance;
    double strip_start_x_location = this->mStripStartXLocation;
    double strip_start_y_location = this->mStripStartYLocation;

    double small_change = this->mSmallChangeForAreaCalculation;
    double preferred_sample_dis = small_change/5.0;

    ///\todo This should probably be returning the nearest node
    c_vector<double, DIM> centroid = zero_vector<double>(DIM);
    unsigned num_nodes_elem = pElement->GetNumNodes();
    for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
    {
        // Find location of current node and add it to the centroid
        centroid += pElement->GetNodeLocation(local_index);
    }
    centroid /= num_nodes_elem;

    int strip_num = round((centroid[0] - strip_start_x_location)/strip_distance);

    double strip_location = strip_start_x_location + strip_distance*strip_num;
    double strip_left = strip_location - strip_width/2;
    double strip_right = strip_location + strip_width/2;

    double sample_area_bottom = 0.0;
    double sample_area_top = 0.0;
    double sample_area_left = 0.0;
    double sample_area_right = 0.0;

    Node<DIM>* p_this_node = pElement->GetNode(localIndex);
    unsigned previous_node_local_index = (num_nodes_elem+localIndex-1)%num_nodes_elem;
    Node<DIM>* p_previous_node = pElement->GetNode(previous_node_local_index);
    unsigned next_node_local_index = (localIndex+1)%num_nodes_elem;
    Node<DIM>* p_next_node = pElement->GetNode(next_node_local_index);

    c_vector<double, DIM> previous_node_location = p_previous_node->rGetLocation();
    c_vector<double, DIM> this_node_location = p_this_node->rGetLocation();
    c_vector<double, DIM> next_node_location = p_next_node->rGetLocation();
    
    double expanded_triangle_box_bottom = -small_change+std::min(std::min(previous_node_location[1],this_node_location[1]),next_node_location[1]);
    double expanded_triangle_box_top = small_change+std::max(std::max(previous_node_location[1],this_node_location[1]),next_node_location[1]);
    double expanded_triangle_box_left = -small_change+std::min(std::min(previous_node_location[0],this_node_location[0]),next_node_location[0]);
    double expanded_triangle_box_right = small_change+std::max(std::max(previous_node_location[0],this_node_location[0]),next_node_location[0]);
    
    bool has_substrate_adhesion_area = true;
    // consider periodicity. Important!!!
    if (fabs(expanded_triangle_box_left-expanded_triangle_box_right) > strip_distance/2)
        has_substrate_adhesion_area = false;

    if (expanded_triangle_box_top<strip_start_y_location)
        has_substrate_adhesion_area = false;
    else
    {
        sample_area_top = expanded_triangle_box_top;
        sample_area_bottom = expanded_triangle_box_bottom>strip_start_y_location? expanded_triangle_box_bottom : strip_start_y_location;
    }

    if ((expanded_triangle_box_left <= strip_left)&&(expanded_triangle_box_right >= strip_left)&&(expanded_triangle_box_right <= strip_right))
    {   
        sample_area_left = strip_left;
        sample_area_right = expanded_triangle_box_right;
    }
    else if ((expanded_triangle_box_left <= strip_left)&&(expanded_triangle_box_right >= strip_right))
    {   
        sample_area_left = strip_left;
        sample_area_right = strip_right;
    }
    else if ((expanded_triangle_box_left >= strip_left)&&(expanded_triangle_box_right <= strip_right))
    {
        sample_area_left = expanded_triangle_box_left;
        sample_area_right = expanded_triangle_box_right;
    }
    else if ((expanded_triangle_box_left >= strip_left)&&(expanded_triangle_box_left <= strip_right)&&(expanded_triangle_box_right >= strip_right))
    {   
        sample_area_left = expanded_triangle_box_left;
        sample_area_right = strip_right;
    }
    else
        has_substrate_adhesion_area = false;

    unsigned num_across = 0;
    unsigned num_up = 0;
    unsigned sample_num = 0;
    double   sample_area = 0.0;

    if (has_substrate_adhesion_area)
    {
        num_across = round((sample_area_right - sample_area_left) / preferred_sample_dis);
        num_up = round((sample_area_top - sample_area_bottom) / preferred_sample_dis);
        sample_num = num_across * num_up;
        sample_area = (sample_area_top - sample_area_bottom)*(sample_area_right - sample_area_left);
    }
    
    // Calculate strip_substrate_adhesion_area_gradient!
    if (has_substrate_adhesion_area && (num_across>0) && (num_up>0))
    {
        c_vector<c_vector<double, DIM>, 3> points;
        points[0] = previous_node_location;
        points[1] = this_node_location;
        points[2] = next_node_location;

        // Calculate initial adhesive area
        double adhesive_sample_num = 0.0;
        for (unsigned i = 0; i <sample_num; i++)
        {
            double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
            double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

            c_vector<double, DIM> point;
            point[0] = x_coord;
            point[1] = y_coord;
            c_vector<bool, 3> point_at_left_of_vector;
            for (unsigned j = 0; j < 3; j++)
            {
                c_vector<double, DIM> vec1 = point - points[j];
                c_vector<double, DIM> vec2 = points[(j+1)%3] - points[j];

                if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                    point_at_left_of_vector[j] = false;
                else
                    point_at_left_of_vector[j] = true;
            }
            if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                adhesive_sample_num += 1.0;
            else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                adhesive_sample_num += -1.0;
        }
        double substrate_adhesion_area = adhesive_sample_num/double(sample_num) * sample_area;

        // Calculate adhesive area *changed* with small displacement of the node along the x or y axis
        for (unsigned j = 0; j<2; j++)
        {
            c_vector<double, DIM> unit_vector_in_small_change_direction;
            unit_vector_in_small_change_direction[0] = cos(j*M_PI/2);
            unit_vector_in_small_change_direction[1] = sin(j*M_PI/2);
            // my new changes for cosistent movement of nodes at edges of the strip.
            if (j==0 && this->mConsiderConsistencyForSSA)
            {
                if ((points[1][j]-strip_right)>-small_change && (points[1][j]-strip_right)<small_change)
                    unit_vector_in_small_change_direction[0] *= -1.0;
            }

            c_vector<c_vector<double, DIM>, 3> new_points = points;
            new_points[1] += small_change*unit_vector_in_small_change_direction;

            double adhesive_sample_num = 0.0;

            for (unsigned i = 0; i<sample_num; i++)
            {
                double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
                double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

                c_vector<double, DIM> point;
                point[0] = x_coord;
                point[1] = y_coord;
                c_vector<bool, 3> point_at_left_of_vector;
                for (unsigned j = 0; j < 3; j++)
                {
                    c_vector<double, DIM> vec1 = point - new_points[j];
                    c_vector<double, DIM> vec2 = new_points[(j+1)%3] - new_points[j];

                    if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                        point_at_left_of_vector[j] = false;
                    else
                        point_at_left_of_vector[j] = true;
                }
                if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                    adhesive_sample_num += 1.0;
                else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                    adhesive_sample_num += -1.0;
            }
            double substrate_adhesion_area_new = adhesive_sample_num/double(sample_num) * sample_area;
            
            strip_substrate_adhesion_area_gradient += (substrate_adhesion_area_new - substrate_adhesion_area)/small_change*unit_vector_in_small_change_direction; 
        }
    }
    return strip_substrate_adhesion_area_gradient;
}

template <unsigned DIM>
c_vector<double, DIM> MyNagaiHondaForceWithStripesAdhesion<DIM>::GetReservoirSubstrateAdhesionAreaGradientOfElementAtNode(AbstractCellPopulation<DIM>& rCellPopulation, VertexElement<DIM,DIM>* pElement, unsigned localIndex)
{
    assert(DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    c_vector<double, DIM> reservoir_substrate_adhesion_area_gradient = zero_vector<double>(DIM);

    // Parameters
    double strip_start_y_location = this->mStripStartYLocation;
    double small_change = this->mSmallChangeForAreaCalculation;
    double preferred_sample_dis = small_change/5.0;

    unsigned num_nodes_elem = pElement->GetNumNodes();

    double sample_area_bottom = 0.0;
    double sample_area_top = 0.0;
    double sample_area_left = 0.0;
    double sample_area_right = 0.0;

    Node<DIM>* p_this_node = pElement->GetNode(localIndex);
    unsigned previous_node_local_index = (num_nodes_elem+localIndex-1)%num_nodes_elem;
    Node<DIM>* p_previous_node = pElement->GetNode(previous_node_local_index);
    unsigned next_node_local_index = (localIndex+1)%num_nodes_elem;
    Node<DIM>* p_next_node = pElement->GetNode(next_node_local_index);

    c_vector<double, DIM> previous_node_location = p_previous_node->rGetLocation();
    c_vector<double, DIM> this_node_location = p_this_node->rGetLocation();
    c_vector<double, DIM> next_node_location = p_next_node->rGetLocation();
    
    // consider periodicity: modify the node location!
    if (dynamic_cast<MyXToroidal2dVertexMesh*>(&p_cell_population->rGetMesh()) != nullptr)
    {
        bool triangle_straddles_left_right_boundary = false;

        c_vector<double, 2> vector1 = this_node_location - previous_node_location;
        if (fabs(vector1[0]) > 0.5*mWidth)
            triangle_straddles_left_right_boundary = true;
        c_vector<double, 2> vector2 = next_node_location - this_node_location;
        if (fabs(vector2[0]) > 0.5*mWidth)
            triangle_straddles_left_right_boundary = true;
        if (triangle_straddles_left_right_boundary)
        {
            if (previous_node_location[0] < mCenterOfWidth)
                previous_node_location[0] += mWidth;
            if (this_node_location[0] < mCenterOfWidth)
                this_node_location[0] += mWidth;
            if (next_node_location[0] < mCenterOfWidth)
                next_node_location[0] += mWidth;
        }
    }
    
    double expanded_triangle_box_bottom = -small_change+std::min(std::min(previous_node_location[1],this_node_location[1]),next_node_location[1]);
    double expanded_triangle_box_top = small_change+std::max(std::max(previous_node_location[1],this_node_location[1]),next_node_location[1]);
    double expanded_triangle_box_left = -small_change+std::min(std::min(previous_node_location[0],this_node_location[0]),next_node_location[0]);
    double expanded_triangle_box_right = small_change+std::max(std::max(previous_node_location[0],this_node_location[0]),next_node_location[0]);

    bool has_substrate_adhesion_area = true;
    if ((expanded_triangle_box_top >= strip_start_y_location)&&(expanded_triangle_box_bottom >= 0.0)&&(expanded_triangle_box_bottom <= strip_start_y_location))
    {
        sample_area_top = strip_start_y_location;
        sample_area_bottom = expanded_triangle_box_bottom;
    }
    else if ((expanded_triangle_box_top <= strip_start_y_location)&&(expanded_triangle_box_bottom >= 0.0)) // note: we get errors here previously.
    {
        sample_area_top = expanded_triangle_box_top;
        sample_area_bottom = expanded_triangle_box_bottom;
    }
    else if ((expanded_triangle_box_top <= strip_start_y_location)&&(expanded_triangle_box_top >= 0.0)&&(expanded_triangle_box_bottom <= 0.0))
    {
        sample_area_top = expanded_triangle_box_top;
        sample_area_bottom = 0.0;
    }     
    else if ((expanded_triangle_box_top >= strip_start_y_location)&&(expanded_triangle_box_bottom <= 0.0))                   
    {
        sample_area_top = strip_start_y_location;
        sample_area_bottom = 0.0;
    }
    else
        has_substrate_adhesion_area = false;

    sample_area_left = expanded_triangle_box_left;
    sample_area_right = expanded_triangle_box_right;

    unsigned num_across = 0;
    unsigned num_up = 0;
    unsigned sample_num = 0;
    double sample_area = 0.0;
    if (has_substrate_adhesion_area)
    {
        num_across = (unsigned)round((sample_area_right - sample_area_left) / preferred_sample_dis);
        num_up = (unsigned)round((sample_area_top - sample_area_bottom) / preferred_sample_dis);
        sample_num = num_across * num_up;
        sample_area = (sample_area_top - sample_area_bottom)*(sample_area_right - sample_area_left);
    }

    // Calculate reservoir_substrate_adhesion_area_gradient!
    if (has_substrate_adhesion_area && (num_across>0) && (num_up>0))
    {
        c_vector<c_vector<double, DIM>, 3> points;
        points[0] = previous_node_location;
        points[1] = this_node_location;
        points[2] = next_node_location;

        // Calculate initial adhesive area
        double adhesive_sample_num = 0.0;
        for (unsigned i = 0; i<sample_num; i++)
        {
            double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
            double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

            c_vector<double, DIM> point;
            point[0] = x_coord;
            point[1] = y_coord;
            c_vector<bool, 3> point_at_left_of_vector;
            for (unsigned j = 0; j < 3; j++)
            {
                c_vector<double, DIM> vec1 = point - points[j];
                c_vector<double, DIM> vec2 = points[(j+1)%3] - points[j];

                if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                    point_at_left_of_vector[j] = false;
                else
                    point_at_left_of_vector[j] = true;
            }
            if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                adhesive_sample_num += 1.0;
            else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                adhesive_sample_num += -1.0;
        }
        double substrate_adhesion_area = adhesive_sample_num/double(sample_num) * sample_area;

        // Calculate adhesive area *changed* with small displacement of the node along the x or y axis
        for (unsigned j = 0; j<2; j++)
        {
            c_vector<double, DIM> unit_vector_in_small_change_direction;
            unit_vector_in_small_change_direction[0] = cos(j*M_PI/2);
            unit_vector_in_small_change_direction[1] = sin(j*M_PI/2);

            c_vector<c_vector<double, DIM>, 3> new_points = points;
            new_points[1] += small_change*unit_vector_in_small_change_direction;

            double adhesive_sample_num = 0.0;

            for (unsigned i = 0; i <sample_num; i++)
            {
                double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
                double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

                c_vector<double, DIM> point;
                point[0] = x_coord;
                point[1] = y_coord;
                c_vector<bool, 3> point_at_left_of_vector;
                for (unsigned j = 0; j < 3; j++)
                {
                    c_vector<double, DIM> vec1 = point - new_points[j];
                    c_vector<double, DIM> vec2 = new_points[(j+1)%3] - new_points[j];

                    if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                        point_at_left_of_vector[j] = false;
                    else
                        point_at_left_of_vector[j] = true;
                }
                if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                    adhesive_sample_num += 1.0;
                else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                    adhesive_sample_num += -1.0;
            }
            double substrate_adhesion_area_new = adhesive_sample_num/double(sample_num) * sample_area;

            reservoir_substrate_adhesion_area_gradient += (substrate_adhesion_area_new - substrate_adhesion_area)/small_change*unit_vector_in_small_change_direction;
        }// end of calculate adhesive area *changed* with small displacement of the node along the x or y axis  
    }// end of statement 'if (has_substrate_adhesion_area)'

    return reservoir_substrate_adhesion_area_gradient;
}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaDeformationEnergyParameter()
{
    return mNagaiHondaDeformationEnergyParameter;
}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaMembraneSurfaceEnergyParameter()
{
    double time_now = SimulationTime::Instance()->GetTime();
    double Ga = mNagaiHondaMembraneSurfaceEnergyParameter;
    if (time_now > mTimeForChanging)
    {
        Ga =  mChangedNagaiHondaMembraneSurfaceEnergyParameter;
    }    

    return Ga;
}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaCellCellAdhesionEnergyParameter()
{
    double time_now = SimulationTime::Instance()->GetTime();
    double Lambda = mNagaiHondaCellCellAdhesionEnergyParameter;
    if (time_now>mTimeForChanging)
    {
        Lambda = mChangedNagaiHondaCellCellAdhesionEnergyParameter;
    }    
        
    return Lambda;
}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaCellBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaDeformationEnergyParameter(double deformationEnergyParameter)
{
    mNagaiHondaDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter)
{
    mNagaiHondaMembraneSurfaceEnergyParameter = membraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mNagaiHondaCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetChangedNagaiHondaMembraneSurfaceEnergyParameter(double changedMembraneSurfaceEnergyParameter)
{
    mChangedNagaiHondaMembraneSurfaceEnergyParameter = changedMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetChangedNagaiHondaCellCellAdhesionEnergyParameter(double changedCellCellAdhesionEnergyParameter)
{
    mChangedNagaiHondaCellCellAdhesionEnergyParameter = changedCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NagaiHondaDeformationEnergyParameter>" << mNagaiHondaDeformationEnergyParameter << "</NagaiHondaDeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaMembraneSurfaceEnergyParameter>" << mNagaiHondaMembraneSurfaceEnergyParameter << "</NagaiHondaMembraneSurfaceEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellCellAdhesionEnergyParameter>" << mNagaiHondaCellCellAdhesionEnergyParameter << "</NagaiHondaCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellBoundaryAdhesionEnergyParameter>" << mNagaiHondaCellBoundaryAdhesionEnergyParameter << "</NagaiHondaCellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MyNagaiHondaForceWithStripesAdhesion<1>;
template class MyNagaiHondaForceWithStripesAdhesion<2>;
template class MyNagaiHondaForceWithStripesAdhesion<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyNagaiHondaForceWithStripesAdhesion)
