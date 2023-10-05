/***************************************************************************
 *            probability_map.cpp
 *
 *  Copyright  2023  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "ariadne_main.hpp"
#include "utility/stopwatch.hpp"

#include "reach_avoid.hpp"

using namespace Ariadne;

void ariadne_main()
{
    RealVariable x("x"),y("y"),theta("theta"),u("u");
    RealConstant vel("vel",3);

    double probability_threshold = 1e-3;

    RealSpace state_spc({x,y,theta});
    RealSpace control_spc({u});

    List<Pair<RealVariable,Real>> state_grid_lengths({{x,0.5_dec},{y,0.5_dec},{theta,pi/4}});
    List<Pair<RealVariable,Real>> control_grid_lengths({{u,pi/4}});

    Vector<Real> real_state_grid_lengths(state_spc.dimension());
    {
        SizeType idx = 0;
        for (auto const& v : state_spc.variables()) {
            for (auto const& gl : state_grid_lengths) {
                if (gl.first == v) {
                    real_state_grid_lengths[idx] = gl.second;
                    ++idx;
                    break;
                }
            }
        }
    }
    Grid state_grid(real_state_grid_lengths);

    GridTreePaving domain_paving(state_grid);
    UpperBoxType state_domain({{0.01_dec,4.99_dec},{0.01_dec,4.99_dec}, {0.01_dec,2*pi-0.01_dec}});
    domain_paving.adjoin_outer_approximation(state_domain,0);

    RealSpace full_spc = join(state_spc,control_spc);

    List<DottedRealAssignment> assignments({dot(x)=vel*cos(theta),
                                                   dot(y)=vel*sin(theta),
                                                   dot(theta)=u,
                                                   dot(u)=0
                                           });

    List<DottedRealAssignment> inverted_assignments;
    for (auto const& a : assignments) {
        inverted_assignments.push_back(DottedRealAssignment(a.lhs,-a.rhs));
    }

    RealVariablesBox initial({0.0_dec<=x<=0.5_dec,0.0_dec<=y<=0.5_dec,0.0_dec<=theta<=pi/4,0.0_dec<=u<=pi/4});

    StepMaximumError max_err=1e-6;
    TaylorPicardIntegrator integrator(max_err);

    VectorFieldEvolver evolver(assignments,integrator);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    evolver.configuration().set_maximum_step_size(0.1);
    evolver.configuration().set_maximum_spacial_error(1e-4);

    CONCLOG_PRINTLN_VAR(initial);
    Real evolution_time(0.1_dec);

    Stopwatch<Milliseconds> sw;
    CONCLOG_PRINTLN("Computing using Interval Arithmetic...")
    auto forward_orbit = evolver.orbit(initial,evolution_time,Semantics::UPPER);

    auto final_box = forward_orbit.final().bounding_box();

    UpperBoxType final_projected_box(3,UpperIntervalType::empty_interval());
    {
    SizeType idx = 0;
    for (auto const& v : final_box.space().variables()) {
        if (state_spc.contains(v)) {
            final_projected_box[idx] = final_box[v];
            ++idx;
        }
    }
    }

    GridTreePaving paving(state_grid);
    paving.adjoin_outer_approximation(final_projected_box,0);
    CONCLOG_PRINTLN(paving.size())

    CONCLOG_PRINTLN(state_grid)
    CONCLOG_PRINTLN(final_projected_box)

    IdentifiedCellFactory::HashTableType state_cells_ids;

    auto default_cell_extent = paving.begin()->root_extent();
    CONCLOG_PRINTLN_VAR(default_cell_extent);
    CONCLOG_PRINTLN_VAR(paving.size());
    for (auto const& c : paving) {
        state_cells_ids.insert(make_pair(word_to_id(c.word(),(default_cell_extent)*paving.dimension()),state_cells_ids.size()));
    }

    auto identified_cell_factory = IdentifiedCellFactory(default_cell_extent,state_cells_ids);

    Map<SizeType,double> interval_probabilities;
    double total_volume = 0;
    for (auto const& cell : paving) {
        auto icell = identified_cell_factory.create(cell);
        auto starting_box = intersection(cell.box(),final_projected_box);
        auto current_volume = starting_box.volume().get_d();
        interval_probabilities.insert(icell.id(),current_volume);
        total_volume += current_volume;
    }

    Set<SizeType> to_remove;
    for (auto& p : interval_probabilities) {
        p.second = p.second/total_volume;
        if (p.second < probability_threshold) to_remove.insert(p.first);
    }
    interval_probabilities.remove_keys(to_remove);
    CONCLOG_PRINTLN("Interval probabilities:")
    for (auto& p : interval_probabilities) {
        CONCLOG_PRINTLN(p.first << ": " << p.second)
    }

    sw.click();
    CONCLOG_PRINTLN("Done in " << sw.elapsed_seconds() << " seconds.")
/*
    sw.restart();

    SizeType accuracy = 4;
    UpperBoxType initial_ubx({{0.0_dec,0.5_dec},{0.0_dec,0.5_dec}, {0.0_dec,pi/4_dec}});

    for (SizeType i=0; i<accuracy;++i) {
        GridTreePaving sampling(state_grid);
        sampling.adjoin_outer_approximation(initial_ubx,i);

        //GridTreePaving intersected_sampling = intersection()

        Map<SizeType,SizeType> occurrencies;
        for (auto const& c : sampling) {
            auto starting_pt = c.box().centre();

            RealBox real_starting_pt(full_spc.dimension());
            for (SizeType i=0; i<state_spc.size();++i) {
                real_starting_pt[i] = {cast_exact(starting_pt[i].get_d()),cast_exact(starting_pt[i].get_d())};
            }
            for (SizeType i=state_spc.size(); i<full_spc.size();++i) {
                real_starting_pt[i] = {cast_exact(initial_ubx[i].lower_bound().get_d()),cast_exact(initial_ubx[i].upper_bound().get_d())};
            }
            RealVariablesBox starting_variables_pt(full_spc,real_starting_pt);

            auto orbit = evolver.orbit(starting_variables_pt,evolution_time,Semantics::UPPER);
            auto final_bx = orbit.final().bounding_box().continuous_set();

            UpperBoxType final_projected_bx(state_spc.dimension());
            for (SizeType i=0; i<state_spc.size();++i) {
                final_projected_bx[i] = {cast_exact(final_bx[i].lower_bound().get_d()),cast_exact(final_bx[i].upper_bound().get_d())};
            }

            GridTreePaving target(state_grid);
            target.adjoin_outer_approximation(final_projected_bx,0);
            for (auto const& tcell : target) {
                auto itcell = identified_cell_factory.create(tcell);
                if (occurrencies.has_key(itcell.id()))
                    occurrencies[itcell.id()]++;
                else
                    occurrencies.insert(itcell.id(),1);
            }
        }
        for (auto const& occ : occurrencies) {
            CONCLOG_PRINTLN(occ.first << ": " << occ.second)
        }
    }
*/
    /*
    RealBox real_starting_box(full_spc.dimension());
    for (SizeType i=0; i<state_spc.size();++i) {
        real_starting_box[i] = {cast_exact(starting_box[i].lower_bound().get_d()),cast_exact(starting_box[i].upper_bound().get_d())};
    }
    for (SizeType i=state_spc.size(); i<full_spc.size();++i) {
        real_starting_box[i] = {cast_exact(final_box.continuous_set()[i].lower_bound().get_d()),cast_exact(final_box.continuous_set()[i].upper_bound().get_d())};
    }
    RealVariablesBox initial(full_spc,real_starting_box);
    */
}
