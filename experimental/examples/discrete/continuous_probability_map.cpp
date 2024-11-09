/***************************************************************************
 *            continuous_probability_map.cpp
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

#include "verification/reach_avoid.hpp"
using namespace Ariadne;

Grid make_grid_from_lengths(RealSpace const& spc, List<Pair<RealVariable,Real>> const& lengths) {
    Vector<Real> result(spc.dimension());
    {
        SizeType idx = 0;
        for (auto const& v : spc.variables()) {
            for (auto const& gl : lengths) {
                if (gl.first == v) {
                    result[idx] = gl.second;
                    ++idx;
                    break;
                }
            }
        }
    }
    return Grid(result);
}

GridTreePaving make_paving(Grid const& grid, UpperBoxType const& domain) {
    GridTreePaving result(grid);
    result.adjoin_outer_approximation(domain,0);
    result.mince(0);

    return result;
}

void ariadne_main()
{
    //! Setup

    RealConstant v("v",3);
    RealVariable x("x"), y("y"), theta("theta"), u("u");

    RealSpace state_spc({x,y,theta});
    RealSpace control_spc({u});

    Real step_size = 0.1_x;

    List<DottedRealAssignment> forward_dynamics = {{dot(x)=v*cos(theta),dot(y)=v*sin(theta),dot(theta)=u,dot(u)=0}};

    FloatDP eps(1e-8_x,DoublePrecision());
    SizeType point_accuracy = 6;
    SizeType preimage_iterations = 1;

    List<Pair<RealVariable,Real>> state_grid_lengths({{x,0.5_dec},{y,0.5_dec},{theta,2*pi/8}});
    List<Pair<RealVariable,Real>> control_grid_lengths({{u,pi/4*10}});

    UpperBoxType state_domain({{0.0_dec+eps,5_dec-eps},{0.0_dec+eps,5_dec-eps}, {-2*pi+eps,2*pi-eps}});
    UpperBoxType control_domain({{-10*pi+eps,10*pi-eps}});

    //! Processing

    RealSpace full_spc = state_spc;
    full_spc.adjoin(control_spc);

    Grid state_grid = make_grid_from_lengths(state_spc,state_grid_lengths);
    Grid control_grid = make_grid_from_lengths(control_spc,control_grid_lengths);

    GridTreePaving state_paving = make_paving(state_grid,state_domain);
    GridTreePaving control_paving = make_paving(control_grid,control_domain);

    IdentifiedCellFactory::HashTableType state_cells_ids;
    auto default_state_cell_extent = state_paving.begin()->root_extent();
    CONCLOG_PRINTLN_VAR(default_state_cell_extent)
    CONCLOG_PRINTLN_VAR(state_paving.size())
    for (auto const& c : state_paving) {
        state_cells_ids.insert(make_pair(word_to_id(c.word(),(default_state_cell_extent)*state_paving.dimension()),state_cells_ids.size()));
    }
    auto state_icell_factory = IdentifiedCellFactory(default_state_cell_extent,state_cells_ids);


    IdentifiedCellFactory::HashTableType control_cells_ids;
    auto default_control_cell_extent = control_paving.begin()->root_extent();
    CONCLOG_PRINTLN_VAR(default_control_cell_extent)
    CONCLOG_PRINTLN_VAR(control_paving.size())
    for (auto const& c : control_paving) {
        control_cells_ids.insert(make_pair(word_to_id(c.word(),(default_control_cell_extent)*control_paving.dimension()),control_cells_ids.size()));
    }
    auto control_icell_factory = IdentifiedCellFactory(default_control_cell_extent,control_cells_ids);

    List<Pair<IdentifiedCell,IdentifiedCell>> state_control_icells;

    for (auto const& state_cell : state_paving)
        for (auto const& control_cell : control_paving)
            state_control_icells.push_back(make_pair(state_icell_factory.create(state_cell), control_icell_factory.create(control_cell)));

    std::atomic<double> sum_xratio = 0;
    std::atomic<double> sum_effective_accuracy = 0;
    std::atomic<SizeType> num_processed = 0;

    TaylorPicardIntegrator integrator(1e-4);
    VectorFieldEvolver forward_evolver(VectorField(forward_dynamics),integrator);
    forward_evolver.configuration().set_maximum_step_size(step_size/10);

    List<DottedRealAssignment> backward_dynamics;
    for (auto const& d : forward_dynamics)
        backward_dynamics.push_back({d.left_hand_side(),-d.right_hand_side()});
    VectorFieldEvolver backward_evolver(VectorField(backward_dynamics),integrator);
    backward_evolver.configuration().set_maximum_step_size(step_size/10);

    BetterThreads::StaticWorkload<Pair<IdentifiedCell,IdentifiedCell>> workload([&](Pair<IdentifiedCell,IdentifiedCell> const& sc_icells){

        Stopwatch<Microseconds> sw;

        Map<SizeType,double> interval_probabilities;

        auto source_icell = sc_icells.first;
        auto source_box = source_icell.cell().box();
        auto control_icell = sc_icells.second;
        auto control_box = control_icell.cell().box();

        auto source_volume = volume(source_box).get_d();

        CONCLOG_PRINTLN_AT(2,"Computing for cell " << source_icell << " using Interval Arithmetic...")

        CONCLOG_RUN_MUTED(auto forward_orbit = forward_evolver.orbit(forward_evolver.enclosure(product(source_box,control_box)),step_size, Semantics::UPPER))

        auto final_forward_orbit_box = forward_orbit.final()[0].euclidean_set().bounding_box();
        auto final_forward_box = source_box;
        for (SizeType i=0; i<final_forward_box.dimension(); ++i)
            final_forward_box[i] = cast_exact_interval(final_forward_orbit_box[i]);

        CONCLOG_PRINTLN_VAR_AT(2,final_forward_box)

        GridTreePaving final_paving(state_grid);
        final_paving.adjoin_outer_approximation(final_forward_box,0);
        final_paving.mince(0);

        if (final_paving.subset(state_paving)) {

            for (auto const& cell : final_paving) {
                auto icell = state_icell_factory.create(cell);
                auto image_intersection_box = intersection(cell.box(),final_forward_box);

                CONCLOG_RUN_MUTED(auto backward_orbit = backward_evolver.orbit(backward_evolver.enclosure(product(image_intersection_box,control_box)),step_size, Semantics::UPPER))

                auto final_backward_orbit_box = backward_orbit.final()[0].euclidean_set().bounding_box();
                auto final_backward_box = source_box;
                for (SizeType i=0; i<final_backward_box.dimension(); ++i)
                    final_backward_box[i] = cast_exact_interval(final_backward_orbit_box[i]);

                auto preimage_intersection_box = intersection(source_box,final_backward_box);

                CONCLOG_PRINTLN_VAR_AT(2,preimage_intersection_box)

                if (definitely(not preimage_intersection_box.is_empty()))
                    interval_probabilities.insert(icell.id(),volume(preimage_intersection_box).get_d()/source_volume);
            }

            CONCLOG_PRINTLN_AT(2,"Interval probabilities:")
            for (auto& p : interval_probabilities) {
                CONCLOG_PRINTLN_AT(2,p.first << ": " << p.second)
            }

            sw.click();
            CONCLOG_PRINTLN_AT(2,"Done in " << sw.elapsed_seconds()*1000 << " ms.")
/*
            SizeType ival_xtime = static_cast<SizeType>(sw.elapsed_seconds()*1000000);

            Vector<Map<SizeType,double>> point_probabilities(point_accuracy+1,Map<SizeType,double>());

            Vector<SizeType> pt_xtime(point_accuracy+1);

            for (SizeType acc=0; acc<=point_accuracy;++acc) {
                CONCLOG_PRINTLN_VAR_AT(2,acc)
                sw.restart();
                GridTreePaving sampling(state_grid);
                sampling.adjoin_outer_approximation(shrink(sc_boxes.first,eps),acc);

                GridTreePaving intersected_sampling = intersection(state_paving,sampling);
                intersected_sampling.mince(acc);

                Map<SizeType,SizeType> occurrencies;
                SizeType total_occurrencies = 0;

                for (auto const& c : sampling) {
                    auto input_box = product(ExactBoxType(c.box().midpoint()), sc_boxes.second);

                    auto final_box = cast_exact_box(apply(forward_dynamics, input_box).bounding_box());

                    GridTreePaving target(state_grid);
                    target.adjoin_outer_approximation(final_box,0);
                    target.mince(0);

                    for (auto const& tcell : target) {
                        auto itcell = identified_cell_factory.create(tcell);
                        if (occurrencies.has_key(itcell.id()))
                            occurrencies[itcell.id()]++;
                        else
                            occurrencies.insert(itcell.id(),1);
                        ++total_occurrencies;
                    }
                }

                for (auto const& occ : occurrencies) {
                    point_probabilities.at(acc).insert(occ.first,static_cast<double>(occ.second)/total_occurrencies);
                }

                sw.click();
                CONCLOG_PRINTLN_AT(2,"Done in " << sw.elapsed_seconds() << " seconds.")
                pt_xtime.at(acc) = static_cast<SizeType>(sw.elapsed_seconds()*1000000);

                for (auto const& pp : point_probabilities.at(acc)) {
                    CONCLOG_PRINTLN_AT(2,pp.first << ": " << pp.second)
                }
            }

            double interval_deviation = 0;
            {
                for (auto const& p : interval_probabilities) {
                    auto val = point_probabilities.at(point_accuracy).find(p.first);
                    if (val == point_probabilities.at(point_accuracy).end())
                        interval_deviation += p.second*p.second;
                    else
                        interval_deviation += pow(p.second-val->second,2);
                }
                interval_deviation = sqrt(interval_deviation/point_probabilities.at(point_accuracy).size());
                CONCLOG_PRINTLN_AT(2,"Interval deviation: " << interval_deviation)
            }

            SizeType reference_accuracy = 0;

            Vector<double> point_deviations(point_accuracy,0.0);
            bool found_reference = false;
            for (SizeType acc=0;acc<point_accuracy;++acc) {
                CONCLOG_PRINTLN_AT(2,"Accuracy " << acc)
                for (auto const& p : point_probabilities.at(point_accuracy)) {
                    auto val = point_probabilities.at(acc).find(p.first);
                    if (val == point_probabilities.at(acc).end())
                        point_deviations.at(acc) += p.second*p.second;
                    else
                        point_deviations.at(acc) += pow(p.second-val->second,2);
                }
                point_deviations.at(acc) = sqrt(point_deviations.at(acc)/point_probabilities.at(point_accuracy).size());
                CONCLOG_PRINTLN_AT(2,"Deviation: " << point_deviations.at(acc))
                if ((not found_reference) and acc > 0 and point_deviations.at(acc) < interval_deviation) {
                    found_reference = true;
                    reference_accuracy = acc-1;
                }
            }

            auto effective_pt_xtime = static_cast<double>(pt_xtime.at(reference_accuracy));
            auto effective_accuracy = static_cast<double>(reference_accuracy);
            if (point_deviations.at(reference_accuracy) > interval_deviation) {
                auto log_ref_pt_dev = log(point_deviations.at(reference_accuracy));
                auto log_next_pt_dev = log(point_deviations.at(reference_accuracy+1));
                auto log_ivl_dev = log(interval_deviation);
                auto fraction = (log_ref_pt_dev-log_ivl_dev)/(log_ref_pt_dev-log_next_pt_dev);
                effective_pt_xtime = static_cast<double>(pt_xtime.at(reference_accuracy))+fraction*(static_cast<double>(pt_xtime.at(reference_accuracy+1)-pt_xtime.at(reference_accuracy)));
                effective_accuracy = static_cast<double>(reference_accuracy)+fraction;
            }

            auto xratio = effective_pt_xtime/static_cast<double>(ival_xtime);

            if (interval_deviation != 0) {
                CONCLOG_PRINTLN_AT(1,"Execution time ratio: " << xratio)
                ++num_processed;
                sum_xratio = sum_xratio + xratio;
                sum_effective_accuracy = sum_effective_accuracy + effective_accuracy;
                CONCLOG_SCOPE_PRINTHOLD("avg ratio: " << sum_xratio/num_processed << ", avg effective accuracy: " << sum_effective_accuracy/num_processed)
            }*/
        } else {
            CONCLOG_PRINTLN_AT(1,"Some targets are outside the domain, skipping")
        }
    });
    
    workload.append(state_control_icells).process();

    CONCLOG_PRINTLN("Avg ratio: " << sum_xratio/num_processed << ", avg effective accuracy: " << sum_effective_accuracy/num_processed)
}