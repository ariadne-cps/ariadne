/***************************************************************************
 *            reachability_graph.cpp
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

#include "reachability_graph.hpp"

#include "geometry/grid_cell.hpp"
#include "geometry/grid_paving.hpp"

#include "conclog/logging.hpp"
#include "conclog/progress_indicator.hpp"

using namespace ConcLog;

namespace Ariadne {

String word_to_id(BinaryWord const& w, SizeType length) {
    std::stringstream ss;
    SizeType effective_size = (w.size()<length? w.size() : length);
    auto size_offset = (w.size()<length? 0 : w.size()-length);
    for (SizeType i=0; i<effective_size; ++i) ss << to_string(w.at(size_offset+i));
    return ss.str();
}

IdentifiedCellFactory::IdentifiedCellFactory(SizeType default_extent, HashTableType const& table) : _default_extent(default_extent), _table(table) { }

IdentifiedCell IdentifiedCellFactory::create(GridCell const& cell) const { return IdentifiedCell(_to_identifier(cell),cell); }

SizeType IdentifiedCellFactory::_to_identifier(GridCell const& cell) const {
    SizeType size_to_use = (cell.root_extent() >= _default_extent ? cell.word().size() - (cell.root_extent() - _default_extent)*cell.dimension() : cell.word().size());
    String expanded_id = word_to_id(cell.word(),size_to_use);

    if (cell.root_extent() < _default_extent) {

        SizeType temporary_root_extent = cell.root_extent();
        while (temporary_root_extent < _default_extent) {
            expanded_id = (*expanded_id.begin() == '1' ? std::string(cell.dimension(),'0') : std::string(cell.dimension(),'1')) + expanded_id;
            ++temporary_root_extent;
        }
    }
    return _table.at(expanded_id);
}

DirectedHashedGraph::DirectedHashedGraph(IdentifiedCellFactory const& vertex_factory, IdentifiedCellFactory const& edge_factory) :
    _vertex_factory(vertex_factory), _edge_factory(edge_factory) { }

DirectedHashedGraph::DirectedHashedGraph(DirectedHashedGraph const& other) :
    _vertex_factory(other._vertex_factory), _edge_factory(other._edge_factory) {

    _map = TransitionScoreMap();
    for (auto const& s : other._map) {
        _map.insert(s.first,Map<IdentifiedCell, Map<IdentifiedCell, TargetScore>>());
        auto& s_it = _map[s.first];
        for (auto const& c : s.second) {
            s_it.insert(c.first,Map<IdentifiedCell, TargetScore>());
            auto& c_it = s_it[c.first];
            for (auto const& d : c.second) {
                c_it.insert(d.first,d.second);
            }
        }
    }
}

SizeType DirectedHashedGraph::num_sources() const { return _map.size(); }

SizeType DirectedHashedGraph::num_transitions() const {
    SizeType result = 0;
    for (auto const& src : _map) {
        for (auto const& trans : src.second) {
            result += trans.second.size();
        }
    }
    return result;
}

IdentifiedCell DirectedHashedGraph::vertex_icell(NCell const& cell) const {
    return _vertex_factory.create(cell);
}

IdentifiedCell DirectedHashedGraph::edge_icell(NCell const& cell) const {
    return _edge_factory.create(cell);
}

SizeType DirectedHashedGraph::vertex_id(NCell const& cell) const {
    return _vertex_factory.create(cell).id();
}

SizeType DirectedHashedGraph::edge_id(NCell const& cell) const {
    return _edge_factory.create(cell).id();
}

SPaving DirectedHashedGraph::destinations_from(NCell const& source_cell) const {
    SPaving result(source_cell.grid());
    auto isrc = _vertex_factory.create(source_cell);
    auto src_ref = _map.find(isrc);
    if (src_ref != _map.end()) {
        for (auto const& trans : src_ref->second) {
            for (auto const& dest : trans.second) {
                result.adjoin(dest.first.cell());
            }
        }
    }
    return result;
}

Map<IdentifiedCell, Map<IdentifiedCell, TargetScore>> const& DirectedHashedGraph::transitions(IdentifiedCell const& cell) const {
    return _map.at(cell);
}

void DirectedHashedGraph::insert_forward(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,TargetScore>> const& destination_cells) {
        auto itrans = _edge_factory.create(transition_cell);
        auto isrc = _vertex_factory.create(source_cell);
        auto src_ref = _map.find(isrc);
        if (src_ref == _map.end()) {
            _map.insert(make_pair(isrc, Map<IdentifiedCell,Map<IdentifiedCell,TargetScore>>()));
            src_ref = _map.find(isrc);
        }

        src_ref->second.insert(make_pair(itrans, Map<IdentifiedCell,TargetScore>()));
        auto trans_ref = src_ref->second.find(itrans);
        for (auto const& dst : destination_cells) {
            auto idst = _vertex_factory.create(dst.first);
            trans_ref->second.insert(make_pair(idst, dst.second));
        }
    }

void DirectedHashedGraph::insert_backward(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,TargetScore>> const& destination_cells) {
    auto itrans = _edge_factory.create(transition_cell);
    for (auto const& src : destination_cells) {
        auto isrc = _vertex_factory.create(src.first);
        auto src_ref = _map.find(isrc);
        if (src_ref == _map.end()) {
            _map.insert(make_pair(isrc, Map<IdentifiedCell,Map<IdentifiedCell,TargetScore>>()));
            src_ref = _map.find(isrc);
        }
        auto trans_ref = src_ref->second.find(itrans);
        if (trans_ref == src_ref->second.end()) {
            src_ref->second.insert(make_pair(itrans, Map<IdentifiedCell,TargetScore>()));
            trans_ref = src_ref->second.find(itrans);
        }
        auto idst = _vertex_factory.create(source_cell);
        trans_ref->second.insert(make_pair(idst, src.second));
    }
}

DirectedHashedGraph::Iterator DirectedHashedGraph::find(IdentifiedCell const& source_icell) {
    return _map.find(source_icell);
}

DirectedHashedGraph::Iterator DirectedHashedGraph::find(NCell const& source_cell) {
    return _map.find(_vertex_factory.create(source_cell));
}

bool DirectedHashedGraph::contains(Iterator const& iterator) const {
    return iterator != _map.end();
}

void DirectedHashedGraph::erase(NCell const& source_cell) {
    auto icell = _vertex_factory.create(source_cell);
    _map.erase(_vertex_factory.create(source_cell));
}

void DirectedHashedGraph::erase(NCell const& source, ECell const& transition, NCell const& destination) {
    auto const& src_ref = _map.find(_vertex_factory.create(source));
    if (src_ref != _map.end()) {
        auto const& trans_ref = src_ref->second.find(_edge_factory.create(transition));
        if (trans_ref != src_ref->second.end()) {
            trans_ref->second.erase(_vertex_factory.create(destination));
        }
    }
}

List<Tuple<IdentifiedCell,IdentifiedCell,IdentifiedCell>> DirectedHashedGraph::deadlock_transitions() const {
    List<Tuple<IdentifiedCell,IdentifiedCell,IdentifiedCell>> result;
    for (auto const& src : _map) {
        for (auto const& trans : src.second) {
            for (auto const& dst : trans.second) {
                if (not _map.contains(dst.first)) {
                    result.push_back({src.first,trans.first,dst.first});
                }
            }
        }
    }
    return result;
}

void DirectedHashedGraph::restrict_sources_to(SPaving const& paving) {
    for (auto it = _map.cbegin(); it != _map.cend(); ) {
        if (not paving.superset(it->first.cell())) _map.erase(it++);
        else ++it;
    }
}

void DirectedHashedGraph::restrict_destinations_to(SPaving const& paving) {
    for (auto src_it = _map.begin(); src_it != _map.end(); ) {
        for (auto trans_it = src_it->second.begin(); trans_it != src_it->second.end(); ) {
            for (auto dest_it = trans_it->second.begin(); dest_it != trans_it->second.end(); ) {
                if (not paving.superset(dest_it->first.cell())) trans_it->second.erase(dest_it++);
                else ++dest_it;
            }
            if (trans_it->second.empty()) src_it->second.erase(trans_it++);
            else ++trans_it;
        }
        if (src_it->second.empty()) _map.erase(src_it++);
        else ++src_it;
    }
}

void DirectedHashedGraph::apply_source_restriction_to(SPaving& paving) const {
    SPaving sources(paving.grid());
    for (auto const& entry : _map) {
        sources.adjoin(entry.first.cell());
    }

    paving.restrict(sources);
}

void DirectedHashedGraph::apply_source_removal_to(SPaving& paving) const {
    for (auto const& entry : _map) {
        paving.remove(entry.first.cell());
    }
}

void DirectedHashedGraph::clear() {
    _map.clear();
}

bool DirectedHashedGraph::is_empty() const {
    return _map.empty();
}

void DirectedHashedGraph::sweep() {
    auto src_it = _map.begin();
    while (src_it != _map.end()) {
        auto trans_it = src_it->second.begin();
        while (trans_it != src_it->second.end()) {
            if (trans_it->second.empty()) src_it->second.erase(trans_it++);
            else ++trans_it;
        }
        if (src_it->second.empty()) _map.erase(src_it++);
        else  ++src_it;
    }
}

OutputStream& operator<<(OutputStream& os, DirectedHashedGraph const& g) {
    for (auto const& src : g._map) {
        os << src.first.id() << src.first.cell().box() << ":[";
        for (auto const& trans : src.second) {
            os << trans.first.id() << "->";
            os << "(";
            for (auto const& dst : trans.second) {
                os << dst.first.id() << "@" << dst.second << ",";
            }
            os.seekp(-1, std::ios_base::end);
            os << ")";
        }
        os << "]\n";
    }
    os << "}";

    return os;
}

ForwardBackwardReachabilityGraph::ForwardBackwardReachabilityGraph(IdentifiedCellFactory const& vertex_factory, IdentifiedCellFactory const& edge_factory) :
    _forward_graph(vertex_factory,edge_factory), _backward_graph(vertex_factory,edge_factory) { }

ForwardBackwardReachabilityGraph::ForwardBackwardReachabilityGraph(ForwardBackwardReachabilityGraph const& other) :
    _forward_graph(other._forward_graph), _backward_graph(other._backward_graph) { }

SizeType ForwardBackwardReachabilityGraph::vertex_id(NCell const& cell) const {
    return _forward_graph.vertex_id(cell);
}

SizeType ForwardBackwardReachabilityGraph::edge_id(NCell const& cell) const {
    return _forward_graph.edge_id(cell);
}

SizeType ForwardBackwardReachabilityGraph::num_transitions() const {
    auto result = _forward_graph.num_transitions();
    ARIADNE_ASSERT_EQUAL(result,_backward_graph.num_transitions())
    return result;
}

SizeType ForwardBackwardReachabilityGraph::num_sources() const {
    return _forward_graph.num_sources();
}

SizeType ForwardBackwardReachabilityGraph::num_destinations() const {
    return _backward_graph.num_sources();
}

void ForwardBackwardReachabilityGraph::insert(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,TargetScore>> const& destination_cells) {
    _forward_graph.insert_forward(source_cell,transition_cell,destination_cells);
    _backward_graph.insert_backward(source_cell,transition_cell,destination_cells);
}

void ForwardBackwardReachabilityGraph::clear() {
    _forward_graph.clear();
    _backward_graph.clear();
}

List<Tuple<IdentifiedCell,IdentifiedCell,IdentifiedCell>> ForwardBackwardReachabilityGraph::deadlock_transitions() const {
    List<Tuple<IdentifiedCell,IdentifiedCell,IdentifiedCell>> result = _forward_graph.deadlock_transitions();
    return result;
}

void ForwardBackwardReachabilityGraph::reduce_to_avoiding(SPaving const& unsafe) {
    CONCLOG_SCOPE_CREATE

    std::deque<NCell> unsafe_cells_queue;
    for (auto const& c : unsafe) unsafe_cells_queue.push_back(c);

    auto indicator_final_original = static_cast<double>(unsafe.size());
    double indicator_current_value = 0;
    ProgressIndicator indicator(indicator_final_original);
    while (not unsafe_cells_queue.empty()) {
        CONCLOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ");
        auto const& u = unsafe_cells_queue.front();
        auto const& bw_unsafe_ref = _backward_graph.find(u);
        if (_backward_graph.contains(bw_unsafe_ref)) {
            for (auto const& bw_unsafe_trans : bw_unsafe_ref->second) {
                std::deque<NCell> sources;
                for (auto const& src_of_unsafe_trans : bw_unsafe_trans.second)
                    sources.push_back(src_of_unsafe_trans.first.cell());

                for (auto const& src_cell : sources) {
                    auto const& fw_src_of_unsafe_trans = _forward_graph.find(src_cell);
                    if (_forward_graph.contains(fw_src_of_unsafe_trans)) {
                        auto const& fw_trans = fw_src_of_unsafe_trans->second.at(bw_unsafe_trans.first);
                        for (auto const& dst_to_prune : fw_trans) {
                            _backward_graph.erase(dst_to_prune.first.cell(),bw_unsafe_trans.first.cell(),src_cell);
                        }
                        fw_src_of_unsafe_trans->second.erase(bw_unsafe_trans.first);
                        if (fw_src_of_unsafe_trans->second.empty()) {
                            unsafe_cells_queue.push_back(src_cell);
                            indicator.update_final(++indicator_final_original);
                            _forward_graph.erase(src_cell);
                        }
                    }
                }
            }
        }
        _backward_graph.erase(u);
        unsafe_cells_queue.pop_front();
        indicator.update_current(++indicator_current_value);
    }

    _backward_graph.sweep();
}

void ForwardBackwardReachabilityGraph::reduce_to_possibly_reaching(SPaving const& goals) {

    SPaving analysed(goals.grid());
    SPaving newly_reached = goals;
    while (not newly_reached.is_empty()) {
        auto destination = *newly_reached.begin();
        analysed.adjoin(destination);

        auto sources = _backward_graph.destinations_from(destination);
        sources.remove(analysed);
        newly_reached.adjoin(sources);
        newly_reached.remove(destination);
    }

    _forward_graph.restrict_destinations_to(analysed);
    _forward_graph.restrict_sources_to(analysed);
    _backward_graph.restrict_destinations_to(analysed);
    _backward_graph.restrict_sources_to(analysed);
}

void ForwardBackwardReachabilityGraph::apply_source_removal_to(SPaving& paving) const {
    _forward_graph.apply_source_removal_to(paving);
}

void ForwardBackwardReachabilityGraph::apply_source_restriction_to(SPaving& paving) const {
    _forward_graph.apply_source_restriction_to(paving);
}

List<Set<IdentifiedCell>> ForwardBackwardReachabilityGraph::sets_equidistant_to_goals(SPaving const& goal) const {

    List<Set<IdentifiedCell>> result;

    auto backward_copy = _backward_graph;

    Set<IdentifiedCell> goal_cells;
    for (auto const& g : goal)
        goal_cells.insert(_backward_graph.vertex_icell(g));
    result.append(goal_cells);

    Set<IdentifiedCell> currently_reached;
    currently_reached.adjoin(goal_cells);

    auto num_sources = _backward_graph.num_sources();

    while (currently_reached.size() < num_sources) {

        SizeType idx = result.size()-1;
        Set<IdentifiedCell> sources;

        for (auto const& g : result[idx]) {
            auto it = backward_copy.find(g);
            for (auto const& t : it->second) {
                for (auto const& s : t.second) {
                    if (not currently_reached.contains(s.first) and not sources.contains(s.first)) {
                        sources.insert(s.first);
                    }
                }
            }
        }

        for (auto const& s : result[idx]) {
            backward_copy.erase(s.cell());
        }

        result.append(sources);
        currently_reached.adjoin(sources);
    }

    return result;
}

Map<IdentifiedCell, Map<IdentifiedCell, TargetScore>> const& ForwardBackwardReachabilityGraph::forward_transitions(IdentifiedCell const& source) const {
    return _forward_graph.transitions(source);
}

Map<IdentifiedCell, Map<IdentifiedCell, TargetScore>> const& ForwardBackwardReachabilityGraph::backward_transitions(IdentifiedCell const& destination) const {
    return _backward_graph.transitions(destination);
}

ReachabilityGraphInterface* ForwardBackwardReachabilityGraph::clone() const {
    return new ForwardBackwardReachabilityGraph(*this);
}

void ForwardBackwardReachabilityGraph::write(std::ostream& os) const {
    os << "Forward:\n" << _forward_graph << "\nBackward:\n" << _backward_graph;
}

BoundedDomainRAG::BoundedDomainRAG(SharedPointer<ReachabilityGraphInterface> graph) :
    _internal(graph) { }

BoundedDomainRAG& BoundedDomainRAG::operator=(BoundedDomainRAG const& other) {
    _internal = other._internal;
    return *this;
}

void BoundedDomainRAG::apply_source_restriction_to(SPaving& paving) const {
    _internal->apply_source_restriction_to(paving);
}

List<Tuple<IdentifiedCell,IdentifiedCell,IdentifiedCell>> BoundedDomainRAG::deadlock_transitions() const {
    return _internal->deadlock_transitions();
}

BoundedDomainRAG::BoundedDomainRAG(BoundedDomainRAG const& other) {
    this->_internal = other._internal;
}

ReachabilityGraphInterface const& BoundedDomainRAG::internal() const {
    return *_internal;
}

bool BoundedDomainRAG::is_empty() const {
    return _internal == nullptr;
}

SizeType BoundedDomainRAG::num_sources() const {
    return _internal->num_sources();
}

SizeType BoundedDomainRAG::num_destinations() const {
    return _internal->num_destinations();
}

AvoidingRAG BoundedDomainRAG::reduce_to_not_reaching(SPaving const& unsafe) const {
    return AvoidingRAG(*this, unsafe);
}

AvoidingRAG::AvoidingRAG(BoundedDomainRAG const& free_graph, SPaving const& unsafe) : _internal(free_graph.internal().clone()), _unsafe(unsafe) {
    _internal->reduce_to_avoiding(unsafe);
}

AvoidingRAG& AvoidingRAG::operator=(AvoidingRAG const& other) {
    _internal = other._internal;
    return *this;
}

AvoidingRAG::AvoidingRAG(AvoidingRAG const& other) {
    this->_internal = other._internal;
}

ReachabilityGraphInterface const& AvoidingRAG::internal() const {
    return *_internal;
}

bool AvoidingRAG::is_empty() const {
    return _internal == nullptr;
}

SizeType AvoidingRAG::num_sources() const {
    return _internal->num_sources();
}

SizeType AvoidingRAG::num_destinations() const {
    return _internal->num_destinations();
}

List<Tuple<IdentifiedCell,IdentifiedCell,IdentifiedCell>> AvoidingRAG::deadlock_transitions() const {
    return _internal->deadlock_transitions();
}

PossiblyReachingRAG AvoidingRAG::reduce_to_possibly_reaching(SPaving const& goals) const {
    return PossiblyReachingRAG(*this, goals);
}

void AvoidingRAG::apply_source_restriction_to(SPaving& paving) const {
    _internal->apply_source_restriction_to(paving);
}

PossiblyReachingRAG::PossiblyReachingRAG(AvoidingRAG const& avoid_graph, SPaving const& goals) : _internal(avoid_graph.internal().clone()), _goals(goals) {
    _internal->reduce_to_possibly_reaching(goals);
}

PossiblyReachingRAG::PossiblyReachingRAG(PossiblyReachingRAG const& other) {
    this->_internal = other._internal;
    this->_goals = other._goals;
}

PossiblyReachingRAG& PossiblyReachingRAG::operator=(PossiblyReachingRAG const& other) {
    _internal = other._internal;
    _goals = other._goals;
    return *this;
}

List<Set<IdentifiedCell>> PossiblyReachingRAG::sets_equidistant_to_goals() const {
    return _internal->sets_equidistant_to_goals(_goals);
}

ReachabilityGraphInterface const& PossiblyReachingRAG::internal() const {
    return *_internal;
}

bool PossiblyReachingRAG::is_empty() const {
    return _internal == nullptr;
}

SizeType PossiblyReachingRAG::num_sources() const {
    return _internal->num_sources();
}

SizeType PossiblyReachingRAG::num_destinations() const {
    return _internal->num_destinations();
}

void PossiblyReachingRAG::apply_source_removal_to(SPaving& paving) const {
    _internal->apply_source_removal_to(paving);
}

List<Tuple<IdentifiedCell,IdentifiedCell,IdentifiedCell>> PossiblyReachingRAG::deadlock_transitions() const {
    return _internal->deadlock_transitions();
}

} // namespace Ariadne