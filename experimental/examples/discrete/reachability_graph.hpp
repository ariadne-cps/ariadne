/***************************************************************************
 *            reachability_graph.hpp
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

#include "ariadne.hpp"

typedef double ProbabilityType;
typedef GridCell NCell;
typedef GridCell ECell;
typedef GridTreePaving SPaving;
typedef GridTreePaving CPaving;

struct IdentifiedCell;

typedef Map<IdentifiedCell,Map<IdentifiedCell,Map<IdentifiedCell,ProbabilityType>>> TransitionProbabilityMap;
typedef Interval<double> BoundType;
typedef Vector<BoundType> BoundsBoxType;

String word_to_id(BinaryWord const& w, SizeType size) {
    std::stringstream ss;
    SizeType effective_size = (w.size()<size? w.size() : size);
    auto size_offset = (w.size()<size? 0 : w.size()-size);
    for (SizeType i=0; i<effective_size; ++i) ss << to_string(w.at(size_offset+i));
    return ss.str();
}

class IdentifiedCellFactory;

struct IdentifiedCell {
    SizeType id() const { return _id; }
    GridCell cell() const { return _cell; }
    friend class IdentifiedCellFactory;
    bool operator<(IdentifiedCell const& other) const { return this->_id < other._id; }
  protected:
    IdentifiedCell(SizeType id, GridCell const& cell) : _id(id), _cell(cell) { }
  private:
    SizeType const _id;
    GridCell const _cell;
};

class IdentifiedCellFactory {
  public:
    typedef std::map<String,SizeType> HashTableType;
  public:
    IdentifiedCellFactory(SizeType default_extent, HashTableType const& table) : _default_extent(default_extent), _table(table) { }
    IdentifiedCell create(GridCell const& cell) const { return IdentifiedCell(_to_identifier(cell),cell); }
  private:
    SizeType _to_identifier(GridCell const& cell) const {
        SizeType size_to_use = (cell.root_extent() >= _default_extent ?
                                cell.word().size() - (cell.root_extent() - _default_extent)*cell.dimension() :
                                cell.word().size());
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
  private:
    SizeType const _default_extent;
    HashTableType const _table;
};

class DirectedHashedGraph {
  public:
    typedef TransitionProbabilityMap::iterator Iterator;
    typedef std::map<String,SizeType> HashTableType;
  public:
    DirectedHashedGraph(IdentifiedCellFactory const& vertex_factory, IdentifiedCellFactory const& edge_factory) :
            _vertex_factory(vertex_factory), _edge_factory(edge_factory) { }

    SizeType num_sources() const { return _map.size(); }

    SizeType num_transitions() const {
        SizeType result = 0;
        for (auto const& src : _map) {
            for (auto const& trans : src.second) {
                result += trans.second.size();
            }
        }
        return result;
    }

    SizeType vertex_id(NCell const& cell) const {
        return _vertex_factory.create(cell).id();
    }

    SizeType edge_id(NCell const& cell) const {
        return _edge_factory.create(cell).id();
    }

    SPaving destinations_from(NCell const& source_cell) const {
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

    //! \brief Insert a forward entry from \a source_cell using \a transition_cell with associated \a destination_cells
    void insert_forward(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,ProbabilityType>> const& destination_cells) {
        auto itrans = _edge_factory.create(transition_cell);
        auto isrc = _vertex_factory.create(source_cell);
        auto src_ref = _map.find(isrc);
        if (src_ref == _map.end()) {
            _map.insert(make_pair(isrc, Map<IdentifiedCell,Map<IdentifiedCell,ProbabilityType>>()));
            src_ref = _map.find(isrc);
        }

        src_ref->second.insert(make_pair(itrans, Map<IdentifiedCell,ProbabilityType>()));
        auto trans_ref = src_ref->second.find(itrans);
        for (auto const& dst : destination_cells) {
            auto idst = _vertex_factory.create(dst.first);
            trans_ref->second.insert(make_pair(idst, dst.second));
        }
    }

    //! \brief Insert backward entries from each of \a destination_cells to \a source_cell using \a transition_cell, hence
    //! hashing on the destination and transition
    void insert_backward(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,ProbabilityType>> const& destination_cells) {
        auto itrans = _edge_factory.create(transition_cell);
        for (auto const& src : destination_cells) {
            auto isrc = _vertex_factory.create(src.first);
            auto src_ref = _map.find(isrc);
            if (src_ref == _map.end()) {
                _map.insert(make_pair(isrc, Map<IdentifiedCell,Map<IdentifiedCell,ProbabilityType>>()));
                src_ref = _map.find(isrc);
            }
            auto trans_ref = src_ref->second.find(itrans);
            if (trans_ref == src_ref->second.end()) {
                src_ref->second.insert(make_pair(itrans, Map<IdentifiedCell,ProbabilityType>()));
                trans_ref = src_ref->second.find(itrans);
            }
            auto idst = _vertex_factory.create(source_cell);
            trans_ref->second.insert(make_pair(idst, src.second));
        }
    }

    Iterator find(NCell const& source_cell) {
        return _map.find(_vertex_factory.create(source_cell));
    }

    bool contains(Iterator const& iterator) const {
        return iterator != _map.end();
    }

    void erase(NCell const& source_cell) {
        _map.erase(_vertex_factory.create(source_cell));
    }

    //! \brief Erase for a given \a source and \a transition its \a destination
    void erase(NCell const& source, ECell const& transition, NCell const& destination) {
        auto const& src_ref = _map.find(_vertex_factory.create(source));
        if (src_ref != _map.end()) {
            auto const& trans_ref = src_ref->second.find(_edge_factory.create(transition));
            if (trans_ref != src_ref->second.end()) {
                trans_ref->second.erase(_vertex_factory.create(destination));
            }
        }
    }

    //! \brief Remove the source cells of the graph that are not in \a paving
    void restrict_sources_to(SPaving const& paving) {
        auto it = _map.begin();
        while (it != _map.end()) {
            if (not paving.superset(it->first.cell())) _map.erase(it++);
            else ++it;
        }
    }

    //! \brief Remove the destination cells of the graph that are not in \a paving
    //! \details Removes transitions/sources when empty
    void restrict_destinations_to(SPaving const& paving) {
        auto src_it = _map.begin();
        while (src_it != _map.end()) {
            auto trans_it = src_it->second.begin();
            while (trans_it != src_it->second.end()) {
                for (auto const& dest : trans_it->second) {
                    if (not paving.superset(dest.first.cell()))
                        trans_it->second.erase(dest.first);
                    if (trans_it->second.empty()) {
                        src_it->second.erase(trans_it++);
                        break;
                    }
                }
                if (not trans_it->second.empty())
                    ++trans_it;
            }
            if (src_it->second.empty()) _map.erase(src_it++);
            else ++src_it;
        }
    }

    //! \brief Remove the source cells of the graph from \a paving
    void apply_source_removal_to(SPaving& paving) const {
        for (auto const& entry : _map) {
            paving.remove(entry.first.cell());
        }
    }

    void clear() { _map.clear(); }

    //! \brief Remove transitions for which the set of destinations is empty
    //! and remove sources for which the set of transitions is empty
    void sweep() {
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

    friend OutputStream& operator<<(OutputStream& os, DirectedHashedGraph const& g) {
        for (auto const& src : g._map) {
            os << src.first.id() << src.first.cell() << ":[";
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
private:
    IdentifiedCellFactory const _vertex_factory;
    IdentifiedCellFactory const _edge_factory;
    TransitionProbabilityMap _map;
};

class ReachabilityGraphInterface {
  public:
    virtual SizeType num_transitions() const = 0;
    virtual SizeType num_sources() const = 0;
    virtual SizeType num_destinations() const = 0;

    virtual SizeType vertex_id(NCell const& cell) = 0;
    virtual SizeType edge_id(NCell const& cell) = 0;

    virtual void insert(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,ProbabilityType>> const& destination_cells) = 0;
    virtual void clear() = 0;

    //! \brief Remove those sources that can reach the \a avoidance paving
    virtual void reduce_to_not_reaching(SPaving const& avoidance) = 0;

    //! \brief Remove those sources that can not reach the \a goal paving
    virtual void reduce_to_possibly_reaching(SPaving const& goal) = 0;

    //! \brief Remove from \a unverified all the sources of the graph
    virtual void apply_source_removal_to(SPaving& unverified) const = 0;

    virtual ReachabilityGraphInterface* clone() const = 0;
    virtual void write(std::ostream& os) const = 0;
    virtual ~ReachabilityGraphInterface() = default;

    friend OutputStream& operator<<(OutputStream& os, ReachabilityGraphInterface const& g) {
        g.write(os);
        return os;
    }
};

class ForwardBackwardReachabilityGraph : public ReachabilityGraphInterface {
  public:
    typedef DirectedHashedGraph::HashTableType HashTableType;
  public:
    ForwardBackwardReachabilityGraph(IdentifiedCellFactory const& vertex_factory, IdentifiedCellFactory const& edge_factory) :
        _vertex_factory(vertex_factory), _edge_factory(edge_factory),
        _forward_graph(vertex_factory,edge_factory), _backward_graph(vertex_factory,edge_factory) { }

    SizeType vertex_id(NCell const& cell) override {
        return _forward_graph.vertex_id(cell);
    }

    SizeType edge_id(NCell const& cell) override {
        return _forward_graph.edge_id(cell);
    }

    SizeType num_transitions() const override {
        auto result = _forward_graph.num_transitions();
        ARIADNE_ASSERT_EQUAL(result,_backward_graph.num_transitions())
        return result;
    }

    SizeType num_sources() const override {
        return _forward_graph.num_sources();
    }

    SizeType num_destinations() const override {
        return _backward_graph.num_sources();
    }

    void insert(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,ProbabilityType>> const& destination_cells) override {
        _forward_graph.insert_forward(source_cell,transition_cell,destination_cells);
        _backward_graph.insert_backward(source_cell,transition_cell,destination_cells);
    }

    void clear() override {
        _forward_graph.clear();
        _backward_graph.clear();
    }

    void reduce_to_not_reaching(SPaving const& unsafe) override {
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

    void reduce_to_possibly_reaching(SPaving const& goals) override {
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

        _forward_graph.restrict_sources_to(analysed);
        _backward_graph.restrict_destinations_to(analysed);
    }

    void apply_source_removal_to(SPaving& paving) const override {
        _forward_graph.apply_source_removal_to(paving);
    }

    ReachabilityGraphInterface* clone() const override {
        return new ForwardBackwardReachabilityGraph(*this);
    }

    void write(std::ostream& os) const override {
        os << "Forward:\n" << _forward_graph << "\nBackward:\n" << _backward_graph;
    }

  private:
    IdentifiedCellFactory const _vertex_factory;
    IdentifiedCellFactory const _edge_factory;

    DirectedHashedGraph _forward_graph;
    DirectedHashedGraph _backward_graph;
};