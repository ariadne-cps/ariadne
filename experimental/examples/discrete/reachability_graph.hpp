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

typedef GridCell NCell;
typedef GridCell ECell;
typedef GridTreePaving SPaving;
typedef GridTreePaving CPaving;

template<class K, class V> using SizeTypeMap = Map<SizeType,Pair<K,V>>;
typedef Interval<double> BoundType;
typedef Vector<BoundType> BoundsBoxType;

String word_to_id(BinaryWord const& w, SizeType size) {
    std::stringstream ss;
    auto size_offset = w.size()-size;
    for (SizeType i=0; i<size; ++i) ss << to_string(w.at(size_offset+i));
    return ss.str();
}

SizeType to_identifier(GridCell const& cell, SizeType default_extent, std::map<String,SizeType> const& hashed) {
    SizeType size_to_use = cell.word().size() - (cell.root_extent() - default_extent)*cell.dimension();
    String expanded_id = word_to_id(cell.word(),size_to_use);
    return hashed.at(expanded_id);
}

class DirectedHashedGraph {
  public:
    typedef SizeTypeMap<NCell,SizeTypeMap<ECell,SPaving>>::iterator Iterator;
    typedef std::map<String,SizeType> HashTableType;
  public:
    DirectedHashedGraph(HashTableType const& hashed_vertices, HashTableType const& hashed_edges, SizeType default_vertex_extent, SizeType default_edge_extent) :
            _hashed_vertices(hashed_vertices), _hashed_edges(hashed_edges),
            _default_vertex_extent(default_vertex_extent), _default_edge_extent(default_edge_extent) { }

    SizeType num_sources() const { return _graph.size(); }

    SizeType num_transitions() const {
        SizeType result = 0;
        for (auto const& src : _graph) {
            for (auto const& ctrl : src.second.second) {
                result += ctrl.second.second.size();
            }
        }
        return result;
    }

    //! \brief Insert a forward entry from \a source_cell using \a transition_cell with associated \a destination_cells
    void insert_forward(NCell const& source_cell, ECell const& transition_cell, SPaving const& destination_cells) {
        auto trans_it = to_identifier(transition_cell, _default_edge_extent, _hashed_edges);
        auto src_id = to_identifier(source_cell, _default_vertex_extent, _hashed_vertices);
        auto src_ref = _graph.find(src_id);
        if (src_ref == _graph.end()) {
            _graph.insert(make_pair(src_id,make_pair(source_cell,Map<SizeType,Pair<ECell,SPaving>>())));
            src_ref = _graph.find(src_id);
        }
        src_ref->second.second.insert(make_pair(trans_it, make_pair(transition_cell, destination_cells)));
    }

    //! \brief Insert backward entries from each of \a destination_cells to \a source_cell using \a transition_cell, hence
    //! hashing on the destination and transition
    void insert_backward(NCell const& source_cell, ECell const& transition_cell, SPaving const& destination_cells) {
        auto ctrl_id = to_identifier(transition_cell, _default_edge_extent, _hashed_edges);
        for (auto const& src : destination_cells) {
            auto src_id = to_identifier(src, _default_vertex_extent, _hashed_vertices);
            auto src_ref = _graph.find(src_id);
            if (src_ref == _graph.end()) {
                _graph.insert(make_pair(src_id,make_pair(src,Map<SizeType,Pair<ECell,SPaving>>())));
                src_ref = _graph.find(src_id);
            }
            auto dst_ref = src_ref->second.second.find(ctrl_id);
            if (dst_ref == src_ref->second.second.end()) {
                src_ref->second.second.insert(make_pair(ctrl_id,make_pair(transition_cell, SPaving(source_cell.grid()))));
                dst_ref = src_ref->second.second.find(ctrl_id);
            }
            dst_ref->second.second.adjoin(source_cell);
        }
    }

    Iterator find(NCell const& source_cell) {
        return _graph.find(to_identifier(source_cell, _default_vertex_extent, _hashed_vertices));
    }

    bool contains(Iterator const& iterator) const {
        return iterator != _graph.end();
    }

    void erase(NCell const& source_cell) {
        _graph.erase(to_identifier(source_cell, _default_vertex_extent, _hashed_vertices));
    }

    //! \brief Erase for a given \a source and \a transition its \a destination
    void erase(NCell const& source, ECell const& transition, NCell const& destination) {
        auto const& src_ref = _graph.find(to_identifier(source, _default_vertex_extent, _hashed_vertices));
        if (src_ref != _graph.end()) {
            auto const& trans_ref = src_ref->second.second.find(to_identifier(transition, _default_edge_extent, _hashed_edges));
            if (trans_ref != src_ref->second.second.end()) {
                trans_ref->second.second.remove(destination);
            }
        }
    }

    //! \brief Subtract the source cells of the graph from \a paving
    void apply_to(SPaving& paving) const {
        for (auto const& entry : _graph) {
            paving.remove(entry.second.first);
        }
    }

    void clear() { _graph.clear(); }

    //! \brief Remove transitions for which the set of destinations is empty
    //! and remove sources for which the set of transitions is empty
    void sweep() {
        auto src_it = _graph.begin();
        while (src_it != _graph.end()) {
            auto trans_it = src_it->second.second.begin();
            while (trans_it != src_it->second.second.end()) {
                if (trans_it->second.second.is_empty()) src_it->second.second.erase(trans_it++);
                else ++trans_it;
            }
            if (src_it->second.second.empty()) _graph.erase(src_it++);
            else  ++src_it;
        }
    }

    friend OutputStream& operator<<(OutputStream& os, DirectedHashedGraph const& g) {
        for (auto const& src : g._graph) {
            os << src.first << src.second.first.box() << ":[";
            for (auto const& ctrl : src.second.second) {
                os << ctrl.first << "->";
                os << "(";
                for (auto const& dst : ctrl.second.second) {
                    os << to_identifier(dst, g._default_vertex_extent, g._hashed_vertices) << ",";
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
    std::map<String,SizeType> const _hashed_vertices;
    std::map<String,SizeType> const _hashed_edges;
    SizeType const _default_vertex_extent;
    SizeType const _default_edge_extent;
    SizeTypeMap<NCell,SizeTypeMap<ECell,SPaving>>  _graph;
};

class ReachabilityGraphInterface {
  public:
    virtual SizeType num_transitions() const = 0;
    virtual SizeType num_sources() const = 0;
    virtual SizeType num_destinations() const = 0;

    virtual void insert(NCell const& source_cell, ECell const& transition_cell, SPaving const& destination_cells) = 0;
    virtual void clear() = 0;
    virtual void refine_to_safety_graph(SPaving const& unsafe, SPaving& unverified) = 0;
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
    ForwardBackwardReachabilityGraph(HashTableType const& hashed_vertices, HashTableType const& hashed_edges, SizeType default_vertex_extent, SizeType default_edge_extent) :
        _forward_graph(hashed_vertices,hashed_edges,default_vertex_extent,default_edge_extent),
        _backward_graph(hashed_vertices,hashed_edges,default_vertex_extent,default_edge_extent) { }

    SizeType num_transitions() const override {
        return _forward_graph.num_transitions();
    }

    SizeType num_sources() const override {
        return _forward_graph.num_sources();
    }

    SizeType num_destinations() const override {
        return _backward_graph.num_sources();
    }

    void insert(NCell const& source_cell, ECell const& transition_cell, SPaving const& destination_cells) override {
        _forward_graph.insert_forward(source_cell,transition_cell,destination_cells);
        _backward_graph.insert_backward(source_cell,transition_cell,destination_cells);
    }

    void clear() override {
        _forward_graph.clear();
        _backward_graph.clear();
    }

    void refine_to_safety_graph(SPaving const& unsafe, SPaving& unverified) override {
        CONCLOG_SCOPE_CREATE

        std::deque<NCell> unsafe_cells_queue;
        for (auto const& c : unsafe) unsafe_cells_queue.push_back(c);

        double indicator_final_original = unsafe.size();
        double indicator_current_value = 0;
        ProgressIndicator indicator(indicator_final_original);
        while (not unsafe_cells_queue.empty()) {
            CONCLOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ");
            auto const& u = unsafe_cells_queue.front();
            auto const& bw_unsafe_ref = _backward_graph.find(u);
            if (_backward_graph.contains(bw_unsafe_ref)) {
                for (auto const& bw_unsafe_trans : bw_unsafe_ref->second.second) {
                    for (auto const& src_of_unsafe_trans : bw_unsafe_trans.second.second) {
                        auto const& fw_src_of_unsafe_trans = _forward_graph.find(src_of_unsafe_trans);
                        if (_forward_graph.contains(fw_src_of_unsafe_trans)) {
                            auto const& fw_trans = fw_src_of_unsafe_trans->second.second.at(bw_unsafe_trans.first);
                            for (auto const& dst_to_prune : fw_trans.second) {
                                _backward_graph.erase(dst_to_prune,bw_unsafe_trans.second.first,src_of_unsafe_trans);
                            }
                            fw_src_of_unsafe_trans->second.second.erase(bw_unsafe_trans.first);
                            if (fw_src_of_unsafe_trans->second.second.empty()) {
                                unsafe_cells_queue.push_back(src_of_unsafe_trans);
                                indicator.update_final(++indicator_final_original);
                                _forward_graph.erase(src_of_unsafe_trans);
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

        _forward_graph.apply_to(unverified);
    }

    ReachabilityGraphInterface* clone() const override {
        return new ForwardBackwardReachabilityGraph(*this);
    }

    void write(std::ostream& os) const override {
        os << "Forward:\n" << _forward_graph << "\nBackward:\n" << _backward_graph;
    }

  private:
    DirectedHashedGraph _forward_graph;
    DirectedHashedGraph _backward_graph;
};