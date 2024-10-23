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

#ifndef ARIADNE_REACHABILITY_GRAPH_HPP
#define ARIADNE_REACHABILITY_GRAPH_HPP

#include "utility/string.hpp"
#include "utility/container.hpp"
#include "algebra/vector.hpp"
#include "geometry/paving_interface.hpp"
#include "geometry/grid_paving.hpp"
#include "geometry/point.hpp"
#include "numeric/floatdp.hpp"

namespace Ariadne {

typedef double AlignmentType;
typedef GridCell NCell;
typedef GridCell ECell;
typedef GridTreePaving SPaving;
typedef GridTreePaving CPaving;
typedef Vector<double> PointType;

struct IdentifiedCell;

String word_to_id(BinaryWord const& w, SizeType length);

typedef Map<IdentifiedCell,Map<IdentifiedCell,Map<IdentifiedCell,AlignmentType>>> TransitionAlignmentMap;

class IdentifiedCellFactory;

struct IdentifiedCell {
    SizeType const& id() const { return _id; }
    GridCell const& cell() const { return _cell; }
    friend class IdentifiedCellFactory;
    bool operator<(IdentifiedCell const& other) const { return this->_id < other._id; }
    friend OutputStream& operator<<(OutputStream& os, IdentifiedCell const& ic) { os << "{" << ic._id << ":" << ic._cell.box() << "}"; return os; }
  protected:
    IdentifiedCell(SizeType id, GridCell const& cell) : _id(id), _cell(cell) { }
  private:
    SizeType _id;
    GridCell _cell;
};

class IdentifiedCellFactory {
  public:
    typedef std::map<String,SizeType> HashTableType;
  public:
    IdentifiedCellFactory(SizeType default_extent, HashTableType const& table);
    IdentifiedCell create(GridCell const& cell) const;
  private:
    SizeType _to_identifier(GridCell const& cell) const;
  private:
    SizeType const _default_extent;
    HashTableType const _table;
};

//! \brief A class for a directed graph with hashed sources, transitions and destinations
class DirectedHashedGraph {
  public:
    typedef TransitionAlignmentMap::iterator Iterator;
    typedef std::map<String,SizeType> HashTableType;
  public:
    DirectedHashedGraph(IdentifiedCellFactory const& vertex_factory, IdentifiedCellFactory const& edge_factory);
    DirectedHashedGraph(DirectedHashedGraph const& other);

    SizeType num_sources() const;

    SizeType num_transitions() const;

    IdentifiedCell vertex_icell(NCell const& cell) const;

    IdentifiedCell edge_icell(NCell const& cell) const;

    SizeType vertex_id(NCell const& cell) const;

    SizeType edge_id(NCell const& cell) const;

    SPaving destinations_from(NCell const& source_cell) const;

    //! \brief The transitions from a given \a cell
    Map<IdentifiedCell, Map<IdentifiedCell, AlignmentType>> const& transitions(IdentifiedCell const& cell) const;

    //! \brief Insert a forward entry from \a source_cell using \a transition_cell with associated \a destination_cells
    void insert_forward(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,AlignmentType>> const& destination_cells);

    //! \brief Insert backward entries from each of \a destination_cells to \a source_cell using \a transition_cell, hence
    //! hashing on the destination and transition
    void insert_backward(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,AlignmentType>> const& destination_cells);

    Iterator find(IdentifiedCell const& source_icell);

    Iterator find(NCell const& source_cell);

    bool contains(Iterator const& iterator) const;

    void erase(NCell const& source_cell);

    //! \brief Erase for a given \a source and \a transition its \a destination
    void erase(NCell const& source, ECell const& transition, NCell const& destination);

    //! \brief Remove the source cells of the graph that are not in \a paving
    void restrict_sources_to(SPaving const& paving);

    //! \brief Remove the destination cells of the graph that are not in \a paving
    //! \details Removes transitions/sources when empty
    void restrict_destinations_to(SPaving const& paving);

    //! \brief Remove the source cells of the graph from \a paving
    void apply_source_removal_to(SPaving& paving) const;

    //! \brief Remove all sources
    void clear();

    //! \brief Check if empty
    bool is_empty() const;

    //! \brief Remove transitions for which the set of destinations is empty
    //! and remove sources for which the set of transitions is empty
    void sweep();

    friend OutputStream& operator<<(OutputStream& os, DirectedHashedGraph const& g);

  private:
    IdentifiedCellFactory const _vertex_factory;
    IdentifiedCellFactory const _edge_factory;
    TransitionAlignmentMap _map;
};

//! \brief Interface for graphs used for discrete reachability under control laws
class ReachabilityGraphInterface {
  public:
    virtual SizeType num_transitions() const = 0;
    virtual SizeType num_sources() const = 0;
    virtual SizeType num_destinations() const = 0;

    virtual SizeType vertex_id(NCell const& cell) const = 0;
    virtual SizeType edge_id(NCell const& cell) const = 0;

    virtual void insert(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,AlignmentType>> const& destination_cells) = 0;
    virtual void clear() = 0;

    //! \brief Find the list of sets of cells having a given distance to the \a goal
    //! \details The position in the list determines the distance (0: within goals)
    virtual List<Set<IdentifiedCell>> sets_equidistant_to_goal(SPaving const& goal) const = 0;

    //! \brief The transitions from a given \a source forward
    virtual Map<IdentifiedCell, Map<IdentifiedCell, AlignmentType>> const& forward_transitions(IdentifiedCell const& source) const = 0;

    //! \brief The transitions from a given \a destination backward
    virtual Map<IdentifiedCell, Map<IdentifiedCell, AlignmentType>> const& backward_transitions(IdentifiedCell const& destination) const = 0;

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

//! \brief Graph having a dual forward and backward structure
class ForwardBackwardReachabilityGraph : public ReachabilityGraphInterface {
  public:
    typedef DirectedHashedGraph::HashTableType HashTableType;
  public:
    ForwardBackwardReachabilityGraph(IdentifiedCellFactory const& vertex_factory, IdentifiedCellFactory const& edge_factory);
    ForwardBackwardReachabilityGraph(ForwardBackwardReachabilityGraph const& other);

    List<Set<IdentifiedCell>> sets_equidistant_to_goal(SPaving const& goal) const override;

    SizeType vertex_id(NCell const& cell) const override;
    SizeType edge_id(NCell const& cell) const override;

    SizeType num_transitions() const override;
    SizeType num_sources() const override;
    SizeType num_destinations() const override;

    Map<IdentifiedCell, Map<IdentifiedCell, AlignmentType>> const& forward_transitions(IdentifiedCell const& source) const override;
    Map<IdentifiedCell, Map<IdentifiedCell, AlignmentType>> const& backward_transitions(IdentifiedCell const& destination) const override;

    void insert(NCell const& source_cell, ECell const& transition_cell, List<Pair<NCell,AlignmentType>> const& destination_cells) override;
    void clear() override;

    void reduce_to_not_reaching(SPaving const& unsafe) override;
    void reduce_to_possibly_reaching(SPaving const& goals) override;
    void apply_source_removal_to(SPaving& paving) const override;

    ReachabilityGraphInterface* clone() const override;

    void write(std::ostream& os) const override;

  private:
    DirectedHashedGraph _forward_graph;
    DirectedHashedGraph _backward_graph;
};

class AvoidingRAG;

//! \brief Graph with no reach or avoid restrictions
class UnconstrainedRAG {
  public:
    UnconstrainedRAG() = default;
    UnconstrainedRAG(UnconstrainedRAG const& other);
    bool is_empty() const;
    UnconstrainedRAG(SharedPointer<ReachabilityGraphInterface> graph);
    AvoidingRAG reduce_to_not_reaching(SPaving const& unsafe) const;
    ReachabilityGraphInterface const& internal() const;
    SizeType num_sources() const;
    SizeType num_destinations() const;
    UnconstrainedRAG& operator=(UnconstrainedRAG const& other);
  private:
    SharedPointer<ReachabilityGraphInterface> _internal;
};

class PossiblyReachingRAG;

//! \brief Graph with avoid restrictions
class AvoidingRAG {
  public:
    friend class UnconstrainedRAG;
    AvoidingRAG() = default;
    AvoidingRAG(AvoidingRAG const& other);
    ReachabilityGraphInterface const& internal() const;
    PossiblyReachingRAG reduce_to_possibly_reaching(SPaving const& goals) const;
    SizeType num_sources() const;
    SizeType num_destinations() const;
    bool is_empty() const;
    AvoidingRAG& operator=(AvoidingRAG const& other);
  protected:
    AvoidingRAG(UnconstrainedRAG const& free_graph, SPaving const& unsafe);
  private:
    SharedPointer<ReachabilityGraphInterface> _internal;
    SPaving const _unsafe;
};

//! \brief Graph with reach and avoid restrictions
class PossiblyReachingRAG {
  public:
    friend class AvoidingRAG;
    PossiblyReachingRAG() = default;
    PossiblyReachingRAG(PossiblyReachingRAG const& other);
    ReachabilityGraphInterface const& internal() const;

    List<Set<IdentifiedCell>> sets_equidistant_to_goal() const;

    SizeType num_sources() const;
    SizeType num_destinations() const;
    bool is_empty() const;
    void apply_source_removal_to(SPaving& paving) const;
    PossiblyReachingRAG& operator=(PossiblyReachingRAG const& other);
  protected:
    PossiblyReachingRAG(AvoidingRAG const& avoid_graph, SPaving const& goals);
  private:
    SharedPointer<ReachabilityGraphInterface> _internal;
    SPaving _goals;
};

} // namespace Ariadne

#endif