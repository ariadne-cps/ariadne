#ifndef ARIADNE_APPROXIMATOR_INTERFACE_H
#define ARIADNE_APPROXIMATOR_INTERFACE_H


template<class ES> class ListSet;

namespace Ariadne {
  
class GridCellListSet;
class BoxListSet;
  
    /*! \brief Interface for approximating enclosure sets on a paving of space.
     *  \ingroup EvaluatorInterfaces \ingroup Approximators
     */
    template<class BS, class ES> class ApproximatorInterface;

    template<class ES> class ApproximatorInterface<Box,ES>
    { 
      typedef Box BS;
      typedef BS BasicSet;
      typedef ES EnclosureSet;
      typedef ListSet<ES> EnclosureSetList;     
      typedef BoxListSet CoverListSet;
      typedef GridCellListSet PartitionListSet;
 
     public:
      /*! \brief Virtual destructor. */
      virtual ~ApproximatorInterface() { }

      /*! \brief Create a dynamically-allocated copy. */
      virtual ApproximatorInterface<BS,ES>* clone() const = 0;


      /*! \brief Computes whether an enclosure set is disjoint from a basic set. */
      virtual tribool disjoint(const ES& es, const BS& bs) const = 0;


      /*! \brief Computes an over-approximation of a set from a rectangle. */
      virtual EnclosureSet enclosure_set(const BasicSet& bs) const = 0;

      /*! \brief Computes a bounding box for a set. */
      virtual BasicSet bounding_box(const EnclosureSet& es) const = 0;


      /*! \brief Computes a bounding box for a list set. */
      virtual BasicSet bounding_box(const EnclosureSetList& es) const = 0;


      /*! \brief Computes a over-approximation of a set on a grid. */
      virtual void adjoin_over_approximation(CoverListSet&, const EnclosureSet& es) const = 0;

      /*! \brief Computes a over-approximation of a set on a grid. */
      virtual void adjoin_inner_approximation(PartitionListSet&, const EnclosureSet& es) const = 0;

      /*! \brief Computes a over-approximation of a set on a grid. */
      virtual void adjoin_outer_approximation(PartitionListSet&, const EnclosureSet& es) const = 0;

      /*! \brief Computes a over-approximation of a set on a grid. */
      virtual void adjoin_outer_approximation(PartitionListSet&, const EnclosureSetList& es) const = 0;


      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;
    };

    template<class Aprx, class ES>
    std::ostream& operator<<(std::ostream& os, const ApproximatorInterface<Aprx,ES>& ai) {
      return ai.write(os); 
    }




  
} // namespace Ariadne


#endif /* ARIADNE_APPROXIMATOR_INTERFACE_H */
