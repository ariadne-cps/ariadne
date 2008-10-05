
#ifndef ARIADNE_LIST_SET_H
#define ARIADNE_LIST_SET_H


#include <vector>

namespace Ariadne {


template<class BS>
class ListSet 
{
 private:
  uint _dimension;
  std::vector<BS> _vector;

 public:
  typedef typename std::vector<BS>::const_iterator const_iterator;

  virtual ~ListSet() { }

  ListSet() { };
  explicit ListSet(uint);
  explicit ListSet(const BS& bs);
  template<class BST> ListSet(const ListSet<BST>& ls);
  template<class Iter> ListSet(Iter first, Iter last);



  /*! \brief Returns the number of basic sets forming this object. */
  uint size() const;

  /*! \brief Accesses the i-th BasicSet. */
  const BS& operator[](uint index) const;
  
  /*! \brief Make the set empty. */
  void clear();
      
  /*! \brief A constant iterator to the beginning of the list of basic sets. */
  const_iterator begin() const;
  
  /*! \brief A constant iterator to the end of the list of basic sets. */
  const_iterator end() const;

  /*! \brief Returns the denotable set's space dimension. */
  uint dimension() const;

  /*! \brief Removes a set from the list and return it. */
  BS pop();

  /*! \brief Adjoins (makes union with) a basic set. */
  void adjoin(const BS& bs) { this->_vector.push_back(bs); }

  /*! \brief Adjoins (makes union with) another list set. */
  void adjoin(const ListSet<BS>& ls);
};      
  

} // namespace Ariadne

#endif /* ARIADNE_LIST_SET_H */
