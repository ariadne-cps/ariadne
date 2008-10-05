#include "numeric.h"              
#include "vector.h"               
#include <iostream>
#include <list>

namespace Ariadne {

typedef Vector<Interval> IntervalVector;
typedef std::list<IntervalVector> ImageList;

typedef IntervalVector (*IntervalFunction)(const IntervalVector&);
typedef tribool (*IntervalPredicate) (const IntervalVector&);


class Node;
class Paving;

typedef Node* Subpaving;  //makes Subpaving alias of
                          //pointer to a Node


class Node {
 private:
  IntervalVector theBox;
  tribool theValue;
  Subpaving leftChild;
  Subpaving rightChild;
  
 public:
  //constructors
  Node(const IntervalVector& v, tribool tb); //initialized
  Node(const Node& n); //copy
  
  //destructor
  ~Node();
  
  //other Methods
  friend const IntervalVector& box(const Subpaving a);
  
  friend bool is_empty(const Subpaving);
  friend bool is_full(const Subpaving);
  friend bool is_leaf(const Subpaving);
  friend bool is_null(const Subpaving);
  
  friend uint size(const Subpaving a);
  friend double volume(const Subpaving a);
  
  friend std::ostream& operator<< (std::ostream&, const Subpaving);
  friend std::ostream& write (std::ostream&, const Node*);
  
  friend tribool subset(const IntervalVector&, const Subpaving);
  
  friend Node* sivia(IntervalPredicate, const IntervalVector& box, double, tribool);
  friend void expand(Node*);
  friend void contract(Node*);
  friend Node* reunite(Node*, Node*, const IntervalVector&);
  
  friend Node* image(IntervalFunction, Node*, double);
  friend void mince(Node*, double);
  friend void flatten(Node*, ImageList&);
  friend void evaluate(IntervalFunction, Node*, ImageList&, IntervalVector&);
  friend Node* regularize(const ImageList&, const IntervalVector&, double);
};

double diam(const IntervalVector&, int&);
uint size(const Subpaving a);
double volume(const Subpaving a);
tribool subset(const IntervalVector&, const Subpaving);

Node* regularize(const ImageList&, const IntervalVector&, double);


class Paving {
 private:
  Subpaving _base;
 public:
  Paving() : _base(new Node(IntervalVector(),false)) { }
  Paving(const Subpaving a) : _base(new Node(*a)) { }
  Paving(const IntervalVector& bx) : _base(new Node(bx,true)) { }
  Paving(const Paving& pv) : _base(new Node(*pv._base)) { }
  Paving& operator=(const Paving& pv) { 
    if(this!=&pv) { delete this->_base; this->_base=new Node(*pv._base); } 
    return *this; }
  ~Paving() { delete this->_base; }
  Subpaving root() { return this->_base; }
  const Subpaving root() const { return this->_base; }
  uint size() const { return Ariadne::size(this->_base); }
  double volume() const { return Ariadne::volume(this->_base); }
  friend std::ostream& operator<< (std::ostream&, const Paving&);
};

inline tribool subset(const IntervalVector& z, const Paving& pv) {
  return subset(z,pv.root()); }

inline std::ostream& operator<< (std::ostream& os, const Paving& pv) {
  return os << pv._base; }

Paving mince(const Paving& pv, double eps);
Paving contract(const Paving& pv);

Paving inner_set(IntervalPredicate, const IntervalVector&, double);
Paving outer_set(IntervalPredicate, const IntervalVector&, double);

Paving image(IntervalFunction, const Paving&, double);

Paving outer_set(const ImageList&, const IntervalVector&, double);

ImageList flatten(const Paving& pv);

} // namespace Ariadne
