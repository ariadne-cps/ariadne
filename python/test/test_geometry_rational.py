from ariadne import *

plt=Polytope([[-1,1],[0,-1],[2,3]])
plh=approx_polyhedron(plt)
plt=approx_polytope(plh)

qpt=QPoint([1,2])
qhs=QHalfspace([1,2],3)
qbx=QBox([[1,2],[3,4]])
qplt=QPolytope([[-1,1],[0,-1],[2,3]])
qplh=QPolyhedron([([-1,0],-1),([0,-1],-2),([1,1],3)])
qplt=QPolytope([[1,1],[0,0],[2,0],[1.75,0.5]])
qaf=QAffineFunction(QMatrix([[1.5,1],[1,1]]),QVector([-1,1]))

print qaf
print repr(qplt), qplt
print repr(qplh), qplh

QPolytope(qbx)
QPolytope(qplt)
QPolytope(qplh)
QPolyhedron(qhs)
QPolyhedron(qbx)
QPolyhedron(qplt)
QPolyhedron(qplh)

qplt.dimension()
qplt.vertices()
qplt.vertices()[0]
qplt.vertex(0)
qplt.empty()
qplt.bounded()
qplt.bounding_box()
qplt.contains(qpt)

qplh.dimension()
qplh.constraint(0)
qplh.empty()
qplh.bounded()
qplh.contains(qpt)

image(qaf,qplt)
preimage(qaf,qplh)

polyhedron(qbx)
polyhedron(qhs)
polyhedron(qplt)
polytope(qbx)
polytope(qplh)

subset(qplt,qplt)
subset(qplt,qplh)
subset(qplh,qplh)
#subset(qplh,qplt)
disjoint(qplt,qplt)
disjoint(qplt,qplh)
disjoint(qplh,qplt)
disjoint(qplh,qplh)
intersection(qplh,qplh)

eps=EpsPlot()
eps.open("geometry-rational.eps",Box([[-5,5],[-5,5]]))
eps.write(qplt)
eps.write(qplh)
eps.close()
