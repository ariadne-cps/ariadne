from ariadne import *

plt=QPolytope([[-1,1],[0,-1],[2,3]])
plh=polyhedron(plt)
plt=QPolytope([[1,1],[0,0],[2,0],[1.75,0.5]])
af=QAffineFunction(QMatrix([[1.5,1],[1,1]]),QVector([-1,1]))

print af
print repr(plt), plt
print repr(plh), plh


plt.dimension()
plh.dimension()
plh.empty()
plt.vertices()
plt.vertices()[0]
plt.vertex(0)

image(af,plt)
preimage(af,plh)
polyhedron(plt)
polytope(plh)
subset(plt,plt)
subset(plt,plh)
subset(plh,plh)
#subset(plh,plt)
disjoint(plt,plt)
disjoint(plt,plh)
disjoint(plh,plt)
disjoint(plh,plh)
intersection(plh,plh)
#assert(equal(plt,plt))
assert(equal(plh,plh))

eps=EpsPlot()
eps.open("geometry-rational.eps",Box([[-5,5],[-5,5]]))
eps.write(plt)
eps.write(plh)
eps.close()
