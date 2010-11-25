from ariadne import *
#from ariadne import tribool,Indeterminate,Interval, Box, split,disjoint,subset

Indeterminate=indeterminate()
def unknown(tb):
    return possibly(tb) and not definitely(tb)

def definitely_not(tb):
    return not possibly(tb)

def not_possibly(tb):
    return not possibly(tb)

def negate(tb):
    if unknown(tb): return Indeterminate
    else: return not tb

def onside(set,box):
    if set.covers(box): return True
    elif set.separated(box): return False
    else: return Indeterminate

def tribool_repr(tb):
    if unknown(tb): return 'Indeterminate'
    if definitely(tb): return 'True'
    return 'False'

tribool.__repr__=tribool_repr

class Set(set):
    def insert(self,element):
        self.add(element)
    def contains(self,element):
        return self.__contains__(element)

def box_str(box):
    result=""
    for i in range(0,len(box)):
        if(i!=0): result+="x"
        result+="["+str(box[i].lower())+","+str(box[i].upper())+"]"
    return result

Box.__str__=box_str


class Node:
    def __init__(self,code,set,image,precessors=Set(),successors=Set(),initial=True,safe=Indeterminate,valid=Indeterminate):
        self.code=code
        self.set=set
        self.image=image
        self.precessors=precessors
        self.successors=successors
        self.initial=initial
        self.safe=safe
        self.valid=valid

    def __str__(self):
        return "Node( code="+str(self.code)+", set="+str(self.set)+", image="+str(self.image)+", initial="+str(self.initial)+", safe="+str(self.safe)+", valid="+str(self.valid)+")"

    def __repr__(self):
        return "\'"+str(self.code)+"\'"

class Validator:
    def __init__(self,system,initial_set,safe_set,bounding_box):
        self.system=system
        self.initial_set=initial_set
        self.safe_set=safe_set
        self.bounding_box=bounding_box
        self.partition=Set()
        image_set=system(bounding_box)
        is_initial=initial_set.overlaps(bounding_box)
        is_safe=onside(safe_set,bounding_box)
        base_node=Node("",bounding_box,image_set,Set(),Set(),is_initial,is_safe,is_safe)
        base_node.successors.insert(base_node)
        base_node.precessors.insert(base_node)
        self.partition.insert(base_node)

    def __str__(self):
        result="Validator\n    initial="+str(self.initial_set)+"\n    safe="+str(self.safe_set)+"\n"
        for node in self.partition:
            result+="  "
            if(node.initial): result+="I"
            else: result+=" "
            if(node.safe): result+="S"
            elif not possibly(node.safe): result+="U"
            else: result+=" "
            if(node.valid): result+="V"
            elif not possibly(node.valid): result+="X"
            else: result+="?"
            result+=node.code+str(node.set)+"->"+str(node.image)
            result+=repr(node.precessors)
            result+="->"
            result+=repr(node)
            result+="->"
            result+=repr(node.successors)
            result+="\n"
        return result

    def split_node(self,node):
        new_nodes=(node,Node(node.code,node.set,node.image,node.precessors.copy(),node.successors.copy(),node.initial,node.safe,node.valid))
        self.partition.insert(new_nodes[1])
        #Test for self-loop
        if new_nodes[0] in new_nodes[0].successors:
            new_nodes[0].successors.insert(new_nodes[1])
            new_nodes[0].precessors.insert(new_nodes[1])
            new_nodes[1].successors.insert(new_nodes[1])
            new_nodes[1].precessors.insert(new_nodes[1])
        for prec_node in new_nodes[1].precessors:
            prec_node.successors.add(new_nodes[1])
        for succ_node in new_nodes[1].successors:
            succ_node.precessors.add(new_nodes[1])
        new_boxes=split(node.set)
        new_nodes[0].set=new_boxes[0]
        new_nodes[1].set=new_boxes[1]
        new_nodes[0].code+='0'
        new_nodes[1].code+='1'

        for new_node in new_nodes:
            #Compute image set
            new_node.image=self.system(new_node.set)

            #Update precessor nodes
            for prec_node in Set(new_node.precessors):
                if disjoint(prec_node.image,new_node.set):
                    self.remove_image(prec_node,new_node)

            #Update sucessor nodes
            for succ_node in Set(new_node.successors):
                if disjoint(new_node.image,succ_node.set):
                    self.remove_image(new_node,succ_node)

            #Update whether node is an initial node
            if new_node.initial==False:
                pass
            else:
                if self.initial_set.separated(new_node.set):
                    new_node.initial=False
                elif self.initial_set.overlaps(new_node.set):
                    new_node.initial=True

            #Update whether node is safe
            if unknown(new_node.safe): #True: #new_node.safe==Indeterminate:
                if self.safe_set.covers(new_node.set):
                    new_node.safe=True
                elif self.safe_set.separated(new_node.set):
                    new_node.safe=False

        return new_nodes

    def remove_image(self,node,image_node):
        #print "remove_image("+node.code+","+image_node.code+")"
        #print node.precessors,"->",node.code,"->",node.successors
        #print image_node.precessors,"->",image_node.code,"->",image_node.successors
        assert(node.successors.contains(image_node))
        assert(image_node.precessors.contains(node))
        node.successors.remove(image_node)
        image_node.precessors.remove(node)
        #print node.precessors,"->",node.code,"->",node.successors
        #print image_node.precessors,"->",image_node.code,"->",image_node.successors
        #print
        #print self
        #print

    #Compute which nodes are safe or unsafe using depth first search
    def compute_safety(self):
        validity={}
        for node in self.partition:
            validity[node]=True
        visited=Set()
        stack=[]
        for node in self.partition:
            if node not in visited:
                self.safety_visit(node,validity,visited,stack)
        print validity
        for node in self.partition:
            node.valid=validity[node]


    #Compute which nodes are safe or unsafe using depth first search
    def safety_visit(self,node,validity,visited,stack):
        visited.insert(node)
        stack.append(node)
        validity[node]=onside(self.safe_set,node.set)
        if validity[node]==False:
            return
        all_successors_invalid=True
        all_successors_valid=True
        for succ in node.successors:
            if succ not in visited:
                self.safety_visit(succ,validity,visited,stack)
            if not definitely(validity[succ]):
                all_successors_valid=False
            if possibly(validity[succ]):
                all_successors_invalid=False
        if all_successors_invalid:
            validity[node]=False
        elif not all_successors_valid:
            validity[node]=Indeterminate
        stack.pop()
        return

    def shortest_paths(self):
        """Returns a map of Node:(Node,Integer) tuples of the form node:(prev,len) where prev is the
           previous node along the path, and len is the lenght. prev is None if len is 0."""




        # Choose node to refine adaptively
    # This method just picks from a number of alternatives
    def find_node_to_subdivide(self):
        self.find_node_to_subdivide1()

    # Choose node to refine adaptively
    # FIXME: This method is not finished
    def find_node_to_subdivide1(self):
        self.initial_nodes=Set()
        for node in self.partition:
            if node.is_initial:
                self.initial_nodes.insert(node)
        for node in self.initial_modes:
            if not possibly(node.is_valid):
                # System is unsafe
                pass
            if not definitely(node.is_valid):
                path=[node]
                for succ in node.successors:
                    if not possibly(node.is_valid):
                        return node
                    if not definitely(node.is_valid):
                        if succ not in path:
                            path.adjoin(succ)
                            break
                return node


def henon(pt):
    a=1.5
    b=0.3125
    x=Interval(pt[0])
    y=Interval(pt[1])
    return Box([a-x*x+b*y,x])

def scaling(pt):
    x=Interval(pt[0])
    y=Interval(pt[1])
    return Box([0.5*x+0.125,0.25*y+0.125])


if __name__=='__main__':
    print dir()
    print dir(tribool)
    z=Box([[0.01,0.1],[-0.1,0.02]])
    print z, henon(z)

    #system=henon
    system=scaling
    initial=Box([[0.125,0.25],[0.375,0.5]])
    safe=Box([[-1.5,+1.25],[-1.75,+1.25]])
    bound=Box([[-2,+2],[-2,+2]])

    validator=Validator(system,initial,safe,bound)
    print validator
    for node in validator.partition: print node,node.precessors,node.successors

    for i in range(0,8):
        nodes=Set(validator.partition)
        for node in nodes:
            subnode=validator.split_node(node)
        print validator

    validator.compute_safety()
    print validator
