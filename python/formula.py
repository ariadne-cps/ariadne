#!/usr/bin/python

##############################################################################
#            formula.py
#
#  Copyright 2009  Pieter Collins <Pieter.Collins@cwi.nl>
##############################################################################

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

import sys
import math
import __builtin__

from ariadne import tribool,indeterminate
Indeterminate=indeterminate()

class Real: float

class UnaryOperator:
    pass

class BinaryOperator:
    pass

class BinaryComparison(BinaryOperator):
    pass

class BinaryArithmetic(BinaryOperator):
    pass

def add(*args):
    res=0
    for arg in args:
        res=res+arg
    return res

def _neg_(arg): return -arg
def _add_(arg1,arg2): return arg1+arg2
def _mul_(arg1,arg2): return arg1*arg2
def _sub_(arg1,arg2): return arg1-arg2
def _div_(arg1,arg2): return arg1/arg2
def _pow_(arg1,arg2): return arg1**arg2

def _and_(arg1,arg2): return arg1 & arg2
def _or_(arg1,arg2): return arg1 | arg2
def _not_(arg):
    if isinstance(arg,bool): return not arg;
    return ~arg

def _eq_(arg1,arg2): return arg1 == arg2
def _ne_(arg1,arg2): return arg1 != arg2

def _sgn_(arg):
    if arg>0: return True
    elif arg<0: return False
    else: return Indeterminate

def _gt_(arg1,arg2): return arg1 > arg2
def _lt_(arg1,arg2): return arg1 < arg2
def _ge_(arg1,arg2): return arg1 >= arg2
def _le_(arg1,arg2): return arg1 <= arg2

def _exp_(arg): return math.exp(arg)
def _log_(arg): return math.log(arg)
def _sin_(arg): return math.sin(arg)
def _cos_(arg): return math.cos(arg)
def _tan_(arg): return math.tan(arg)

class Expression:
    def simplify(self,subexpr={}):
        subexpr[self]=self
        return self

    def __neg__(self):
        return UnaryExpression(_neg_,self)

    def __add__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_add_,self,other)

    def __radd__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_add_,other,self)

    def __sub__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_sub_,self,other)

    def __rsub__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_sub_,other,self)

    def __mul__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_mul_,self,other)

    def __rmul__(self,other):
        if(isinstance(other,float)): other=Constant(str(other),'Real',other)
        return BinaryExpression(_mul_,other,self)

    def __div__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_div_,self,other)

    def __rdiv__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_div_,other,self)

    def __pow__(self,other):
        assert isinstance(other,int)
        return BinaryExpression(_pow_,self,Constant(other))

    def __and__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_and_,self,other)

    def __or__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_or_,self,other)

    def __invert__(self):
        return UnaryExpression(_not_,self)

    def __gt__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_gt_,self,other)

    def __lt__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_lt_,self,other)

    def __ge__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_ge_,self,other)

    def __le__(self,other):
        if(not isinstance(other,Expression)): other=Constant(other)
        return BinaryExpression(_le_,self,other)

    def __str__(self):
        res=""
        for arg in self.args: res=res+str(arg)+","
        return str(self.op.__name__)+"("+res[:-1]+")"

    def __repr__(self):
        return str(self)


def sgn(expr):
    if isinstance(expr,Expression): return UnaryExpression(_sgn_,expr)


def exp(expr):
    if isinstance(expr,Expression): return UnaryExpression(_exp_,expr)

def log(expr):
    if isinstance(expr,Expression): return UnaryExpression(_log_,expr)

def sin(expr):
    if isinstance(expr,Expression): return UnaryExpression(_sin_,expr)

def cos(expr):
    if isinstance(expr,Expression): return UnaryExpression(_cos_,expr)

def tan(expr):
    if isinstance(expr,Expression): return UnaryExpression(_tan_,expr)


class UnaryExpression(Expression):
    def __init__(self,op,arg):
        self.op=op
        self.args=(arg,)
        for arg in self.args:
            assert(isinstance(arg,Expression))

    def variables(self):
        return self.args[0].variables()

    def evaluate(self,val):
        return self.op(self.args[0].evaluate(val))

    def substitute(self,val):
        return self.op(self.args[0].substitute(val))

    def simplify(self,subexpr={}):
        arg=self.args[0]
        if arg not in subexpr:
            sarg=arg.simplify()
            subexpr[arg]=sarg
        else:
            sarg=subexpr[arg]
        sself=UnaryExpression(self.op,sarg)
        if sarg.op=="cnst":
            sself=Constant(str(self),'Real',self.op(sarg.value))
        subexpr[self]=sself
        return sself

    def __str__(self):
        opstr={_neg_:'neg', _exp_:'exp', _log_:'log', _sin_:'sin', _cos_:'cos', _tan_:'&tan', _sgn_:'sgn'}
        if self.op==_not_:
            return '~'+repr(self.args[0])+''
        elif self.op==_neg_:
            return '-'+repr(self.args[0])+''
        else:
            return opstr[self.op]+'('+repr(self.args[0])+')'


class BinaryExpression(Expression):
    def __init__(self,op,arg1,arg2):
        self.op=op
        self.args=(arg1,arg2)
        for arg in self.args:
            assert(isinstance(arg,Expression))

    def variables(self):
        return (self.args[0].variables()).union(self.args[1].variables())

    def evaluate(self,val):
        return self.op(self.args[0].evaluate(val),self.args[1].evaluate(val))

    def substitute(self,val):
        return self.op(self.args[0].substitute(val),self.args[1].substitute(val))

    def simplify(self,subexpr={}):
        op=self.op
        args=self.args
        for arg in self.args:
            if arg not in subexpr:
                subexpr[arg]=arg.simplify(subexpr)
        is_const=True
        sargs=[]
        for arg in self.args:
            sargs.append(subexpr[arg])
        sself=BinaryExpression(op,sargs[0],sargs[1])
        if sargs[0].op=="cnst" and sargs[1].op=="cnst":
            val=op(sargs[0].value,sargs[1].value)
            sself=Constant(val)
            subexpr[self]=sself
            return sself
        if sargs[0].op=="cnst":
            if (self.op==_and_ and sargs[0].value==True) \
                    or (self.op==_or_ and sargs[0].value==False):
                sself=sargs[1]
            elif (self.op==_or_ and sargs[0].value==True) \
                    or (self.op==_and_ and sargs[0].value==False):
                sself=sargs[0]
        elif sargs[1].op=="cnst":
            if (self.op==_and_ and sargs[1].value==True) \
                    or (self.op==_or_ and sargs[1].value==False):
                sself=sargs[0]
            elif (self.op==_or_ and sargs[1].value==True) \
                    or (self.op==_and_ and sargs[1].value==False):
                sself=sargs[1]
        subexpr[self]=sself
        return sself

    def __str__(self):
        opstr={_add_:' + ', _sub_:' - ', _mul_:'*', _div_:'/', \
                   _and_:' & ', _or_:' | ', \
                   _eq_:'==', _ne_:'!=', _gt_:'>',_lt_:'<', _ge_:'>=', _le_:'<='}
        if self.op==_pow_: return 'pow('+repr(self.args[0])+','+repr(self.args[1])+')'
        return repr(self.args[0])+''+opstr[self.op]+''+repr(self.args[1])


class Variable(Expression):
    def __init__(self,name,type):
        assert(isinstance(name,str))
        self.op="var"
        self.args=(name,)
        self.name=name
        self.type=str(type)
    def variables(self):
        return set([self])
    def evaluate(self,val):
        assert(self in val)
        return val[self]
    def substitute(self,val):
        if self in val:
            c=val[self]
            return Constant(str(c),self.type,c)
        else:
            return self
    def __hash__(self):
        return hash(self.name)
    def __str__(self):
        #return self.name
        return self.name+":"+str(self.type)
    def __repr__(self):
        return self.name
    def __eq__(self,expression):
        if isinstance(expression,str): expression=Constant(expression,'String',expression)
        return BinaryExpression(_eq_,self,expression)

class BooleanVariable(Variable):
    def __init__(self,name):
        Variable.__init__(self,name,'Bool')

class TriboolVariable(Variable):
    def __init__(self,name):
        Variable.__init__(self,name,'Tribool')

class StringVariable(Variable):
    def __init__(self,name):
        Variable.__init__(self,name,'String')

class IntegerVariable(Variable):
    def __init__(self,name):
        Variable.__init__(self,name,'Integer')

class RealVariable(Variable):
    def __init__(self,name):
        Variable.__init__(self,name,'Real')


class Constant(Expression):
    __types={ type(True):"Bool", type(Indeterminate): "Tribool", type("str"):"String", type(1): 'Integer', type(1.0): 'Real' }

    def __init__(self,name,type=None,value=None):
        if value==None:
            self.name=str(name)
            self.value=name
        else:
            assert(isinstance(name,str))
            self.name=name
            self.value=value
        if type==None:
            self.type=self.__types[__builtin__.type(self.value)]
        else:
            self.type=str(type)
        self.op="cnst"
        self.args=(value,)
    def variables(self):
        return set()
    def evaluate(self,val):
        return self.value
    def substitute(self,val):
        return self
    def __str__(self):
        return str(self.value)+":Constant"+self.type
    def __repr__(self):
        return str(self.name)

class BooleanConstant(Constant):
    def __init__(self,name,value):
        Constant.__init__(self,name,'Boolean',bool(value))

class RealConstant(Constant):
    def __init__(self,name,value):
        Constant.__init__(self,name,'Real',value)


def evaluate(expr,val):
    return expr.evaluate(val)


if __name__=='__main__':
    x=RealVariable('x')
    v=RealVariable('v')
    t=RealVariable('t')
    g=Constant('g','Real',9.8)
    h=x+v*t-0.5*g*t*t
    print h
    print evaluate(h,{x:1.5,v:0.25,t:0.01})
    print evaluate(x*exp(g*t),{x:1.5,v:0.25,t:0.01})
    print

    print pow(x,3)
    print exp(x+3)
    print

    s=((2.0+t)*x*2.0*(t+5.0))
    print s
    ss=s.substitute({t:1.0})
    print ss
    sss=ss.simplify()
    print sss
    print

    a=Variable('a','Bool')
    b=Variable('b','Bool')
    c=Variable('c','Bool')

    e=(a&b)|~c
    print e
    es=e.substitute({b:True,c:True})
    print es
    es=es.simplify()
    print es
