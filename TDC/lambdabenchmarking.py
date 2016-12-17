import dis
import ast
import os
from time import *
from math import *
from random import *
from timeit import timeit
from statistics import *

arr=[random() for r in range(100)]

def dothing1(num,arg1):
    return exp(1/(num+arg1))*sqrt(num%arg1)/max(num,arg1)*((1.0052346+arg1/100)**num)

def dothing2(num,arg1):
    return sum([dothing1(num*arr[j],arg1) for j in range(100)])

def main():
#   actual benchmarking
    seed(42)
    arg=random();
    tim=time()
    arr1=[]
    for i in range(100000):arr1.append(timeit("dothing2(i/100,arg)",setup="from __main__ import dothing1,dothing2"))

    print(str(time()-tim));arr1.sort();
    print("min: "+str(arr1[:100])+", max: "+str(arr1[-100:])+", av: "+str(mean(arr1))+", av: "+str(mean(arr1)))
    tim=time()
    return
    func=lambda num:sum([exp(1/(num*arr[j]+arg))*sqrt(num*arr[j] % arg)/max(num*arr[j],arg)*((1.0052346+arg/100)**num*arr[j]) for j in range(100)])

    print(str(time()-tim));tim=time()
    for i in range(100000):j=func(i/100)
    
    print(str(time()-tim));tim=time()
    
    return
#
    l1=lambda x:x**2
    l2=lambda y:y*l1(y)
    print(l2(5))
    print("dis'es")
    dis.dis(l1)
    print("---------dis2-----------")
    dis.dis(l2)
    print("--------showcode---------")
    dis.show_code(l2)
#    print("--------disco---------")
#    dis.disco(l2)
    print("--------getinstructions---------")
    i=dis.get_instructions(l2)
    print(i)
    thiscode=open(os.path.basename(__file__))
    thissource=thiscode.read()
    thiscode.close()
    
    # we pass <string> as a recognizable identifier since
    # we have no filename
    ast_root = compile(source, "<string>", "exec", ast.PyCF_ONLY_AST)
    walk_function_defs(ast_root)
    print("source2")
    ast_root = compile(source2, "<string>", "exec", ast.PyCF_ONLY_AST)
    walk_function_defs(ast_root)
    print("now this program")
    ast_root2=compile(thissource, "<string>", "exec", ast.PyCF_ONLY_AST)
    walk_function_defs(ast_root2)


source = """
def foo():
    bar()

def bar():
    baz()

def baz():
    print("hello!")
"""

source2="""
l1=lambda x:x**2
l2=lambda y:y*l1(y)
"""

def get_method_name_for_call(call_obj):
    for fieldname, value in ast.iter_fields(call_obj):
        if fieldname == "func":
            print(value)
            return value.id

def walk_method_calls(node, cur_func):
    if not node:
        return

    for cur_node in ast.iter_child_nodes(node):
        if type(cur_node) == ast.Call:
            method_called = get_method_name_for_call(cur_node)
            print("Found a call to "+method_called+" in body of "+cur_func+".")
        walk_method_calls(cur_node, cur_func)


def walk_function_defs(node):
    if not node:
        return

    for cur_node in ast.iter_child_nodes(node):
        if type(cur_node) == ast.FunctionDef:
            walk_method_calls(cur_node, cur_node.name)


if __name__=="__main__":main()
