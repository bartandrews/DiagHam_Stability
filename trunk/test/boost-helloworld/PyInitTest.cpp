#include <boost/python.hpp> //boost libraries
#include <iostream>
#include <Python.h> //python libraries

using namespace boost::python;

int main()
{
    Py_Initialize();

    object main_module = import("__main__");
    object main_namespace = main_module.attr("__dict__");

    std::cout << "Hello!" << std::endl;

    return 0;
}
