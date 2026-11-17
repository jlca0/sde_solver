#include <vector>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include "solver.hpp" 

namespace py = pybind11;

PYBIND11_MODULE(sde, m) {
    m.doc() = "Simulación estocástica en C++";

    // Camino browniano
    py::class_<BrownianPath>(m, "BrownianPath")
        .def(py::init<double, int>(), py::arg("T"), py::arg("pasos"))
        .def_readwrite("camino", &BrownianPath::camino);

    // Puente browniano
    py::class_<BrownianBridge>(m, "BrownianBridge")
        .def(py::init<double, const BrownianPath&>(), py::arg("T"), py::arg("W"))
        .def_readwrite("camino", &BrownianBridge::camino);

    // Camino geométrico exacto
    py::class_<ExactGeomPath>(m, "ExactGeomPath")
        .def(py::init<double, double, double, double, const BrownianPath&>(),
             py::arg("X0"), py::arg("T"), py::arg("mu"), py::arg("sigma"), py::arg("W"))
        .def("get_trayectoria", &ExactGeomPath::get_trayectoria);

    // Euler-Maruyama
    py::class_<EulerMaruyama>(m, "EulerMaruyama")
        .def(py::init<int>(), py::arg("pasos"))
        .def("get_trayectoria", &EulerMaruyama::get_trayectoria)
        .def("solve", &EulerMaruyama::solve, 
             py::arg("X0"), py::arg("T"), py::arg("mu"), py::arg("sigma"), py::arg("W"));

    // Milstein
    py::class_<Milstein>(m, "Milstein")
        .def(py::init<int>(), py::arg("pasos"))
        .def("get_trayectoria", &Milstein::get_trayectoria)
        .def("solve", &Milstein::solve, 
             py::arg("X0"), py::arg("T"), py::arg("mu"), py::arg("sigma"), py::arg("sigma_x"), py::arg("W"));
}