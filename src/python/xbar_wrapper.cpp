#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include "core/crossbar_simulator.h"
#include "core/simulation_settings.h"

namespace py = pybind11;

// Helper function to convert numpy array to vector of vectors
template <typename T>
std::vector<std::vector<T>> numpy_to_vector2d(py::array_t<T> array) {
    auto buf = array.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Input must be a 2D array");
    }
    
    size_t rows = buf.shape[0];
    size_t cols = buf.shape[1];
    T* ptr = static_cast<T*>(buf.ptr);
    
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols));
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            result[i][j] = ptr[i * cols + j];
        }
    }
    
    return result;
}

// Helper function to convert vector of vectors to numpy array
template <typename T>
py::array_t<T> vector2d_to_numpy(const std::vector<std::vector<T>>& vec) {
    if (vec.empty()) {
        std::vector<size_t> shape = {0, 0};
        return py::array_t<T>(shape);
    }
    
    size_t rows = vec.size();
    size_t cols = vec[0].size();
    
    std::vector<size_t> shape = {rows, cols};
    py::array_t<T> result(shape);
    auto buf = result.request();
    T* ptr = static_cast<T*>(buf.ptr);
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            ptr[i * cols + j] = vec[i][j];
        }
    }
    
    return result;
}

// Helper function to convert numpy array to Eigen vector
Eigen::VectorXf numpy_to_eigen_vector(py::array_t<float> array) {
    auto buf = array.request();
    if (buf.ndim != 1) {
        throw std::runtime_error("Input must be a 1D array");
    }
    
    size_t size = buf.shape[0];
    float* ptr = static_cast<float*>(buf.ptr);
    
    Eigen::VectorXf result(size);
    for (size_t i = 0; i < size; i++) {
        result(i) = ptr[i];
    }
    
    return result;
}

class PyXbarSimulator {
private:
    CrossbarSimulator simulator;

public:
    PyXbarSimulator(int M, int N) : simulator(M, N) {}
    
    void set_weights(py::array_t<bool> weights) {
        auto weights_vec = numpy_to_vector2d<bool>(weights);
        simulator.SetRRAM(weights_vec);
    }
    
    py::array_t<float> run_inference(py::array_t<float> input, float dt = simulation_time_step) {
        if (input.ndim() != 1) {
            throw std::runtime_error("Input must be a 1D array");
        }
        
        auto buf = input.request();
        float* ptr = static_cast<float*>(buf.ptr);
        size_t M = simulator.M;
        size_t N = simulator.N;
        
        // Set up input voltages
        Eigen::VectorXf Vappwl1 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappwl2 = Eigen::VectorXf::Zero(M);
        Eigen::VectorXf Vappbl1 = Eigen::VectorXf::Zero(N);
        Eigen::VectorXf Vappbl2 = Eigen::VectorXf::Zero(N);
        
        // Set access transistors based on input
        for (size_t j = 0; j < M; j++) {
            if (j < buf.shape[0] && ptr[j] != 0) {
                Vappwl1(j) = voltage_pulse_height;
                simulator.access_transistors[j] = std::vector<bool>(N, true);
            } else {
                simulator.access_transistors[j] = std::vector<bool>(N, false);
            }
        }
        
        // Initial voltage guess
        Eigen::VectorXf Vguess = Eigen::VectorXf::Zero(2*M*N);
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                Vguess(i*N + j) = Vappwl1(i);
            }
        }
        
        // Apply voltage and get output currents
        auto currents = simulator.ApplyVoltage(Vguess, Vappwl1, Vappwl2, Vappbl1, Vappbl2, dt);
        
        // Calculate MAC outputs (sum of currents per column)
        std::vector<float> mac_outputs(N, 0.0f);
        for (size_t n = 0; n < N; n++) {
            for (size_t m = 0; m < M; m++) {
                mac_outputs[n] += currents[m][n];
            }
        }
        
        // Convert to numpy array
        std::vector<size_t> shape = {N};
        py::array_t<float> result(shape);
        auto result_buf = result.request();
        float* result_ptr = static_cast<float*>(result_buf.ptr);
        
        for (size_t i = 0; i < N; i++) {
            result_ptr[i] = mac_outputs[i];
        }
        
        return result;
    }
    
    py::array_t<float> run_multiple_inferences(py::array_t<float> input, int num_inferences, float dt = simulation_time_step) {
        if (input.ndim() != 1) {
            throw std::runtime_error("Input must be a 1D array");
        }
        
        size_t N = simulator.N;
        std::vector<std::vector<float>> mac_outputs_all_iterations;
        
        for (int iteration = 0; iteration < num_inferences; iteration++) {
            auto mac_outputs = run_inference(input, dt);
            
            // Convert 1D numpy array to vector
            auto buf = mac_outputs.request();
            float* ptr = static_cast<float*>(buf.ptr);
            std::vector<float> mac_outputs_vec(ptr, ptr + N);
            
            mac_outputs_all_iterations.push_back(mac_outputs_vec);
        }
        
        return vector2d_to_numpy<float>(mac_outputs_all_iterations);
    }
    
    // Expose additional parameters and methods
    void set_parasitic_resistances(float Rswl1, float Rswl2, float Rsbl1, float Rsbl2, float Rwl, float Rbl) {
        simulator.Rswl1 = Rswl1;
        simulator.Rswl2 = Rswl2;
        simulator.Rsbl1 = Rsbl1;
        simulator.Rsbl2 = Rsbl2;
        simulator.Rwl = Rwl;
        simulator.Rbl = Rbl;
        
        // Recompute partial_G_ABCD
        simulator.partial_G_ABCD = PartiallyPrecomputeG_ABCD(
            simulator.M, simulator.N, 
            simulator.Rswl1, simulator.Rswl2, 
            simulator.Rsbl1, simulator.Rsbl2, 
            simulator.Rwl, simulator.Rbl
        );
    }
};

PYBIND11_MODULE(xbar_simulator, m) {
    m.doc() = "Python bindings for the memristor crossbar simulator";
    
    py::class_<PyXbarSimulator>(m, "CrossbarSimulator")
        .def(py::init<int, int>(), py::arg("M"), py::arg("N"))
        .def("set_weights", &PyXbarSimulator::set_weights, 
             "Set the weights of the crossbar (boolean matrix)")
        .def("run_inference", &PyXbarSimulator::run_inference, 
             py::arg("input"), py::arg("dt") = simulation_time_step,
             "Run a single inference with the given input vector")
        .def("run_multiple_inferences", &PyXbarSimulator::run_multiple_inferences,
             py::arg("input"), py::arg("num_inferences"), py::arg("dt") = simulation_time_step,
             "Run multiple inferences with the same input vector")
        .def("set_parasitic_resistances", &PyXbarSimulator::set_parasitic_resistances,
             py::arg("Rswl1"), py::arg("Rswl2"), py::arg("Rsbl1"), py::arg("Rsbl2"), 
             py::arg("Rwl"), py::arg("Rbl"),
             "Set the parasitic resistances of the crossbar");
    
    // Expose simulation settings
    m.attr("voltage_pulse_height") = voltage_pulse_height;
    m.attr("simulation_time_step") = simulation_time_step;
} 