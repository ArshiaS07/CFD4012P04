#pragma once
#include <iostream>
#include <sstream>
#include <vector>
#include <functional>
#include <execution>
#include <algorithm>
#include <unordered_map>
#include <any>
#include <fstream>
#include <omp.h>

#include <chrono>


struct Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<float> duration;
    std::string message;

    Timer(const std::string& msg) : message(msg) {
        start = std::chrono::high_resolution_clock::now();
    }

    ~Timer() {
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        float ms = duration.count() * 1000.0f;
        std::cout << message << " took " << ms << "ms\n";
    }
};

class WaveEquation1D {
public:

    WaveEquation1D();

    WaveEquation1D(std::pair<double, double> xLimits, size_t nodeCounts
        , double deltaT, double deltaX, double c, std::function<double(double)> initialConditionFunction);
    
    void createProblem();

    void updateU(const size_t, const std::string);

    void initializeDefaultSchemes();

    //Accessors

    double getValue(double _x, double _t);
    // Returns std::pair object with firsr = x_min and second = x_max
    std::pair<double, double> getXLimits() const;
    // Returns total node count
    size_t getNodeCount() const;
    // Returns delta t (time step)
    double getDeltaT() const;
    // Returns delta x (x step)
    double getDeltaX() const;
    // Returns the left boundary value
    double getLeftBoundaryValue() const;
    // Returns the right boundary value
    double getRightBoundaryValue() const;
    // Returns wave speed c (coeff. of partial(u)/partial(x))
    double getC() const;
    // Returns initial condition as a function of x
    std::function<double(double)> getInitialConditionFunctionOfX() const;
    // Saves u(x,t) for all x values at current time step in proper format to be plotted by gnuplot
    void saveToFileGnuplot(const std::string&);

    //Set as std::make_pair(x_min, x_max)
    void setXLimits(std::pair<double, double>);

    void setRightBoundaryValue(double);

    void setLeftBoundaryValue(double);
    //Set node count
    void setNodeCount(size_t);
    //Set delta t (time step)
    void setDeltaT(double);
    //Set delta x (x step)
    void setDeltaX(double);
    //Set wave speed c (coeff. of partial(u)/partial(x))
    void setC(double);
    //Set initial condition as a function of x
    void setInitialConditionFunctionOfX(std::function<double(double)>);
    //void setInitialCondition(double(*)(double));

private:
    std::pair<double, double> m_xLimits;
    double m_deltaX;
    double m_deltaT;
    double C;

    std::function<double(double)> m_initialConditionFunctionOfX;
    double m_leftBoundaryValue;
    double m_rightBoundaryValue;
    size_t m_nodeCount;

    std::vector<double> m_uNextValues;
    std::vector<double> m_uPrevValues;
    std::vector<double> m_uCurrValues;
    std::vector<double> m_xValues;
    std::vector<double> m_tValues;

    size_t m_currentTimeStep = 0;

    //double (*m_initialCondition)(double _x);


    // Allocates enough memory for all vectors
    void allocateMemory();  
    // Sets m_xValues with proper x values
    void discretizeDomain();
    // Sets m_uCurrValues.at(i) with m_initialConditionFunctionOfX(m_xValues.at(i))
    void setInitialValues();

    // Unorderes map to store all schemes.
    //The key is method names as strings.
    //The value is a function taking current u values returning next u values.
    std::unordered_map<std::string, std::function<std::vector<double>(std::vector<double>)> > m_schemes;


    //Upwind scheme. Takes Current u values and returns next u values.
    inline std::vector<double> scheme_upwind(const std::vector<double>);
    //Lax scheme. Takes Current u values and returns next u values.
    inline std::vector<double> scheme_lax(const std::vector<double>);
    //Lax-Wendroff scheme. Takes Current u values and returns next u values.
    inline std::vector<double> scheme_lax_wendroff(const std::vector<double>);

    inline std::vector<double> scheme_leapfrog(const std::vector<double>);

    void freeMemory();
};