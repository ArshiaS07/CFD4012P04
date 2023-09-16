// CFD4012P04.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include "WaveEquation1D.h"
#include <exprtk.hpp>


///*

int main() {
    
    // Create WaveEquation1D object
    WaveEquation1D waveEquation;
    // Loop to accept user input
    while (true)
    {
        // Prompt user for simulation parameters
        std::cout << "Enter simulation parameters (leftBoundaryValue rightBoundaryValue nodeCount domainStart domainEnd simulationTime timestep wavespeed): \n";
        double leftBoundaryValue, rightBoundaryValue, domainStart, domainEnd, simulationTime, timeStep, c;
        size_t nodeCount;
        std::cin >> leftBoundaryValue >> rightBoundaryValue >> nodeCount >> domainStart >> domainEnd >> simulationTime >> timeStep >> c;

        // Prompt user for initial condition function expression
        std::cout << "Enter initial condition function expression (e.g. sin(x), x^2, etc.): ";
        // Create function string
        std::string expression_string;
        // Remove white spaces and store the function entered by user in expression_string .
        std::getline(std::cin >> std::ws, expression_string);
        // Create the function from string

        double x;
        exprtk::symbol_table <double> symbol_table;
        symbol_table.add_variable("x", x);
        symbol_table.add_constants();
        exprtk::expression<double> expression;
        expression.register_symbol_table(symbol_table);
        exprtk::parser<double> parser;
        if (!parser.compile(expression_string, expression)) {
            std::cerr << "Error: Failed to compile expression '" << expression_string << "'." << std::endl;
        }

        std::function<double(double)> initialFunc = [&](double _x) {
            x = _x;
            return expression.value();
        };

        waveEquation.setInitialConditionFunctionOfX(initialFunc);

        // Set simulation parameters
        waveEquation.setLeftBoundaryValue(leftBoundaryValue);
        waveEquation.setRightBoundaryValue(rightBoundaryValue);
        waveEquation.setNodeCount(nodeCount);
        waveEquation.setXLimits(std::make_pair(domainStart, domainEnd));
        waveEquation.setDeltaT(timeStep);
        waveEquation.setC(c);


        // Prompt user for numerical scheme and number of time steps
        std::cout << "Enter numerical scheme (upwind/leapfrog/lax-wendroff/lax): ";
        std::string numericalScheme;
        std::cin >> std::ws >> numericalScheme;

        // Prompt user for output file name
        std::cout << "Enter output file name (press enter for default name): ";
        std::string outputFileName;
        std::getline(std::cin >> std::ws, outputFileName);
        if (outputFileName.empty()) {
            outputFileName = "output.txt";
        }

        waveEquation.initializeDefaultSchemes();
        // Solve the wave equation
        waveEquation.updateU(size_t(simulationTime/waveEquation.getDeltaT()), numericalScheme);

        // Output results to file
        waveEquation.saveToFileGnuplot(outputFileName);

        std::cout << "Simulation completed. Output saved to " << outputFileName << "." << std::endl;

        // Prompt user to repeat simulation or exit
        std::cout << "Do you want to run another simulation? (Y/N): ";
        char repeat;
        std::cin >> repeat;
        if (repeat != 'Y' && repeat != 'y') {
            break;
        }
    }
    

    return 0;
}

//*/
/*
int main() {
    // Create a new instance of the WaveEquation1D class
    WaveEquation1D waveEquation;

    // Set the wave speed and initial condition
    waveEquation.setC(1.0);
    waveEquation.setInitialConditionFunctionOfX([](double x) { return std::sin(2.0 * M_PI * x); });
    // Set the domain limits, number of nodes, and time and space step sizes
    waveEquation.setXLimits(std::make_pair(0.0, 1.0));
    waveEquation.setNodeCount(100);
    waveEquation.setDeltaT(0.01);
    waveEquation.setDeltaX(0.01);

    

    // Set the left and right boundary values
    waveEquation.setLeftBoundaryValue(0.0);
    waveEquation.setRightBoundaryValue(0.0);

    //
    waveEquation.saveToFileGnuplot("atTime0.dat");

    waveEquation.initializeDefaultSchemes();

    // Solve the equation using the leapfrog scheme for 100 time steps
    waveEquation.updateU(10, "upwind");
    waveEquation.saveToFileGnuplot("atTime1s.dat");
    // Retrieve the solution at a point (x, t)
    double x = 0.1;
    double t = 1.0;
    double u = waveEquation.getValue(x, t);

    std::cout << "The solution at (x, t) = (" << x << ", " << t << ") is " << u << std::endl;

    return 0;
}
*/
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
