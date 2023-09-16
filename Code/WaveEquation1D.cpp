#include "WaveEquation1D.h"

WaveEquation1D::WaveEquation1D() {

}

WaveEquation1D::WaveEquation1D(std::pair<double, double> xLimits, size_t nodeCount
    , double deltaT, double deltaX, double c, std::function<double(double)> initialConditionFunction)
    : m_xLimits(xLimits), m_nodeCount(nodeCount), m_deltaT(deltaT), m_deltaX(deltaX)
    , C(c), m_initialConditionFunctionOfX(initialConditionFunction)
{
    allocateMemory();
    discretizeDomain();
    setInitialValues();
    initializeDefaultSchemes();
}

void WaveEquation1D::allocateMemory() {
    //double* p_uNexValues  = new double [m_nodeCount];
    //double* p_uPrevValues = new double [m_nodeCount];
    //double* p_uCurrValues = new double [m_nodeCount];
    //double* p_xValues     = new double [m_nodeCount];
    //double* p_tValues     = new double [m_nodeCount];

    m_uNextValues.assign(m_nodeCount, 0);
    m_uPrevValues.assign(m_nodeCount, 0);
    m_uCurrValues.assign(m_nodeCount, 0);
    m_xValues.assign(m_nodeCount, 0);
    m_tValues.assign(m_nodeCount, 0);


}

void WaveEquation1D::freeMemory() {
    //delete [] p_uNexValues; 
    //delete [] p_uPrevValues;
    //delete [] p_uCurrValues;
    //delete [] p_xValues;
    //delete [] p_tValues;

    m_uNextValues.clear();
    m_uPrevValues.clear();
    m_uCurrValues.clear();
    m_xValues.clear();
    m_tValues.clear();
}

void WaveEquation1D::discretizeDomain() {
    //Discretize by x_i = a + i * deltaX
    //This function sets m_deltaX as well as m_xValues
    m_deltaX = (m_xLimits.second - m_xLimits.first) / (double(m_nodeCount) - 1);

    for (int i = 0; i < m_nodeCount; i++) {
        m_xValues.at(i) = m_xLimits.first + i * m_deltaX;
    }
}

void WaveEquation1D::setInitialValues() {
    if (!m_initialConditionFunctionOfX) {
        throw std::runtime_error("No initial condition function set.");
    }

    for (size_t i = 0; i < m_nodeCount; ++i) {
        double x = m_xValues.at(i);
        m_uCurrValues.at(i) = m_initialConditionFunctionOfX(x);
    }
}

void WaveEquation1D::initializeDefaultSchemes() {
    m_schemes["upwind"] = [&](std::vector<double> _uCurrValues) -> std::vector<double> {
        return scheme_upwind(_uCurrValues);
    };
    m_schemes["lax"] = [&](std::vector<double> _uCurrValues) -> std::vector<double> {
        return scheme_lax(_uCurrValues);
    };
    m_schemes["lax-wendroff"] = [&](std::vector<double> _uCurrValues) -> std::vector<double> {
        return scheme_lax_wendroff(_uCurrValues);
    };
    m_schemes["leapfrog"] = [&](std::vector<double> _uCurrValues) -> std::vector<double> {
        return scheme_leapfrog(_uCurrValues);
    };

}

void WaveEquation1D::createProblem() {
    allocateMemory();
    discretizeDomain();
    setInitialValues();
    initializeDefaultSchemes();
}

void WaveEquation1D::updateU(const size_t _update_time_in_timestep, const std::string _numerical_scheme) {

    Timer time("updateU()");

    // n_new = n_current + _tSteps
    std::cout << "Looking for scheme: " << _numerical_scheme << std::endl;
    auto schemeIter = m_schemes.find(_numerical_scheme);
    if (schemeIter == m_schemes.end()) {
        throw std::runtime_error("Scheme not found: " + _numerical_scheme);
    }
    else { std::cout << "method formula found\n";}

    for (size_t tSteps = 0; tSteps < _update_time_in_timestep; tSteps++) {
        // Update all nodes between using the specified numerical scheme
        m_uNextValues = schemeIter->second(m_uCurrValues);
        // Update first node (m_uNextValues.at(0)) using left boundary condition
        m_uNextValues.at(0) = m_leftBoundaryValue;
        // Update last node (m_uNextValues.at(m_nodeCount - 1)) using right boundary condition
        m_uNextValues.at(m_nodeCount - 1) = m_rightBoundaryValue;
        // Move data
        std::copy(std::execution::par, m_uCurrValues.begin(), m_uCurrValues.end(), m_uPrevValues.begin());
        std::copy(std::execution::par, m_uNextValues.begin(), m_uNextValues.end(), m_uCurrValues.begin());
        m_currentTimeStep++;
    }
    
};

double WaveEquation1D::getValue(double _x, double _t) {
    return m_uCurrValues.at((_x - m_xLimits.first) / (m_deltaX));
}

void WaveEquation1D::saveToFileGnuplot(const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    // Write the output to the file
    for (size_t i = 0; i < m_nodeCount; i++) {
        outputFile << m_xValues[i] << ' ' << m_uCurrValues[i] << '\n';
    }

    // Close the file
    outputFile.close();
}

//Default schemes
inline std::vector<double> WaveEquation1D::scheme_upwind(const std::vector<double> _uCurrValues) {
    std::vector<double> uNextValues(m_nodeCount, 0);
    for (size_t i = 1; i < m_nodeCount - 1; i++) {
        uNextValues.at(i) = _uCurrValues.at(i) - C * m_deltaT * (_uCurrValues.at(i) - _uCurrValues.at(i - 1)) / m_deltaX;
    }
    return uNextValues;
}
inline std::vector<double> WaveEquation1D::scheme_lax(const std::vector<double> _uCurrValues) {
    std::vector<double> uNextValues(m_nodeCount, 0);
    double average_u;
        for (size_t i = 1; i < m_nodeCount - 1; i++) {
        average_u = _uCurrValues.at(i + 1) + _uCurrValues.at(i - 1);
        average_u /= 2;
        uNextValues.at(i) = _uCurrValues.at(i) - C * m_deltaT * (_uCurrValues.at(i) - _uCurrValues.at(i - 1)) / m_deltaX;
    }
    return uNextValues;
}
inline std::vector<double> WaveEquation1D::scheme_lax_wendroff(const std::vector<double> _uCurrValues) {
    std::vector<double> uNextValues(m_nodeCount, 0);
    double first_order_derivative_centeral;
    double second_order_derivative_centeral;
    for (size_t i = 1; i < m_nodeCount - 1; i++) {
        first_order_derivative_centeral = (_uCurrValues.at(i + 1) - _uCurrValues.at(i - 1)) / (2 * m_deltaX);
        second_order_derivative_centeral = (_uCurrValues.at(i + 1) - 2 * _uCurrValues.at(i) + _uCurrValues.at(i - 1)) / (m_deltaX * m_deltaX);
        uNextValues.at(i) = _uCurrValues.at(i) - C * m_deltaT * first_order_derivative_centeral + C * C * m_deltaT * m_deltaT
            * second_order_derivative_centeral / 4;
    }
    return uNextValues;
}
inline std::vector<double> WaveEquation1D::scheme_leapfrog(const std::vector<double> _uCurrValues) {
    std::vector<double> uNextValues(m_nodeCount, 0);
    // Compute the next values using the leapfrog scheme
    for (size_t i = 1; i < m_nodeCount - 1; i++) {
        uNextValues.at(i) = _uCurrValues.at(i) + m_deltaT * C * C / (m_deltaX * m_deltaX)
            * (_uCurrValues.at(i + 1) - 2.0 * _uCurrValues.at(i) + _uCurrValues.at(i - 1));
    }
    return uNextValues;
}


std::pair<double, double> WaveEquation1D::getXLimits() const {
    return m_xLimits;
}

size_t WaveEquation1D::getNodeCount() const {
    return m_nodeCount;
}

double WaveEquation1D::getDeltaT() const {
    return m_deltaT;
}

double WaveEquation1D::getDeltaX() const {
    return m_deltaX;
}

double WaveEquation1D::getLeftBoundaryValue() const {
    return m_leftBoundaryValue;
};

double WaveEquation1D::getRightBoundaryValue() const {
    return m_rightBoundaryValue;
};

double WaveEquation1D::getC() const {
    return C;
}

std::function<double(double)> WaveEquation1D::getInitialConditionFunctionOfX() const {
    return m_initialConditionFunctionOfX;
}

void WaveEquation1D::setXLimits(std::pair<double, double> xLimits) {
    m_xLimits = xLimits;
    discretizeDomain();
    setInitialValues();
}

void WaveEquation1D::setRightBoundaryValue(double RightBoundaryValue) {
    m_rightBoundaryValue = RightBoundaryValue;
};

void WaveEquation1D::setLeftBoundaryValue(double LeftBoundaryValue) {
    m_leftBoundaryValue = LeftBoundaryValue;
};

void WaveEquation1D::setNodeCount(size_t nodeCount) {
    m_nodeCount = nodeCount;
    allocateMemory();
    discretizeDomain();
    setInitialValues();
}

void WaveEquation1D::setDeltaT(double deltaT) {
    m_deltaT = deltaT;
}

void WaveEquation1D::setDeltaX(double deltaX) {
    m_deltaX = deltaX;
    discretizeDomain();
    setInitialValues();
}

void WaveEquation1D::setC(double c) {
    C = c;
}

void WaveEquation1D::setInitialConditionFunctionOfX(std::function<double(double)> initialConditionFunction) {
    m_initialConditionFunctionOfX = initialConditionFunction;
    setInitialValues();
}