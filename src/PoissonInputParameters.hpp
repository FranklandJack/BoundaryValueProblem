#ifndef PoissonInputParameters_hpp
#define PoissonInputParameters_hpp
#include <iostream>
#include <iomanip>
/**
 *\file 
 *\class PoissonInputParameters
 *\brief Class for easily handling input parameters of Poisson differential equation solver.
 *
 * This class essentially just holds some values and has an operator to easily output
 * them to a stream.
 */
class PoissonInputParameters 
{

public:

    /**
     *\enum to hold the three solution methods we will consider when solving the Poisson equation
     */
    enum SolutionMethod
    {
        Jacobi,
        GaussSeidel,
        SOR
    };

    /**
     *\enum to hold the problem we are going to solve
     */
    enum ProblemSolved
    {
    	Electro,
    	Magneto
    };

    /// Solution method.
    SolutionMethod solutionMethod;

    /// Problem to be solved.
    ProblemSolved problem;

	/// Spatial discretisation step.
    double spaceStep;

    /// Permittivity in the Poisson equation.
    double permittivity;

    /// \phi_0 initial value of potential.
    double initialValue;

    /// Maximum magnitude of initial noise.
    double noise;

    /// Precision of the final answer in terms of convergence.
    double precision;

    /// Range of x-values lattice domain.
    int xRange;

    /// Range of y-values lattice domain.
    int yRange;

    /// Range of z-values lattice domain.
    int zRange;

    /// Name of output directory to save any output into.
    std::string outputName;

    /// Successive over relaxation parameter.
    double sorParameter;

    /** 
	 *\brief operator<< overload for outputting the results.
	 *\param out std::ostream reference that is the stream being outputted to.
	 *\param params constant PoissonInputParameters instance to be output.
	 *\return std::ostream reference so the operator can be chained.
	 *
	 * Results will be output in a formatted table for easy viewing in the command line or a file.
	 */
    friend std::ostream& operator<<(std::ostream& out, const PoissonInputParameters& params);

};
#endif /* PoissonInputParameters_hpp */