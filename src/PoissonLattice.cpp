#include "PoissonLattice.hpp"

PoissonLattice::PoissonLattice(int xRange, int yRange, int zRange, double permativity, double dx): m_xRange(xRange),
																								   m_yRange(yRange),
																								   m_zRange(zRange),
																								   m_permativity(permativity),
																								   m_dx(dx),
																								   m_chargeDensity(xRange * yRange * zRange, 0.0),
																								   m_potential((xRange+1)*(yRange+1)*(zRange+1), 0.0)
{
	
}

void PoissonLattice::initialise(double initialValue, double noise, std::default_random_engine &generator)
{
	// Create the uniform distribution for generating the random numbers.
	std::uniform_real_distribution<double> noiseDistribution(-noise, noise);
	for(int k = 0; k < m_zRange; ++k)
	{
		for(int j = 0; j < m_yRange; ++j)
		{
			for(int i = 0; i < m_xRange; ++i )
			{
				(*this)(i,j,k) = initialValue + noiseDistribution(generator);
			}
		}
	}

}

double& PoissonLattice::operator()(int i, int j, int k)
{
	// Increment each index so that we take into account the halo of boundary conditions.
	i+=1;
	j+=1;
	k+=1;

	return m_potential[i + j*m_xRange + k*m_xRange*m_yRange];

}

const double& PoissonLattice::operator()(int i, int j, int k) const
{
	// Increment each index so that we take into account the halo of boundary conditions.
	i+=1;
	j+=1;
	k+=1;

	return m_potential[i + j*m_xRange + k*m_xRange*m_yRange];


}

double PoissonLattice::getChargeDensity(int i, int j, int k) const
{
	return m_chargeDensity[i + j*m_xRange + k*m_xRange*m_yRange];
}


void PoissonLattice::setChargeDensity(int i, int j, int k, double charge)
{
	m_chargeDensity[i + j*m_xRange + k*m_xRange*m_yRange] = charge;
}

double PoissonLattice::nextValueJacobi(int i, int j, int k) const
{
	return (((*this)(i+1,j,k) + (*this)(i-1,j,k) 
		         + (*this)(i,j+1,k) + (*this)(i,j-1,k)
		         + (*this)(i,j,k+1) + (*this)(i,j,k-1) + (pow(m_dx,2)/m_permativity) * getChargeDensity(i,j,k))/6.0);

}


double jacobiUpdate(PoissonLattice &currentLattice, PoissonLattice &updatedLattice)
{
	double convergenceMeasure = 0;

	for(int i = 0; i < currentLattice.m_xRange; ++i)
	{
		for(int j = 0; j < currentLattice.m_yRange; ++j)
		{
			for(int k = 0; k < currentLattice.m_zRange; ++k)
			{
				updatedLattice(i,j,k) = currentLattice.nextValueJacobi(i,j,k);

				convergenceMeasure += abs(updatedLattice(i,j,k)-currentLattice(i,j,k));

			}
		}
	}

	return convergenceMeasure;

}

double gaussSeidelUpdate(PoissonLattice &lattice)
{
	double convergenceMeasure = 0;

	double updatedValue = 0;
	double currentValue = 0;

	for(int i = 0; i < lattice.m_xRange; ++i)
	{
		for(int j = 0; j < lattice.m_yRange; ++j)
		{
			for(int k = 0; k < lattice.m_zRange; ++k)
			{

				updatedValue = lattice.nextValueJacobi(i,j,k);

				currentValue = lattice(i,j,k);

				convergenceMeasure += abs(updatedValue-currentValue);

				lattice(i,j,k) = updatedValue;

			}
		}
	}

	return convergenceMeasure;

}

double sorUpdate(double sorParameter, PoissonLattice &lattice)
{
	// Local variable to measure how much the lattice changes on this update.
	double convergenceMeasure = 0;

	// Local variable to hold the standard GS update.
	double updatedGSValue = 0;

	// Local variable to hold current value at given lattice site.
	double currentValue = 0;

	// Local variable to hold the final update value under the SOR update.
	double updatedSORValue = 0;

	for(int i = 0; i < lattice.m_xRange; ++i)
	{
		for(int j = 0; j < lattice.m_yRange; ++j)
		{
			for(int k = 0; k < lattice.m_zRange; ++k)
			{
				currentValue = lattice(i,j,k);

				updatedGSValue = lattice.nextValueJacobi(i,j,k);

				updatedSORValue = (1-sorParameter) * currentValue + sorParameter * updatedGSValue;

				lattice(i,j,k) = updatedSORValue;

				convergenceMeasure += abs(updatedSORValue-currentValue);

			}
		}
	}

	return convergenceMeasure;

}


double latticeDifference(PoissonLattice &lattice1, PoissonLattice &lattice2)
{
	double deltaPhi = 0;
	for(int k = 0; k < lattice1.m_zRange; ++k)
	{
		for(int j = 0; j < lattice1.m_yRange; ++j)
		{
			for(int i = 0; i < lattice1.m_xRange; ++i )
			{
				deltaPhi += abs(lattice1(i,j,k) - lattice2(i,j,k));
			}
		}
	}

	return deltaPhi;
}

void PoissonLattice::printPotential(std::ostream &out)
{
	for(int k = 0; k < m_zRange; ++k)
	{
		for(int j = 0; j < m_yRange; ++j)
		{
			for(int i = 0; i < m_xRange; ++i )
			{
				out << i << ' ' << j << ' ' << k << ' ' << (*this)(i,j,k) <<'\n';
			}

			out << '\n';
		}

		out << '\n';
	}

}

std::array<double,3> PoissonLattice::electricField(int i, int j, int k) const
{
	 std::array<double,3> electricField = {-((*this)(i+1,j,k)-(*this)(i-1,j,k))/(2*m_dx),
								-((*this)(i,j+1,k)-(*this)(i,j-1,k))/(2*m_dx),
								-((*this)(i,j,k+1)-(*this)(i,j,k-1))/(2*m_dx)};

	return electricField;
}

std::array<double,3> PoissonLattice::magneticField(int i, int j, int k) const
{
	std::array<double,3> magneticField = {((*this)(i,j+1,k)-(*this)(i,j-1,k))/(2*m_dx),
							-((*this)(i+1,j,k)-(*this)(i-1,j,k))/(2*m_dx),
							0};

	return magneticField;

}

void PoissonLattice::printElectricField(std::ostream& out) const
{
	std::array<double,3> electricFieldTemp;
	for(int k = 0; k < m_zRange; ++k)
	{
		for(int j = 0; j < m_yRange; ++j)
		{
			for(int i = 0; i < m_xRange; ++i )
			{
				electricFieldTemp = electricField(i, j, k);
				out << i << ' ' << j << ' ' << k << ' ' << 
				electricFieldTemp[0] << ' ' << electricFieldTemp[1] << ' ' << electricFieldTemp[2] << '\n';
			}

			out << '\n';
		}

		out << '\n';
	}

}

std::ostream& operator<<(std::ostream &out, const PoissonLattice &lattice)
{
	// Calculate the centre of the lattice utilising integer division.
	double xCentre = lattice.m_xRange/2;
	double yCentre = lattice.m_yRange/2;
	double zCentre = lattice.m_zRange/2;


	std::array<double,3> electricFieldTemp;
	std::array<double,3> magneticFieldTemp;
	for(int k = 0; k < lattice.m_zRange; ++k)
	{
		for(int j = 0; j < lattice.m_yRange; ++j)
		{
			for(int i = 0; i < lattice.m_xRange; ++i )
			{
				// Calculate distance of current point from centre of lattice.
				double xDistance = xCentre - i;
				double yDistance = yCentre - j;
				double zDistance = zCentre - k;


				electricFieldTemp = lattice.electricField(i, j, k);
				magneticFieldTemp = lattice.magneticField(i,j,k);

				out << i << ' ' << j << ' ' << k << ' ' << 
				sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance) << ' ' << lattice(i,j,k) << 
				' ' << electricFieldTemp[0] << ' ' << electricFieldTemp[1] << ' ' << electricFieldTemp[2] << 
				' ' << magneticFieldTemp[0] << ' ' << magneticFieldTemp[1] << ' ' << magneticFieldTemp[2] << '\n';
			}

			out << '\n';
		}

		out << '\n';
	}

	return out;

}