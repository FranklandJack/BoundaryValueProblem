#include "PoissonLattice.hpp"

PoissonLattice::PoissonLattice(int xRange, int yRange, int zRange, double permativity, double dx): m_xRange(xRange),
																								   m_yRange(yRange),
																								   m_zRange(zRange),
																								   m_permativity(permativity),
																								   m_dx(dx),
																								   m_chargeDensity(xRange * yRange * zRange, 0.0),
																								   m_potential(xRange*yRange*zRange, 0.0)
{

}

void PoissonLattice::initialise(double initialValue, double noise, std::default_random_engine &generator)
{
	// Create the uniform distribution for generating the random numbers.
	std::uniform_real_distribution<double> noiseDistribution(-noise, noise);
	// Only update from 1 to range-1 since boundaries should be fixed by initial conditions.
	for(int k = 1; k < m_zRange-1; ++k)
	{
		for(int j = 1; j < m_yRange-1; ++j)
		{
			for(int i = 1; i < m_xRange-1; ++i )
			{
				(*this)(i,j,k) = initialValue + noiseDistribution(generator);
			}
		}
	}

}

double& PoissonLattice::operator()(int i, int j, int k)
{
	return m_potential[i + j*m_xRange + k*m_xRange*m_yRange];
}

const double& PoissonLattice::operator()(int i, int j, int k) const
{
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
		         + (*this)(i,j,k+1) + (*this)(i,j,k-1) + (std::pow(m_dx,2)/m_permativity) * getChargeDensity(i,j,k))/6.0);

}


double jacobiUpdate(PoissonLattice &currentLattice, PoissonLattice &updatedLattice)
{
	double convergenceMeasure = 0;

	// Only need to check from 1 to range-1 since bounaries are fixed.
	for(int i = 1; i < currentLattice.m_xRange-1; ++i)
	{
		for(int j = 1; j < currentLattice.m_yRange-1; ++j)
		{
			for(int k = 1; k < currentLattice.m_zRange-1; ++k)
			{
				updatedLattice(i,j,k) = currentLattice.nextValueJacobi(i,j,k);

				convergenceMeasure += std::abs(updatedLattice(i,j,k)-currentLattice(i,j,k));

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

	// Only update from 1 to range-1 since boundaries should be fixed by initial conditions.
	for(int i = 1; i < lattice.m_xRange-1; ++i)
	{
		for(int j = 1; j < lattice.m_yRange-1; ++j)
		{
			for(int k = 1; k < lattice.m_zRange-1; ++k)
			{

				updatedValue = lattice.nextValueJacobi(i,j,k);

				currentValue = lattice(i,j,k);

				convergenceMeasure += std::abs(updatedValue-currentValue);

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

	// Only update from 1 to range-1 since boundaries should be fixed by initial conditions.
	for(int i = 1; i < lattice.m_xRange-1; ++i)
	{
		for(int j = 1; j < lattice.m_yRange-1; ++j)
		{
			for(int k = 1; k < lattice.m_zRange-1; ++k)
			{
				currentValue = lattice(i,j,k);

				updatedGSValue = lattice.nextValueJacobi(i,j,k);

				updatedSORValue = (1-sorParameter) * currentValue + sorParameter * updatedGSValue;

				lattice(i,j,k) = updatedSORValue;

				convergenceMeasure += std::abs(updatedSORValue-currentValue);

			}
		}
	}

	return convergenceMeasure;

}





std::array<double,3> PoissonLattice::electricField(int i, int j, int k) const
{
	 std::array<double,3> electricField = {-((*this)(i+1,j,k)-(*this)(i-1,j,k))/(2*m_dx),
								-((*this)(i,j+1,k)-(*this)(i,j-1,k))/(2*m_dx),
								-((*this)(i,j,k+1)-(*this)(i,j,k-1))/(2*m_dx)};

	return electricField;
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
				double radialDistance = std::sqrt(xDistance*xDistance+yDistance*yDistance+zDistance*zDistance);

				double fieldStrength;



				if(i==0 || j==0 || k==0 || i==lattice.m_xRange-1 || j==lattice.m_yRange-1 || k==lattice.m_zRange-1)
				{
					electricFieldTemp = std::array<double,3>{0,0,0};
				}
				else
				{
					electricFieldTemp = lattice.electricField(i, j, k);
				}

				fieldStrength = std::sqrt(electricFieldTemp[0]*electricFieldTemp[0] + electricFieldTemp[1]*electricFieldTemp[1] + electricFieldTemp[2]*electricFieldTemp[2]);

				out << i << ' ' << j << ' ' << k << ' ' <<
				radialDistance << ' ' << lattice(i,j,k) <<
				' ' << electricFieldTemp[0] << ' ' << electricFieldTemp[1] << ' ' << electricFieldTemp[2] <<
				' ' << fieldStrength << ' ' << '\n';
			}

			out << '\n';
		}

		out << '\n';
	}

	return out;

}


 void PoissonLattice::setPointChargeDist()
 {
	     // Utilise integer division to find the centre of the box.
	     int xCentre = m_xRange/2;
	     int yCentre = m_yRange/2;
	     int zCentre = m_zRange/2;

			 // Point charge magnitude.
	     double deltaCharge = 1;

			 // Set the charge.
			 setChargeDensity(xCentre, yCentre, zCentre, deltaCharge);
 }
