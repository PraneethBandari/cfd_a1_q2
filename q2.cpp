#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
	// Define problem parameters
	double L = 0.1;          // Wall thickness (m)
	int N = 50;              // Number of nodes
	double dx = L / (N - 1); // Step size (m)
	double T_hotfluid = 300; // Temperature of hot fluid (B0C)
	double T_right = 50;     // Temperature at the right boundary (B0C)
	double k = 0.1;          // Thermal conductivity (W/m-K)
	double h = 25;           // Heat transfer coefficient (W/mB2-K)
	int q = 200;             // Heat flux at the left boundary (W/mB2)
	double accuracy = 1e-6;  // Convergence criterion

	// Initialize x positions for each node
	vector<double> x;
	for (int i = 0; i < N; i++) {
		double d = (L * i) / (N - 1);
		x.push_back(d);
	}

	// Analytical temperature distribution (for comparison)
	vector<double> T_analytical;
	double T_left = 298.076;  // Precomputed temperature at left boundary (B0C)
	for (int i = 0; i < N; i++) {
		double T = T_left + (T_right - T_left) * x[i] / L;
		T_analytical.push_back(T);
		cout << "Analytical temperature at node " << i + 1 << " is " << T_analytical[i] << "\n";
	}

	// Initialize FDM temperature array with an initial guess
	vector<double> Temp(N, (T_right + T_hotfluid) / 2);
	vector<double> T_new(N);

	// Set the right boundary condition (fixed temperature)
	T_new[N - 1] = T_right;

	double max_dt = 1; // Maximum temperature difference (for convergence check)
	int n = 0;         // Iteration counter

	// Main iterative loop
	while (max_dt > accuracy) {
		max_dt = 0;

		// Apply boundary condition at the left face (x = 0)
		T_new[0] = ((k / dx) * Temp[1] + h * T_hotfluid + q) / (k / dx + h);

		// Update interior nodes using the FDM approximation
		for (int i = 1; i < N - 1; i++) {
			T_new[i] = (Temp[i + 1] + Temp[i - 1]) * 0.5;
		}

		// Check for convergence and update temperatures
		for (int i = 0; i < N; i++) {
			double dt = fabs(Temp[i] - T_new[i]);
			if (dt > max_dt) {
				max_dt = dt;
			}
			Temp[i] = T_new[i];
		}

		// Increment iteration counter
		n++;
	}

	// Output the number of iterations and final temperatures
	cout << "Number of iterations: " << n << "\n";
	for (int i = 0; i < N; i++) {
		cout << "Node " << i + 1 << ": Analytical temp = " << T_analytical[i] << " B0C, FDM temp = " << Temp[i] << " B0C\n";
	}
	ofstream myfile("fdmq2sol.txt");

	if (myfile.is_open()) {

		for (int i = 0; i < N; i++) {
			myfile<< i + 1 << ":" << Temp[i] << endl;
		}
		myfile.close();
		cout << "Temperature values have been written to fdmq2sol.txt" << endl;
	} else {
		cout << "Error opening file!" << endl;
	}
	ofstream outfile("fdmq2analyticalsol.txt");

	if (outfile.is_open()) {

		for (int i = 0; i < N; i++) {
			outfile << i + 1 << ":" << T_analytical[i] << endl;
		}
		outfile.close();
		cout << "Temperature values have been written to fdmq2analyticalsol.txt" << endl;
	} else {
		cout << "Error opening file!" << endl;
	}


	return 0;
}
