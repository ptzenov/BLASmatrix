#include <matrix.hpp>
#include <iostream>


int main() {

	
	dat::Matrix<float> M2(2,2,0); // row major
	dat::Matrix<float> M3(2,2,1); // col major 
	for(int i= 0; i< 2; i ++)
		for(int j = 0; j<2;j++){
			M2(i,j) = i+j*2;
			M3(i,j) = i+j*2; 
	}


	std::cout << M2 << "\n";
	std::cout << M3 << "\n"; 
	
	std::cout << M2*M2 << "\n"; 
	std::cout << M3*M3 << "\n";
	std::cout << M2*M3 << "\n";
	std::cout << M3*M2 << "\n";

}

