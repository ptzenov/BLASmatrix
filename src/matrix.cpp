/*
 * matrix.cpp
 *
 * Author: petz
 *
 */

#include <matrix.hpp>

using namespace dat;

std::ostream& operator<<(std::ostream& out, complex double nr) {
	out << "(" << creal(nr) << ";" << cimag(nr) << ")";
	return out;
}

std::ostream& operator<<(std::ostream& out, complex float nr) {
	out << "(" << creal(nr) << ";" << cimag(nr) << ")";
	return out;
}

template<typename _Tp>
void Matrix<_Tp>::registerview(unsigned int viewid, Matrix<_Tp>& view) {

	subscribers.insert(std::pair<unsigned int, Matrix<_Tp> >(viewid, view));

}

/*
 * Remove a view with a unique view_id from the list of subscribers. 
 * This method shall be called upon view destruction.
 */
template<typename _Tp>
void Matrix<_Tp>::removeview(unsigned int view_id) {
	subscribers.erase(view_id);
}

template<typename _Tp>
unsigned int Matrix<_Tp>::getDim_i() const {
	return _dimI;
}

template<typename _Tp>
unsigned int Matrix<_Tp>::getDim_j() const {
	return _dimJ;
}

template<typename _Tp> 
unsigned int Matrix<_Tp>::getLDim() const {
	return _ldim;
}

template<typename _Tp>
_Tp** const Matrix<_Tp>::getMtxData() const {
	return _data;
}

/**
 * the standard constructor - takes a row and column dimension.
 */
template<typename _Tp>
void Matrix<_Tp>::init(unsigned int dim_i, unsigned int dim_j, unsigned int ldim) {

	_dimI = dim_i;
	_dimJ = dim_j;

	assert(ldim <= 1); 

	_ldim = ldim; // leading dimension (0 for rows and 1 for cols) 

	// contiguous memory allocation  	
	_Tp* tmp_data = (_Tp*) calloc( _dimI * _dimJ, sizeof(_Tp));
	_data = (_Tp**) calloc((_ldim == 0 ? _dimI:_dimJ), sizeof(_Tp*));

	for (int i = 0; i < (_ldim == 0 ? _dimI:_dimJ) ; i++)
		_data[i] = tmp_data + i * (_ldim==0?_dimJ:_dimI);
	
	_bytes = _dimI * _dimJ * sizeof(_Tp);
	subscribers.clear();
}

template<typename _Tp>
Matrix<_Tp>::Matrix(unsigned int dim_i, unsigned int dim_j, unsigned int ldim) {
	init(dim_i, dim_j,ldim);
}

template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator()(unsigned int i1, unsigned int i2,
		unsigned int j1, unsigned int j2){
	//ToDo
}

/**
 *	Copy constructor
 */
template<typename _Tp>
Matrix<_Tp>::Matrix(const Matrix<_Tp>& arg) {

	init(arg._dimI, arg._dimJ, arg._ldim);

	switch (gettype<_Tp>()) {
		case FLT:
			cblas_scopy(getDim_i() * getDim_j(), (float*) arg.getMtxData()[0], 1,
					(float*) this->getMtxData()[0], 1);
			break;

		case DBL:
			cblas_dcopy(getDim_i() * getDim_j(), (double*) arg.getMtxData()[0], 1,
					(double*) this->getMtxData()[0], 1);
			break;
		case CPLXFLT:
			cblas_ccopy(getDim_i() * getDim_j(), (void*) arg.getMtxData()[0], 1,
					(void*) this->getMtxData()[0], 1);
			break;
		case CPLXDBL:
			cblas_zcopy(getDim_i() * getDim_j(), (void*) arg.getMtxData()[0], 1,
					(void*) this->getMtxData()[0], 1);
			break;
		default:
			throw std::domain_error(
					"Unsupported data type for matrix copy operation");
	}

}
/*
 * Destructor
 */
template<typename _Tp>
Matrix<_Tp>::~Matrix() {
	free(_data); // free the whole data!
}

// make the function call operator retrive the i,j th data element
template<typename _Tp>
_Tp& Matrix<_Tp>::operator()(unsigned int i, unsigned int j) const {

	if (i >= this->getDim_i() || i < 0) {
		DATASTRUCT_OUT_ERR(
				"Array index out of bounds. Cannot retrieve Matrix element. ",
				__FILE__, __LINE__)
			throw std::out_of_range("row index out of bounds.");

	}

	if (j >= this->getDim_j() || j < 0) {
		DATASTRUCT_OUT_ERR(
				"Array index out of bounds. Cannot retrieve Matrix element. ",
				__FILE__, __LINE__)
			throw std::out_of_range("column index out of bounds.");
	}

	if (_ldim == 0)
		return _data[i][j];
	else
		return _data[j][i];
}
/* 
 * overload the assignment operator!
 *
 */
template<typename _Tp>
void Matrix<_Tp>::operator=(const Matrix<_Tp>& arg) {

	if (this->getDim_i() != arg.getDim_i()
			|| this->getDim_j() != arg.getDim_j()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	if (_ldim != arg.getLDim()){
		DATASTRUCT_OUT_WRN("Matrix leading dimensions do not agree." 
				"Enforce slow copy", __FILE__,__LINE__)
			for(int i =0; i<_dimI;i++)
				for (int j = 0; j < _dimJ; j++)
					(*this)(i,j) = arg.getMtxData()[j][i];

	}else{
		switch (gettype<_Tp>()) {
			case FLT:
				cblas_scopy(getDim_i()*getDim_j(), (float*) arg.getMtxData()[0], 1,
						(float*) this->getMtxData()[0], 1);
				break;

			case DBL:
				cblas_dcopy(getDim_i()*getDim_j(), (double*) arg.getMtxData()[0], 1,
						(double*) this->getMtxData()[0], 1);
				break;
			case CPLXFLT:
				cblas_ccopy(getDim_i()*getDim_j(), (void*) arg.getMtxData()[0], 1,
						(void*) this->getMtxData()[0], 1);
				break;
			case CPLXDBL:
				cblas_zcopy(getDim_i()*getDim_j(), (void*) arg.getMtxData()[0], 1,
						(void*) this->getMtxData()[0], 1);
				break;
			default:
				throw std::domain_error(
						"Unsupported data type for matrix copy operation");
		}

	}
}

// overload the multiplication by matrix operator...
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator*(const Matrix<_Tp>& arg) const {

	if (_dimJ != arg.getDim_i()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	int L = arg.getDim_j();
	Matrix<_Tp> res(_dimI, L,_ldim);
	if (_ldim != arg.getLDim()){
		DATASTRUCT_OUT_WRN("Matrix leading dimensions do not agree." 
				" Enforce slow multiplication!",__FILE__,__LINE__)
			for(int i =0; i<_dimI;i++)
				for (int j = 0; j < arg.getDim_j(); j++)
					for (int k=0; k <L; k++)
						res(i,j) += (*this)(i,k)*arg(k,j);
		return res; 
	}

	complex float alpha_f = 1.0;
	complex double alpha_d = 1.0;
	complex float beta_f = 0.;
	complex double beta_d = 0.;

	switch (gettype<_Tp>()) {
		case FLT:
			cblas_sgemm((_ldim==0?CblasRowMajor:CblasColMajor),
					CblasNoTrans, CblasNoTrans, this->getDim_i(),
					arg.getDim_j(), this->getDim_j(), 1.f,
					(float*) *this->getMtxData(), this->getDim_i(),
					(float*) *arg.getMtxData(), arg.getDim_i(), 0.f,
					(float*) *res.getMtxData(), this->getDim_i());
			break;
		case DBL:
			cblas_dgemm((_ldim==0?CblasRowMajor:CblasColMajor),
					CblasNoTrans, CblasNoTrans, this->getDim_i(),
					arg.getDim_j(), this->getDim_j(), 1.,
					(double*) *this->getMtxData(), this->getDim_i(),
					(double*) *arg.getMtxData(), arg.getDim_i(), 0.,
					(double*) *res.getMtxData(), this->getDim_i());
			break;
		case CPLXFLT:
			cblas_cgemm((_ldim==0?CblasRowMajor:CblasColMajor),
					CblasNoTrans, CblasNoTrans, this->getDim_i(),
					arg.getDim_j(), this->getDim_j(), (void*) &alpha_f,
					(void*) *this->getMtxData(), this->getDim_i(),
					(void*) *arg.getMtxData(), arg.getDim_i(), (void*) &beta_f,
					(void*) *res.getMtxData(), this->getDim_i());
			break;

		case CPLXDBL:
			cblas_zgemm((_ldim==0?CblasRowMajor:CblasColMajor),
					CblasNoTrans, CblasNoTrans, this->getDim_i(),
					arg.getDim_j(), this->getDim_j(), (void*) &alpha_d,
					(void*) *this->getMtxData(), this->getDim_i(),
					(void*) *arg.getMtxData(), arg.getDim_i(), (void*) &beta_d,
					(void*) *res.getMtxData(), this->getDim_i());
			break;
		default:
			throw std::domain_error(
					"Unsupported data type for matrix-matrix mutiplication!");

	}

	return res;
}

//overload the addition operator
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator+(const Matrix<_Tp>& arg) const {

	if (this->getDim_i() != arg.getDim_i()
			|| this->getDim_j() != arg.getDim_j()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	Matrix<_Tp> res(_dimI,_dimJ,_ldim);

	if (_ldim != arg.getLDim()){
		DATASTRUCT_OUT_WRN("Matrix leading dimensions do not agree. Enforce slow addition",__FILE__,__LINE__)
			for(int i =0; i<_dimI;i++)
				for (int j = 0; j < _dimJ; j++)
					res(i,j) = (*this)(i,j)+arg(i,j);
		return res; 
	}






	complex double alpha_d = 1.;
	complex float alpha_f = 1.;

	switch (gettype<_Tp>()) {
		case FLT:
			cblas_saxpy(this->getDim_j()*this->getDim_i(), 1.0, (float*) this->getMtxData()[0],
					1, (float*) res.getMtxData()[0], 1);

			break;
		case DBL:

			cblas_daxpy(this->getDim_j()*this->getDim_i(), 1.0, (double*) this->getMtxData()[0],
					1, (double*) res.getMtxData()[0], 1);
			break;
		case CPLXFLT:

			cblas_caxpy(this->getDim_j()*this->getDim_i(), (void*) &alpha_f,
					(void*) this->getMtxData()[0], 1,
					(void*) res.getMtxData()[0], 1);
			break;
		case CPLXDBL:
			cblas_zaxpy(this->getDim_j()*this->getDim_i(), (void*) &alpha_d,
					(void*) this->getMtxData()[0], 1,
					(void*) res.getMtxData()[0], 1);
			break;
		default:
			throw std::domain_error(
					"Unsupported data type for matrix-matrix addition!");

	}
	return res;
}

//overload the subtraction operator
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator-(const Matrix<_Tp>& arg) {

	if (this->getDim_i() != arg.getDim_i()
			|| this->getDim_j() != arg.getDim_j()) {
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	Matrix<_Tp> res = ((_Tp) -1.) * arg;
	return (*this) + res;

}

template class dat::Matrix<float>;
template class dat::Matrix< complex float>;
template class dat::Matrix<double>;
template class dat::Matrix<complex double>;
