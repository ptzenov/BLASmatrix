#ifndef _MATRIX_
#define _MATRIX_

#include <utils.hpp>

/***************************************************************
 * Usual 2D matrix of size dimI x dimJ
 ***************************************************************/
namespace dat {

	template<typename _Tp>
		class Matrix {

			/*
			 * A generic matrix class implement fast (via BLAS) matrix-matrix and matrix-vector 
			 * operations with convenient matlab-style syntax.  
			 *
			 * 
			 * Implements the observer pattern for quick memory management
			 *
			 */

			private:
				// list of subscribers 
				std::map<unsigned int, Matrix<_Tp> > subscribers;

			protected:
				unsigned int _dimI, _dimJ, _bytes;
				unsigned int _ldim; // store the leading dimension.
				void init(unsigned int dim_i, unsigned int dim_j,unsigned int ldim);

				_Tp** _data;
				// cant initialize an empty matrix 
				Matrix(){
					;
				}
				/*
				 *
				 *
				 *
				 *  as every matrix M can be the parent of more than 1 views,
				 *  for memory management issues the parent matrix needs to keep track 
				 *  of all "active views" on itself. This is needed, for instance, if the 
				 *  parent matrix is freed. In such cases all correspondingly registered
				 *  views have to be freed as well since otherwise  they will refer
				 *  to invalid memory location.
				 */

				void registerview(unsigned int viewid, Matrix<_Tp>& view);

				// remove a view with a unique view_id from the list 
				// of subscribers. this method shall be called upon view destruction.
				void removeview(unsigned int view_id);

			public:

				unsigned int getDim_i() const;
				unsigned int getDim_j() const;
				unsigned int getLDim() const; 

				unsigned int getBytesize() const {
					return _bytes;
				}

				_Tp** const getMtxData() const;
				bool is_square() const {
					return _dimI == _dimJ;
				}


				// initialize an empty complex Nr matrix with leading dimension ldim 
				// ldim = 0 -> C++ style memory allocation (row-major) 
				// ldim = 1 -> Fortran style memry allocaiton (column-major)
				Matrix(unsigned int dim_i, unsigned int dim_j, unsigned int ldim);

				// copy constructor
				Matrix(const Matrix<_Tp>& arg);
				// destructor
				virtual ~Matrix();

				/**
				 * Overload the assignment operator !
				 * Extremely important in order to be able to correctly execute 
				 * arithmetic operations on matrices with matlab-like syntax!!!
				 */
				void operator=(const Matrix<_Tp>& arg);

				/**
				 * Return a reference (implicit pointer) to the i,j-th data element
				 * so that the user can actually overwrite the corresponding entry
				 */
				_Tp& operator()(unsigned int i, unsigned int j) const;

				/**
				 * Overload the function call operator ->  indices
				 *
				 *  @param i1 - the row index of the upper left element  of the	submatrix
				 *  @param j1 - the col index of the upper left element  of the	submatrix
				 *  @param i2 - the row index of the lower right element of the submatrix;
				 *  @param j2 - the col index of the lower right element of the submatrix;
				 *  @return - a submatrix VIEW-object of the current matrix 
				 *  with indices (i1,j1)-->(i2,j2)!
				 *
				 */
				Matrix<_Tp> operator()(unsigned int i1, unsigned int i2, unsigned int j1,
						unsigned int j2);

				// overload the multiplication by matrix operator...
				Matrix<_Tp> operator*(const Matrix<_Tp>& arg) const;

				//overload the addition operator
				Matrix<_Tp> operator+(const Matrix<_Tp>& arg) const;

				//overload the subtraction operator
				Matrix<_Tp> operator-(const Matrix<_Tp>& arg);

				//overload the print operator
				friend std::ostream& operator<<(std::ostream& os,
						const dat::Matrix<_Tp>& obj) {

					int dimROW = obj.getDim_i();
					int dimCOL = obj.getDim_j();

					for (int i = 0; i < dimROW; i++) {
						for (int j = 0; j < dimCOL; j++)
							os << obj(i,j)<< " ";
						os << "\n";
					}

					return os;
				}

				/**
				 * overload the matrix * number operator
				 */
				friend Matrix<_Tp> operator*(_Tp lhs, const Matrix<_Tp>& rhs) {
					dat::Matrix<_Tp> res = rhs;
					int M = rhs.getDim_i();
					int N = rhs.getDim_j();
					int type = gettype<_Tp>();

					switch (type) {
						case FLT:
							cblas_sscal(N * M, explicit_cast(float,lhs), 
									(float*) *res._data,1);
							break;
						case DBL:
							cblas_dscal(N * M, explicit_cast(double,lhs), 
									(double*) *res._data,1);
							break;
						case CPLXFLT:
							cblas_cscal(N * M, (void*) &lhs, 
									(void*) (*res._data), 1);
							break;
						case CPLXDBL:
							cblas_zscal(N * M, (void*) &lhs,
									(void*) (*res._data), 1);
							break;
						default:
							throw std::domain_error("Unsupported matrix scale operation");
					}

					return res;
				}

				friend Matrix<_Tp> operator*(Matrix<_Tp> lhs, _Tp rhs) {
					return rhs * lhs;
				}

		}
	;

	template<typename _Tp>
		dat::Matrix<_Tp> eye(unsigned int M) {

			dat::Matrix<_Tp> res(M, M,0);
			for (int i = 0; i < M; i++)
				res(i, i) = 1.0;

			return res;
		}
}
#endif
