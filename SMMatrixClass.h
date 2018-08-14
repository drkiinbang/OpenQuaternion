/*
* Copyright (c) 2001-2001, KI IN Bang
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS `AS IS'
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

namespace SMATICS_MATRIX
{
	/**@enum defined errors
	*/
	enum ErrorName
	{
		R_OR_C_IS_LESS_THAN_ZERO,
		SIZE_NOT_SAME,
		COL_ROW_NOT_SAME,
		OUT_OF_RANGE,
		NOT_SQUARE_MATRIX,
		SCALAR_VALUE,
		SINGULAR_MATRIX,
		EMPTY_DATA,
		IMPROPER_DIMENSION,
		NOT_DEFINED
	};

	/**@class MatrixError
	*@brief CSMMatrix Error Handling
	*/
	class CSMMatrixError
	{
	public:
		CSMMatrixError(ErrorName e)
		{
			error = e;
		}
		const char* WhichError()
		{
			if (error == SIZE_NOT_SAME)
			{
				const char* reason = "Two Matrixs are NOT same size!"; return reason;
			}
			else if (error == R_OR_C_IS_LESS_THAN_ZERO)
			{
				const char* reason = "Row or Column is not bigger than zero"; return reason;
			}
			else if (error == OUT_OF_RANGE)
			{
				const char* reason = "Request element out of range"; return reason;
			}
			else if (error == COL_ROW_NOT_SAME)
			{
				const char* reason = "Matrix Mutiple is impossible. Because Column of A matrix is not equal to Row of B matrix"; return reason;
			}
			else if (error == NOT_SQUARE_MATRIX)
			{
				const char* reason = "Not square matrix. So Inverse impossible"; return reason;
			}
			else if (error == SCALAR_VALUE)
			{
				const char* reason = "Matrix rows one and Matrix cols is one"; return reason;
			}
			else if (error == SINGULAR_MATRIX)
			{
				const char* reason = "Singular matrix"; return reason;
			}
			else if (error == NOT_DEFINED)
			{
				const char* reason = "No name error"; return reason;
			}
			else
			{
				const char* reason = "I don't know reason!"; return reason;
			}
		}

	private:
		ErrorName error;

	};

	/**
	*@class CSMMatrix Class
	*@brief CSMMatrix handling class.
	*/
	template<typename T>
	class CSMMatrix
	{
	public:

		/**constructor*/
		CSMMatrix()
		{
			R = 1;
			C = 1;
			Data = new T[R*C];
			Data[0] = (T)0;
		}
		
		/**copy constructor*/
		CSMMatrix(const CSMMatrix<T> &copy)
		{
			R = copy.R;
			C = copy.C;
			Data = new T[R*C];

			for (unsigned int i = 0;i<(R*C);i++)
			{
				Data[i] = copy.Data[i];
			}
		}
		
		/**constructor*/
		CSMMatrix(const unsigned int rows, const unsigned int cols = 1, const T init_val = T(0))
		{
			if ((rows <= 0) || (cols <= 0))
			{
				Data = nullptr;
			}
			else
			{
				R = rows;
				C = cols;

				Data = new T[R*C];

				for (unsigned int i = 0;i<(R*C);i++)
				{
					Data[i] = T(init_val);
				}
			}

		}

		/**desturctor*/
		virtual ~CSMMatrix()
		{
			if (Data != nullptr)
			{
				delete[] Data;
				Data = nullptr;
			}
		}

		/**Add rows to existing matrix*/
		void AddRow(const CSMMatrix<T> &copy)
		{
			if (C != copy.C)
				throw CSMMatrixError(SIZE_NOT_SAME);

			CSMMatrix<T> temp;
			temp = *this;

			unsigned int rows, cols; //Original size
			rows = R;
			cols = C;

			Resize(R + copy.R, C, 0);

			unsigned int i, j;

			for (i = 0; i<rows; i++)
			{
				for (j = 0; j<cols; j++)
				{
					operator()(i, j) = temp(i, j);
				}
			}

			for (i = 0; i<copy.R; i++)
			{
				for (j = 0; j<copy.C; j++)
				{
					operator()(rows + i, j) = copy(i, j);
				}
			}
		}

		/**Add columns to an existing matrix*/
		void AddCol(const CSMMatrix<T> &copy)
		{
			if (R != copy.R)
				throw CSMMatrixError(SIZE_NOT_SAME);

			CSMMatrix<T> temp;
			temp = *this;

			unsigned int rows, cols; //Original size
			rows = R;
			cols = C;

			Resize(R, C + copy.C, 0); //Resizing

			unsigned int i, j;

			for (i = 0; i<rows; i++)
			{
				for (j = 0; j<cols; j++)
				{
					operator()(i, j) = temp(i, j);
				}
			}

			for (i = 0; i<copy.R; i++)
			{
				for (j = 0; j<copy.C; j++)
				{
					operator()(i, cols + j) = copy(i, j);
				}
			}
		}

		/**bind two matrices<br>
		*|A|add|B|=|A||0|<br>
		*          |0||B| = C.
		*/
		void AddRowCol(const CSMMatrix<T> &copy)
		{
			CSMMatrix<T> temp;
			temp = *this;

			unsigned int rows, cols; //Original size
			rows = R;
			cols = C;

			Resize(R + copy.R, C + copy.C, 0);

			unsigned int i, j;

			for (i = 0; i<rows; i++)
			{
				for (j = 0; j<cols; j++)
				{
					operator()(i, j) = temp(i, j);
				}
			}

			for (i = 0; i<copy.R; i++)
			{
				for (j = 0; j<copy.C; j++)
				{
					operator()(rows + i, cols + j) = copy(i, j);
				}
			}

		}

		/**change specific part with given matrix<br>
		*if the size of an existing matrix is not available, resize the matrix.
		*/
		void Insert(const unsigned int row, const unsigned int col, const CSMMatrix<T> &copy)
		{
			if ((R<(row + copy.R)) || (C<(col + copy.C)))
			{
				//throw CSMMatrixError(OUT_OF_RANGE);
				CSMMatrix<T> temp = *this;
				unsigned int newrow = R, newcol = C;
				if (R<(row + copy.R)) newrow = row + copy.R;
				if (C<(col + copy.C)) newcol = col + copy.C;

				Resize(newrow, newcol, 0.);
				Insert(0, 0, temp);

				unsigned int i, j;

				for (i = 0; i<copy.R; i++)
				{
					for (j = 0; j<copy.C; j++)
					{
						operator()(row + i, col + j) = copy(i, j);
					}
				}
			}
			else
			{
				unsigned int i, j;

				for (i = 0; i<copy.R; i++)
				{
					for (j = 0; j<copy.C; j++)
					{
						operator()(row + i, col + j) = copy(i, j);
					}
				}
			}
		}

		/**using given matrix and position accumulate existing matrix<br>
		*if the existing matrix size is not proper, throw error.
		*/
		void InsertPlus(const unsigned int row, const unsigned int col, const CSMMatrix<T> &copy)
		{
			if ((R<(row + copy.R)) || (C<(col + copy.C)))
			{
				throw CSMMatrixError(OUT_OF_RANGE);
			}
			else
			{
				unsigned int i, j;

				for (i = 0; i<copy.R; i++)
				{
					for (j = 0; j<copy.C; j++)
					{
						operator()(row + i, col + j) += copy(i, j);
					}
				}
			}
		}

		/**get a subset from an existing matrix*/
		CSMMatrix<T> Subset(const unsigned int row, const unsigned int col, const unsigned int size_R, const unsigned int size_C) const
		{
			CSMMatrix<T> ret;
			if ((R<(row + size_R)) || (C<(col + size_C)))
			{
				throw CSMMatrixError(OUT_OF_RANGE);
			}
			else
			{
				ret.Resize(size_R, size_C);
				unsigned int i, j;

				for (i = 0; i<size_R; i++)
				{
					for (j = 0; j<size_C; j++)
					{
						ret(i, j) = operator()(row + i, col + j);
					}
				}
			}
			return ret;
		}

		/**Get a largest absolute value.*/
		double MaxabsElement() const
		{
			double val = fabs(Data[0]);
			for (unsigned int i = 0; i<R; i++)
			{
				for (unsigned int j = 0; j<C; j++)
				{
					if (val<(fabs(Data[i*C + j])))
						val = fabs(Data[i*C + j]);
				}
			}
			return val;
		}

		/**Get a largest value.*/
		double MaxElement() const
		{
			double val = Data[0];
			for (unsigned int i = 0; i<R; i++)
			{
				for (unsigned int j = 0; j<C; j++)
				{
					if (val<Data[i*C + j])
						val = Data[i*C + j];
				}
			}
			return val;
		}

		/**Get a smallest absolute value.*/
		double MinabsElement() const
		{
			double val = fabs(Data[0]);
			for (unsigned int i = 0; i<R; i++)
			{
				for (unsigned int j = 0; j<C; j++)
				{
					if (val>(fabs(Data[i*C + j])))
						val = fabs(Data[i*C + j]);
				}
			}
			return val;
		}

		/**Get a smallest value.*/
		double MinElement() const
		{
			double val = Data[0];
			for (unsigned int i = 0; i<R; i++)
			{
				for (unsigned int j = 0; j<C; j++)
				{
					if (val>Data[i*C + j])
						val = Data[i*C + j];
				}
			}
			return val;
		}

		/**Assignment operator*/
		void operator = (const CSMMatrix<T> &copy)
		{
			this->Resize(copy.R, copy.C, 0);

			for (unsigned int i = 0;i<(R*C);i++)
			{
				Data[i] = copy.Data[i];
			}
		}

		/**+ operator*/
		CSMMatrix<T> operator + (const CSMMatrix<T> & mat) const//Matrix + Matrix
		{
			if ((R != mat.R) || (C != mat.C))
				throw CSMMatrixError(SIZE_NOT_SAME);

			CSMMatrix<T> result(R, C, T(0));
			for (unsigned int i = 0; i<(R*C); i++)
				result.Data[i] = Data[i] + mat.Data[i];

			return result;
		}
		
		/**+= operator*/
		void operator += (const CSMMatrix<T> & mat) //Matrix += Matrix
		{
			if ((R != mat.R) || (C != mat.C))
				throw CSMMatrixError(SIZE_NOT_SAME);

			for (unsigned int i = 0; i<(R*C); i++)
				Data[i] = Data[i] + mat.Data[i];
		}
		
		/**+ operator*/
		CSMMatrix<T> operator + (const T & val) const//Matrix + val
		{
			CSMMatrix<T> result(R, C, T(0));
			for (unsigned int i = 0; i<(R*C); i++)
				result.Data[i] = Data[i] + val;

			return result;
		}
		
		/**+= operator*/
		void operator += (const T & val) //Matrix += val
		{
			for (unsigned int i = 0; i<(R*C); i++)
				Data[i] = Data[i] + val;
		}

		/**- operator*/
		CSMMatrix<T> operator - (const CSMMatrix<T> & mat) const//Matrix - Matrix
		{
			if ((R != mat.R) || (C != mat.C))
				throw CSMMatrixError(SIZE_NOT_SAME);

			CSMMatrix<T> result(R, C, T(0));
			for (unsigned int i = 0; i<(R*C); i++)
				result.Data[i] = Data[i] - mat.Data[i];

			return result;
		}

		/**-= operator*/
		void operator -= (const CSMMatrix<T> & mat) //Matrix -= Matrix
		{
			if ((R != mat.R) || (C != mat.C))
				throw CSMMatrixError(SIZE_NOT_SAME);

			for (unsigned int i = 0; i<(R*C); i++)
				Data[i] = Data[i] - mat.Data[i];
		}

		/**- operator*/
		CSMMatrix<T> operator - (const T & val) const//Matrix - val
		{
			CSMMatrix<T> result(R, C, T(0));
			for (unsigned int i = 0; i<(R*C); i++)
				result.Data[i] = Data[i] - val;

			return result;
		}

		/**-= operator*/
		void operator -= (const T & val) //Matrix -= val
		{
			for (unsigned int i = 0; i<(R*C); i++)
				Data[i] = Data[i] - val;
		}

		/**scalar * operator*/
		CSMMatrix<T> operator * (const CSMMatrix<T> & mat) const//Matrix * Matrix
		{
			if ((R != mat.R) || (C != mat.C))
				throw CSMMatrixError(SIZE_NOT_SAME);

			CSMMatrix<T> result(R, C, T(0));
			for (unsigned int i = 0; i<(R*C); i++)
				result.Data[i] = Data[i] * mat.Data[i];

			return result;
		}

		/**scalar *= operator*/
		void operator *= (const CSMMatrix<T> & mat) //Matrix *= Matrix
		{
			if ((R != mat.R) || (C != mat.C))
				throw CSMMatrixError(SIZE_NOT_SAME);

			for (unsigned int i = 0; i<(R*C); i++)
				Data[i] = Data[i] * mat.Data[i];
		}

		/**scalar * operator*/
		CSMMatrix<T> operator * (const T & val) const//Matrix * val
		{
			CSMMatrix<T> result(R, C, T(0));
			for (unsigned int i = 0; i<(R*C); i++)
				result.Data[i] = Data[i] * val;

			return result;
		}

		/**scalar *= operator*/
		void operator *= (const T & val) //Matrix *= val
		{
			for (unsigned int i = 0; i<(R*C); i++)
				Data[i] = Data[i] * val;
		}

		/**scalar / operator*/
		CSMMatrix<T> operator / (const CSMMatrix<T> & mat) const//Matrix / Matrix
		{
			if ((R != mat.R) || (C != mat.C))
				throw CSMMatrixError(SIZE_NOT_SAME);

			CSMMatrix<T> result(R, C, T(0));
			for (unsigned int i = 0; i<(R*C); i++)
				result.Data[i] = Data[i] / mat.Data[i];

			return result;
		}

		/**scalar /= operator*/
		void operator /= (const CSMMatrix<T> & mat) //Matrix /= Matrix
		{
			if ((R != mat.R) || (C != mat.C))
				throw CSMMatrixError(SIZE_NOT_SAME);

			for (unsigned int i = 0; i<(R*C); i++)
				Data[i] = Data[i] / mat.Data[i];
		}

		/**scalar / operator*/
		CSMMatrix<T> operator / (const T & val) const//Matrix / val
		{
			CSMMatrix<T> result(R, C, T(0));
			for (unsigned int i = 0; i<(R*C); i++)
				result.Data[i] = Data[i] / val;

			return result;
		}

		/**scalar /= operator*/
		void operator /= (const T & val) //Matrix / val
		{
			for (unsigned int i = 0; i<(R*C); i++)
				Data[i] = Data[i] / val;
		}

		/**matrix multiplication*/
		CSMMatrix<T> operator % (const CSMMatrix<T> & mat) const//Matrix % Matrix
		{
			if (C != mat.R)
				throw CSMMatrixError(COL_ROW_NOT_SAME);

			T sum;
			CSMMatrix<T> result(R, mat.C, T(0));

			CSMMatrix<T> mat_copy = mat;
			for (unsigned int i = 0;i<R; i++)
			{
				for (unsigned int j = 0; j<mat.C; j++)
				{
					sum = T(0);
					for (unsigned int k = 0; k<C; k++)
					{
						sum += operator()(i, k)*mat_copy(k, j);
					}
					result(i, j) = sum;
				}
			}
			return result;
		}

		/**minus matrix*/
		CSMMatrix<T> operator - () const
		{
			CSMMatrix<T> result(R, C, T(0));
			for (unsigned int i = 0; i<(R*C); i++)
			{
				result.Data[i] = -Data[i];
			}

			return result;
		}

		/**() operator: access a element of matrix*/
		T& operator () (const unsigned int row, const unsigned int col = 0) const
		{
			unsigned int i;
			i = FindIndex(row, col);
			return Data[i];
		}

		CSMMatrix<T> cross(const CSMMatrix<T> & mat) const
		{
			if (R != mat.R || C != 1 || mat.C != 1)
				throw CSMMatrixError(IMPROPER_DIMENSION);

			CSMMatrix<T> result(R, C);
			for (unsigned int i = 0;i<R; i++)
			{
				unsigned int idx1, idx2;
				idx1 = fmod(i + 1, R);
				idx2 = fmod(i + 2, R);

				result(i, 0) = Data[idx1] * mat.Data[idx2] - mat.Data[idx1] * Data[idx2];
			}

			return result;
		}

		T dot(const CSMMatrix<T> & mat) const
		{
			if (R != mat.R || C != 1 || mat.C != 1)
				throw CSMMatrixError(IMPROPER_DIMENSION);

			T sum = 0;

			for (unsigned int i = 0; i<(R*C); i++)
			{
				sum += Data[i] * mat.Data[i];
			}

			return sum;
		}

	public:
		T* getDataHandle() const { return Data; }

		/**make identity matrix*/
		void Identity(const T & init_val = T(1))
		{
			if (C != R)
				throw CSMMatrixError(NOT_SQUARE_MATRIX);
			if (C<2)
				throw CSMMatrixError(SCALAR_VALUE);

			for (unsigned int i = 0; i<R; i++)
			{
				for (unsigned int j = 0; j<C; j++)
				{
					if (i == j) operator()(i, j) = init_val;
					else operator()(i, j) = T(0);
				}
			}
		};

		/**Resize a matrix*/
		void Resize(const unsigned int rows, const unsigned int cols = 1, const T & init_val = T(0))
		{
			if ((C*R) != 0)
				delete[] Data;

			R = rows;
			C = cols;

			if ((rows <= 0) || (cols <= 0))
			{
				Data = nullptr;
				//throw CSMMatrixError(R_OR_C_IS_LESS_THAN_ZERO);
			}
			else
			{
				Data = new T[R*C];
				for (unsigned int i = 0; i<(R*C); i++)
					Data[i] = init_val;
			}
		}
		
		/**Get number of rows*/
		unsigned int GetRows() const { return R; }

		/**Get number of cols*/
		unsigned int GetCols() const { return C; }

		/**make a transpose matrix*/
		CSMMatrix<T> Transpose() const
		{
			CSMMatrix<T> result(C, R, T(0));
			for (unsigned int i = 0; i<R; i++)
			{
				for (unsigned int j = 0; j<C; j++)
				{
					result(j, i) = operator()(i, j);
				}
			}

			return result;
		};

		/**make a inverse matrix*/
		CSMMatrix<T> Inverse() const
		{
			if (C != R)
				throw CSMMatrixError(NOT_SQUARE_MATRIX);
			if (C<1)
				throw CSMMatrixError(SCALAR_VALUE);

			CSMMatrix<T> copy(R, C, 0);

			if (C == 1)
			{
				copy.Data[0] = 1.0 / Data[0];
				return copy;
			}
			else
			{
				for (unsigned int i = 0; i<(R*C); i++)
					copy.Data[i] = Data[i];

				CSMMatrix<unsigned int> LUD = copy.LUDcompose();

				CSMMatrix<T> mat_inv = copy.LUInvert(LUD);

				return mat_inv;
			}

		}

		/**Get Max value*/
		T GetMax() const
		{
			return MaxabsElement();
		};
		
		/**Get square sum*/
		double GetSumOfSquares() const
		{
			double sum = 0.0;
			unsigned int maxNum = R*C;
			for (unsigned int i = 0; i < maxNum; ++i)
				sum += Data[i] * Data[i];

			return sum;
		}

		/**Get data handle */
		const T* pData() const { return Data; }

	private:

		/**LU Decompose*/
		CSMMatrix<unsigned int> LUDcompose();

		/**LU Invert*/
		CSMMatrix<T> LUInvert(CSMMatrix<unsigned int> &perm);

		/**LU solve*/
		CSMMatrix<T> LUSolve(CSMMatrix<unsigned int> & perm, CSMMatrix<T> & b);

		/**get index of Data for (row, col)*/
		unsigned int FindIndex(const unsigned int row, const unsigned int col) const
		{
			if ((R <= row) || (C <= col))
				throw CSMMatrixError(OUT_OF_RANGE);

			unsigned int Index = (row*(C)) + col;
			return Index;
		}

	private:
		/**<number of rows*/
		unsigned int R;
		/**<number of cols*/
		unsigned int C;

		T* Data;


	}; //Matrix Class

	template<typename T>
	CSMMatrix<unsigned int> CSMMatrix<T>::LUDcompose()
	{
		// LU decomposition
		unsigned int i, j, k, k2, t;
		T p, temp;
		// permutation matrix
		CSMMatrix<unsigned int> perm(R, 1);

		// initialize permutation
		for (i = 0; i < R; ++i)
			perm(i, 0) = i;

		for (k = 0; k < (R - 1); ++k)
		{
			p = T(0);

			for (i = k; i < R; ++i)
			{
				temp = (Data[FindIndex(i, k)] < 0 ? (-Data[FindIndex(i, k)]) : Data[FindIndex(i, k)]);

				if (temp > p)
				{
					p = temp;
					k2 = i;
				}
			}

			if (p == T(0))
				throw CSMMatrixError(SINGULAR_MATRIX);

			// exchange rows
			t = perm(k, 0);
			perm(k, 0) = perm(k2, 0);
			perm(k2, 0) = t;

			for (i = 0; i < R; ++i)
			{
				temp = Data[FindIndex(k, i)];
				Data[FindIndex(k, i)] = Data[FindIndex(k2, i)];
				Data[FindIndex(k2, i)] = temp;
			}

			for (i = k + 1; i < R; ++i)
			{
				Data[FindIndex(i, k)] /= Data[FindIndex(k, k)];

				for (j = k + 1; j < R; ++j)
					Data[FindIndex(i, j)] -= Data[FindIndex(i, k)] * Data[FindIndex(k, j)];
			}
		}
		return perm;
	};

	template <typename T>
	CSMMatrix<T> CSMMatrix<T>::LUInvert(CSMMatrix<unsigned int> &perm)
	{
		unsigned int i, j;
		CSMMatrix<T> p(R, 1);
		CSMMatrix<T> result(R, R);

		for (j = 0; j < R; ++j)
		{
			for (i = 0; i < R; ++i)
				p(i, 0) = T(0);

			p(j, 0) = T(1);

			p = LUSolve(perm, p);

			for (i = 0; i < R; ++i)
				result(i, j) = p(i, 0);
		}

		return result;
	};

	template <typename T>
	CSMMatrix<T> CSMMatrix<T>::LUSolve(CSMMatrix<unsigned int> & perm, CSMMatrix<T> & b)
	{
		unsigned int i, j, j2;
		T sum, u;
		CSMMatrix<T> y(R, 1), x(R, 1);

		for (i = 0; i < R; ++i)
		{
			sum = T(0);
			j2 = 0;

			for (j = 1; j <= i; ++j)
			{
				sum += Data[FindIndex(i, j2)] * y(j2, 0);
				++j2;
			}
			y(i, 0) = b(perm(i, 0), 0) - sum;
		}

		i = R - 1;

		while (1)
		{
			sum = T(0);
			u = Data[FindIndex(i, i)];

			for (j = i + 1; j < R; ++j)
				sum += Data[FindIndex(i, j)] * x(j, 0);

			x(i, 0) = (y(i, 0) - sum) / u;

			if (i == 0)
				break;

			--i;
		}

		return x;
	};
}