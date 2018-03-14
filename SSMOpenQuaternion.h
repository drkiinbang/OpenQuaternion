/*
 * Copyright (c) 2007-2017, SSMatics
 * All rights reserved.
 * SPATIAL&SPATIAL-MATICS for EarthOnAir
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
 * 
 * Any questions are welcome. earthliveonair@gmail.com
 */

#pragma once

#include <math.h>
#include <iostream>

#include "SMRotationMatrix.h"

namespace SSMATICS_EOA
{
	#define W_THRESHOLD 1.0e-99

	enum CSMQuaternionError
	{
		_DIVIDED_ZERO_, //Divided by zero
		_NOT_POINT_ //the quaternion does not mean a normal 3D point (w is not zero).
	};

	//
	// Info about quaternion
	// Q1Q2 is not equal Q2Q1. Not commutative.
	// (Q1Q2)Q3 = Q1(Q2Q3)
	// (Q1+Q2)Q3 = Q1Q3 + Q2Q3
	// a*(Q1+Q2) = a*Q1 + a*Q2
	// (a*Q1)Q2 = Q1(a*Q2)
	//

	class CSMQuaternion
	{
	public:

		//
		//Constructors
		//
		CSMQuaternion() : w_threshold(W_THRESHOLD)
		{
			w = x = y = z = 0.0;
		}

		CSMQuaternion(const double w, const double x, const double y, const double z) : w_threshold(W_THRESHOLD)
		{
			this->w = w;
			this->x = x;
			this->y = y;
			this->z = z;
		}

		CSMQuaternion(const double x, const double y, const double z) : w_threshold(W_THRESHOLD)
		{
			this->w = 0.0;
			this->x = x;
			this->y = y;
			this->z = z;
		}
			
		CSMQuaternion(const CSMMatrix<double>& R) : w_threshold(W_THRESHOLD)
		{
			EulerToQuaternion(R);
		}

		CSMQuaternion(const CSMQuaternion& V, const double th) : w_threshold(W_THRESHOLD)
		{
			w = cos(th*0.5);
			x = V.x*sin(th*0.5);
			y = V.y*sin(th*0.5);
			z = V.z*sin(th*0.5);
		}

		CSMQuaternion(const CSMQuaternion& copy) : w_threshold(W_THRESHOLD)
		{
			w = copy.w;
			x = copy.x;
			y = copy.y;
			z = copy.z;
		}

		//
		//Operators
		//

		void operator = (const CSMQuaternion& copy)
		{
			w = copy.w;
			x = copy.x;
			y = copy.y;
			z = copy.z;
		}

		CSMQuaternion operator % (const CSMQuaternion& q) const
		{
			//v(x,y,z)
			//p*q = ( p.w*q.w - p.v*q.v, p.w*q.v + p.v*q.w + p.vXq.v )  
			CSMQuaternion ret;
			
			ret.w = w*q.w - ( x*q.x + y*q.y + z*q.z );
			ret.x = w*q.x + x*q.w + ( y*q.z - z*q.y );
			ret.y = w*q.y + y*q.w + ( z*q.x - x*q.z );
			ret.z = w*q.z + z*q.w + ( x*q.y - y*q.x );
			
			return ret;
		}

		CSMQuaternion& operator %= (const CSMQuaternion& q)
		{
			*this = this->operator%(q);
			return *this;
		}

		CSMQuaternion operator * (const double& val) const
		{
			CSMQuaternion ret;
			
			ret.w = this->w*val;
			ret.x = this->x*val;
			ret.y = this->y*val;
			ret.z = this->z*val;
			
			return ret;
		}

		CSMQuaternion& operator *= (const double& val)
		{
			if (val == 0)
				throw CSMQuaternionError(_DIVIDED_ZERO_);

			this->w *= val;
			this->x *= val;
			this->y *= val;
			this->z *= val;

			return *this;
		}

		double operator * (const CSMQuaternion& val) const
		{
			return (this->w*val.w + this->x*val.x + this->y*val.y + this->z*val.z);
		}

		CSMQuaternion operator / (const double& val) const
		{
			if(val == 0)
				throw CSMQuaternionError(_DIVIDED_ZERO_);

			CSMQuaternion ret;
			
			ret.w = this->w/val;
			ret.x = this->x/val;
			ret.y = this->y/val;
			ret.z = this->z/val;
			
			return ret;
		}

		CSMQuaternion& operator /= (const double& val)
		{
			if (val == 0)
				throw CSMQuaternionError(_DIVIDED_ZERO_);

			this->w /= val;
			this->x /= val;
			this->y /= val;
			this->z /= val;

			return *this;
		}

		CSMQuaternion operator + (const CSMQuaternion& q) const
		{
			//v(x,y,z)
			//p+q = ( p.w+q.w, p.x+q.x, p.y+q.y, p.z+q.z)  
			CSMQuaternion ret;
			
			ret.w = this->w+q.w;
			ret.x = this->x+q.x;
			ret.y = this->y+q.y;
			ret.z = this->z+q.z;
			
			return ret;
		}

		CSMQuaternion& operator += (const CSMQuaternion& q)
		{
			//*this = this->operator + (q);
			this->w += q.w;
			this->x += q.x;
			this->y += q.y;
			this->z += q.z;

			return *this;
		}

		CSMQuaternion operator - (const CSMQuaternion& q) const
		{
			CSMQuaternion ret;
			
			ret.w = this->w-q.w;
			ret.x = this->x-q.x;
			ret.y = this->y-q.y;
			ret.z = this->z-q.z;
			
			return ret;
		}

		CSMQuaternion& operator -= (const CSMQuaternion& q)
		{
			this->w -= q.w;
			this->x -= q.x;
			this->y -= q.y;
			this->z -= q.z;

			return *this;
		}

		CSMQuaternion operator - () const
		{
			CSMQuaternion ret = *this * (-1.0);
			return ret;
		}

		//
		//public functions
		//

		CSMQuaternion Conjugate() const
		{
			CSMQuaternion ret(*this);

			ret.x *= -1.0;
			ret.y *= -1.0;
			ret.z *= -1.0;

			return ret;
		}

		CSMQuaternion Rotate(const CSMQuaternion P) const
		{
			if(P.w != 0)
				throw CSMQuaternionError(_NOT_POINT_);

			return Rotate(P.x, P.y, P.z);
		}

		CSMQuaternion Rotate(const double ax, const double ay, const double az) const
		{
			double norm = Norm();

			if(norm == 0)
				throw CSMQuaternionError(_DIVIDED_ZERO_);
			
			double w_, x_, y_, z_;

			w_ = w/norm;
			x_ = x/norm;
			y_ = y/norm;
			z_ = z/norm;

			CSMQuaternion ret;

			//
			//R % P % (R.Conjugate())
			//(w,x,y,z) % (0,ax,ay,az) % (w,-x,-y,-z)
			//
			ret.w = 0;
			ret.x = ax*(w_*w_ + x_*x_ - y_*y_ - z_*z_) - 2*ay*w_*z_ + 2*ay*x_*y_ + 2*az*w_*y_ + 2*az*x_*z_;
			ret.y = ay*(w_*w_ - x_*x_ + y_*y_ - z_*z_) + 2*ax*w_*z_ + 2*ax*x_*y_ - 2*az*w_*x_ + 2*az*y_*z_;
			ret.z = az*(w_*w_ - x_*x_ - y_*y_ + z_*z_) - 2*ax*w_*y_ + 2*ax*x_*z_ + 2*ay*w_*x_ + 2*ay*y_*z_;

			return ret;
		}

		CSMQuaternion Rotate_dw(const double& ax, const double& ay, const double& az) const
		{
			double norm = Norm();

			if(norm == 0)
				throw CSMQuaternionError(_DIVIDED_ZERO_);
			
			double w_, x_, y_, z_;

			w_ = w/norm;
			x_ = x/norm;
			y_ = y/norm;
			z_ = z/norm;

			CSMQuaternion ret;

			//
			//(w,x,y,z) % (0,ax,ay,az) % (w,-x,-y,-z)
			//
			ret.w = 0;
			ret.x = ax*(2*w_) - 2*ay*z_ + 2*az*y_;
			ret.y = ay*(2*w_) + 2*ax*z_ - 2*az*x_;
			ret.z = az*(2*w_) - 2*ax*y_ + 2*ay*x_;

			return ret;
		}

		CSMQuaternion Rotate_dx(const double& ax, const double& ay, const double& az) const
		{
			double norm = Norm();

			if(norm == 0)
				throw CSMQuaternionError(_DIVIDED_ZERO_);
			
			double w_, x_, y_, z_;

			w_ = w/norm;
			x_ = x/norm;
			y_ = y/norm;
			z_ = z/norm;

			CSMQuaternion ret;

			//
			//(w,x,y,z) % (0,ax,ay,az) % (w,-x,-y,-z)
			//
			ret.w = 0;
			ret.x = ax*(2*x_)  - 2*ay*y_ + 2*az*z_;
			ret.y = ay*(-2*x_) + 2*ax*y_ - 2*az*w_;
			ret.z = az*(-2*x_) + 2*ax*z_ + 2*ay*w_;

			return ret;
		}

		CSMQuaternion Rotate_dy(const double& ax, const double& ay, const double& az) const
		{
			double norm = Norm();

			if(norm == 0)
				throw CSMQuaternionError(_DIVIDED_ZERO_);
			
			double w_, x_, y_, z_;

			w_ = w/norm;
			x_ = x/norm;
			y_ = y/norm;
			z_ = z/norm;

			CSMQuaternion ret;

			//
			//(w,x,y,z) % (0,ax,ay,az) % (w,-x,-y,-z)
			//
			ret.w = 0;
			ret.x = ax*(-2*y_) + 2*ay*x_ + 2*az*w_;
			ret.y = ay*(2*y_ ) + 2*ax*x_ + 2*az*z_;
			ret.z = az*(-2*y_) - 2*ax*w_ + 2*ay*z_;

			return ret;
		}

		CSMQuaternion Rotate_dz(const double& ax, const double& ay, const double& az) const
		{
			double norm = Norm();

			if(norm == 0)
				throw CSMQuaternionError(_DIVIDED_ZERO_);
			
			double w_, x_, y_, z_;

			w_ = w/norm;
			x_ = x/norm;
			y_ = y/norm;
			z_ = z/norm;

			CSMQuaternion ret;

			//
			//(w,x,y,z) % (0,ax,ay,az) % (w,-x,-y,-z)
			//
			ret.w = 0;
			ret.x = ax*(-2*z_) - 2*ay*w_ + 2*az*x_;
			ret.y = ay*(-2*z_) + 2*ax*w_ + 2*az*y_;
			ret.z = az*(+2*z_) + 2*ax*x_ + 2*ay*y_;

			return ret;
		}

		double Norm2() const
		{
			return ( w*w + x*x + y*y + z*z);
		}

		double Norm() const
		{
			return sqrt(Norm2());
		}

		void Normalize()
		{
			double norm = Norm();
			if (0.0 == norm)
			{
				throw CSMQuaternionError(_DIVIDED_ZERO_);
			}
			else
			{
				if(1.0 != norm) *this /= norm;
			}
		}		

		CSMMatrix<double> QuaternionToEuler() const
		{
			CSMMatrix<double> R(3,3);
			double norm = this->Norm();
			return UnitQuaternionToEuler(this->w / norm, this->x / norm, this->y / norm, this->z / norm);
		}

		static CSMMatrix<double> QuaternionToEuler(const double w, const double x, const double y, const double z)
		{
			CSMMatrix<double> R(3, 3);
			
			try
			{				
				double norm = sqrt(w*w + x*x + y*y + z*z);
				return UnitQuaternionToEuler(w / norm, x / norm, y / norm, z / norm);
			}
			catch (...)
			{
				std::string msg("Error: EulerToQuaternion(const CSMMatrix<double>& R) in SMQuaternion.h");
				std::cerr << msg;
				throw std::runtime_error(msg);
			}

			return R;
		}

		static CSMMatrix<double> UnitQuaternionToEuler(const double q0, const double q1, const double q2, const double q3)
		{
			CSMMatrix<double> R(3, 3);

			R(0, 0) = 1.0 - 2.0 * q2*q2 - 2 * q3*q3;	R(0, 1) = 2.0 * q1*q2 - 2.0 * q0*q3;		R(0, 2) = 2.0 * q1*q3 + 2.0 * q0*q2;
			R(1, 0) = 2.0 * q1*q2 + 2.0 * q0*q3;		R(1, 1) = 1.0 - 2.0 * q1*q1 - 2.0 * q3*q3;	R(1, 2) = 2.0 * q2*q3 - 2.0 * q0*q1;
			R(2, 0) = 2.0 * q1*q3 - 2.0 * q0*q2;		R(2, 1) = 2.0 * q2*q3 + 2.0 * q0*q1;		R(2, 2) = 1.0 - 2.0 * q1*q1 - 2.0 * q2*q2;

			return R;
		}

		bool EulerToQuaternion(const CSMMatrix<double>& R)
		{
			try
			{
				w = 0.5*sqrt(R(0, 0) + R(1, 1) + R(2, 2) + 1.0);

				if (fabs(w) > fabs(w_threshold))
				{
					x = (0.25 / w)*(R(2, 1) - R(1, 2));
					y = (0.25 / w)*(R(0, 2) - R(2, 0));
					z = (0.25 / w)*(R(1, 0) - R(0, 1));
				}
				else
				{
					x = 0.5*sqrt(R(0, 0) - R(1, 1) - R(2, 2) + 1.0);
					y = 0.5*sqrt(-R(0, 0) + R(1, 1) - R(2, 2) + 1.0);
					z = 0.5*sqrt(-R(0, 0) - R(1, 1) + R(2, 2) + 1.0);

					if (fabs(x) >= fabs(y) && fabs(x) >= fabs(z) && fabs(x) > fabs(w_threshold))
					{
						y = (0.25 / x)*(R(0, 1) + R(1, 0));
						z = (0.25 / x)*(R(2, 0) + R(0, 2));
						w = (0.25 / x)*(R(2, 1) - R(1, 2));
					}
					else if (fabs(y) >= fabs(x) && fabs(y) >= fabs(z) && fabs(y) > fabs(w_threshold))
					{
						x = (0.25 / y)*(R(0, 1) + R(1, 0));
						z = (0.25 / y)*(R(1, 2) + R(2, 1));
						w = (0.25 / y)*(R(0, 2) - R(2, 0));
					}
					else if (fabs(z)>fabs(x) && fabs(z)>fabs(y) && fabs(z) > fabs(w_threshold))
					{
						x = (0.25 / z)*(R(0, 2) + R(2, 0));
						y = (0.25 / z)*(R(1, 2) + R(2, 1));
						w = (0.25 / z)*(R(1, 0) - R(0, 1));
					}
					else
						return false;
				}

				return true;
			}
			catch(...)
			{
				std::string msg("Error: EulerToQuaternion(const CSMMatrix<double>& R) in SMQuaternion.h");
				std::cerr << msg;
				throw std::runtime_error(msg);
			}
		}

		//spherical interpolation
		static CSMQuaternion Interpolation(const CSMQuaternion& qa, const CSMQuaternion& qb, const double ta, const double tb, const double t, const double th = 0.95)
		{
			double dt = tb - ta;
			if (dt == 0. || t == ta)
				return qa;
			else if (t == tb)
				return qb;

			return Interpolation(qa, qb, (t-ta)/(dt), th);
		}

		//Quaternion interpolation including slerp and lerp
		static CSMQuaternion Interpolation(const CSMQuaternion& qa, const CSMQuaternion& qb0, const double u, const double th = 0.95)
		{
			// 0 < |u| < 1.0

			CSMQuaternion qb2;
			double dot = qa * qb0;
			
			CSMQuaternion retval;
			//	
			//dot is cos(theta),
			//if (dot < 0), q1 and q2 are more than 90 degrees apart,
			//so we can invert one to reduce spinning
			//
			CSMQuaternion qb;
			if(fabs(dot) < th)// spherical interpolation (slerp)
			{
				if (dot < 0.0)
				{
					dot = dot * -1.0;
					qb = qb0 * -1.0;
				}

				double angle = acos(dot);
				retval = (qa*sin(angle*(1.0 - u)) + qb*sin(angle*u)) / sin(angle);
			}
			else//linear interpolartion (lerp)
			{
				retval = (qa*(1.0 - u) + qb*u);
			}

			retval.Normalize();
			return retval;
		}

		//slerp interpolation
		static CSMQuaternion InterpolationSlerp(const CSMQuaternion& qa, const CSMQuaternion& qb0, const double dot0, const double u)
		{
			// 0 < |u| < 1.0

			//acos: 0.0 to 1.0 (therefore, -pi/2 to +pi/2)
			//if dot(qa*qb) is less than zero, then make it positive.
			double dot = dot0;
			CSMQuaternion qb = qb0;
			if (dot0 < 0.0)
			{
				dot = dot * -1.0;
				qb = qb * -1.0;;
			}
			double angle = acos(dot);
			CSMQuaternion retval = (qa*sin(angle*(1.0 - u)) + qb*sin(angle*u)) / sin(angle);
			retval.Normalize();
			return retval;
		}

		//lerp interpolation
		static CSMQuaternion InterpolationLerp(const CSMQuaternion& qa, const CSMQuaternion& qb, const double u)
		{
			CSMQuaternion retval = (qa*(1.0 - u) + qb*u);
			retval.Normalize();
			return retval;
		}

		// it does not check > 90deg (refer to quaternion.h)
		static CSMQuaternion InterpolationNoInvert(const CSMQuaternion& qa, const CSMQuaternion& qb, const double u, const double th = 0.95)
		{
			double dot = qa * qb;
			CSMQuaternion retval;
			if (fabs(dot) < th)
			{
				double angle = acos(dot);
				retval =  (qa*sin(angle*(1.0 - u)) + qb*sin(angle*u)) / sin(angle);
			}
			else  // if the angle is small, use linear interpolation								
				retval = (qa*(1.0 - u) + qb*u);

			retval.Normalize();
			return retval;
		}

		//
		//static functions
		//

		//Error propagation for quaternion derived from Euler angles
		static CSMQuaternion EulerToQuaternion_Variance(const double omega, const double phi, const double kappa, const double VarO, const double VarP, const double VarK)
		{
			CSMMatrix<double> RmatVar = GetRmatVarfromEulerAngleVar(omega,phi,kappa,VarO,VarP,VarK);

			double w, x, y, z;

			_Rmat_ Rotation(omega, phi, kappa);
			CSMMatrix<double> R = Rotation.rmatrix();
			CSMQuaternion Var_q;

			w = 0.5*sqrt(R(0,0) + R(1,1) + R(2,2) + 1.0);

			if(sqrt(R(0,0) + R(1,1) + R(2,2) + 1.0) > 0.0)
			{
				CSMMatrix<double> G(1,3,0);
				CSMMatrix<double> V(3,3,0);

				V(0,0) = RmatVar(0,0);
				V(1,1) = RmatVar(1,1);
				V(2,2) = RmatVar(2,2);

				G(0,0) = 0.5*0.5/sqrt(R(0,0) + R(1,1) + R(2,2) + 1.0);
				G(0,1) = 0.5*0.5/sqrt(R(0,0) + R(1,1) + R(2,2) + 1.0);
				G(0,2) = 0.5*0.5/sqrt(R(0,0) + R(1,1) + R(2,2) + 1.0);

				Var_q.w = (G%V%G.Transpose())(0,0);
			}
			
			if( fabs(w) > fabs(W_THRESHOLD) ) 
			{
				CSMMatrix<double> G(1,3,0);
				CSMMatrix<double> V(3,3,0);

				//x = (0.25/w)*(R(2,1) - R(1,2));
				V(0,0) = Var_q.w;
				V(1,1) = RmatVar(2,1);
				V(2,2) = RmatVar(1,2);

				G(0,0) = -0.25/w/w*(R(2,1) - R(1,2));
				G(0,1) = 0.25/w;
				G(0,2) = -0.25/w;

				Var_q.x = (G%V%G.Transpose())(0,0);

				//y = (0.25/w)*(R(0,2) - R(2,0));
				V(0,0) = Var_q.w;
				V(1,1) = RmatVar(0,2);
				V(2,2) = RmatVar(2,0);

				G(0,0) = -0.25/w/w*(R(0,2) - R(2,0));
				G(0,1) = 0.25/w;
				G(0,2) = -0.25/w;

				Var_q.y = (G%V%G.Transpose())(0,0);

				//z = (0.25/w)*(R(1,0) - R(0,1));
				V(0,0) = Var_q.w;
				V(1,1) = RmatVar(1,0);
				V(2,2) = RmatVar(0,1);

				G(0,0) = -0.25/w/w*(R(1,0) - R(0,1));
				G(0,1) = 0.25/w;
				G(0,2) = -0.25/w;

				Var_q.z = (G%V%G.Transpose())(0,0);
			}
			else
			{
				x = 0.5*sqrt(R(0,0) - R(1,1) - R(2,2) + 1.0);
				y = 0.5*sqrt(-R(0,0) + R(1,1) - R(2,2) + 1.0);
				z = 0.5*sqrt(-R(0,0) - R(1,1) + R(2,2) + 1.0);

				if(fabs(x)>=fabs(y) && fabs(x)>=fabs(z) && fabs(x) > fabs(W_THRESHOLD))
				{
					if(sqrt(R(0,0) - R(1,1) - R(2,2) + 1.0) > 0.0)
					{
						CSMMatrix<double> G(1,3,0);
						CSMMatrix<double> V(3,3,0);

						V(0,0) = RmatVar(0,0);
						V(1,1) = RmatVar(1,1);
						V(2,2) = RmatVar(2,2);

						G(0,0) = 0.5*0.5/sqrt(R(0,0) - R(1,1) - R(2,2) + 1.0);
						G(0,1) = -0.5*0.5/sqrt(R(0,0) - R(1,1) - R(2,2) + 1.0);
						G(0,2) = -0.5*0.5/sqrt(R(0,0) - R(1,1) - R(2,2) + 1.0);

						Var_q.x = (G%V%G.Transpose())(0,0);
					}

					CSMMatrix<double> G(1,3,0);
					CSMMatrix<double> V(3,3,0);

					//y = (0.25/x)*(R(0,1) + R(1,0));
					V(0,0) = Var_q.x;
					V(1,1) = RmatVar(0,1);
					V(2,2) = RmatVar(1,0);

					G(0,0) = -0.25/x/x*(R(0,1) + R(1,0));
					G(0,1) = 0.25/x;
					G(0,2) = 0.25/x;

					Var_q.y = (G%V%G.Transpose())(0,0);

					//z = (0.25/x)*(R(2,0) + R(0,2));
					V(0,0) = Var_q.x;
					V(1,1) = RmatVar(2,0);
					V(2,2) = RmatVar(0,2);

					G(0,0) = -0.25/x/x*(R(2,0) + R(0,2));
					G(0,1) = 0.25/x;
					G(0,2) = 0.25/x;

					Var_q.z = (G%V%G.Transpose())(0,0);

					//w = (0.25/x)*(R(2,1) - R(1,2));
					V(0,0) = Var_q.x;
					V(1,1) = RmatVar(2,1);
					V(2,2) = RmatVar(1,2);

					G(0,0) = -0.25/x/x*(R(2,1) - R(1,2));
					G(0,1) = 0.25/x;
					G(0,2) = -0.25/x;

					Var_q.w = (G%V%G.Transpose())(0,0);
				}
				else if(fabs(y)>=fabs(x) && fabs(y)>=fabs(z) && fabs(y) > fabs(W_THRESHOLD))
				{
					if(sqrt(-R(0,0) + R(1,1) - R(2,2) + 1.0) > 0.0)
					{
						CSMMatrix<double> G(1,3,0);
						CSMMatrix<double> V(3,3,0);

						V(0,0) = RmatVar(0,0);
						V(1,1) = RmatVar(1,1);
						V(2,2) = RmatVar(2,2);

						G(0,0) = -0.5*0.5/sqrt(-R(0,0) + R(1,1) - R(2,2) + 1.0);
						G(0,1) = 0.5*0.5/sqrt(-R(0,0) + R(1,1) - R(2,2) + 1.0);
						G(0,2) = -0.5*0.5/sqrt(-R(0,0) + R(1,1) - R(2,2) + 1.0);

						Var_q.y = (G%V%G.Transpose())(0,0);
					}

					CSMMatrix<double> G(1,3,0);
					CSMMatrix<double> V(3,3,0);

					//x = (0.25/y)*(R(0,1) + R(1,0));
					V(0,0) = Var_q.y;
					V(1,1) = RmatVar(0,1);
					V(2,2) = RmatVar(1,0);

					G(0,0) = -0.25/y/y*(R(0,1) + R(1,0));
					G(0,1) = 0.25/y;
					G(0,2) = 0.25/y;

					Var_q.x = (G%V%G.Transpose())(0,0);

					//z = (0.25/y)*(R(1,2) + R(2,1));
					V(0,0) = Var_q.y;
					V(1,1) = RmatVar(1,2);
					V(2,2) = RmatVar(2,1);

					G(0,0) = -0.25/y/y*(R(1,2) + R(2,1));
					G(0,1) = 0.25/y;
					G(0,2) = 0.25/y;

					Var_q.z = (G%V%G.Transpose())(0,0);

					//w = (0.25/y)*(R(0,2) - R(2,0));
					V(0,0) = Var_q.y;
					V(1,1) = RmatVar(0,2);
					V(2,2) = RmatVar(2,0);

					G(0,0) = -0.25/y/y*(R(0,2) - R(2,0));
					G(0,1) = 0.25/y;
					G(0,2) = -0.25/y;

					Var_q.w = (G%V%G.Transpose())(0,0);
				}
				else if(fabs(z)>fabs(x) && fabs(z)>fabs(y) && fabs(z) > fabs(W_THRESHOLD))
				{
					if(sqrt(-R(0,0) - R(1,1) + R(2,2) + 1.0) > 0.0)
					{
						CSMMatrix<double> G(1,3,0);
						CSMMatrix<double> V(3,3,0);

						V(0,0) = RmatVar(0,0);
						V(1,1) = RmatVar(1,1);
						V(2,2) = RmatVar(2,2);

						G(0,0) = -0.5*0.5/sqrt(-R(0,0) - R(1,1) + R(2,2) + 1.0);
						G(0,1) = -0.5*0.5/sqrt(-R(0,0) - R(1,1) + R(2,2) + 1.0);
						G(0,2) = 0.5*0.5/sqrt(-R(0,0) - R(1,1) + R(2,2) + 1.0);

						Var_q.z = (G%V%G.Transpose())(0,0);
					}

					CSMMatrix<double> G(1,3,0);
					CSMMatrix<double> V(3,3,0);

					//x = (0.25/z)*(R(0,2) + R(2,0));
					V(0,0) = Var_q.z;
					V(1,1) = RmatVar(0,2);
					V(2,2) = RmatVar(2,0);

					G(0,0) = -0.25/z/z*(R(0,2) + R(2,0));
					G(0,1) = 0.25/z;
					G(0,2) = 0.25/z;

					Var_q.x = (G%V%G.Transpose())(0,0);

					//y = (0.25/z)*(R(1,2) + R(2,1));
					V(0,0) = Var_q.z;
					V(1,1) = RmatVar(1,2);
					V(2,2) = RmatVar(2,1);

					G(0,0) = -0.25/z/z*(R(1,2) + R(2,1));
					G(0,1) = 0.25/z;
					G(0,2) = 0.25/z;

					Var_q.y = (G%V%G.Transpose())(0,0);

					//w = (0.25/z)*(R(1,0) - R(0,1));
					V(0,0) = Var_q.z;
					V(1,1) = RmatVar(1,0);
					V(2,2) = RmatVar(0,1);

					G(0,0) = -0.25/z/z*(R(1,0) - R(0,1));
					G(0,1) = 0.25/z;
					G(0,2) = -0.25/z;

					Var_q.w = (G%V%G.Transpose())(0,0);
				}
			}

			return Var_q;
		}

		//Error propagation for direction consines derived from Euler angles
		static CSMMatrix<double> GetRmatVarfromEulerAngleVar(double omega, double phi, double kappa, const double VarO, const double VarP, const double VarK)
		{
			if(omega < 1.0e20)
				omega = 1.0e20;

			if(phi < 1.0e20)
				phi = 1.0e20;

			if(kappa < 1.0e20)
				kappa = 1.0e20;

			CSMMatrix<double> Varmat(3,3,0);
			CSMMatrix<double> G(1,3);
			CSMMatrix<double> V(3,3,0);
			V(0,0) = VarO;
			V(1,1) = VarP;
			V(2,2) = VarK;

			//Rmatrix(0,0) = cos(phi)*cos(kappa);
			G(0,0) = 0;
			G(0,1) = -sin(phi)*cos(kappa);
			G(0,2) = -cos(phi)*sin(kappa);
			Varmat(0,0) = (G%V%G.Transpose())(0,0);

			//Rmatrix(1,0) = sin(omega)*sin(phi)*cos(kappa) + cos(omega)*sin(kappa);
			G(0,0) = cos(omega)*sin(phi)*cos(kappa) - sin(omega)*sin(kappa);
			G(0,1) = sin(omega)*cos(phi)*cos(kappa);
			G(0,2) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
			Varmat(1,0) = (G%V%G.Transpose())(0,0);

			//Rmatrix(2,0) = -cos(omega)*sin(phi)*cos(kappa) + sin(omega)*sin(kappa);
			G(0,0) = sin(omega)*sin(phi)*cos(kappa) - cos(omega)*sin(kappa);
			G(0,1) = -cos(omega)*cos(phi)*cos(kappa);
			G(0,2) = -cos(omega)*cos(phi)*cos(kappa);
			Varmat(2,0) = (G%V%G.Transpose())(0,0);

			//Rmatrix(0,1) = -cos(phi)*sin(kappa);
			G(0,0) = 0;
			G(0,1) = sin(phi)*sin(kappa);
			G(0,2) = -cos(phi)*cos(kappa);
			Varmat(0,1) = (G%V%G.Transpose())(0,0);

			//Rmatrix(1,1) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
			G(0,0) = -cos(omega)*sin(phi)*sin(kappa) - sin(omega)*cos(kappa);
			G(0,1) = -sin(omega)*cos(phi)*sin(kappa);
			G(0,2) = -sin(omega)*sin(phi)*cos(kappa) - cos(omega)*sin(kappa);
			Varmat(1,1) = (G%V%G.Transpose())(0,0);

			//Rmatrix(2,1) = cos(omega)*sin(phi)*sin(kappa) + sin(omega)*cos(kappa);
			G(0,0) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
			G(0,1) = cos(omega)*cos(phi)*sin(kappa);
			G(0,2) = cos(omega)*sin(phi)*cos(kappa) - sin(omega)*sin(kappa);
			Varmat(2,1) = (G%V%G.Transpose())(0,0);

			//Rmatrix(0,2) = sin(phi);
			G(0,0) = 0;
			G(0,1) = cos(phi);
			G(0,2) = 0;
			Varmat(0,2) = (G%V%G.Transpose())(0,0);

			//Rmatrix(1,2) = -sin(omega)*cos(phi);
			G(0,0) = -cos(omega)*cos(phi);
			G(0,1) = sin(omega)*sin(phi);
			G(0,2) = 0;
			Varmat(1,2) = (G%V%G.Transpose())(0,0);

			//Rmatrix(2,2) = cos(omega)*cos(phi);
			G(0,0) = -sin(omega)*cos(phi);
			G(0,1) = -cos(omega)*sin(phi);
			G(0,2) = 0;
			Varmat(2,2) = (G%V%G.Transpose())(0,0);

			return Varmat;
		}

		//3D Transformation
		static CSMQuaternion Transformation3D(CSMQuaternion P, CSMQuaternion R, CSMQuaternion T, double S)
		{
			return ( (R % P % (R.Conjugate()))*S + T );
		}
		
		// Get a quaternion rotating about an arbitrary axis (x, y, z)
		static CSMQuaternion RotationAboutArbitraryAixs(const double x, const double y, const double z, const double ang_rad)
		{
			//Rotation axis
			CSMQuaternion R(0.0, x, y, z);
			R.Normalize();
			//Rotation angle
			double half_ang = ang_rad * 0.5;
			R.w = cos(half_ang);
			R.x *= sin(half_ang);
			R.y *= sin(half_ang);
			R.z *= sin(half_ang);

			return R;
		}

		// Get a rotated quaternion derived using rotation axis and angle
		static CSMQuaternion RotateAboutArbitraryAixs(const double x, const double y, const double z, const double ang_rad, const CSMQuaternion& P)
		{
			//Rotation axis
			CSMQuaternion R = RotationAboutArbitraryAixs(x, y, z, ang_rad);
			//return a rotated quaternion
			return R.Rotate(P);
		}

		// Reflected quaternion (vector)
		static CSMQuaternion Reflection(const CSMQuaternion& P, const CSMQuaternion& normalVec)
		{
			CSMQuaternion P0 = P; P0.Normalize();
			CSMQuaternion normalVec0 = normalVec; normalVec0.Normalize();
			return (normalVec0 % P0 % normalVec0);
		}

	public:
		double w, x, y, z;

		const double w_threshold;
	};
}