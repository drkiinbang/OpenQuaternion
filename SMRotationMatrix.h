/*
 * Copyright (c) 1999-2001, KI IN Bang
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
 */

//////////////////////////
//SMRotationMatrix.h
//made by BbaraB
///////////////////////////////////////////////////////

#ifndef ___SMRotationMatrix___BBARAB___
#define ___SMRotationMatrix___BBARAB___

#include "SMMatrixClass.h"

#include <math.h>
#include <sstream>
#include <string>

#define GET_PHI asin(1.0)*2

#define CSMMatrix SMATICS_MATRIX::CSMMatrix

namespace SSMATICS_EOA
{

static CSMMatrix<double> getM1M2M3(double omega, double phi, double kappa)
{
	double cosO = cos(omega),	sinO = sin(omega);
	double cosP = cos(phi),		sinP = sin(phi);
	double cosK = cos(kappa),	sinK = sin(kappa);

	CSMMatrix<double> M1M2M3(3, 3);
	M1M2M3(0,0) = cosP*cosK;						M1M2M3(0,1) = cosP*sinK;						M1M2M3(0,2) = -sinP;
	M1M2M3(1,0) = sinO*sinP*cosK - cosO*sinK;		M1M2M3(1,1) = sinO*sinP*sinK + cosO*cosK;		M1M2M3(1,2) = sinO*cosP;
	M1M2M3(2,0) = cosO*sinP*cosK + sinO*sinK;		M1M2M3(2,1) = cosO*sinP*sinK - sinO*cosK;		M1M2M3(2,2) = cosO*cosP;

	return M1M2M3;
}

//It is not completed yet.
static void extractEulerAnglesM1M2M3(const CSMMatrix<double>& M1M2M3, double& O, double& P, double& K)
{
	const CSMMatrix<double>& M = M1M2M3;
	double phi = asin(-M(0,2));

	double omega_1, phi_1, kappa_1;
	double omega_2, phi_2, kappa_2;

	// two possible values for phi
	phi_1 = phi;
	if(phi < 0) phi_2 = -M_PI - phi;
	else phi_2 = M_PI - phi;

	double sin_omega_1 = M1M2M3(1,2)/cos(phi_1);
	double cos_omega_1 = M1M2M3(2,2)/cos(phi_1);

	double sin_omega_2 = M1M2M3(1,2)/cos(phi_2);
	double cos_omega_2 = M1M2M3(2,2)/cos(phi_2);

	omega_1 = atan2(sin_omega_1, cos_omega_1);
	omega_2 = atan2(sin_omega_2, cos_omega_2);
}

#define _Mmat_ CRotationcoeff
#define _Rmat_ CRotationcoeff2
#define _Rmat_R1R2R3_ CRmat_R1R2R3
#define _Rmat_R3R2R1_ CRmat_R3R2R1
#define _Rmat_Y3P1R2_ CRmat_Y3P1R2

class CRotationcoeff //Ground->Photo
{
public:
	CSMMatrix<double> Mmatrix;//Rotation Matrix(Ground->Photo)
	
public:
	CRotationcoeff()
	{
	}
	
	CRotationcoeff(double omega,double phi,double kappa)
	{
		ReMake(omega, phi, kappa);
	}

	CRotationcoeff(const double(&rot)[3][3])
	{
		Mmatrix.Resize(3, 3);
		for (unsigned int r=0; r < 3; ++r)
		{
			for (unsigned int c = 0; c < 3; ++c)
			{
				Mmatrix(r, c) = rot[r][c];
			}
		}
	}

	void ExtractRotation(double &omega, double &phi, double &kappa)
	{
		ExtractRotation(Mmatrix, omega, phi, kappa);
	}

	static void ExtractRotation(CSMMatrix<double> Mmatrix, double &omega, double &phi, double &kappa)
	{
		double omega_1, phi_1, kappa_1;
		double omega_2, phi_2, kappa_2;
		
		phi = asin(Mmatrix(2,0));
		
		// two possible values for phi
		phi_1 = phi;
		if(phi < 0) phi_2 = -M_PI - phi;
		else phi_2 = GET_PHI - phi;
		
		omega_1 = atan2(-Mmatrix(2,1)/cos(phi_1), Mmatrix(2,2)/cos(phi_1));
		kappa_1 = atan2(-Mmatrix(1,0)/cos(phi_1), Mmatrix(0,0)/cos(phi_1));

		omega_2 = atan2(-Mmatrix(2,1)/cos(phi_2), Mmatrix(2,2)/cos(phi_2));
		kappa_2 = atan2(-Mmatrix(1,0)/cos(phi_2), Mmatrix(0,0)/cos(phi_2));
		
		CRotationcoeff temp_1(omega_1,  phi_1, kappa_1);
		CRotationcoeff temp_2(omega_2, phi_2, kappa_2);

		double sum_1=0, sum_2=0;
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				sum_1 += fabs(Mmatrix(i, j) - temp_1.Mmatrix(i, j));
				sum_2 += fabs(Mmatrix(i, j) - temp_2.Mmatrix(i, j));
			}
		}
		
		if(sum_1 < sum_2)
		{
			phi = phi_1;
			omega = omega_1;
			kappa = kappa_1;
		}
		else
		{
			phi = phi_2;
			omega = omega_2;
			kappa = kappa_2;
		}
	}
	
	virtual ~CRotationcoeff(){}

	void ReMake(double omega,double phi,double kappa)
	{
		// Mmatrix = M(0,0,kappa)%M(0,phi,0)%M(omega,0,0)
		Mmatrix.Resize(3,3);
		
		Mmatrix(0,0) = cos(phi)*cos(kappa);
		Mmatrix(0,1) = sin(omega)*sin(phi)*cos(kappa) + cos(omega)*sin(kappa);
		Mmatrix(0,2) = -cos(omega)*sin(phi)*cos(kappa) + sin(omega)*sin(kappa);
		
		Mmatrix(1,0) = -cos(phi)*sin(kappa);
		Mmatrix(1,1) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
		Mmatrix(1,2) = cos(omega)*sin(phi)*sin(kappa) + sin(omega)*cos(kappa);
		
		Mmatrix(2,0) = sin(phi);
		Mmatrix(2,1) = -sin(omega)*cos(phi);
		Mmatrix(2,2) = cos(omega)*cos(phi);
	}

	bool GetOPK(double &O, double &P, double &K)
	{
		if( true == GetOPK(this->Mmatrix, O, P, K) )
			return true;
		else
			return false;
	}

	static bool GetOPK(CSMMatrix<double> Mmatrix, double &O, double &P, double &K)
	{
		if(Mmatrix(2,2) == 0.0) return false;
		O = atan2(-Mmatrix(2,1), Mmatrix(2,2));

		if(Mmatrix(0,0) == 0.0) return false;
		K = atan2(-Mmatrix(1,0), Mmatrix(0,0));

		//if(Mmatrix(2,2) == 0.0) return false;
		P = asin(Mmatrix(2,0));
		double cos_phi_1 = Mmatrix(0,0)/cos(K);
		double cos_phi_2 = -Mmatrix(1,0)/sin(K);
		double cos_phi_3 = -Mmatrix(2,1)/sin(O);
		double cos_phi_4 = Mmatrix(2,2)/cos(O);

		double cos_phi = (cos_phi_1+cos_phi_2+cos_phi_3+cos_phi_4)/4.0;
		double sin_phi = Mmatrix(2,0);
		P = atan2(sin_phi, cos_phi);

		return true;
	}
	
	static CSMMatrix<double> Partial_dMdO(const double omega, const double phi, const double kappa)
	{
		CSMMatrix<double> dMdO(3,3);

		dMdO(0,0) = 0;
		dMdO(1,0) = 0;
		dMdO(2,0) = 0;
		
		dMdO(0,1) = cos(omega)*sin(phi)*cos(kappa) - sin(omega)*sin(kappa);
		dMdO(1,1) = -cos(omega)*sin(phi)*sin(kappa) - sin(omega)*cos(kappa);
		dMdO(2,1) = -cos(omega)*cos(phi);
		
		dMdO(0,2) = sin(omega)*sin(phi)*cos(kappa) + cos(omega)*sin(kappa);
		dMdO(1,2) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
		dMdO(2,2) = -sin(omega)*cos(phi);

		return dMdO;
	}
	
	static CSMMatrix<double> Partial_dMdP(const double omega, const double phi, const double kappa)
	{
		CSMMatrix<double> dMdP(3, 3);

		dMdP(0,0) = -sin(phi)*cos(kappa);
		dMdP(1,0) = sin(phi)*sin(kappa);
		dMdP(2,0) = cos(phi);
		
		dMdP(0,1) = sin(omega)*cos(phi)*cos(kappa);
		dMdP(1,1) = -sin(omega)*cos(phi)*sin(kappa);
		dMdP(2,1) = sin(omega)*sin(phi);
		
		dMdP(0,2) = -cos(omega)*cos(phi)*cos(kappa);
		dMdP(1,2) = cos(omega)*cos(phi)*sin(kappa);
		dMdP(2,2) = -cos(omega)*sin(phi);

		return dMdP;
	}
	
	static CSMMatrix<double> Partial_dMdK(const double omega, const double phi, const double kappa)
	{
		CSMMatrix<double> dMdK(3, 3);

		dMdK(0,0) = -cos(phi)*sin(kappa);
		dMdK(1,0) = -cos(phi)*cos(kappa);
		dMdK(2,0) = 0.;
		
		dMdK(0,1) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
		dMdK(1,1) = -sin(omega)*sin(phi)*cos(kappa) - cos(omega)*sin(kappa);
		dMdK(2,1) = 0.;
		
		dMdK(0,2) = cos(omega)*sin(phi)*sin(kappa) + sin(omega)*cos(kappa);
		dMdK(1,2) = cos(omega)*sin(phi)*cos(kappa) - sin(omega)*sin(kappa);
		dMdK(2,2) = 0.;

		return dMdK;
	}

	static void GetPartial(const double omega, const double phi, const double kappa, CSMMatrix<double> &dRdO, CSMMatrix<double> &dRdP, CSMMatrix<double> &dRdK)
	{
		dRdO = Partial_dMdO(omega, phi, kappa);
		dRdP = Partial_dMdP(omega, phi, kappa);
		dRdK = Partial_dMdK(omega, phi, kappa);

	}

	CRotationcoeff operator % (const CRotationcoeff& copy)
	{
		CSMMatrix<double> temp = this->Mmatrix % copy.Mmatrix;
		CRotationcoeff retval;
		retval.Mmatrix = temp;
		return retval;
	}

	CSMMatrix<double> operator % (const CSMMatrix<double>& copy)
	{
		CSMMatrix<double> retval;

		try
		{
			retval = this->Mmatrix % copy;
		}
		catch (...)
		{
			std::string msg("Error: this->Mmatrix % copy in SMRotationMatrix.h");
			std::cerr << msg << std::endl;
			throw std::runtime_error(msg);
		}

		return retval;
	}
};

class CRotationcoeff2 //Photo->Ground
{
private:
	CSMMatrix<double> Rmatrix;//Inverse Rotation Matrix(Photo->Ground)
public:
	CSMMatrix<double>& rmatrix() { return Rmatrix; }
	
	virtual ~CRotationcoeff2(){}
	
	CRotationcoeff2()
	{
	}
	
	CRotationcoeff2(double omega,double phi,double kappa)
	{
		ReMake(omega, phi, kappa);	
	}

	CRotationcoeff2(const double(&rot)[3][3])
	{
		Rmatrix.Resize(3, 3);
		for (unsigned int r = 0; r < 3; ++r)
		{
			for (unsigned int c = 0; c < 3; ++c)
			{
				Rmatrix(r, c) = rot[r][c];
			}
		}
	}

	void ExtractRotation(double &omega, double &phi, double &kappa)
	{
		ExtractRotation(Rmatrix, omega, phi, kappa);
	}

	static void ExtractRotation(CSMMatrix<double> Rmatrix, double &omega, double &phi, double &kappa)
	{
		double omega_1, phi_1, kappa_1;
		double omega_2, phi_2, kappa_2;
		
		phi = asin(Rmatrix(0,2));
		
		if(phi < 0)
		{ 
			phi_1 = phi;
			phi_2 = -GET_PHI - phi;
		}
		else
		{
			phi_1 = phi;
			phi_2 = GET_PHI - phi;
		}
		
		omega_1 = atan2(-Rmatrix(1,2)/cos(phi_1), Rmatrix(2,2)/cos(phi_1));
		kappa_1 = atan2(-Rmatrix(0,1)/cos(phi_1), Rmatrix(0,0)/cos(phi_1));

		omega_2 = atan2(-Rmatrix(1,2)/cos(phi_2), Rmatrix(2,2)/cos(phi_2));
		kappa_2 = atan2(-Rmatrix(0,1)/cos(phi_2), Rmatrix(0,0)/cos(phi_2));
		
		CRotationcoeff2 temp_1(omega_1,  phi_1, kappa_1);
		CRotationcoeff2 temp_2(omega_2, phi_2, kappa_2);

		double sum_1=0, sum_2=0;
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				sum_1 += fabs(Rmatrix(i, j) - temp_1.Rmatrix(i, j));
				sum_2 += fabs(Rmatrix(i, j) - temp_2.Rmatrix(i, j));
			}
		}
		
		if(sum_1 < sum_2)
		{
			phi = phi_1;
			omega = omega_1;
			kappa = kappa_1;
		}
		else
		{
			phi = phi_2;
			omega = omega_2;
			kappa = kappa_2;
		}
	}
	
	void ReMake(double omega,double phi,double kappa)
	{
		// Rmatrix = R(omega,0,0)%R(0,phi,0)%R(0,0,kappa)
		
		Rmatrix.Resize(3,3);
		Rmatrix(0,0) = cos(phi)*cos(kappa);
		Rmatrix(0,1) = -cos(phi)*sin(kappa);
		Rmatrix(0,2) = sin(phi);
		
		Rmatrix(1,0) = sin(omega)*sin(phi)*cos(kappa) + cos(omega)*sin(kappa);
		Rmatrix(1,1) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
		Rmatrix(1,2) = -sin(omega)*cos(phi);
		
		Rmatrix(2,0) = -cos(omega)*sin(phi)*cos(kappa) + sin(omega)*sin(kappa);
		Rmatrix(2,1) = cos(omega)*sin(phi)*sin(kappa) + sin(omega)*cos(kappa);
		Rmatrix(2,2) = cos(omega)*cos(phi);
	}

	static CSMMatrix<double> GetRotationMat(double omega, double phi, double kappa)
	{
		CSMMatrix<double> RotationMat(3,3);
		RotationMat(0,0) = cos(phi)*cos(kappa);
		RotationMat(1,0) = sin(omega)*sin(phi)*cos(kappa) + cos(omega)*sin(kappa);
		RotationMat(2,0) = -cos(omega)*sin(phi)*cos(kappa) + sin(omega)*sin(kappa);
		
		RotationMat(0,1) = -cos(phi)*sin(kappa);
		RotationMat(1,1) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
		RotationMat(2,1) = cos(omega)*sin(phi)*sin(kappa) + sin(omega)*cos(kappa);
		
		RotationMat(0,2) = sin(phi);
		RotationMat(1,2) = -sin(omega)*cos(phi);
		RotationMat(2,2) = cos(omega)*cos(phi);

		return RotationMat;
	}

	// Rempove GetOPK() in _Rmat_, wrong compatation.Instead of this, use GetOPK in _Mmat_ or ExtractRotation() function.
	/*
	bool GetOPK(double &O, double &P, double &K)
	{
		...
	}

	static bool GetOPK(CSMMatrix<double> Rmatrix, double &O, double &P, double &K)
	{
		...
	}
	*/
	
	static CSMMatrix<double> Partial_dRdO(const double omega, const double phi, const double kappa)
	{
		CSMMatrix<double> dRdO(3, 3, 0.0);

		dRdO(0,0) = 0;
		dRdO(0,1) = 0;
		dRdO(0,2) = 0;
		
		dRdO(1,0) = cos(omega)*sin(phi)*cos(kappa) - sin(omega)*sin(kappa);
		dRdO(1,1) = -cos(omega)*sin(phi)*sin(kappa) - sin(omega)*cos(kappa);
		dRdO(1,2) = -cos(omega)*cos(phi);
		
		dRdO(2,0) = sin(omega)*sin(phi)*cos(kappa) + cos(omega)*sin(kappa);
		dRdO(2,1) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
		dRdO(2,2) = -sin(omega)*cos(phi);

		return dRdO;
	}
	
	static CSMMatrix<double> Partial_dRdP(const double omega, const double phi, const double kappa)
	{
		CSMMatrix<double> dRdP(3, 3, 0.0);

		dRdP(0,0) = -sin(phi)*cos(kappa);
		dRdP(0,1) = sin(phi)*sin(kappa);
		dRdP(0,2) = cos(phi);
		
		dRdP(1,0) = sin(omega)*cos(phi)*cos(kappa);
		dRdP(1,1) = -sin(omega)*cos(phi)*sin(kappa);
		dRdP(1,2) = sin(omega)*sin(phi);
		
		dRdP(2,0) = -cos(omega)*cos(phi)*cos(kappa);
		dRdP(2,1) = cos(omega)*cos(phi)*sin(kappa);
		dRdP(2,2) = -cos(omega)*sin(phi);

		return dRdP;
	}
	
	static CSMMatrix<double> Partial_dRdK(const double omega, const double phi, const double kappa)
	{
		CSMMatrix<double> dRdK(3, 3, 0.0);

		dRdK(0,0) = -cos(phi)*sin(kappa);
		dRdK(0,1) = -cos(phi)*cos(kappa);
		dRdK(0,2) = 0.;
		
		dRdK(1,0) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
		dRdK(1,1) = -sin(omega)*sin(phi)*cos(kappa) - cos(omega)*sin(kappa);
		dRdK(1,2) = 0.;
		
		dRdK(2,0) = cos(omega)*sin(phi)*sin(kappa) + sin(omega)*cos(kappa);
		dRdK(2,1) = cos(omega)*sin(phi)*cos(kappa) - sin(omega)*sin(kappa);
		dRdK(2,2) = 0.;

		return dRdK;
	}

	static void GetPartial(const double omega, const double phi, const double kappa, CSMMatrix<double> &dRdO, CSMMatrix<double> &dRdP, CSMMatrix<double> &dRdK)
	{
		dRdO = Partial_dRdO(omega, phi, kappa);
		dRdP = Partial_dRdP(omega, phi, kappa);
		dRdK = Partial_dRdK(omega, phi, kappa);
	}

	CRotationcoeff2 operator % (const CRotationcoeff2& copy)
	{
		CSMMatrix<double> temp = this->Rmatrix % copy.Rmatrix;
		CRotationcoeff2 retval;
		retval.Rmatrix = temp;
		return retval;
	}

	CSMMatrix<double> operator % (const CSMMatrix<double>& copy)
	{
		CSMMatrix<double> retval;

		try
		{
			retval = this->Rmatrix % copy;
		}
		catch (...)
		{
			std::string msg("Error: this->Rmatrix % copy in SMRotationMatrix.h");
			std::cerr << msg << std::endl;
			throw std::runtime_error(msg);
		}

		return retval;
	}
};

class CRmat_R1R2R3
{
public:
	CSMMatrix<double> Rmatrix;//Inverse Rotation Matrix(Photo->Ground)
	
	virtual ~CRmat_R1R2R3(){}
	
	CRmat_R1R2R3()
	{
	}
	
	CRmat_R1R2R3(double omega,double phi,double kappa)
	{
		ReMake(omega, phi, kappa);	
	}

	void ReMake(double omega,double phi,double kappa)
	{
		_Rmat_ R1(omega, 0.0, 0.0);
		_Rmat_ R2(0.0, phi, 0.0);
		_Rmat_ R3(0.0, 0.0, kappa);

		Rmatrix = R1.rmatrix() % R2.rmatrix() % R3.rmatrix();
	}
};

class CRmat_R3R2R1
{
public:
	CSMMatrix<double> Rmatrix;//Inverse Rotation Matrix(Photo->Ground)
	
	virtual ~CRmat_R3R2R1(){}
	
	CRmat_R3R2R1()
	{
	}
	
	CRmat_R3R2R1(double omega,double phi,double kappa)
	{
		ReMake(omega, phi, kappa);	
	}

	void ReMake(double omega,double phi,double kappa)
	{
		_Rmat_ R1(omega, 0.0, 0.0);
		_Rmat_ R2(0.0, phi, 0.0);
		_Rmat_ R3(0.0, 0.0, kappa);

		Rmatrix = R3.rmatrix() % R2.rmatrix() % R1.rmatrix();
	}
};

class CRmat_Y3P1R2 //1st rotation: roll (2nd axis), 2nd rotation: pitch (1st axis), 3rd rotation: yaw(3rd axis)
{
public:
	CSMMatrix<double> Rmatrix;//Inverse Rotation Matrix(Photo->Ground)
	
	virtual ~CRmat_Y3P1R2(){}
	
	CRmat_Y3P1R2()
	{
	}
	
	CRmat_Y3P1R2(double omega,double phi,double kappa)
	{
		ReMake(omega, phi, kappa);	
	}

	void ReMake(double omega,double phi,double kappa)
	{
		_Rmat_ Roll(0, omega, 0);
		_Rmat_ Pitch(phi, 0, 0);
		_Rmat_ Yaw(0, 0, kappa);
		
		Rmatrix = Yaw.rmatrix() % Pitch.rmatrix() % Roll.rmatrix();
	}
};

//Forward rotation matrix Differential (omega)
inline CSMMatrix<double> makedRdO(CSMMatrix<double> R)
{
	CSMMatrix<double> M(3,3,0.0);
	M(0,0) = 0.0;		M(0,1) = 0.0;		M(0,2) = 0.0;
	M(1,0) = -R(2,0);	M(1,1) = -R(2,1);	M(1,2) = -R(2,2);
	M(2,0) = R(1,0);	M(2,1) = R(1,1);	M(2,2) = R(1,2);
	return M;
};

//Forward rotation matrix Differential (phi)
inline CSMMatrix<double> makedRdP(double o,double p,double k)
{
	CSMMatrix<double> M(3,3,0.0);
	M(0,0) = -cos(k)*sin(p);		M(0,1) = sin(p)*sin(k);			M(0,2) = cos(p);
	M(1,0) = cos(k)*cos(p)*sin(o);	M(1,1) = -sin(k)*cos(p)*sin(o);	M(1,2) = sin(p)*sin(o);
	M(2,0) = -cos(k)*cos(o)*cos(p);	M(2,1) = cos(o)*cos(p)*sin(k);	M(2,2) = -sin(p)*cos(o);
	return M;
};

//Forward rotation matrix Differential (kappa)
inline CSMMatrix<double> makedRdK(double o,double p,double k)
{
	CSMMatrix<double> M(3,3,0.0);
	M(0,0) = -sin(k)*cos(p);						M(0,1) = -cos(p)*cos(k);						M(0,2) = 0.0;
	M(1,0) = cos(o)*cos(k)-sin(k)*sin(p)*sin(o);	M(1,1) = -sin(k)*cos(o)-cos(k)*sin(p)*sin(o);	M(1,2) = 0.0;
	M(2,0) = cos(k)*sin(o)+sin(k)*cos(o)*sin(p);	M(2,1) = -sin(o)*sin(k)+cos(o)*sin(p)*cos(k);	M(2,2) = 0.0;
	return M;
};

}//namespace
#endif