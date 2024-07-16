////////////////////////////////////////////////////////////////////////////////////
/// Liutex method C code.
/// This code calculates liutex when given the veloctiy gradient tensor 
/// for 3D flow fields.
/// 
/// Author:  Oscar Alvarez
/// 
/// email: oscar.alvarez@uta.edu
/// University of Texas at Arlington
/// Department of Mathematics
/// Center for Numerical Simulation and Modeling (CNSM)
/// Arlington, Texas, United States
////////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>

/* ORF this was needed, and compile with -lm */
//#define abs fabs
void liutex(double velocity_gradient_tensor[3][3], double r[3])
{
  /* velocity_gradient_tensor = [  du/dx  du/dy  du/dz
                                   dv/dx  dv/dy  dv/dz
                                   dw/dx  dw/dy  dw/dz  ]
  */
  // r = liutex vector r.

  /// Extract the eigenvalues from the velocity gradient tensor. ///
  /// This part of the method finds the eigenvalues of the velocity gradient tensor.
  /// If you have better methods that do this for you, you can use them.

  double a[3][3];

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      a[i][j] = velocity_gradient_tensor[i][j];
    }
  }
  /// ! Cubic Formula
  /// ! Reference: Numerical Recipes in FORTRAN 77, Second Edition
  /// ! 5.6 Quadratic and Cubic Equations
  /// ! Page 179
  /// !----- Cardano's solution of the cubic equation -----
  /// !---------------------------------------------------------------------
  /// 
  /// ! cubic equation
  /// ! x**3 + aa * x**2 + bb * x + cc = 0
  /// 
  /// ! coefficients of characteristic equation for the velocity gradient tensor.

  double aa = -( a[0][0] + a[1][1] + a[2][2] );

  /// If you can find a method for multiplying matrices, you can use it.
  /// Squaring the velocity gradient tensor (i.e., tt = a^2 = a*a).
  double tt[3][3];
// ORF initialize tt to zero!
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        tt[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        tt[i][j] += a[i][k] * a[k][j];
      }
    }
  }

  double bb = -0.5 * ( tt[0][0] + tt[1][1] + tt[2][2] - pow(a[0][0] + a[1][1] + a[2][2], 2) );

  double cc = -(  a[0][0] * (a[1][1]*a[2][2]-a[1][2]*a[2][1])
                - a[0][1] * (a[1][0]*a[2][2]-a[1][2]*a[2][0])
                + a[0][2] * (a[1][0]*a[2][1]-a[1][1]*a[2][0]) );

  /// delta is the discriminant of characteristic equation for the velocity graidient tensor.
  double qq = (aa * aa - 3.0 * bb) / 9.0;
  double rr = (2.0 * pow(aa, 3) - 9.0 * aa * bb + 27.0 * cc) / 54.0;

  double delta = pow(qq, 3) - pow(rr, 2);

//  printf("bb = %f cc = %f qq = %f rr = %f\n",bb,cc,qq,rr);
//  printf("delta = %f\n",delta);
    
  /// If the discriminant is less than 0 (delta < 0) then the velocity gradient tensor has one
  /// real eigenvalue and two complex conjugate eigenvalues and thus Liutex exists; 
  /// Else, the velocity gradient tensor has three real eigenvalues and Liutex is equal to 0.
  if (delta < 0.0)
  {
    double R;

    double sign;
    if (rr < 0)
    {
      sign = -1.0;
    }
    else if (rr > 0)
    {
      sign = 1.0;
    }
    else
    {
      sign = 0.0;
    }
    
    double aaaa = -sign * pow(abs(rr) + sqrt(-delta), 1.0/3.0);

    double bbbb;
    if (aaaa == 0.0)
    {
      bbbb = 0.0;
    }
    else
    {
      bbbb = qq / aaaa;
    }

    /// eigenvalues 1 and 2 are the complex conjugate eigenvalues of the velocity gradient tensor.
    /// eig3r (eigenvalue 3) is the real eigenvalue of the velocity gradient tensor.
    double eig3r = aaaa + bbbb - aa/3.0;
//    ORF commented this out
//    printf("real eigenvalue: %f \n", eig3r);

    /// Calculating the real right eigenvalue.
    double delta1 = (a[0][0] - eig3r) * (a[1][1] - eig3r) - a[1][0]*a[0][1];
    double delta2 = (a[1][1] - eig3r) * (a[2][2] - eig3r) - a[1][2]*a[2][1];
    double delta3 = (a[0][0] - eig3r) * (a[2][2] - eig3r) - a[0][2]*a[2][0];

    if (abs(delta1) >= abs(delta2) && abs(delta1) >= abs(delta3))
    {
      r[0] = (-(a[1][1]-eig3r)*a[0][2] +         a[0][1]*a[1][2]) / delta1;
      r[1] = (         a[1][0]*a[0][2] - (a[0][0]-eig3r)*a[1][2]) / delta1;
      r[2] = 1.0;
    }
    else if (abs(delta2) >= abs(delta1) && abs(delta2) >= abs(delta3))
    {
      r[0] = 1.0;
      r[1] = (-(a[2][2]-eig3r)*a[1][0] +         a[1][2]*a[2][0])/delta2;
      r[2] = (         a[2][1]*a[1][0] - (a[1][1]-eig3r)*a[2][0])/delta2;
    }
    else if (abs(delta3) >= abs(delta1) && abs(delta3) >= abs(delta2))
    {
      r[0] = (-(a[2][2]-eig3r)*a[0][1] +         a[0][2]*a[2][1])/delta3;
      r[1] = 1.0;
      r[2] = (         a[2][0]*a[0][1] - (a[0][0]-eig3r)*a[2][1])/delta3;
    }
    else
    {
      printf("ERROR: delta1, delta2, delta3\n");
      return;
    }

    double r_norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    r[0] = r[0] / r_norm;
    r[1] = r[1] / r_norm;
    r[2] = r[2] / r_norm;

    /// Calculate rotation matrix r which rotates unit vector u to unit vector v.
    /// z0 = 0, 0, 1; 
    /// r = liutex_vec
    /// qqq = transformation matrix

    double z0[3] = { 0.0, 0.0, 1.0 };
    double temp_vec[3];
    double transformation_matrix[3][3];

	/* ORF had to change from eps to epsilon, eps must be taken */
    double epsilon = 1.0e-9;

    temp_vec[0] = z0[1] * r[2] - z0[2] * r[1];
    temp_vec[1] = z0[2] * r[0] - z0[0] * r[2];
    temp_vec[2] = z0[0] * r[1] - z0[1] * r[0];

    aa = sqrt( temp_vec[0]*temp_vec[0] + temp_vec[1]*temp_vec[1] + temp_vec[2]*temp_vec[2] );
        
    if (aa < epsilon)
    {
      transformation_matrix[0][0] = 1.0;
      transformation_matrix[1][1] = 1.0;
      transformation_matrix[2][2] = 1.0;
    }
    else
    {
      temp_vec[0] = temp_vec[0] / aa;
      temp_vec[1] = temp_vec[1] / aa;
      temp_vec[2] = temp_vec[2] / aa;

      double t = z0[0] * r[0] + z0[1] * r[1] + z0[2] * r[2];
            
      if (t > 1.0)
      {
        t = 1.0;
      }

      if (t < -1.0)
      {
        t = -1.0;
      }
            
      double alpha = acos(t);
        
      double c = cos(alpha);
      double s = sin(alpha);
        
      transformation_matrix[0][0] = temp_vec[0] * temp_vec[0] * (1.0 - c) + c;
      transformation_matrix[0][1] = temp_vec[0] * temp_vec[1] * (1.0 - c) - temp_vec[2] * s;
      transformation_matrix[0][2] = temp_vec[0] * temp_vec[2] * (1.0 - c) + temp_vec[1] * s;
        
      transformation_matrix[1][0] = temp_vec[1] * temp_vec[0] * (1.0 - c) + temp_vec[2] * s;
      transformation_matrix[1][1] = temp_vec[1] * temp_vec[1] * (1.0 - c) + c;
      transformation_matrix[1][2] = temp_vec[1] * temp_vec[2] * (1.0 - c) - temp_vec[0] * s;
        
      transformation_matrix[2][0] = temp_vec[2] * temp_vec[0] * (1.0 - c) - temp_vec[1] * s;
      transformation_matrix[2][1] = temp_vec[2] * temp_vec[1] * (1.0 - c) + temp_vec[0] * s;
      transformation_matrix[2][2] = temp_vec[2] * temp_vec[2] * (1.0 - c) + c;
    }

    /// If you can find a method for multiplying matrices, you can use it.
    /// vg = transpose(transformation_matrix) * velocity_gradient_tensor * transformation_matrix.
        
    /// Transpose the transformation matrix.
    double transposed_transformation_matrix[3][3];
        
    for (int i = 0; i < 3 ; i++)
    {
      for (int j = 0; j < 3 ; j++)
      {
        transposed_transformation_matrix[i][j] = transformation_matrix[j][i];
      }

    }

    /// matrix_product_1 = transpose_transformation_matrix * velocity_gradient_tensor.
    double matrix_product_1[3][3];

    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        matrix_product_1[i][j] = 0.0;

        for (int k = 0; k < 3; k++)
        {
          matrix_product_1[i][j] += transposed_transformation_matrix[i][k] * a[k][j];
        }
      }
    }

    /// vgt = velocity_gradient_tensor * transformation_matrix.
    double vgt[3][3];

    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        vgt[i][j] = 0.0;

        for (int k = 0; k < 3; k++)
        {
          vgt[i][j] += a[i][k] * transformation_matrix[k][j];
        }
      }
    }

    double alpha = 0.5 * sqrt( pow(vgt[1][1] - vgt[0][0], 2) + pow(vgt[1][0] + vgt[0][1], 2) );
    double beta  = 0.5 * (vgt[1][0] - vgt[0][1]);

    if (beta*beta > alpha*alpha)
    {
      if(beta > 0.0)
      {
        R = 2.0 * (beta - alpha);
        r[0] = R * r[0];
        r[1] = R * r[1];
        r[2] = R * r[2];
      }
      else
      {
        R = 2.0 * (beta + alpha);
        r[0] = R * r[0];
        r[1] = R * r[1];
        r[2] = R * r[2];
      }
    }
    else
    {
      r[0] = 0.0;
      r[1] = 0.0;
      r[2] = 0.0;
    }
  }
  else
  {
    /// All eigenvalue are real so Liutex = 0.
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;
  }

}

