 /*
    This file is part of spectralHMM.

    spectralHMM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    spectralHMM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with spectralHMM.  If not, see <http://www.gnu.org/licenses/>.
  */

package edu.berkeley.spectralHMM.matrix;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Random;

public class LinearSolver {
	
	//solve square system Ax = B using Gaussian elimination with partial pivoting
	//Results will be overwritten into B, and A will be converted into row-echelon form
	//NOTE: Both A & B WILL be changed by the procedure, so make a copy before passing them in.
	//The desired precision and rounding mode can be passed through the MathContext object mc
	public static void geppLinearSolve(BigDecimal[][] A, BigDecimal[] B, MathContext mc)	{
		int n = A.length;
		assert((A[0].length == n) && (B.length == n));
		for (int i = 0; i < n; i++)	{
			//pivot
			int besti = i;
			for (int j = i; j < n; j++)	{
				if (A[j][i].abs(mc).compareTo(A[besti][i].abs(mc)) > 0)	
					besti = j;
			}
		
			//swap rows besti and i
			for (int j = i; j < n; j++)	{
				BigDecimal temp = A[besti][j];
				A[besti][j] = A[i][j];
				A[i][j] = temp;
			}
			BigDecimal temp = B[besti];
			B[besti] = B[i];
			B[i] = temp;
		
			//scale i-th row
			BigDecimal div = A[i][i];
			assert(!div.unscaledValue().equals(BigInteger.ZERO));
			for (int j = i; j < n; j++)	{
				A[i][j] = A[i][j].divide(div, mc);
			}
			B[i] = B[i].divide(div, mc);
		
			//subtract i^th row from all rows
			for (int j = i+1; j < n; j++)	{
				BigDecimal m = A[j][i];
				if (m.unscaledValue().equals(BigInteger.ZERO))	continue;
				for (int k = i+1; k < n; k++)	{
					if (A[i][k].unscaledValue().equals(BigInteger.ZERO))	continue;
					A[j][k] = A[j][k].subtract(A[i][k].multiply(m, mc), mc);
				}
				B[j] = B[j].subtract(B[i].multiply(m, mc), mc);
			}
		}
	
		//back substitution
		for (int i = n-1; i >= 0; i--)	{
			for (int j = i+1; j < n; j++)	{
				B[i] = B[i].subtract(A[i][j].multiply(B[j], mc), mc);
			}
		}
	}
	

	
	//solve square system Ax = B using Gaussian elimination with partial pivoting
	//where A is a banded matrix with the specified bandwidth.
	//Results will be overwritten into B, and A will be converted into row-echelon form
	//NOTE: Both A & B WILL be changed by the procedure, so make a copy before passing them in.
	//The desired precision and rounding mode can be passed through the MathContext object mc		
	public static void geppLinearSolveBanded(BigDecimal[][] A, BigDecimal[] B, int bandwidth, MathContext mc) {
		int n = A.length;
		assert(A[0].length == n & B.length == n);
		for (int i = 0; i < n; i++)	{
			//pivot
			int besti = i;
			for (int j = i; j <= Math.min(i + bandwidth, n-1); j++)	{
				if (A[j][i].abs(mc).compareTo(A[besti][i].abs(mc)) > 0)	
					besti = j;
			}
		
			//swap rows besti and i
			for (int j = i; j <= Math.min(i + 2*bandwidth, n-1); j++)	{
				BigDecimal temp = A[besti][j];
				A[besti][j] = A[i][j];
				A[i][j] = temp;
			}
			BigDecimal temp = B[besti];
			B[besti] = B[i];
			B[i] = temp;
		
			//scale i-th row
			BigDecimal div = A[i][i];
			assert(!div.unscaledValue().equals(BigInteger.ZERO));
			for (int j = i; j <= Math.min(i + 2*bandwidth, n-1); j++)	{
				A[i][j] = A[i][j].divide(div, mc);
			}
			B[i] = B[i].divide(div, mc);
		
			//subtract i^th row from all rows
			for (int j = i+1; j <= Math.min(i + bandwidth, n-1); j++)	{
				BigDecimal m = A[j][i];
				if (m.unscaledValue().equals(BigInteger.ZERO))	continue;
				for (int k = i+1; k <= Math.min(j + 2*bandwidth, n-1); k++)	{
					if (A[i][k].unscaledValue().equals(BigInteger.ZERO))	continue;
					A[j][k] = A[j][k].subtract(A[i][k].multiply(m, mc), mc);
				}
				B[j] = B[j].subtract(B[i].multiply(m, mc), mc);
			}
		}
	
		//back substitution
		for (int i = n-1; i >= 0; i--)	{
			for (int j = i+1; j <= Math.min(i + 2*bandwidth, n-1); j++)	{
				B[i] = B[i].subtract(A[i][j].multiply(B[j], mc), mc);
			}
		}
	}
	
	
	//tests the function, and shows example usage
	public static void main(String[] args)	{
		//100 digits of precision.
		int precision = 150;
		MathContext mc = new MathContext(precision, RoundingMode.HALF_EVEN);
		
		//generate a random square 10x10 linear system
		int dim = 300;
		BigDecimal[][] A = new BigDecimal[dim][dim];
		BigDecimal[][] A2 = new BigDecimal[dim][dim];
		BigDecimal[][] Acopy = new BigDecimal[dim][dim];		
		BigDecimal[] B = new BigDecimal[dim];
		BigDecimal[] B2 = new BigDecimal[dim];
		BigDecimal[] Bcopy = new BigDecimal[dim];
		int band = 4;
		
		Random rnd = new Random();
		for (int i = 0; i < dim; i++)	{
			for (int j = 0; j < dim; j++)	{
				double x = 0.;
				if (Math.abs(i - j) <= band)	{
					x = rnd.nextDouble();
				}
				Acopy[i][j] = BigDecimal.valueOf(x);
				//System.out.print(A[i][j] + "\t");
			}
			double x = rnd.nextDouble();
			Bcopy[i] = BigDecimal.valueOf(x);
			//System.out.println(B[i]);
		}
		//System.out.println();		
		
		long start = System.currentTimeMillis();
		for (int runs = 0; runs < 1200; runs++)	{
			for (int i = 0; i < dim; i++)	{
				for (int j = 0; j < dim; j++)	{
					A[i][j] = Acopy[i][j];
					//A2[i][j] = Acopy[i][j];
				}
				B[i] = Bcopy[i];
				//B2[i] = Bcopy[i];
			}
			
			//Solve A*x = B, where the output x will be stored in B
			LinearSolver.geppLinearSolve(A, B, mc);
		}
		long end = System.currentTimeMillis();
		System.out.printf("Time taken for full solver = %d ms\n", (end - start));

		start = System.currentTimeMillis();
		for (int runs = 0; runs < 1200; runs++)	{
			for (int i = 0; i < dim; i++)	{
				for (int j = 0; j < dim; j++)	{
					A2[i][j] = Acopy[i][j];
				}
				B2[i] = Bcopy[i];
			}
			
			//Solve A*x = B, where the output x will be stored in B
			//LinearSolver.geppLinearSolve(A, B, mc);
			LinearSolver.geppLinearSolveBanded(A2, B2, band, mc);
		}
		end = System.currentTimeMillis();
		System.out.printf("Time taken for banded solver = %d ms\n", end - start);

		
		System.out.printf("Worst error between banded and full solver\n");
		BigDecimal worstErr = BigDecimal.ZERO.setScale(precision);
		for (int i = 0; i < dim; i++)	{
			worstErr = worstErr.max(B[i].subtract(B2[i], mc).abs(mc));
		}
		System.out.println(worstErr);
		
		System.out.println("Worst relative error between A*x and b, where x is the result of solving A*x = b, upto 100 digits of precision");
		BigDecimal worstRelErr = BigDecimal.ZERO.setScale(precision);
		//check if the ouput is correct by comparing if A*x == Bcopy
		for (int i = 0; i < dim; i++)	{
			BigDecimal tot = BigDecimal.ZERO;
			for (int j = 0; j < dim; j++)	{
				tot = tot.add(Acopy[i][j].multiply(B[j], mc), mc);
			}
			//print relative error
			BigDecimal tmp = tot.subtract(Bcopy[i], mc).divide(tot.min(Bcopy[i]), mc);
			worstRelErr = worstRelErr.max(tmp); 
		}
		System.out.println(worstRelErr);
		
	}


	// solves x*A=b
	public static void geppLinearSolveTransposed (BigDecimal[][] A, BigDecimal[] b, MathContext mc) {
		BigDecimal[][] tA= MatrixPower.transpose (A);
		// and solve the transposed problem
		geppLinearSolve(tA, b, mc);
	}


};
