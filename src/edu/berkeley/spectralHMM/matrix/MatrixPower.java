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
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;

public class MatrixPower {
	
	public static BigDecimal[][] multiplyMatricesBanded (BigDecimal[][] A, BigDecimal[][] B, int bwA, int bwB, MathContext mc) {
		int n = A.length;
		int m = A[0].length;
		int p = B[0].length;
		assert(m == B.length);
		
		// initialize (we assume them squared)
		BigDecimal[][] result = new BigDecimal[n][p];
		for (int i=0; i<result.length; i++) {
			Arrays.fill(result[i], BigDecimal.ZERO);
		}
		
		// compute		
		for (int i=0; i<result.length; i++) {
			for (int j=Math.max(0, i-bwA-bwB); j<=Math.min(p-1, i+bwA+bwB); j++) {
				// multiply row of A with the column of B, recognizing the bandwidth
				for (int k=Math.max(0, Math.max(i-bwA, j-bwB)); k <= Math.min (m-1, Math.min(i+bwA, j+bwB)); k++) {
					result[i][j] = result[i][j].add (A[i][k].multiply(B[k][j], mc), mc);
				}
			}
		}
		
		// return it
		return result;
	}
	

	//returns A*v, where A has bandwidth bw
	public static BigDecimal[] matrixVectorBanded (BigDecimal[][] A, BigDecimal[] v, int bw, MathContext mc)	{
		int n = A.length, m = A[0].length;
		assert(m == v.length);
		BigDecimal[] ret = new BigDecimal[n];
		for (int i = 0; i < n; i++)	{
			ret[i] = BigDecimal.ZERO.setScale(mc.getPrecision());
			for (int j = Math.max(i - bw, 0); j <= Math.min(i + bw, m-1); j++)	{
				ret[i] = ret[i].add(A[i][j].multiply(v[j], mc), mc);
			}
		}
		return ret;
	}
	
	// return v*A
	public static BigDecimal[] vectorMatrixBanded (BigDecimal[] v, BigDecimal[][] A, int bw, MathContext mc)	{
		BigDecimal[][] newA = transpose (A);
		return matrixVectorBanded (newA, v, bw, mc);
	}

	public static BigDecimal[][] getIdentityMatrixBigDecimal(int n) {
		BigDecimal[][] id = new BigDecimal[n][n];
		for (int i = 0; i < n; i++)	{
			Arrays.fill(id[i], BigDecimal.ZERO);
			id[i][i] = BigDecimal.ONE;
		}
		return id;
	}
	
	public static BigDecimal[][] powerMatrixBanded (BigDecimal[][] A, int bwA, int p, MathContext mc) {
		assert(A.length == A[0].length);
		
		if (p == 0)	{	//return identity matrix
			return getIdentityMatrixBigDecimal(A.length);
		}
		
		if (p == 1)	{
			return A;
		}
		
		BigDecimal[][] sq = powerMatrixBanded (A, bwA, p / 2, mc);
		BigDecimal[][] ret = multiplyMatricesBanded(sq, sq, (p/2)*bwA, (p/2)*bwA, mc);
		if (p % 2 == 1)	{
			ret = multiplyMatricesBanded(A, ret, bwA, p*bwA, mc);
		}
		
		return ret;
	}
	
	public static void main(String[] args){
		BigDecimal[][] A = getIdentityMatrixBigDecimal(10);
		for (int i = 0; i < A.length; i++)	{
			if (i > 0)	A[i][i-1] = new BigDecimal("-1");
			if (i < A.length - 1)	A[i][i+1] = new BigDecimal("-1");
		}
		BigDecimal[][] B = getIdentityMatrixBigDecimal(10);
		for (int i = 0; i < B.length; i++)	{			
			if (i > 0)	B[i][i-1] = new BigDecimal("-2");
			if (i < B.length - 1)	B[i][i+1] = new BigDecimal("-2");
		}
		
		MathContext mc = new MathContext(10, RoundingMode.HALF_EVEN);
		BigDecimal[][] C = multiplyMatricesBanded(A, B, 1, 1, mc);
		for (int i = 0; i < C.length; i++)	{
			for (int j = 0; j < C.length; j++)	{
				System.out.print(C[i][j]+ "\t");
			}
			System.out.println();
		}
		
	}

	public static BigDecimal[][] transpose(BigDecimal[][] w) {
		int n = w.length, m = w[0].length;
		BigDecimal[][] result = new BigDecimal[m][n];
		
		for (int i = 0; i < n; i++)	{
			for (int j = 0; j < m; j++)	{
				result[j][i] = w[i][j];
			}
		}
		
		return result;
	}

	// the bandwidth is the full matrix dimension
	public static BigDecimal[] vectorMatrix (BigDecimal[] v, BigDecimal[][] A, MathContext mc) {
		int bandWidth = Math.max(A.length, A[0].length);
		return vectorMatrixBanded(v, A, bandWidth, mc);
	}


	// take the full matrix dimensions as bandwith
	public static BigDecimal[][] multiplyMatrices(BigDecimal[][] A,	BigDecimal[][] B, MathContext mc) {
		int bwA = Math.max(A.length, A[0].length);
		int bwB = Math.max(B.length, B[0].length);
		return multiplyMatricesBanded (A, B, bwA, bwB, mc);
	}

}
