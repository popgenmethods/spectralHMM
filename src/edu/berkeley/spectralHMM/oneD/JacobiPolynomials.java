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

package edu.berkeley.spectralHMM.oneD;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.nevec.rjm.BigDecimalMath;

import edu.berkeley.spectralHMM.datatypes.Pair;
import edu.berkeley.spectralHMM.matrix.LinearSolver;
import edu.berkeley.spectralHMM.matrix.MatrixPower;

public class JacobiPolynomials {

	// so you can't create an object of it
	private JacobiPolynomials() {
		;
	}
	
	// evaluate Jacobi polynomial at x of given DEGREE n, even though this may be the (n-2)-th polynomial in 
	// the sequence of eigenpolynomials for some alpha, beta  
	static public BigDecimal evaluate(BigDecimal alpha, BigDecimal beta, int n, BigDecimal x, MathContext mc)	{
		PolynomialMetaArguments metaArgs = new PolynomialMetaArguments(alpha, beta, mc);
		TreeMap<Pair<Integer, BigDecimal>, BigDecimal> t = findTreeMap (cachePoly, metaArgs);
		assert(t != null);
		Pair<Integer, BigDecimal> cacheEntry = new Pair<Integer, BigDecimal>(n, x);
		BigDecimal value = t.get(cacheEntry);
		if (value != null)	{
			return value;
		}
		
		
		int offset = calculateOffset(alpha, beta);
				
		if (n == 0)	{
			assert (offset == 0);
			value = BigDecimal.ONE.setScale(mc.getPrecision());
		}
		else if (n == 1)	{
			assert ((offset == 0) || (offset == 1));
			// polynomial one  is ((a + b) x - a)
			value = x.multiply(alpha.add(beta, mc), mc).subtract(alpha, mc);
		}
		else if (n == 2 && (offset == 1 || offset == 2)) {
			// polynomial two is 1/2 b (1 + b) + (1 + b) (1 + a + b) (-1 + x) + 1/2 (1 + a + b) (2 + a + b) (-1 + x)^2
			// TODO cache them
			BigDecimal t1 = HALF.multiply(beta, mc).multiply(BigDecimal.ONE.add(beta, mc), mc);
			BigDecimal f2 = BigDecimal.ONE.add(beta, mc).multiply(BigDecimal.ONE.add(alpha, mc).add(beta, mc), mc);
			BigDecimal f3 = HALF.multiply(BigDecimal.ONE.add(alpha, mc).add(beta, mc),mc).multiply(TWO.add(alpha, mc).add(beta, mc),mc);
			
			BigDecimal xm1 = x.subtract(BigDecimal.ONE, mc);
			BigDecimal squared = xm1.multiply(xm1, mc);
			// now put them together
			value = t1.add(f2.multiply(xm1, mc), mc).add(f3.multiply(squared, mc), mc);
		}
		else if (n == 3 && offset == 2) {
			// polynomial three is 2 x (1 - 3 x + 2 x^2)
			value = TWO.multiply(x, mc).multiply(BigDecimal.ONE.subtract((new BigDecimal("3")).multiply(x, mc), mc).add(TWO.multiply(x, mc).multiply(x, mc) , mc), mc);
		}
		else {
			assert(n >= offset + 2);
			// now we have the first two filled, so kick off the recursion
			// but be aware of the offset
			BigDecimal pn1 = evaluate(alpha, beta, n - 1, x, mc);
			BigDecimal pn2 = evaluate(alpha, beta, n - 2, x, mc);
			value = jacobiRecursionHelper(alpha, beta, n, x, pn1, pn2, mc);
		}
		
		// memoize it
		t.put(cacheEntry, value);
		return value;
	}
	
	// evaluate the derivative of the Jacobi polynomial at x of given DEGREE n
	static public BigDecimal evaluateDerivative (BigDecimal alpha, BigDecimal beta, int n, BigDecimal x, MathContext mc) {
		//try to find it in old lists
		PolynomialMetaArguments metaArgs = new PolynomialMetaArguments(alpha, beta, mc);
		TreeMap<Pair<Integer, BigDecimal>, BigDecimal> t = findTreeMap (cachePolyDerivative, metaArgs);
		assert(t != null);
		Pair<Integer, BigDecimal> cacheEntry = new Pair<Integer, BigDecimal>(n, x);
		BigDecimal value = t.get(cacheEntry);
		if (value != null)	{
			return value;
		}
		
		int offset = calculateOffset(alpha, beta);
		
		// now get the derivative right
		if (n == 0)	{
			value = BigDecimal.ZERO;
		}
		else if (n == 1)	{
			// polynomial one  is ((a + b) x - a), so derivative is (a + b)
			value = alpha.add(beta, mc);
		}
		else if (n == 2 && (offset == 1 || offset == 2)) {
			// polynomial two is 1/2 b (1 + b) + (1 + b) (1 + a + b) (-1 + x) + 1/2 (1 + a + b) (2 + a + b) (-1 + x)^2
			// so the derivative is (1 + b) (1 + a + b) + (1 + a + b) (2 + a + b) (-1 + x)
			BigDecimal t1 = BigDecimal.ONE.add(beta, mc).multiply(BigDecimal.ONE.add(alpha, mc).add(beta, mc), mc);
			BigDecimal f2 = BigDecimal.ONE.add(alpha, mc).add(beta, mc).multiply(TWO.add(alpha, mc).add(beta, mc),mc);
			
			// get the xs
			BigDecimal xm1 = x.subtract(BigDecimal.ONE, mc);

			// now put them together
			value = t1.add(f2.multiply(xm1, mc), mc);
		}
		else if (n == 3 && offset == 2) {
			// polynomial three is 2 x (1 - 3 x + 2 x^2)
			// so derivative is 
			value = TWO.subtract(new BigDecimal("12").multiply(x, mc),mc).add(new BigDecimal("12").multiply(x, mc).multiply(x, mc), mc);
		}
		else {
			assert(n >= offset + 2);
			// now we have the first two filled, so kick off the recursion
			// but be aware of the offset
			BigDecimal Pn1 = evaluate (alpha, beta, n - 1, x, mc);
			BigDecimal dPn1 = evaluateDerivative (alpha, beta, n - 1, x, mc);
			BigDecimal dPn2 = evaluateDerivative (alpha, beta, n - 2, x, mc);
			value = jacobiDerivativeRecursionHelper (alpha, beta, n, x, Pn1, dPn1, dPn2, mc);
		}
		
		// memoize it
		t.put(cacheEntry, value);
		return value;
	}
	
	// evaluate n-th Jacobi polynomial squared length, which can even be a degree n+2 polynomial 
	static public BigDecimal evaluateSquaredLength(BigDecimal alpha, BigDecimal beta, int n, MathContext mc)	{
		// is value already computed
		PolynomialMetaArguments metaArgs = new PolynomialMetaArguments(alpha, beta, mc);
		HashMap<Integer, BigDecimal> t = findHashMap (cacheSL, metaArgs);
		BigDecimal value = t.get(n);
		if (value != null)	{
			return value;
		}

		int offset = calculateOffset(alpha, beta);
		
		// here it goes		
		// fill it
		if (n == 0)	{
			assert (offset == 0);
			value = betaFunction (alpha, beta, mc);
		}
		else if (n == 1)	{
			assert ((offset == 0) || (offset == 1));
			value = alpha.multiply(beta, mc).divide(alpha.add(beta, mc).add(BigDecimal.ONE, mc), mc).multiply(betaFunction(alpha, beta, mc), mc);
		}
		else if (offset == 2 && n == 2) {
			value = BigDecimal.ONE.divide(new BigDecimal("6"), mc);
		}
		else {
			assert ((((offset == 0) || (offset == 1)) && (n == 2)) || (n > 2));
			// get the factor
			// (2n + alpha + beta - 3)/(2n + alpha + beta -1) * ((n + alpha - 1)*(n + beta - 1))/(n*(n + alpha + beta - 2))
			BigDecimal bN = new BigDecimal(n);
			BigDecimal TWOnab = TWO.multiply(bN, mc).add(alpha, mc).add(beta, mc);
			BigDecimal f1 = TWOnab.subtract(new BigDecimal("3"), mc).divide(TWOnab.subtract(BigDecimal.ONE, mc), mc);
			BigDecimal n2 = bN.add(alpha, mc).subtract(BigDecimal.ONE, mc).multiply(bN.add(beta, mc).subtract(BigDecimal.ONE, mc), mc);
			BigDecimal d2 = bN.multiply(bN.add(alpha, mc).add(beta, mc).subtract(TWO, mc), mc);
			BigDecimal factor = f1.multiply(n2, mc).divide(d2, mc);
			// set it
			value = factor.multiply (evaluateSquaredLength (alpha, beta, n-1, mc), mc);

		}
		
		// memoize it
		t.put(n, value);
		return value;
	}
	
	public static int calculateOffset(BigDecimal alpha, BigDecimal beta) {
		// calculate the offset
		int offset = 0;
		if (alpha.unscaledValue() == BigInteger.ZERO) {
			offset += 1;
		}
		if (beta.unscaledValue() == BigInteger.ZERO) {
			offset += 1;
		}
		// set it
		assert ((offset >= 0) && (offset <= 2));
		return offset;
	}

	// the G-factor
	public static BigDecimal G (BigDecimal alpha, BigDecimal beta, int n, int k, MathContext mc) {
		// should not be requested, and we would get in trouble (actually we might not)
		// and we will request at least 0,0 and 0,1 in the future
		assert (n>0 || k >= 0);
		
		PolynomialMetaArguments metaArgs = new PolynomialMetaArguments(alpha, beta, mc);
		HashMap<Pair<Integer, Integer>, BigDecimal> t = findHashMap (cacheG, metaArgs);
		Pair<Integer, Integer> recEntry = new Pair<Integer, Integer>(n, k);
		BigDecimal retValue = t.get(recEntry);
		if (retValue != null)	{
			return retValue;
		}
		
		// locals
		BigDecimal t1,t2;
		BigDecimal TWO = new BigDecimal("2");
		BigDecimal N = new BigDecimal(n);
		BigDecimal TWOnab = TWO.multiply(N, mc).add(alpha, mc).add(beta, mc);
		
		retValue = BigDecimal.ZERO;
		
		if (k == n+1) {
			if (n == 0) {
				retValue = BigDecimal.ONE.divide(alpha.add(beta, mc), mc);
				// return 1.0/(alpha+beta);
			} else {
				t1 = N.add(alpha, mc).add(beta, mc).subtract(BigDecimal.ONE, mc);
				retValue = t1.multiply(N.add (BigDecimal.ONE, mc), mc).divide(TWOnab, mc).divide(TWOnab.subtract(BigDecimal.ONE, mc), mc);
				// return (n + alpha + beta - 1) * (n + 1) /((2*n + alpha + beta) * (2*n + alpha + beta - 1));
			}
		} else if (k == n) {
			if (n == 0) {
				retValue = alpha.divide (alpha.add(beta, mc), mc);
				// return alpha/(alpha+beta);
			} else {
				t1 = beta.multiply(beta, mc).subtract(alpha.multiply(alpha, mc), mc).subtract(TWO.multiply(beta.subtract(alpha, mc), mc), mc);
				t2 = TWO.multiply(TWOnab, mc).multiply(TWOnab.subtract(TWO, mc), mc);
				retValue = (new BigDecimal("0.5")).subtract(t1.divide(t2, mc), mc);
				// return (0.5 - (beta*beta - alpha*alpha - 2*(beta - alpha)) /(2 * (2*n + alpha + beta) * (2*n + alpha + beta - 2)));
			}
		} else if (k == n-1) {
			t1 = N.add(beta, mc).subtract(BigDecimal.ONE, mc);
			t2 = N.add(alpha, mc).subtract(BigDecimal.ONE, mc);
			retValue = t1.multiply(t2, mc).divide(TWOnab.subtract(BigDecimal.ONE, mc), mc).divide(TWOnab.subtract(TWO, mc), mc); 
			// return (n + beta - 1) * (n + alpha - 1) /((2*n + alpha + beta - 1) * (2*n + alpha + beta - 2));
		} else {
			//not allowed
			assert (false);
		}

		t.put(new Pair<Integer, Integer>(n, k), retValue);
		return retValue;
	}
	
	//multiply nRescale / nrSteps * G * (nRescale - 1) / (nrSteps - 1) * G * ...
	public static BigDecimal[][] coefficientMatrixRescaled(BigDecimal alpha, BigDecimal beta, int nRescale, int nrSteps, int maxM, MathContext mc, boolean isL, boolean rescale) {
		if (rescale)	{
			int dim = maxM + 1;
			BigDecimal[][] ret = MatrixPower.getIdentityMatrixBigDecimal(dim);
			BigDecimal[][] v = coefficientMatrix(alpha, beta, 1, maxM, mc, isL);
			for (int i = 0; i < nrSteps; i++)	{
				BigDecimal[][] foo = new BigDecimal[dim][dim];
				BigDecimal factor = (new BigDecimal(nRescale - i)).setScale(mc.getPrecision()).divide((new BigDecimal(nrSteps - i)).setScale(mc.getPrecision()), mc);
				for (int k = 0; k < dim; k++)	{
					for (int p = Math.max(0, k - 1); p < Math.min(dim, k + 2); p++)	{
						foo[k][p] = v[k][p].multiply(factor, mc);
					}
				}
				ret = MatrixPower.multiplyMatricesBanded(ret, foo, i, 1, mc);
			}
			return ret;
		}
		else	{
			return coefficientMatrix(alpha, beta, nrSteps, maxM, mc, isL);
		}
		
	}

	public static BigDecimal[][] coefficientMatrix (BigDecimal alpha, BigDecimal beta, int nrSteps, int maxM, MathContext mc, boolean isL) {
		assert(nrSteps >= 0);		

		// get the right cache
		MatrixMetaArguments metaArgs = new MatrixMetaArguments (alpha, beta, maxM, mc);
		HashMap<Integer, BigDecimal[][]> t;
		if (isL) {
			t = findHashMap (cacheMatrixL, metaArgs);
		}
		else {
			t = findHashMap (cacheMatrixM, metaArgs);
		}
		assert(t != null);
		
		// get it from cache, if it is there
		BigDecimal[][] result = t.get(nrSteps);
		if (result != null) {
			return result;
		}

		// for 0 steps we just have the identity
		if (nrSteps == 0) {
			// the identity
			result = MatrixPower.getIdentityMatrixBigDecimal(maxM+1);
		}
		// for 1 step, we have the one step matrix
		else if (nrSteps == 1)	{
			// build the matrix
			BigDecimal[][] opMat = null;
			opMat = getThreeTermRecurrenceMatrix (alpha, beta, maxM, mc);
			if (!isL) {
				// modify it to be right for M
				for (int i=0; i<=maxM; i++) {
					for (int j = Math.max(0, i-1); j <= Math.min(maxM, i+1); j++)	{
						opMat[i][j] = BigDecimal.ZERO.subtract(opMat[i][j], mc);
					}
					opMat[i][i] = BigDecimal.ONE.add(opMat[i][i], mc);
				}
			}
			
			// and remember the matrix
			result = opMat;
		}
		// or it's higher than that
		else {
			// whats the highest power of two in nrSteps
			int hp2 = highestPowerOfTwo(nrSteps);
			
			// is the nr of steps actually the highest power of two
			if (hp2 == nrSteps) {
				// get the half-step matrix
				BigDecimal[][] mat = coefficientMatrix (alpha, beta, nrSteps/2, maxM, mc, isL);
				// and square
				result = MatrixPower.multiplyMatricesBanded (mat, mat, nrSteps/2, nrSteps/2, mc);
			}
			else {
				// get the different components
				BigDecimal[][] mat1 = coefficientMatrix (alpha, beta, hp2, maxM, mc, isL);
				BigDecimal[][] mat2 = coefficientMatrix (alpha, beta, nrSteps - hp2, maxM, mc, isL);
				// and multiply it
				result = MatrixPower.multiplyMatricesBanded (mat1, mat2, hp2, nrSteps - hp2, mc);
			}
		}
		
		// store the result
		t.put (nrSteps, result);

		
		//and return it
		return result;
	}
	
	private static BigDecimal[][] getThreeTermRecurrenceMatrix (BigDecimal alpha, BigDecimal beta, int maxM, MathContext mc) {

		BigDecimal[][] result = new BigDecimal[maxM+1][maxM+1];
		
		// I think we need an offset here
		int offset = calculateOffset (alpha, beta);
//		int offset = 0;
		
		// fill it
		for (int i=0; i<=maxM; i++) {
			Arrays.fill (result[i], BigDecimal.ZERO);
			
			// put the three factors at the right place
			for (int j = Math.max(0, i-1); j <= Math.min(maxM, i+1); j++)	{
				result[i][j] = G (alpha, beta, offset + i, offset + j, mc);
			}
		}
		
		return result;
	}

	public static int highestPowerOfTwo (int number) {
		int result = 1;
		while (result <= number) {
			result *= 2;
		}
		return result/2;
	}
	
	// the L (start at startIdx and go in nrSteps steps to endIdx)
	public static BigDecimal L (BigDecimal alpha, BigDecimal beta, int startIdx, int endIdx, int nrSteps, int offset, MathContext mc) {
		// safety first
		assert (startIdx >= offset && endIdx >= offset && nrSteps >= 0);

		PolynomialMetaArguments metaArgs = new PolynomialMetaArguments(alpha, beta, mc);
		HashMap<RecurrenceEntry, BigDecimal> t = findHashMap (cacheL, metaArgs);
		assert(t != null);
		RecurrenceEntry recEntry = new RecurrenceEntry(startIdx, endIdx, nrSteps);
		BigDecimal result = t.get(recEntry);
		if (result != null) {
			return result;
		}
		
		// the return value
		result = BigDecimal.ZERO;
		// how long can the path be?
		if (nrSteps == 0) {		
			// no path at all
			assert (startIdx == endIdx);
			
			result = BigDecimal.ONE;
		}
		else if (nrSteps == 1) {
			// path of length one, this one should also be doable
			assert (Math.abs (startIdx-endIdx) <= 1);

			// G should do the rest
			// kappa can not be less then zero
			result = G (alpha, beta, startIdx, endIdx, mc);
		}
		else {
			// add multiplications along all paths of a general length from some general mu to some general (but reachable) kappa
			// we do this recursively
			// this might not be the most efficient for large l, but since we are only interested in l <= 4, we should be fine
			// also less error prone
			result = BigDecimal.ZERO;
			
			// so first step can be down, stay or up
			
			// down, but only if kappa stays reachable and we are not to small (mu > 0)
			if (startIdx > offset && Math.abs(startIdx-1 - endIdx) <= nrSteps-1) {
				// first step down, reach same kappa, one step less
				result = result.add(G(alpha, beta, startIdx, startIdx-1, mc).multiply(L(alpha, beta, startIdx-1, endIdx, nrSteps-1, offset, mc), mc), mc);
			}

			// stay, but only if kappa stays reachable
			if (Math.abs(startIdx - endIdx) <= nrSteps-1) {
				// first step stay, reach same kappa, one step less
				result = result.add(G(alpha, beta, startIdx, startIdx, mc).multiply(L(alpha, beta, startIdx, endIdx, nrSteps-1, offset, mc), mc), mc);
			}

			// up, but only if kappa stays reachable
			if (Math.abs(startIdx+1 - endIdx) <= nrSteps-1) {
				// first step up, reach same kappa, one step less
				result = result.add(G(alpha, beta, startIdx, startIdx+1, mc).multiply(L(alpha, beta, startIdx+1, endIdx, nrSteps-1, offset, mc), mc), mc);
			}
			
		}
		
		// give it away now
		t.put(recEntry, result);
		
		if (startIdx == 10 && nrSteps == 6) {
			System.out.println (endIdx + "\t" + result);
		}

		
		return result;
	}

	// the M (start at startIdx and go in nrSteps steps to endIdx)
	public static BigDecimal M (BigDecimal alpha, BigDecimal beta, int startIdx, int endIdx, int nrSteps, int offset, MathContext mc) {
		// safety first
		assert (startIdx >= offset && endIdx >= offset && nrSteps >= 0);
		
		PolynomialMetaArguments metaArgs = new PolynomialMetaArguments(alpha, beta, mc);
		HashMap<RecurrenceEntry, BigDecimal> t = findHashMap (cacheM, metaArgs);
		RecurrenceEntry recEntry = new RecurrenceEntry(startIdx, endIdx, nrSteps);
		BigDecimal result = t.get(recEntry);
		if (result != null)	{
			return result;
		}
		
		// the return value
		result = BigDecimal.ZERO;
		// how long can the path be?
		if (nrSteps == 0) {		
			// no path at all
			assert (startIdx == endIdx);
			
			result = BigDecimal.ONE;
		}
		else if (nrSteps == 1) {
			// path of length one, this one should also be doable
			assert (Math.abs (startIdx-endIdx) <= 1);

			// G should do the rest
			// kappa can not be less then zero
			if (startIdx == endIdx) {
				result = BigDecimal.ONE.subtract(G (alpha, beta, startIdx, endIdx, mc), mc);
			}
			else {
				result = BigDecimal.ZERO.subtract(G (alpha, beta, startIdx, endIdx, mc), mc);
			}
		}
		else {
			// add multiplications along all paths of a general length from some general mu to some general (but reachable) kappa
			// we do this recursively
			// this might not be the most efficient for large l, but since we are only interested in l <= 4, we should be fine
			// also less error prone
			result = BigDecimal.ZERO;
			
			// so first step can be down, stay or up
			
			// down, but only if kappa stays reachable and we are not to small (mu > 0)
			if (startIdx > offset && Math.abs(startIdx-1 - endIdx) <= nrSteps-1) {
				// first step down, reach same kappa, one step less
				result = result.subtract(G(alpha, beta, startIdx, startIdx-1, mc).multiply(M(alpha, beta, startIdx-1, endIdx, nrSteps-1, offset, mc), mc), mc);
			}

			// stay, but only if kappa stays reachable
			if (Math.abs(startIdx - endIdx) <= nrSteps-1) {
				// first step stay, reach same kappa, one step less
				result = result.add(BigDecimal.ONE.subtract(G(alpha, beta, startIdx, startIdx, mc), mc).multiply(M(alpha, beta, startIdx, endIdx, nrSteps-1, offset, mc), mc), mc);
			}

			// up, but only if kappa stays reachable
			if (Math.abs(startIdx+1 - endIdx) <= nrSteps-1) {
				// first step up, reach same kappa, one step less
				result = result.subtract(G(alpha, beta, startIdx, startIdx+1, mc).multiply(M(alpha, beta, startIdx+1, endIdx, nrSteps-1, offset, mc), mc), mc);
			}
			
		}
		
		// give it away now
		t.put(recEntry, result);
		return result;
	}
	
	
	// the lambda 
	public static BigDecimal lambda (int n, BigDecimal alpha, BigDecimal beta, MathContext mc) {
		BigDecimal bn = new BigDecimal(n);
		return (bn.multiply(bn.add(alpha,mc).add(beta,mc).subtract(BigDecimal.ONE,mc), mc).multiply(new BigDecimal("0.5"),mc));
	}

	// the beta function
	public static BigDecimal betaFunction (BigDecimal x, BigDecimal y, MathContext mc) {
		// gamma(x)*gamma(y)/gamma(x+y)
		return BigDecimalMath.Gamma(x, mc).multiply(BigDecimalMath.Gamma(y, mc), mc).divide(BigDecimalMath.Gamma(x.add(y, mc), mc), mc);
	}

	//converts sum_{j=0}^{n-1} coeffs[j]*R_j(x) to sum_{j=0}^{n-1} output[j]*x^j for the given alpha, beta parameters
	public static BigDecimal[] convertJacobiToMonomialSequence (BigDecimal alpha, BigDecimal beta, BigDecimal[] coeffs, MathContext mc) {
		int n = coeffs.length;
		
		if (n == 1)	{	//if constant, nothing to do
			BigDecimal[] output = new BigDecimal[n];
			output[0] = coeffs[0].multiply(evaluate(alpha, beta, 0, BigDecimal.ZERO, mc), mc);
			return output;
		}
		
		//we will eventually solve Ax = b
		BigDecimal[][] A = new BigDecimal[n][n];
		BigDecimal[] b = new BigDecimal[n];
		
		//we have at least a degree 1 polynomial
		//evaluate at a grid of n points and solve the linear system
		//BigDecimal[] grid = new BigDecimal[n];
		BigDecimal nminus1 = new BigDecimal(""+(n-1));
		for (int i = 0; i < n; i++){
			BigDecimal gridpt = BigDecimal.ONE.divide(nminus1, mc);
			gridpt = gridpt.multiply(new BigDecimal(""+i), mc);
			b[i] = BigDecimal.ZERO;
			for (int j = 0; j < n; j++)	{
				b[i] = b[i].add(coeffs[j].multiply(evaluate(alpha, beta, j, gridpt, mc), mc), mc);
				A[i][j] = gridpt.pow(j, mc);
			}
		}
		
		LinearSolver.geppLinearSolve(A, b, mc);
		
		return b;
	}
	
	//converts sum_{j=0}^{n-1} coeffs[j]*x^j to sum_{j=0}^{n-1} output[j]*R_j(x) for the given alpha, beta parameters
	public static BigDecimal[] convertMonomialToJacobiSequence (BigDecimal alpha, BigDecimal beta, BigDecimal[] coeffs, MathContext mc) {
		int n = coeffs.length;
		
		if (n == 1)	{	//if constant, nothing to do
			BigDecimal[] output = new BigDecimal[n];
			output[0] = coeffs[0].divide(evaluate(alpha, beta, 0, BigDecimal.ZERO, mc), mc);
			return output;
		}
		
		//we will eventually solve Ax = b
		BigDecimal[][] A = new BigDecimal[n][n];
		BigDecimal[] b = new BigDecimal[n];
		
		//we have at least a degree 1 polynomial
		//evaluate at a grid of n points and solve the linear system
		//BigDecimal[] grid = new BigDecimal[n];
		BigDecimal nminus1 = new BigDecimal(""+(n-1));
		for (int i = 0; i < n; i++){
			BigDecimal gridpt = BigDecimal.ONE.divide(nminus1, mc);
			gridpt = gridpt.multiply(new BigDecimal(""+i), mc);
			b[i] = BigDecimal.ZERO;
			for (int j = 0; j < n; j++)	{
				b[i] = b[i].add(coeffs[j].multiply(gridpt.pow(j, mc), mc), mc);
				A[i][j] = evaluate(alpha, beta, j, gridpt, mc);
			}
		}
		
		LinearSolver.geppLinearSolve(A, b, mc);
		
		return b;
	}

	//Inverts a formal series containing the Jacobi polynomials.
	//Specifically, takes p(x) = \sum_{j=0}^{n-1} coeffs[j]*R_j(x) and finds a q(x) = \sum_{j=0}^{n-1} output[j]*R_j(x) 
	//such that p(x)*q(x) = 1 + O(x^n) (i.e. the first non-zero monomial after the constant term 1 is of degree n)
	public static BigDecimal[] oldInvertJacobiSequence (BigDecimal alpha, BigDecimal beta, BigDecimal[] coeffs, MathContext mc) {
		//convert to monomial sequence
		BigDecimal[] monomCoeffs = convertJacobiToMonomialSequence(alpha, beta, coeffs, mc);
		
		int n = coeffs.length;
		
		assert (monomCoeffs[0].compareTo(BigDecimal.ZERO) != 0);	//coefficient of the constant term in the formal series being inverted can't be 0!
		
		//invert the sequence by a simple dynamic programming
		BigDecimal[] invCoeffs = new BigDecimal[n];
		invCoeffs[0] = BigDecimal.ONE.divide(monomCoeffs[0], mc);
		for (int i = 1; i < n; i++)	{
			BigDecimal sum = BigDecimal.ZERO;
			for (int j = 0; j < i; j++)	{
				sum = sum.subtract(invCoeffs[j].multiply(monomCoeffs[i - j], mc), mc);
			}
			invCoeffs[i] = sum.divide(monomCoeffs[0], mc);
		}
		
		//convert back to Jacobi sequence
		BigDecimal[] output = convertMonomialToJacobiSequence(alpha, beta, invCoeffs, mc);
		return output;
	}
	
	
	//Inverts a formal series containing the Jacobi polynomials.
	//Specifically, takes p(x) = \sum_{j=0}^{n-1} coeffs[j]*R_j(x) and finds a q(x) = \sum_{j=0}^{n-1} output[j]*R_j(x) 
	//such that p(x)*q(x) = 1 + O(x^n) (i.e. the first non-zero monomial after the constant term 1 is of degree n)
	public static BigDecimal[] invertJacobiSequence (BigDecimal alpha, BigDecimal beta, BigDecimal[] coeffs, MathContext mc) {
		
		int n = coeffs.length;
		
		//we will eventually solve Ax = b
		BigDecimal[][] A = new BigDecimal[n][n];
		BigDecimal[] b = new BigDecimal[n];
		
		// evaluate at a grid of n points and solve the linear system
		// the actual points don' matter, as long as they are in the interval [0,1]
		// n + 1 for now, cause we don't want guys on the boundary yet
		BigDecimal gridStep = BigDecimal.ONE.divide(new BigDecimal(""+(n+1)), mc);
		// iterate over rows
		for (int rowIdx = 0; rowIdx < n; rowIdx++){
			BigDecimal gridPoint = gridStep.multiply(new BigDecimal(""+rowIdx+1), mc);
			
			// and fill the system
			b[rowIdx] = BigDecimal.ZERO;
			// iterate over columns (or the jacobi series expansion)
			for (int colIdx = 0; colIdx < n; colIdx++)	{
				// but value into A
				A[rowIdx][colIdx] = evaluate(alpha, beta, colIdx, gridPoint, mc);

				// sum in b
				b[rowIdx] = b[rowIdx].add(coeffs[colIdx].multiply(A[rowIdx][colIdx], mc), mc);
			}
			// and invert b to get the right thing
			b[rowIdx] = BigDecimal.ONE.divide(b[rowIdx], mc);
		}
		
		// solve it (result is stored in b)
		LinearSolver.geppLinearSolve(A, b, mc);
		
		return b;
	}
	
	
	// evaluate \sum_{j=0}^{n-1} coeffs[j]*R_j(x)
	public static BigDecimal evaluateLinearCombination (BigDecimal alpha, BigDecimal beta, BigDecimal[] coeffs, BigDecimal x, MathContext mc) {
		BigDecimal ret = BigDecimal.ZERO;
		int n = coeffs.length;
		for (int i = 0; i < n; i++)	{
			ret = ret.add(coeffs[i].multiply(evaluate(alpha, beta, i, x, mc), mc), mc);
		}
		return ret;
	}
	
	// the recursion
	private static BigDecimal jacobiRecursionHelper(BigDecimal alpha, BigDecimal beta, int n, BigDecimal x, BigDecimal P_n_1, BigDecimal P_n_2, MathContext mc) {
		assert (n>1);
		// return (x * P_n_1 - MyJacobi::G(alpha,beta,n-1,n-1) * P_n_1 - MyJacobi::G(alpha,beta,n-1,n-2) * P_n_2) / MyJacobi::G(alpha,beta,n-1,n);
		return x.multiply(P_n_1, mc).subtract(G(alpha,beta,n-1,n-1, mc).multiply(P_n_1, mc), mc).subtract(G(alpha,beta,n-1,n-2,mc).multiply(P_n_2, mc), mc).divide(G(alpha,beta,n-1,n,mc), mc);
	}

	// the recursion for the derivative
	private static BigDecimal jacobiDerivativeRecursionHelper(BigDecimal alpha, BigDecimal beta, int n, BigDecimal x, BigDecimal Pn1, BigDecimal dPn1, BigDecimal dPn2, MathContext mc) {
		assert (n>1);
		// return (x * dPn1 + Pn1 - MyJacobi::G(alpha,beta,n-1,n-1) * dPn1 - MyJacobi::G(alpha,beta,n-1,n-2) * dPn2) / MyJacobi::G(alpha,beta,n-1,n);
		return x.multiply(dPn1, mc).add(Pn1, mc).subtract(G(alpha,beta,n-1,n-1, mc).multiply(dPn1, mc), mc).subtract(G(alpha,beta,n-1,n-2,mc).multiply(dPn2, mc), mc).divide(G(alpha,beta,n-1,n,mc), mc);
	}
	
	
	public static void printCacheContents (HashMap<Pair<BigDecimal, BigDecimal>, TreeMap<RecurrenceEntry, BigDecimal> > cache) {
		int cnt = 0;
		for (Pair<BigDecimal, BigDecimal> pr: cache.keySet()) {
			for (RecurrenceEntry ie: cache.get(pr).keySet())	{
				System.out.println(pr.first() + "\t" + pr.second() + "\t" + ie.startIdx + "\t" + ie.endIdx + "\t" + ie.nrSteps + "\t" + cache.get(pr).get(ie));
				cnt++;
			}
		}
		System.out.println("# Entries = " + cnt);
	}

	private static <MapKeyType, TreeKeyType, TreeValueType> TreeMap<TreeKeyType, TreeValueType > findTreeMap (Map<MapKeyType, TreeMap<TreeKeyType, TreeValueType >> cache, MapKeyType key) {
		for (Entry<MapKeyType, TreeMap<TreeKeyType, TreeValueType>> currEntry : cache.entrySet()) {
			if (currEntry.getKey().equals(key))	{
				return currEntry.getValue();
			}
		}
		
		TreeMap<TreeKeyType, TreeValueType > treeMap = new TreeMap<TreeKeyType, TreeValueType >();
		cache.put(key, treeMap);
		return treeMap;
	}

	private static <MapKeyType, HashKeyType, HashValueType> HashMap<HashKeyType, HashValueType> findHashMap (Map<MapKeyType, HashMap<HashKeyType, HashValueType >> cache, MapKeyType key) {
		for (Entry<MapKeyType, HashMap<HashKeyType, HashValueType>> currEntry : cache.entrySet()) {
			if (currEntry.getKey().equals(key))	{
				return currEntry.getValue();
			}
		}
		
		HashMap<HashKeyType, HashValueType > hashMap = new HashMap<HashKeyType, HashValueType >();
		cache.put(key, hashMap);
		return hashMap;
	}

	// test it
	public static void main(String[] args) {
		
		// some parameters
		int numPoly = 20;
//		int maxExp = 10;
		int precision = 50;
		int scale = 50;
		
		// set precision and some math context
		MathContext mc = new MathContext(precision, RoundingMode.HALF_EVEN);

		BigDecimal alpha = new BigDecimal("1.1");
//		BigDecimal alpha = BigDecimal.ZERO;
		alpha = alpha.setScale(scale);
		BigDecimal beta = new BigDecimal("0.3");
//		BigDecimal beta = BigDecimal.ZERO;
		beta = beta.setScale(scale);
		
		// make a grid
		int resolution = 50;
		BigDecimal bigResolution = new BigDecimal(resolution); 
		BigDecimal[] theGrid = new BigDecimal[resolution + 1];
		// beginning
		theGrid[0] = BigDecimal.ZERO;
		// the middle
		for (int i=1; i< resolution; i++) {
			theGrid[i] = new BigDecimal(i).divide(bigResolution, mc);
		}
		// end
		theGrid[resolution] = BigDecimal.ONE;
		
//		for (BigDecimal value : theGrid) {
//			System.out.println(value);
//		}
//		System.out.println();

		
//		// what sayeth the poly ab length
//		for (int n=0; n<=numPoly; n++) {
//			for (BigDecimal value : theGrid) {
//				System.out.print (JacobiPolynomials.evaluate(alpha, beta, n, value, mc) + "\t");
//			}
//			System.out.println();
//		}
		
//		// what sayeth the poly ab length
//		for (int n=0; n<=numPoly; n++) {
//			System.out.println (JacobiPolynomials.evaluateSquaredLength(alpha, beta, n, mc) + "\t");
//		}

//		for (int n=0; n<=numPoly; n++) {
////			System.out.println(n);
//			for (BigDecimal value : theGrid) {
////				System.out.println(value);
//				BigDecimal jac = JacobiPolynomials.evaluate(alpha, beta, n, value, mc);
//				for (int exp = 0; exp <= 5; exp++) {
////					System.out.println("this is k now " + k);
//					// left hand sides
//					BigDecimal lhsY = value.pow(exp, mc).multiply(jac, mc);
//					BigDecimal lhs1mY = BigDecimal.ONE.subtract(value, mc).pow(exp, mc).multiply(jac, mc);
////					BigDecimal lhsCombined = lhsY.multiply(lhs1mY, mc);
//					
//					// right hand sides
//					BigDecimal rhsY = BigDecimal.ZERO;
//					BigDecimal rhs1mY = BigDecimal.ZERO;
//					for (int l = Math.max(0, n-exp); l <= n+exp; l++){
////						System.out.println("this is l now " + l);
//						rhsY = rhsY.add(JacobiPolynomials.L (alpha, beta, n, l, exp, mc).multiply(JacobiPolynomials.evaluate(alpha, beta, l, value, mc), mc), mc);
//						rhs1mY = rhs1mY.add(JacobiPolynomials.M (alpha, beta, n, l, exp, mc).multiply(JacobiPolynomials.evaluate(alpha, beta, l, value, mc), mc), mc);
//					}
//					
//					// print something
//					System.out.println ("n = " + n + "\tx = " + value + "\tk = " + exp);
//					System.out.println (lhsY +"\tvs.\t" + rhsY + "\t[" + lhsY.subtract(rhsY,mc) + "]");
//					System.out.println (lhs1mY +"\tvs.\t" + rhs1mY + "\t[" + lhs1mY.subtract(rhs1mY,mc) + "]");
//					
//				}
//			}				
//		}


		// first poly degree
		for (int j=0; j<=numPoly; j++) {
			// go through grid
			for (BigDecimal value : theGrid) {
				BigDecimal jac = JacobiPolynomials.evaluate(alpha, beta, j, value, mc);
				System.out.print (jac + "\t");
				
				
//				// go through first exponenet
//				for (int expY = 0; expY <= maxExp; expY++) {
//					// go thorugh second exponent
//					for (int exp1mY = 0; exp1mY <= maxExp; exp1mY++) {
//						// left hand side
//						BigDecimal Yt1mY = value.pow(expY, mc).multiply(BigDecimal.ONE.subtract(value, mc).pow(exp1mY, mc), mc);
//						BigDecimal lhs = Yt1mY.multiply(jac, mc);
//						
//						// right hand side
//						BigDecimal rhs1 = BigDecimal.ZERO;
//						BigDecimal rhs2 = BigDecimal.ZERO;
//						// first: how long can it stretch in total
//						for (int l=Math.max (0, j - (expY + exp1mY)); l<= j + (expY + exp1mY); l++) {
//							// calculate the factor for the l-th jacobi polynomial
//							BigDecimal factor = BigDecimal.ZERO;
//							// what are the possible u=intermediate one
//							for (int mu = Math.max(Math.max(0, j-expY), l - exp1mY); mu <= Math.min(j+expY, l+exp1mY); mu++){
//								factor = factor.add (JacobiPolynomials.L (alpha, beta, j, mu, expY, mc).multiply(JacobiPolynomials.M (alpha, beta, mu, l, exp1mY, mc), mc), mc);
//							}
//							// add it to the right hand side
//							rhs1 = rhs1.add(factor.multiply(JacobiPolynomials.evaluate(alpha, beta, l, value, mc), mc), mc);
//							
//							// get the new factor
//							factor = SelectionHMM.N(j, l, expY + exp1mY, expY, alpha, beta, mc).divide(JacobiPolynomials.evaluateSquaredLength (alpha, beta, l, mc), mc);
//							rhs2 = rhs2.add (factor.multiply(JacobiPolynomials.evaluate(alpha, beta, l, value, mc), mc), mc);						
//						}
//							
//						// print something
////						System.out.println ("j = " + j + "\tx = " + value + "\ty^ = " + expY + "\t(1-y)^" + exp1mY);
////						System.out.println (lhs +"\tvs.\t" + rhs1 + "\tvs.\t" + rhs2);
//						System.out.println ("[" + lhs.subtract(rhs1,mc) + "]\t["+ lhs.subtract(rhs2,mc) + "]");
				
			}
			System.out.println();
			
			for (BigDecimal value : theGrid) {
				BigDecimal jac = JacobiPolynomials.evaluateDerivative (alpha, beta, j, value, mc);
				System.out.print (jac + "\t");
			}
			System.out.println();
		}
	}

	// debug
	public static int numberOfLMatrixPowers () {
		for (MatrixMetaArguments key : cacheMatrixL.keySet()) {
			return cacheMatrixL.get(key).size();
		}
		return 0;
	}
	
	public static void clearCaches () {
		cacheL.clear();
		cacheM.clear();
		cacheG.clear();
		cachePoly.clear();
		cachePolyDerivative.clear();
		cacheSL.clear();
		cacheMatrixL.clear();
		cacheMatrixM.clear();
	}
	// containers and stuff
	
	final private static BigDecimal HALF = new BigDecimal("0.5");
	final private static BigDecimal TWO = new BigDecimal("2");	
	
	private static HashMap<PolynomialMetaArguments, HashMap<RecurrenceEntry, BigDecimal> > cacheL = new HashMap<PolynomialMetaArguments, HashMap<RecurrenceEntry, BigDecimal> >(); 
	private static HashMap<PolynomialMetaArguments, HashMap<RecurrenceEntry, BigDecimal> > cacheM = new HashMap<PolynomialMetaArguments, HashMap<RecurrenceEntry, BigDecimal> >(); 
	private static HashMap<PolynomialMetaArguments, HashMap<Pair<Integer, Integer>, BigDecimal> > cacheG = new HashMap<PolynomialMetaArguments, HashMap<Pair<Integer, Integer>, BigDecimal> >();
	// have to use TreeMap for the inner map from (n,x) -> evaluation, because we don't want to use == on BigDecimal (which would be the case if we used HashMap)
	private static HashMap<PolynomialMetaArguments, TreeMap<Pair<Integer, BigDecimal>, BigDecimal> > cachePoly = new HashMap<PolynomialMetaArguments, TreeMap<Pair<Integer, BigDecimal>, BigDecimal> >();
	private static HashMap<PolynomialMetaArguments, TreeMap<Pair<Integer, BigDecimal>, BigDecimal> > cachePolyDerivative = new HashMap<PolynomialMetaArguments, TreeMap<Pair<Integer, BigDecimal>, BigDecimal> >();
	private static HashMap<PolynomialMetaArguments, HashMap<Integer, BigDecimal> > cacheSL = new HashMap<PolynomialMetaArguments, HashMap<Integer, BigDecimal> >();

	private static HashMap<MatrixMetaArguments, HashMap<Integer, BigDecimal[][]> > cacheMatrixL = new HashMap<MatrixMetaArguments, HashMap<Integer, BigDecimal[][]> >(); 
	private static HashMap<MatrixMetaArguments, HashMap<Integer, BigDecimal[][]> > cacheMatrixM = new HashMap<MatrixMetaArguments, HashMap<Integer, BigDecimal[][]> >();

	
	private static class RecurrenceEntry implements Comparable<RecurrenceEntry> {
		int startIdx, endIdx, nrSteps; 
		
		public RecurrenceEntry(int startIdx, int endIdx, int nrSteps) {
			this.startIdx = startIdx;
			this.endIdx = endIdx;
			this.nrSteps = nrSteps;
		}

		public boolean equals(Object o) {
			assert (o instanceof RecurrenceEntry);
			RecurrenceEntry newO = ((RecurrenceEntry) o);
			return (startIdx == newO.startIdx) && (endIdx == newO.endIdx) && (nrSteps == newO.nrSteps);
		}
		
		public int compareTo(RecurrenceEntry o){
			if (o == null)	return 1;
			int r = startIdx - o.startIdx;
			if (r != 0)	return r;
			r = endIdx - o.endIdx;
			if (r != 0)	return r;
			r = nrSteps - o.nrSteps;
			return r;
		}
		
		public int hashCode()	{
			return (((startIdx * 0x1f1f1f1f) ^ endIdx) * 0x1f1f1f1f) ^ nrSteps; 
		}
	}

	// look up class for polynomials
	private static class PolynomialMetaArguments {
		// variables for look up
		BigDecimal alpha, beta;
		MathContext mc; 
		
		public PolynomialMetaArguments(BigDecimal alpha, BigDecimal beta, MathContext mc) {
			this.alpha = alpha;
			this.beta = beta;
			this.mc = mc;
		}

		public boolean equals (Object o) {
			assert (o instanceof PolynomialMetaArguments);
			PolynomialMetaArguments newO = ((PolynomialMetaArguments) o);
			return (alpha.compareTo(newO.alpha) == 0) && (beta.compareTo(newO.beta) == 0) && (mc == newO.mc);
		}
	}

	// look up class for polynomials
	private static class MatrixMetaArguments {
		// variables for look up
		BigDecimal alpha, beta;
		MathContext mc;
		int maxM;
		
		public MatrixMetaArguments (BigDecimal alpha, BigDecimal beta, int maxM, MathContext mc) {
			this.alpha = alpha;
			this.beta = beta;
			this.maxM = maxM;
			this.mc = mc;
		}

		public boolean equals (Object o) {
			assert (o instanceof MatrixMetaArguments);
			MatrixMetaArguments newO = ((MatrixMetaArguments) o);
			return (alpha.compareTo(newO.alpha) == 0) && (beta.compareTo(newO.beta) == 0) && (maxM == newO.maxM) && (mc == newO.mc);
		}
	}

	
}
