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
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.nevec.rjm.BigDecimalMath;

import edu.berkeley.spectralHMM.datatypes.Pair;
import edu.berkeley.spectralHMM.matrix.EigenSystem;
import edu.berkeley.spectralHMM.matrix.TDFEigenSolverBanded;
import edu.berkeley.spectralHMM.matrix.TDFEigenSolverBanded.EigenRefineError;
import edu.berkeley.spectralHMM.matrix.TDFEigenSolverBanded.EigenRefineException;

public class SelectionEigenfunktion {
	
	private SelectionEigenfunktion () {
		;
	}
	
	public static BigDecimal evaluate (BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int maxM, int n, BigDecimal y, MathContext mc) throws EigenRefineException, EigenRefineError {
		int offset = JacobiPolynomials.calculateOffset (alpha, beta);
		// check sizes (+1, cause that's our notation convention)
		assert (maxM < matrixCutoff);
		assert (n < matrixCutoff);

		EigenfunctionMetaArguments metaArgs = new EigenfunctionMetaArguments(alpha, beta, hetF, homF, matrixCutoff, maxM, mc);
		TreeMap<Pair<Integer, BigDecimal>, BigDecimal> t = findTreeMap (cacheEigenfunctions, metaArgs);
		Pair<Integer, BigDecimal> cacheEntry = new Pair<Integer, BigDecimal>(n, y);
		BigDecimal eigenFunction = t.get(cacheEntry);
		if (eigenFunction != null)	{
			return eigenFunction;
		}
		
		// if eigenfunction not found in the cache, check if the eigensystem for the w coefficients was solved earlier
		EigenSystem eigenSystem = findEigenSystem (metaArgs);
		
		assert (matrixCutoff == eigenSystem.eigenVectors.length);
		
		// eigensystem exists, now compute the eigenfunction using it
		eigenFunction = BigDecimal.ZERO.setScale(mc.getPrecision());
		
		// TODO memoize that
		//using heterozygous and homozygous fitnesses
		BigDecimal expHalfWeight = BigDecimalMath.exp (meanFitness(hetF, homF, y, mc).divide (mTWO, mc), mc);
		
		// and do the sum over the coefficients
		for (int m=0; m <= maxM; m++) {
			eigenFunction = eigenFunction.add (eigenSystem.eigenVectors[n][m].multiply(JacobiPolynomials.evaluate(alpha, beta, m + offset, y, mc), mc), mc);
		}
		// and multiply the the half exponential
		eigenFunction = eigenFunction.multiply(expHalfWeight, mc);
		
		t.put(new Pair<Integer, BigDecimal>(n, y), eigenFunction);
		return eigenFunction;
	}
	
	
	static public BigDecimal evaluateSquaredLength (BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int maxM, int n, MathContext mc) throws EigenRefineException, EigenRefineError {
		int offset = JacobiPolynomials.calculateOffset (alpha, beta);
		// check sizes (+1, cause that's our notation convention)
		assert (maxM < matrixCutoff);
		assert (n < matrixCutoff);

		EigenfunctionMetaArguments metaArgs = new EigenfunctionMetaArguments(alpha, beta, hetF, homF, matrixCutoff, maxM, mc);
		HashMap<Integer, BigDecimal> t = findHashMap (cacheSquaredLength, metaArgs);
		BigDecimal squaredEigenLength = t.get(n);
		if (squaredEigenLength != null)	{
			return squaredEigenLength ;
		}
		
		EigenSystem eigenSystem = findEigenSystem(metaArgs);
		
		assert(matrixCutoff == eigenSystem.eigenVectors.length);
		
		squaredEigenLength = BigDecimal.ZERO;
		// TODO Should this be summing till maxM?
		for (int m=0; m <= maxM; m++) {
			squaredEigenLength = squaredEigenLength.add(eigenSystem.eigenVectors[n][m].multiply(eigenSystem.eigenVectors[n][m], mc).multiply(JacobiPolynomials.evaluateSquaredLength(alpha, beta, m + offset, mc), mc), mc);
		}

		t.put(n, squaredEigenLength);
		return squaredEigenLength;
	}
	
	
	// the final matrix
	public static BigDecimal[][] buildMatrix (BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int systemCutoff, MathContext mc) {
		
		int offset = JacobiPolynomials.calculateOffset (alpha, beta);
		// what's the offset
		assert ((offset >= 0) && (offset <= 2));
		
		// build up the matrix	
		BigDecimal[][] H = new BigDecimal[systemCutoff][systemCutoff];
		
		// fill her up
		for (int i=0; i<H.length; i++) {
			Arrays.fill(H[i], BigDecimal.ZERO);
		}
		// be aware of offset
		// main index
		for (int m=0; m<H.length; m++) {
			// diagonal and off-diagonal elements
			for (int d=-4; d<= 4; d++) {
				// this index stuff should yield the right thing
				if ((m+d >= 0) && (m+d < systemCutoff)) {
					H[m+d][m] = b (d, m + offset, alpha, beta, hetF, homF, offset, mc);
				}
			}
			
			// the diagonal needs the eigenvalue
			H[m][m] = H[m][m].add (JacobiPolynomials.lambda (m + offset, alpha, beta, mc), mc);
		}

		// give it away
		return H;
	}

	// the b^(i)_m
	public static BigDecimal b (int i, int m, BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int offset, MathContext mc) {
		// we need a return container
		BigDecimal b = BigDecimal.ZERO;
		
		// now do want we need to do
		// make the sum
		for (int l = Math.abs(i); l <= 4; l++) {
			b = b.add(q(l, alpha, beta, hetF, homF, mc).multiply(JacobiPolynomials.L(alpha, beta, m, m+i, l, offset, mc), mc), mc);
		}
		
		//give it away
		return b;
	}

	/// these are the coefficients in the Q-polynomial 
	public static BigDecimal q (int l, BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, MathContext mc) {
		// safety first
		assert (l >= 0);
		
		// some locals
		BigDecimal t1, t2,t3;
		
		// then give a q
		BigDecimal q = BigDecimal.ZERO;
		switch (l) {
			case 0:
				q = alpha.multiply(hetF, mc);
				// q0 = alpha*hetF;
				break;
			case 1:
				t1 = BigDecimal.ONE.add(alpha, mc).multiply(homF, mc);
				t2 = (new BigDecimal("3")).multiply(alpha, mc).add(beta, mc).add(TWO, mc);
				t3 = TWO.multiply(hetF, mc).subtract(t2, mc);
				q = hetF.multiply(t3, mc).add(t1, mc);
				// q1 = hetF*(2*hetF - (3*alpha + beta + 2)) + (1 + alpha)*homF;
				break;
			case 2:
				t1 = BigDecimal.ONE.add(alpha, mc).add(beta, mc);
				t2 = t1.add(TWO.multiply(homF, mc), mc).subtract((new BigDecimal("5")).multiply(hetF, mc), mc);
				q = TWO.multiply(hetF, mc).multiply(t2, mc).subtract(homF.multiply(t1, mc), mc);
				// q2 = 2*hetF*(1 + alpha + beta + 2*homF - 5*hetF) - homF*(1 + alpha + beta);
				break;
			case 3:
				t1 = (new BigDecimal("16")).multiply(hetF, mc).multiply(hetF, mc);
				t2 = (new BigDecimal("12")).multiply(homF, mc).multiply(hetF, mc);
				t3 = TWO.multiply(homF, mc).multiply(homF, mc);
				q = t1.subtract(t2, mc).add(t3, mc);
				// q3 = 16*hetF*hetF - 12*homF*hetF + 2*homF*homF
				break;
			case 4:
				t1 = homF.subtract(TWO.multiply(hetF, mc), mc);
				q = (new BigDecimal("-2")).multiply(t1, mc).multiply(t1, mc);
				// q4 = -2*(homF - 2*hetF)*(homF - 2*hetF);
				break;
			default:
				// this should not happen
				assert (false);
		}
		// give it away now
		return q;
	}

	
	// the mean fitness
	public static BigDecimal meanFitness (BigDecimal hetF, BigDecimal homF, BigDecimal x, MathContext mc) {
		// return 2*homF*x*x + 4*hetF*x*(1-x);
		BigDecimal zeroOne = (new BigDecimal("4")).multiply(hetF, mc).multiply(x, mc).multiply(BigDecimal.ONE.subtract(x, mc), mc);
		BigDecimal oneOne = TWO.multiply(homF, mc).multiply(x, mc).multiply(x, mc);
		return oneOne.add(zeroOne, mc);
	}

	// the derivative of the mean fitness
	public static BigDecimal meanFitnessDerivative (BigDecimal hetF, BigDecimal homF, BigDecimal x, MathContext mc) {
		// function is 2*homF*x*x + 4*hetF*x*(1-x);
		// so the derivative is 4*homF*x + 4*hetF*(1-2*x)
		BigDecimal zeroOne = (new BigDecimal("4")).multiply(hetF, mc).multiply(BigDecimal.ONE.subtract(TWO.multiply(x, mc), mc), mc);
		BigDecimal oneOne = (new BigDecimal("4")).multiply(homF, mc).multiply(x, mc);
		return oneOne.add(zeroOne, mc);
	}
	
	// the constant C
	public static BigDecimal getC (BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int maxM, MathContext mc) throws EigenRefineException, EigenRefineError {
		
		// get the eigensystem
		// TODO we don't need maxM for this eigensystem actually
		EigenfunctionMetaArguments metaArgs = new EigenfunctionMetaArguments (alpha, beta, hetF, homF, matrixCutoff, maxM, mc);
		EigenSystem theSystem = findEigenSystem (metaArgs);
		
		// initialize return value
		BigDecimal result = BigDecimal.ZERO;
		// initialize the running value of the "multinomial"
		BigDecimal multi = BigDecimal.ONE;
		
		// run
		for (int m=0; m<=maxM; m++) {
			// add the current stuff

			// calcutate the summand
			BigDecimal summand = theSystem.eigenVectors[0][m].multiply(multi, mc);
			
			// even or odd
			if (m % 2 == 0) {
				// even, add it
				result = result.add(summand, mc);
			}
			else {
				// even, subtract it
				result = result.subtract(summand, mc);
			}
			
			// update the multinomial-factor
			BigDecimal bigM = new BigDecimal(m);
			multi = multi.multiply(bigM.add(alpha, mc), mc).divide(bigM.add(BigDecimal.ONE, mc), mc);
		}
		
		// give it away now
		return result;
	}
	
	
	
	private static <MapKeyType, HashKeyType, HashValueType> HashMap<HashKeyType, HashValueType > findHashMap (Map<MapKeyType, HashMap<HashKeyType, HashValueType >> cache, MapKeyType key) {
		for (Entry<MapKeyType, HashMap<HashKeyType, HashValueType > > currEntry: cache.entrySet())	{
			if (currEntry.getKey().equals(key))	{
				return currEntry.getValue();
			}
		}
	
		HashMap<HashKeyType, HashValueType > hashMap = new HashMap<HashKeyType, HashValueType >();
		cache.put(key, hashMap);
		return hashMap;
	}
	
	private static <MapKeyType, TreeKeyType, TreeValueType> TreeMap<TreeKeyType, TreeValueType > findTreeMap (Map<MapKeyType, TreeMap<TreeKeyType, TreeValueType >> cache, MapKeyType key) {
		for (Entry<MapKeyType, TreeMap<TreeKeyType, TreeValueType > > currEntry: cache.entrySet())	{
			if (currEntry.getKey().equals(key))	{
				return currEntry.getValue();
			}
		}
	
		TreeMap<TreeKeyType, TreeValueType > treeMap = new TreeMap<TreeKeyType, TreeValueType >();
		cache.put(key, treeMap);
		return treeMap;
	}
	
	private static EigenSystem findEigenSystem(EigenfunctionMetaArguments metaArgs) throws EigenRefineException, EigenRefineError{
		for (Entry<EigenfunctionMetaArguments, EigenSystem> mp: cacheEigensystem.entrySet())	{
			if (mp.getKey().equals(metaArgs))	{
				return mp.getValue();
			}
		}

		// eigensystem not already solved, so do it now
		EigenSystem eigenSystem;

		// first let us create a coefficient matrix
		// for solving the system appropriately later (especially making 0 an eigenvalue), we need some extra precision
		int extraPrecision = Math.max(0, (int)(3. * Math.log10(1. * metaArgs.matrixCutoff / 100.) + 1)) + 5;
		System.out.println("# Old precision = " + metaArgs.mc.getPrecision());
		System.out.println("# extraPrecision for M matrix = " + extraPrecision);
		// now build the actual matrix
		BigDecimal[][] M = buildMatrix (metaArgs.alpha, metaArgs.beta, metaArgs.hetF, metaArgs.homF, metaArgs.matrixCutoff, new MathContext (metaArgs.mc.getPrecision() + extraPrecision, metaArgs.mc.getRoundingMode()));
		
		// get some eigensystem
		int bw = 4;
		//parameterized by hetF and homF
		if (metaArgs.hetF.compareTo(BigDecimal.ZERO) == 0 && metaArgs.homF.compareTo(BigDecimal.ZERO) == 0)	{	// sigma and sigma * h are both 0
			bw = 0;
		}
		else if (new BigDecimal("2").multiply(metaArgs.hetF).compareTo(metaArgs.homF) == 0)	{	// h = 1/2
			bw = 2;
		}
		
		System.out.printf("# Using bw = %d for coefficient matrix\n", bw);
		
		// we need to know the weighting for the inner product in the hilbert space
		System.out.println("# Get weights of the Jacobi polynomials.");
		int offset = JacobiPolynomials.calculateOffset (metaArgs.alpha, metaArgs.beta);
		BigDecimal jacobiWeights[] = new BigDecimal[M.length];
		for (int k=0; k<jacobiWeights.length; k++) {
			jacobiWeights[k] = JacobiPolynomials.evaluateSquaredLength (metaArgs.alpha, metaArgs.beta, k+offset, metaArgs.mc);
		}
		// solve the eigensystem (with extra precision for the refinement)
		eigenSystem = TDFEigenSolverBanded.solve (M, jacobiWeights, bw, metaArgs.mc.getPrecision(), extraPrecision, offset == 0);		
		
		// put it into cache
		cacheEigensystem.put(metaArgs, eigenSystem);
		return eigenSystem;
	}
	
	public static BigDecimal getEigenValue(BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int maxM, int n, MathContext mc) throws EigenRefineException, EigenRefineError {
		EigenfunctionMetaArguments metaArgs = new EigenfunctionMetaArguments(alpha, beta, hetF, homF, matrixCutoff, maxM, mc);
		EigenSystem eigenSystem = findEigenSystem(metaArgs);

		return eigenSystem.eigenValues[n];
	}

	public static BigDecimal getEigenVectorEntry(BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int maxM, int n, int m, MathContext mc) throws EigenRefineException, EigenRefineError {
		EigenfunctionMetaArguments metaArgs = new EigenfunctionMetaArguments(alpha, beta, hetF, homF, matrixCutoff, maxM, mc);
		EigenSystem eigenSystem = findEigenSystem(metaArgs);
		
		return eigenSystem.eigenVectors[n][m];
	}
	
	public static BigDecimal[] getEigenVectorCopy(BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int maxM, int n, MathContext mc) throws EigenRefineException, EigenRefineError {
		EigenfunctionMetaArguments metaArgs = new EigenfunctionMetaArguments(alpha, beta, hetF, homF, matrixCutoff, maxM, mc);
		EigenSystem eigenSystem = findEigenSystem(metaArgs);
		
		return Arrays.copyOf (eigenSystem.eigenVectors[n], eigenSystem.eigenVectors[n].length);
	}

	public static BigDecimal[][] getEigenVectorMatrixCopy(BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int numRows, int numCols, MathContext mc) throws EigenRefineException, EigenRefineError {
		// maxM is numCols-1
		EigenfunctionMetaArguments metaArgs = new EigenfunctionMetaArguments(alpha, beta, hetF, homF, matrixCutoff, numCols-1, mc);
		EigenSystem eigenSystem = findEigenSystem(metaArgs);
		
		BigDecimal[][] result = new BigDecimal[numRows][numCols];
		for (int i=0; i<numRows; i++) {
			for (int j=0; j<numCols; j++)	{
				result[i][j] = eigenSystem.eigenVectors[i][j];
			}
		}
		
		return result;
	}

	public static void main (String[] args) {
		try {
			// some parameters
			int precision = 80;
			int scale = precision;
			// set precision and some math context
			MathContext mc = new MathContext(precision, RoundingMode.HALF_EVEN);
	
			// get population genetic parameters
			BigDecimal alpha = new BigDecimal("0.02");
			BigDecimal beta = new BigDecimal("0.02");
			BigDecimal sigma = new BigDecimal("2");
			BigDecimal h = new BigDecimal("0.5");
			BigDecimal hetF = new BigDecimal("1");
			BigDecimal homF = new BigDecimal("2");
			h = h.setScale(scale);
			System.out.println("# sigma: " + sigma + "\t h: " + h);
			
			
	
			
			
			
	//		# matrixCutoff = 500
	//		# maxM = 495
	//		# maxN = 490
	
			// set some method parameters
			int matrixCutoff = 150;
			int maxM = 140;
			
			int maxN = 2;
			
			// make a grid
			int resolution = 40;
			BigDecimal bigResolution = new BigDecimal(resolution); 
			BigDecimal[] theGrid = new BigDecimal[resolution + 1];
			// beginning
			theGrid[0] = BigDecimal.ZERO.setScale(scale);
			// the middle
			for (int i=1; i< resolution; i++) {
				theGrid[i] = (new BigDecimal(i).divide(bigResolution, mc)).setScale(scale);
			}
			// end
			theGrid[resolution] = BigDecimal.ONE.setScale(scale);
			
			// print eigenfunctions
			
			for (int n = 0; n<=maxN; n++) {
				for (int y=0; y<theGrid.length; y++) {
					BigDecimal value = SelectionEigenfunktion.evaluate (alpha, beta, hetF, homF, matrixCutoff, maxM, n, theGrid[y], mc);
					System.out.print (value + "\t");
				}
				System.out.println();
				// some derivative
				for (int y=0; y<theGrid.length; y++) {
					BigDecimal value = SelectionEigenfunktion.evaluateDerivative (alpha, beta, hetF, homF, matrixCutoff, maxM, n, theGrid[y], mc);
					System.out.print (value + "\t");
				}
				System.out.println();
			}

	
	//		BigDecimal C = BigDecimalSelectionEigenfunktion.getC (alpha, beta, sigma, h, matrixCutoff, maxM, mc);
	//		System.out.println(C);
	//		
	//		// get second eigenvector
	//		BigDecimal[] blub = BigDecimalSelectionEigenfunktion.getEigenVectorCopy	(alpha, beta, sigma, h, matrixCutoff, maxM, 2, mc);
	//		for (BigDecimal entry : blub) {
	//			System.out.print (entry + "\t");
	//		}
	//		System.out.println();
	//		
	//		// second eigenvalue
	//		System.out.println (BigDecimalSelectionEigenfunktion.getEigenValue (alpha, beta, sigma, h, matrixCutoff, maxM, 2, mc));
	//		
	//		
	//		// get the coefficients for the 0-th eigenfunction
	//		BigDecimal[] w0 = BigDecimalSelectionEigenfunktion.getEigenVectorCopy(alpha, beta, sigma, h, matrixCutoff, maxM, 0, mc);
	//		System.out.print ("# orig: ");
	//		for (BigDecimal value : w0) {
	//			System.out.print (value.divide(C, mc).doubleValue() + "\t");
	//		}
	//		System.out.println();
	//
	//		// print exp(\bar\sigma/2)
	//		for (int y=0; y<theGrid.length; y++) {
	//			System.out.print (BigDecimalMath.exp (meanFitness(sigma, h, theGrid[y], mc).divide(new BigDecimal("2"), mc)) + "\t");
	//		}
	//		System.out.println();
	//		
	//		
	//		// get the jacobi polynomial approximation to it
	//		for (int y=0; y<theGrid.length; y++) {
	//			// calculate for this y
	//			BigDecimal value = BigDecimal.ZERO;
	//			for (int m=0; m<=maxM; m++) {
	//				value = value.add(w0[m].multiply(BigDecimalJacobiPolynomials.evaluate(alpha, beta, m, theGrid[y], mc), mc), mc);
	//			}
	//			value = value.divide(C, mc);
	//			System.out.print (value + "\t");
	//		}
	//		System.out.println();
	//
	//		
	//		// print exp( - \bar\sigma/2)
	//		for (int y=0; y<theGrid.length; y++) {
	//			System.out.print (BigDecimalMath.exp (meanFitness(sigma, h, theGrid[y], mc).divide(new BigDecimal("-2"), mc)) + "\t");
	//		}
	//		System.out.println();
	//		
	//		// get the jacobi polynomial approximation to it
	//		// truncate the w
	//		BigDecimal[] truncW = Arrays.copyOf(w0, maxM+1);
	//		// get the invers
	//		BigDecimal[] v = BigDecimalJacobiPolynomials.invertJacobiSequence(alpha, beta, truncW, mc);
	//		
	//		// print it
	//		System.out.println ("# C = "+ C);
	//		System.out.print ("# inv: ");
	//		for (BigDecimal value : v) {
	//			System.out.print (value.multiply(C, mc).doubleValue() + "\t");
	//		}
	//		System.out.println();
	//		for (int y=0; y<theGrid.length; y++) {
	//			// calculate for this y
	//			BigDecimal value = BigDecimal.ZERO;
	//			for (int m=0; m<v.length; m++) {
	//				value = value.add(v[m].multiply(BigDecimalJacobiPolynomials.evaluate(alpha, beta, m, theGrid[y], mc), mc), mc);
	//			}
	//			value = value.multiply (C, mc);
	//			System.out.print (value + "\t");
	//		}
	//		System.out.println();
		}
		catch (EigenRefineException e) {
			System.out.println ("# [EXCEPTION (maybe recoverable)] " + e.getMessage());
		} catch (EigenRefineError e) {
			System.out.println ("# [ERROR (not recoverable)] " + e.getMessage());
		}
	}

	
	public static void clearCaches () {
		cacheEigenfunctions.clear();
		cacheEigenfunctionDerivatives.clear();
		cacheEigensystem.clear();
		cacheSquaredLength.clear();
	}
	// containers
	
	private static final BigDecimal mTWO = new BigDecimal("-2");
	private static final BigDecimal TWO = new BigDecimal("2");

	// tree map, cause no good hashing for BigDecimal
	private static HashMap<EigenfunctionMetaArguments, TreeMap<Pair<Integer, BigDecimal>, BigDecimal> > cacheEigenfunctions = new HashMap<EigenfunctionMetaArguments, TreeMap<Pair<Integer, BigDecimal>, BigDecimal> >();
	private static HashMap<EigenfunctionMetaArguments, TreeMap<Pair<Integer, BigDecimal>, BigDecimal> > cacheEigenfunctionDerivatives = new HashMap<EigenfunctionMetaArguments, TreeMap<Pair<Integer, BigDecimal>, BigDecimal> >();
	private static HashMap<EigenfunctionMetaArguments, EigenSystem> cacheEigensystem = new HashMap<EigenfunctionMetaArguments, EigenSystem>();
	private static HashMap<EigenfunctionMetaArguments, HashMap<Integer, BigDecimal> > cacheSquaredLength = new HashMap<EigenfunctionMetaArguments, HashMap<Integer, BigDecimal> >();
	

	// look up class for polynomials
	private static class EigenfunctionMetaArguments {
		// variables for look up
		public BigDecimal alpha, beta;
		//hetF is heterozygous fitness
		//homF is homozygous fitness
		public BigDecimal hetF, homF;	
		public int matrixCutoff, maxM;
		
		MathContext mc; 
		
		public EigenfunctionMetaArguments(BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int maxM, MathContext mc) {			
			this.alpha = alpha;
			this.beta = beta;
			this.hetF = hetF;
			this.homF = homF;
			this.matrixCutoff = matrixCutoff;
			this.maxM = maxM;
			this.mc = mc;
		}

		public boolean equals(Object o) {
			assert (o instanceof EigenfunctionMetaArguments);
			EigenfunctionMetaArguments newO = ((EigenfunctionMetaArguments) o);
			return (alpha.compareTo(newO.alpha) == 0) && (beta.compareTo(newO.beta) == 0) && (hetF.compareTo(newO.hetF) == 0) && (homF.compareTo(newO.homF) == 0) && (matrixCutoff == newO.matrixCutoff) && (maxM == newO.maxM) && (mc == newO.mc);
		}
	}


	public static BigDecimal evaluateDerivative (BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int maxM, int n, BigDecimal y, MathContext mc) throws EigenRefineException, EigenRefineError {
		int offset = JacobiPolynomials.calculateOffset (alpha, beta);
		// check sizes (+1, cause that's our notation convention)
		assert (maxM < matrixCutoff);
		assert (n < matrixCutoff);

		EigenfunctionMetaArguments metaArgs = new EigenfunctionMetaArguments(alpha, beta, hetF, homF, matrixCutoff, maxM, mc);
		TreeMap<Pair<Integer, BigDecimal>, BigDecimal> t = findTreeMap (cacheEigenfunctionDerivatives, metaArgs);
		Pair<Integer, BigDecimal> cacheEntry = new Pair<Integer, BigDecimal>(n, y);
		BigDecimal eigenFunctionDerivative = t.get(cacheEntry);
		if (eigenFunctionDerivative != null)	{
			return eigenFunctionDerivative;
		}
		
		// if eigenfunction not found in the cache, check if the eigensystem for the w coefficients was solved earlier
		EigenSystem eigenSystem = findEigenSystem (metaArgs);
		
		assert (matrixCutoff == eigenSystem.eigenVectors.length);
		
		// gets a bit more complicated due to product rule
		
		// we do need the exponential weight
		BigDecimal expHalfWeight = BigDecimalMath.exp (meanFitness(hetF, homF, y, mc).divide (mTWO, mc), mc);

		// also some derivate factor
		BigDecimal dFactor = expHalfWeight.divide(mTWO,mc).multiply(meanFitnessDerivative(hetF, homF, y, mc), mc);
		
		// also the sum of jacobis
		BigDecimal weightedJacobiSum = BigDecimal.ZERO;
		// and the derivatives
		BigDecimal weightedJacobiDerivativeSum = BigDecimal.ZERO;

		// do the sum over the coefficients
		for (int m=0; m <= maxM; m++) {
			weightedJacobiSum = weightedJacobiSum.add (eigenSystem.eigenVectors[n][m].multiply(JacobiPolynomials.evaluate (alpha, beta, m + offset, y, mc), mc), mc);
			weightedJacobiDerivativeSum = weightedJacobiDerivativeSum.add (eigenSystem.eigenVectors[n][m].multiply(JacobiPolynomials.evaluateDerivative (alpha, beta, m + offset, y, mc), mc), mc);
		}
		
		// then combine everything
		eigenFunctionDerivative = dFactor.multiply(weightedJacobiSum, mc).add(expHalfWeight.multiply(weightedJacobiDerivativeSum, mc), mc);
				
		t.put(new Pair<Integer, BigDecimal>(n, y), eigenFunctionDerivative);
		return eigenFunctionDerivative;
	}


}
