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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

import org.netlib.lapack.DGEEV;
import org.netlib.util.intW;

public class TDFEigenSolverBanded {

	public static EigenSystem solve(BigDecimal[][] matrix, BigDecimal[] weights, int bandwidth, int precision, int extaPrecision, boolean firstEValBeZero) throws EigenRefineException, EigenRefineError {
		if (bandwidth == 0)	{
			return diagonalMatrixEigensolution(matrix);
		}
		// get some dimension
		int dim = matrix.length;
		assert (dim == weights.length);
		// also some MathContext
		MathContext mc = new MathContext(precision + 5, RoundingMode.HALF_EVEN);
		
		// now find some eigenvalues using jLAPACK (and time it)
		System.out.println("# Start jLApack.");
		long start = System.currentTimeMillis();
		EigenSystem mySystem = solveJLAPackEigs(matrix);
		long end = System.currentTimeMillis();		
		System.out.printf("# jLAPACK took %d ms\n", (end-start));
		
		// now we have some seed eigenvalues and vectors
		// lets do something

		// copy results we got from jLApack (and normalize it on the way)
		BigDecimal[] eigenValues = new BigDecimal[dim];
		BigDecimal[][] eigenVectors = new BigDecimal[dim][dim];
		for (int i = 0; i < dim; i++) {
			// copy eigenvalue
			eigenValues[i] = mySystem.eigenValues[i];
			
			// calculate norm of vec
			BigDecimal norm = getNorm (mySystem.eigenVectors[i], mc);
			
			// copy normalized vec
			for (int j = 0; j < dim; j++) {
				eigenVectors[i][j] = mySystem.eigenVectors[i][j].divide(norm, mc);
			}
		}
		
		// refine the eigenvectors by inverse Arnoldi iteration (also time it)
		start = System.currentTimeMillis();
		EigenSystem refinedESystem = refineEvecs (matrix, weights, eigenVectors, eigenValues, bandwidth, precision, extaPrecision, firstEValBeZero);
		end = System.currentTimeMillis();
		System.out.printf("# Refining eigenvectors took %d ms\n", (end-start));
		
		System.out.println("# Done solving eigensystem");
		return refinedESystem;
	}

	private static EigenSystem diagonalMatrixEigensolution(BigDecimal[][] matrix) {
		int n = matrix.length;
		
		BigDecimal[][] eigenVectors = new BigDecimal[n][n];
		BigDecimal[] eigenValues = new BigDecimal[n];
		
		for (int i = 0; i < n; i++)	{
			eigenValues[i] = matrix[i][i];
			Arrays.fill(eigenVectors[i], BigDecimal.ZERO);
			eigenVectors[i][i] = BigDecimal.ONE;
		}
		return new EigenSystem(eigenVectors, eigenValues);
	}

	// calls JLapack for solving initial system
	public static EigenSystem solveJLAPackEigs (BigDecimal[][] H) throws EigenRefineException {

		int dim = H.length;
		// assert squaredness
		for (BigDecimal[] row : H) {
			assert (row.length == dim);
		}

		// we should make a copy of the original matrix, cause finding eigensystem might mess up the matrix
		// also convert into double
		double[][] workMatrix = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				workMatrix[i][j] = H[i][j].doubleValue();
			}
		}

		// use jlapack stuff for basic eigendecomposition

		// the sandboxes for the jlapack eigenstuff
		double[] realValues = new double[dim];
		double[] imaginaryValues = new double[dim];
		double[][] leftVectors = new double[dim][dim];
		double[][] rightVectors = new double[dim][dim];
		double[] sandBox = new double[5 * dim];
		intW returnInt = new intW(0);
		// calculate some eigenstuff (hopefully the correct call)
		DGEEV.DGEEV("N", "V", dim, workMatrix, realValues, imaginaryValues, leftVectors, rightVectors, sandBox, sandBox.length, returnInt);

		// put them into the new container, fill class variable
		double[] eigenValues = new double[dim];
		// also check that all imaginary parts are zero, that is we only
		// want to have positive eigenvalues
		for (int i = 0; i < dim; i++) {
			// copy (might actually be too much, but ok)
			eigenValues[i] = realValues[i];
			
			// we allow them to be a bit imaginary, but only if the real part is very small
			// APARENTLY WE HAVE TO ALWAYS ALLOW A SMALL IMAGINARY PART cause of close difference
			// the orthogonalization should take care of this later
			if (Math.abs (imaginaryValues[i]) > jLApackPrecision.doubleValue()) {
				throw new EigenRefineException("After jLApack, eigenvalue " + i + " is too complex.");
			}
			// it's ok that two eigenvalues are the same,
			// cause they will be treated differently anyways
		}
		// seems like the eigenvalues are not ordered
		
		// I think the eigenvector corresponding to eigenvalues i is now in
		// eigenVector[:][i]
		// so transpose the matrix, cause i want to
		// fill class variable
		double[][] eigenVectors = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				eigenVectors[i][j] = rightVectors[j][i];
			}
		}

		//convert to BigDecimal
		BigDecimal[][] bdEigenVectors = new BigDecimal[dim][dim];
		BigDecimal[] bdEigenValues = new BigDecimal[dim]; 
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				bdEigenVectors[i][j] = BigDecimal.valueOf(eigenVectors[i][j]);
			}
			bdEigenValues[i] = BigDecimal.valueOf(eigenValues[i]);
		}
		
		
		// everything said and done
		// this also sorts them 
		return new EigenSystem(bdEigenVectors, bdEigenValues);
	}
	

	// refine eigenvectors up to given precision 
	private static EigenSystem refineEvecs(BigDecimal[][] matrix, BigDecimal[] weights, BigDecimal[][] eigenVectors, BigDecimal[] eigenValues, int bandwidth, int precision, int extraPrecision, boolean firstEValBeZero) throws EigenRefineException, EigenRefineError {
		// get a MathContext and an epsilon that are precise enough
		MathContext globalMC = new MathContext(precision + extraPrecision, RoundingMode.HALF_EVEN);
		// here is the epsilon, but leave some room
		BigDecimal epsilon = BigDecimal.ONE.scaleByPowerOfTen (- precision + 3);
		
		// remember the eigenvalues that are close together in one group
		ArrayList<TreeSet<Integer>> closeToGroups = new ArrayList<TreeSet<Integer>>();
		
		// WE HAVE TO ASSUME THAT THE EIGENVECTORS ARE ORDERED
		
		boolean secondEVSmall = false;
		// this is a bit ugly hacking:
		if ((eigenValues[0].abs(globalMC).compareTo(jLApackPrecision) > 0) && firstEValBeZero) {
			throw new EigenRefineError("The given seed for the 0-th eigenvalue is NOT smaller than " + jLApackPrecision);
		}
		if (firstEValBeZero) {
			if (! (eigenValues[1].abs(globalMC).compareTo(jLApackPrecision) < 0) ) {
				// only the first one is very small
				secondEVSmall = false;
				// so we can seed it with zero
				eigenValues[0] = BigDecimal.ZERO;
			}
			else {
				// first two are very small
				secondEVSmall = true;
				System.out.println ("# [WARNING] eigenvalue 1 smaller than " + jLApackPrecision + " in absolute value: " + eigenValues[1]);
				// seed first one to zero
				eigenValues[0] = BigDecimal.ZERO;
				// and second one to precision
				eigenValues[1] = jLApackPrecision;
			}
		}

		assert (jLApackPrecision.compareTo(groupingPreciscion) < 0);
		// we need a tmp group for the loop
		TreeSet<Integer> tmpGroup = new TreeSet<Integer>();
		// loop over rest
		for (int eigIdx=0; eigIdx<eigenVectors.length-1; eigIdx++) {
						
			// make sure that all the other eigenvalues are bigger
			if ((eigIdx > 1) && (eigenValues[eigIdx].abs(globalMC).compareTo(groupingPreciscion) < 0)) {
				throw new EigenRefineError("Seed for eigenvalue " + eigIdx + " is also smaller than " + groupingPreciscion);
			}
			
			// remember the eigenvalues that are very close to their upper neighbour
			if (!(eigenValues[eigIdx+1].subtract (eigenValues[eigIdx], globalMC).compareTo (groupingPreciscion) > 0)) {
				// also, it can only be eigenvalue number one, if we are not in the second eigenvalue small mode
				if ((eigIdx == 1) && secondEVSmall && firstEValBeZero) {
					throw new EigenRefineError("Eigenvalue 1 is close to zero, but also to the neighbour above.");
				}

				// and remember it (the lower one of the two)
				tmpGroup.add(new Integer(eigIdx));
			}
			else {
				// the current one is NOT close to its upper neighbor
				// put it into the current group anyway
				tmpGroup.add(new Integer(eigIdx));
				
				// we have to end the current group and add it to the list of groups
				// but only if it is big enough
				if (tmpGroup.size() > 1) {
					closeToGroups.add (tmpGroup);
				}
				
				// and start a new group
				tmpGroup = new TreeSet<Integer>(); 
			}
		}

		
		// one set to rule them all
		TreeSet<Integer> closeToFullSet = new TreeSet<Integer>();
		// go through and save (maybe have a look)
		System.out.println("# After jLApack, the following eigenvalues are close to one another:");
		for (TreeSet<Integer> thisGroup : closeToGroups) {
			System.out.print("# ");
			for (Integer thisInt : thisGroup) {
				// add it to full set
				closeToFullSet.add (thisInt);
				// but also show it
				System.out.print(thisInt + ", ");
			}
			System.out.println();
		}		
		
		// remember seed for eigenvector 0
		BigDecimal[] zeroSeed = new BigDecimal[eigenVectors[0].length];
		// fill it
		for (int i=0; i<zeroSeed.length; i++) {
			zeroSeed[i] = eigenVectors[0][i];
		}
		
		
		// now start the real refinement
		int n = eigenVectors.length;
		int maxRestarts = 10;
		int maxIters = 40;	//increase this if worried that won't converge to within precision even in these many iterations
		int iterReset = 20; // if we reset iterations, go back by this much
		int totIters = 0;	//track the total number of iterations of refinement across all eigenvalues
		int totItersSq = 0;	//track the sum of squares of iterations, to calculate the variance later
		// vector to calculate the squaredlengthes as we go
		BigDecimal[] squaredLength = new BigDecimal[n];
		// initialize it to something negative
		for (int i=0; i<squaredLength.length; i++) {
			squaredLength[i] = mOne;
		}
		
		System.out.println("# Starting refinement");
		
		// attention, now in reversed order
		for (int eigIdx = n-1; eigIdx >= 0; eigIdx--) {

			// we need a local mc, because we might have to adjust it in some situations
			// also reset it with each iteration
			MathContext localMC = globalMC;
			int thisMaxIters = maxIters;
			
			// we might have to re-seed some eigenvalues
			if (closeToFullSet.contains(new Integer(eigIdx))) {

				// if zero it's special, cause we can reseed it against everything else
				if (eigIdx == 0) {
					// first get the set of all indices exceot 0
					TreeSet<Integer> allButOne = new TreeSet<Integer>();
					for (int d=1; d<eigenValues.length; d++) {
						allButOne.add(new Integer(d));
					}
					
					// now start from the given seed and make the eigenvector orthogonal to everything else
					eigenVectors[0] = orthogonalizeToHigherVectorGroup(eigenVectors, eigenVectors[0], 0, allButOne, weights, globalMC);
					// also normalize the vector
					BigDecimal norm = getNorm(eigenVectors[0], globalMC);
					for (int d=0; d<eigenVectors[0].length; d++) {
						eigenVectors[0][d] = eigenVectors[0][d].divide(norm, globalMC);
					}
					
					// then set the value accordingly
					BigDecimal[] vectorResult = MatrixPower.matrixVectorBanded(matrix, eigenVectors[0], bandwidth, globalMC);
					
					// get the fraction of the entries
					// WE CANNOT CHECK ALL, SINCE SOME OF THEM MIGHT BE VERY SMALL
					// so just take the maximum in the original vector
					int maxIdx = 0;
					BigDecimal maxValue = eigenVectors[0][0].abs();
					// see though
					for (int d=1; d<eigenVectors[0].length; d++) {
						// is it the new maximum?
						if (eigenVectors[0][d].abs().compareTo(maxValue) > 0) {
							// yes, update
							maxValue = eigenVectors[0][d].abs();
							maxIdx = d;
						}
					}
					// at the end, take the quotient at the max position 
					BigDecimal newEigenValue = vectorResult[maxIdx].divide(eigenVectors[0][maxIdx], globalMC);
					System.out.println ("# newEigenValue: " + newEigenValue);
					
					// and update it
					eigenValues[0] = newEigenValue;

				}
				// otherwise just standard
				else {
					// re-seed it to be orthogonal to the other guys in his group
					eigenVectors[eigIdx] = orthogonalizeToHigherVectorGroup (eigenVectors, eigenVectors[eigIdx], eigIdx, findHisGroup (eigIdx, closeToGroups), weights, localMC);
				}
				
			}
			
			// is it converged yet?
			boolean converged = false;
			// the error
			BigDecimal maxErr = null;
			// storage for vector from previous iteration
			BigDecimal[] prevVec = null;
			// just get it
			BigDecimal smallestError = BigDecimal.ONE;
			int restarts = 0;
			
			// iterate until conversion
			for (int iter = 0; iter < thisMaxIters; iter++) {

				
				// debugging
				if (eigIdx < 4) {
					System.out.print("# (" + eigIdx + "," + iter + ") \t" + eigenValues[eigIdx]);
					for (int l=0; l<6; l++){
						BigDecimal value = eigenVectors[eigIdx][l];
						if (value.compareTo(BigDecimal.ZERO) < 0)
							System.out.print ("\t" + value);
						else
							System.out.print ("\t " + value);
					}
					System.out.println ();
				}

				
				// how good are we now?
				if (eigIdx == 0 || (eigIdx == 1 && secondEVSmall)) {
					// in this case the error is the paralleleness
					// TODO: maybe do this parallelness business for all of them
					// but we need at least one previous iteration
					if (iter > 0) {
						// check whether direction to previous iteraton changed much
						maxErr = normVectorDifference (eigenVectors[eigIdx], prevVec, localMC);
						// for relative error we should normalize this by the norm of one of the vectors
						// but this should be one
						if (getNorm (eigenVectors[eigIdx], localMC).subtract (BigDecimal.ONE, localMC).abs(localMC).compareTo(epsilon) > 0) {
							throw new EigenRefineException("During refinement, the norm of eigenvector " + eigIdx + " is not 1: " + getNorm (eigenVectors[eigIdx], localMC));
						}
													
						// and remember error (this list works as a queue, put it in front)
						smallestError = smallestError.min (maxErr);
						
						// if we get close to the end, we might have to increase the precision
						if ((iter > maxIters - 3) && (restarts < maxRestarts)) {
							// see if you can update the precision somehow
							
							// first get the magnitude of the smallest error
							int magnitude = getMagnitude (smallestError, localMC);
							
							// and make a new math context (with some wiggling room)
							int newPrecision = localMC.getPrecision() + magnitude + precision + 5;
							System.out.println ("# current precision: " + localMC.getPrecision() + ", magnitude of the error: " + magnitude);
							System.out.println ("# new precision we would like: " + newPrecision);
							
							// make the new context
							// but don't be ridiculous, yeha, another magic constant
							if (magnitude < 0) {
								localMC = new MathContext (newPrecision, localMC.getRoundingMode());
								
								// also clear the old errors, so we start kind of afresh
								smallestError = BigDecimal.ONE;
								// and reset your iterations
								iter -= iterReset;
								// remember the restarting
								restarts += 1;
							}
						}
						
						
					}
					else {
						// not done yet
						maxErr = BigDecimal.ONE;
					}
				}
				else {
					// in the non-small cases, we take the absolute error
					maxErr = normAxMinusLambdaX (matrix, eigenValues[eigIdx], eigenVectors[eigIdx], bandwidth, localMC);
				}

				// show some error
				if (eigIdx < 4) {
					System.out.println ("# maxErr: " + maxErr);
				}
				
				// if we are as good as we want to be, we can stop iterating
				if (maxErr.compareTo(epsilon) <= 0) {
					// remember the number of iterations
					totIters += iter;
					totItersSq += iter*iter;
					// and say that you are done
					converged = true;
					break;
				}
				
				// if we are not good enough, let's try to become better
				// now lets try to improve the current eigenvalue/vector pair
				BigDecimal currValue = eigenValues[eigIdx];
				// we need to copy everything, so we do not loose it
				BigDecimal[] copyVec = new BigDecimal[n];
				for (int i = 0; i < n; i++) {
					copyVec[i] = eigenVectors[eigIdx][i];
				}
				
				// now build the matrix for the linear system (matrix - currValue * Id)
				BigDecimal[][] sysMatrix = new BigDecimal[n][n];
				// fill her up
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n; j++) {
						// copy the value from matrix
						sysMatrix[i][j] = matrix[i][j];
					}
					// and on the diagonal we have - currValue
					sysMatrix[i][i] = sysMatrix[i][i].subtract(currValue, localMC);
				}

				
				// and solve it
				LinearSolver.geppLinearSolveBanded (sysMatrix, copyVec, bandwidth, localMC);
				
				// update the eigenvalue [ALWAYS]
				// we need some inner products to update the eigenvalue
				BigDecimal innProdNumerator = innerProduct(copyVec, eigenVectors[eigIdx], localMC);
				BigDecimal innProdDenominator = innerProduct(eigenVectors[eigIdx], eigenVectors[eigIdx], localMC);
				
				//calculate new value
				BigDecimal newValue = currValue.add(innProdDenominator.divide(innProdNumerator, localMC), localMC);

				// and update it
				eigenValues[eigIdx] = newValue;
				
				// update the eigenvector (always)
				// but first get norm of new copyVec
				BigDecimal norm = getNorm (copyVec, localMC);
				// and normalize it
				for (int i = 0; i < n; i++) {
					copyVec[i] = copyVec[i].divide(norm, localMC);
				}

				
				// save old one
				prevVec = eigenVectors[eigIdx];
				// update vector pair
				eigenVectors[eigIdx] = copyVec;
				
				// and got to next round of checking goodness of fit
			}
			
			
			//if didn't converge in the given number of iterations, it sucks
			if (!converged) {
				// we should actually leave this run
				System.out.println ("# [WARNING] the eigen refinement did not converge. Perhaps your cutoff and recision are not high enough.");
				// that's what we do now
				throw new EigenRefineException (String.format("eigIdx = %d\t DID NOT converge in %d iterations\t maxErr = %s", eigIdx, thisMaxIters, maxErr.toEngineeringString()));
			}
		}
		
		
		// set scales right [since they might have been increased at some point]
		for (int i=0; i<eigenVectors[0].length; i++) {
			eigenVectors[0][i] = eigenVectors[0][i].setScale(extraPrecision + precision, RoundingMode.HALF_EVEN);
		}
		if (secondEVSmall) {
			for (int i=0; i<eigenVectors[1].length; i++) {
				eigenVectors[1][i] = eigenVectors[1][i].setScale(extraPrecision + precision, RoundingMode.HALF_EVEN);
			}
		}
		System.out.println ("# Eigenvalue 0: " + eigenValues[0]);
		System.out.println ("# Eigenvalue 1: " + eigenValues[1]);
		
		
		// first two eigenvalues are allowed to have switched order, but only if the second one is small
		if ((eigenValues[0].compareTo(eigenValues[1]) > 0) && !secondEVSmall) {
			throw new EigenRefineException("Eigenvalue 1 is smaller than 0, although they were far apart berfore refinement.");
		}
		
		// now we have to sort the eigenvectors according to increasing eigenvalue
		// constructor of eigensystem does this (at the end, when we return the values)
		EigenSystem resultSystem = new EigenSystem(eigenVectors, eigenValues);
		// just so that the remaining checks work
		// we should actually also be able to modify them
		// YES, WE CAN
		eigenValues = resultSystem.eigenValues;
		eigenVectors = resultSystem.eigenVectors;
		
		// now look at the result
		if (eigenValues[1].abs(globalMC).compareTo(BigDecimal.ONE) < 0) {
			secondEVSmall = true;
			System.out.println ("# [WARNING] the second eigenvalue is smaller than 1, so maybe the answer is bogus.");
		}
		else {
			secondEVSmall = false;
		}
		// check first two for matrix equation
		// first we check the first
		BigDecimal diff0 = normAxMinusLambdaX(matrix, eigenValues[0], eigenVectors[0], bandwidth, globalMC);
		// norm of vector should be one
		if (getNorm (eigenVectors[0], globalMC).subtract (BigDecimal.ONE, globalMC).abs(globalMC).compareTo(epsilon) > 0) {
			throw new EigenRefineException("After refinement, the norm of eigenvector 0 is not 1: " + getNorm (eigenVectors[0], globalMC));
		}
		// first eigenvalue very small, so just take absolute thing
		if (diff0.compareTo(epsilon) > 0) {
			throw new EigenRefineException ("Eigenvalue 0 does not satisfy eigenequation: " + eigenValues[0]);
		}
		
		// then the second one, but only if we have to
		if (secondEVSmall) {
			BigDecimal diff1 = normAxMinusLambdaX(matrix, eigenValues[1], eigenVectors[1], bandwidth, globalMC);
			// norm of vector should be one
			if (getNorm (eigenVectors[1], globalMC).subtract (BigDecimal.ONE, globalMC).abs(globalMC).compareTo(epsilon) > 0) {
				throw new EigenRefineException("After refinement, the norm of eigenvalue 1 is not 1: " + getNorm (eigenVectors[0], globalMC));
			}
			// second one should not be too small
			if (eigenValues[1].abs().compareTo(epsilon.multiply(BigDecimal.ONE.scaleByPowerOfTen(5),globalMC)) < 0) {
				// too small for comfort
				throw new EigenRefineException("The second eigenvalue is indistinguishable from zero. Try increasing cutoff or precision.");
			}
			// second eigenvalue very small, so just take absolute thing
			if (diff1.compareTo(epsilon) > 0) {
				throw new EigenRefineException ("Eigenvalue 1 does not satisfy eigenequation: " + eigenValues[1]);
			}
		}
		
		// the vectors all have to be pairwise orthogonal with respect to given weight vector
		// but since we know that only those guys that where close to their above partner might be problematic
		// we only check those
		for (TreeSet<Integer> currentGroup : closeToGroups) {
			//no go through and assert that the inner products are small enough
			System.out.println("# Close to group:");
			for (Integer I : currentGroup) {
				for (Integer J : currentGroup) {
					if (I < J) {
						// get the final difference of the eigenvalues
						BigDecimal diff = eigenValues[I.intValue()].subtract(eigenValues[J.intValue()], globalMC).abs();
						// and the inner product
						BigDecimal innerProdAbs = weightedInnerProduct (eigenVectors[I.intValue()], eigenVectors[J.intValue()], weights, globalMC).abs();

						// some debugging
						System.out.println("# difference between " + I.intValue() + " and " + J.intValue() + " is: " + diff);
						System.out.println("# inner product (absolute value) <" + I.intValue() + ", " + J.intValue() + "> is: " + innerProdAbs);

						// see whether they are as orthogonal as our method can make them
						if (diff.compareTo(BigDecimal.ZERO) == 0) {
							// if the two eigencalues are exactly equal, we don't thtrow an exception, but maybe a warning
							System.out.println ("# [WARNING] eigenvalue " + I.intValue() + " and " + J.intValue() + " are exactly equal, so you might not trust your results.");
						}
						else {
							// difference is not zero, so usual stuff
							if (innerProdAbs.compareTo (epsilon.divide (diff, globalMC)) > 0) {
								// somehow it is still too  big
								throw new EigenRefineException ("Too big <" + I.intValue() + "," + J.intValue() + "> " + innerProdAbs);
							}
						}
					}
				}
			}
			

			// then orthogonalify all vectors against the ones above
			for (Integer I : currentGroup) {
				// orthogonalize this one
				eigenVectors[I.intValue()] = orthogonalizeToHigherVectorGroup (eigenVectors, eigenVectors[I.intValue()], I.intValue(), currentGroup, weights, globalMC); 			
			}
			
			// check at least the guys in this group against the rest
			int min = currentGroup.first();
			int max = currentGroup.last();
			for (int l=0; l<eigenVectors.length; l++) {
				// only if it is not the same vector
				// we assume them consecutive
				if (l < min || max<l) {
					for (Integer I : currentGroup) {
						BigDecimal orthFirst = weightedInnerProduct (eigenVectors[I.intValue()], eigenVectors[l], weights, globalMC);
						if (orthFirst.abs(globalMC).compareTo (epsilon) > 0) {
							throw new EigenRefineException ("Too big <" + I.intValue() + "," + l + "> " + orthFirst);
						}
					}
				}
			}			
		}
		
		
		// see whether the pairwise inner products are actually what we expect
		System.out.printf("# Total # of refinement iterations for all %d eigenvalues = %d\n", n, totIters);
		System.out.printf("# Standard deviation of # of refinement iterations = %f\n", Math.sqrt((1.*totItersSq - 1.*totIters*totIters/n)/n));
		
		// and return it
		return resultSystem;
	}

	
	private static TreeSet<Integer> findHisGroup(int value, ArrayList<TreeSet<Integer>> groups) {
		Integer intObject = new Integer(value);
		// find the group that value belongs to 
		for (TreeSet<Integer> thisGroup : groups) {
			if (thisGroup.contains(intObject)) {
				// found it
				return thisGroup;
			}
		}
		// otherwise return null
		return null;
	}

	// function to orthogonalize the vectors
	private static BigDecimal[] orthogonalizeToHigherVectorGroup (BigDecimal[][] basisVectors, BigDecimal[] vec, int vecIdx, TreeSet<Integer> vectorGroup, BigDecimal[] weights, MathContext localMC) throws EigenRefineException {
		
		// something to orthogonalize it against
		// the result
		BigDecimal[] resultVec = new BigDecimal[vec.length];
		// copy old vector
		for (int i=0; i<resultVec.length; i++) {
			resultVec[i] = vec[i];
		}
		
		// and now go through all vectors that we should be orthogonal to and subtract the right thing
		for (Integer groupMember : vectorGroup) {
			int toOrthAgainstIdx = groupMember.intValue(); 
			// the index should be larger, we only orthogonalize against higher guys
			if (toOrthAgainstIdx > vecIdx) {
				// get the requisite inner products
				BigDecimal ab = weightedInnerProduct (vec, basisVectors[toOrthAgainstIdx], weights, localMC);
				
				BigDecimal bb  = weightedInnerProduct (basisVectors[toOrthAgainstIdx], basisVectors[toOrthAgainstIdx], weights, localMC);
				if (bb.compareTo(BigDecimal.ZERO) == 0) {
					throw new EigenRefineException("The length of basisVector " + toOrthAgainstIdx + " is 0.");
				}

				// get the factor
				BigDecimal factor = ab.divide(bb, localMC);
				
				// then adjust the result vector
				for (int i=0; i<resultVec.length; i++) {
					resultVec[i] = resultVec[i].subtract(basisVectors[toOrthAgainstIdx][i].multiply(factor, localMC), localMC);
				}
				// this should be done
			}
		}
		
		// give it away now
		return resultVec;
	}
	
	private static BigDecimal normVectorDifference(BigDecimal[] newVec, BigDecimal[] oldVec, MathContext mc) {
		assert (newVec.length == oldVec.length);
		
		// the first component should point in the right direction
		if (newVec[0].unscaledValue().signum() != oldVec[0].unscaledValue().signum()) {
			// flip oldVec
			for (int l=0; l<newVec.length; l++) {
				oldVec[l] = BigDecimal.ZERO.subtract (oldVec[l]);
			}
		}
		
		// make a new container for the difference
		BigDecimal[] vectorDiff = new BigDecimal[newVec.length];
		for (int i=0; i<newVec.length; i++) {
			vectorDiff[i] = newVec[i].subtract(oldVec[i], mc).abs(mc);
		}
		
		// give it away now (the norm of the vector difference)
		return getNorm (vectorDiff, mc);
	}

	
	private static BigDecimal normAxMinusLambdaX (BigDecimal[][] A, BigDecimal lambda, BigDecimal[] w, int bandwidth, MathContext mc) {
		//MathContext extraMc = new MathContext(mc.getPrecision() + 10, RoundingMode.HALF_EVEN);
		
		assert (w.length == A.length);
		// now do the matrix multiplication
		BigDecimal[] vectorResult = MatrixPower.matrixVectorBanded(A, w, bandwidth, mc);

		// now we only want the distance to the scalar multiplication to be small
		// in every component
		for (int i = 0; i < vectorResult.length; i++) {
			// get the scalar from the normal multiplication
			vectorResult[i] = vectorResult[i].subtract (lambda.multiply (w[i], mc), mc);
		}

		// return the norm of the difference vector
		return getNorm (vectorResult, mc);
	}

	private static BigDecimal getNorm(BigDecimal[] vector, MathContext mc) {
		// store it here
		BigDecimal norm = BigDecimal.ZERO;
		
		// go through and do stuff
		for (int j = 0; j < vector.length; j++) {
			norm = norm.max(vector[j].abs(mc));
		}
		
		// return stuff
		return norm;
	}
	
	private static BigDecimal innerProduct(BigDecimal[] a, BigDecimal[] b, MathContext mc) {
		assert (a.length == b.length);
		BigDecimal ret = BigDecimal.ZERO;
		for (int i = 0; i < a.length; i++) {
			ret = ret.add(a[i].multiply(b[i], mc), mc);
		}
		return ret;
	}

	private static BigDecimal weightedInnerProduct (BigDecimal[] a, BigDecimal[] b, BigDecimal[] weights, MathContext mc) {
		assert (a.length == b.length);
		assert (a.length == weights.length);
		BigDecimal ret = BigDecimal.ZERO;
		for (int i = 0; i < a.length; i++) {
			ret = ret.add(a[i].multiply(b[i], mc).multiply(weights[i], mc), mc);
		}
		return ret;
	}

	// get a rough magnitude of the big decimal at hand
	private static int getMagnitude (BigDecimal value, MathContext mc) {
		boolean smallerThanOne = false;
		// get absolute value
		BigDecimal tmpVal = value.abs (mc);
		// smaller or larger than one
		if (value.compareTo(BigDecimal.ONE) < 0) {
			smallerThanOne = true;
			tmpVal = BigDecimal.ONE.divide (tmpVal, mc);
		}

		int magnitude = 0;
		// now iterate to get the approximate magnitude
		while (tmpVal.compareTo(BigDecimal.ONE) > 0) {
			tmpVal = tmpVal.scaleByPowerOfTen (-1);
			magnitude++;
		}
		
		// return it in the appropriate way
		if (smallerThanOne) return - magnitude;
		else return magnitude;
	}

	
	// just an exception class
	// this one might be recoverable with more precision
	@SuppressWarnings("serial")
	public static class EigenRefineException extends Exception {
		public EigenRefineException(String err) {
			super (err);
		}
	};
	
	// and here goes the real errors
	// the ones that don't seem to be recoverable even with higher precision
	@SuppressWarnings("serial")
	public static class EigenRefineError extends Exception {
		public EigenRefineError(String err) {
			super (err);
		}
	};
	
	
	
	public static BigDecimal groupingPreciscion = BigDecimal.ONE.scaleByPowerOfTen(-4);
	public static BigDecimal jLApackPrecision = BigDecimal.ONE.scaleByPowerOfTen(-5);
	static public BigDecimal ZeroPointOne = BigDecimal.ONE.scaleByPowerOfTen(-1);
	static public BigDecimal mOne = BigDecimal.ZERO.subtract(BigDecimal.ONE);
}
