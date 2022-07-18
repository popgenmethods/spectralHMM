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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintStream;
import java.io.Reader;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Formatter;

import org.nevec.rjm.BigDecimalMath;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.SimpleJSAP;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.stringparsers.FileStringParser;

import edu.berkeley.spectralHMM.algorithms.Combinatorics;
import edu.berkeley.spectralHMM.matrix.MatrixPower;
import edu.berkeley.spectralHMM.matrix.TDFEigenSolverBanded.EigenRefineError;
import edu.berkeley.spectralHMM.matrix.TDFEigenSolverBanded.EigenRefineException;

public class SelectionHMM {

	public static BigDecimal run (BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, InitialConditionEnum initCondition, BigDecimal initFrequency, int matrixCutoff, int maxM, int maxN, TimedSample sample, Boolean ignoreBinomials, boolean condOnLastSegregating, MathContext mc) throws EigenRefineException, EigenRefineError	{	

		// get the sample
		ArrayList<BigDecimal> times = sample.times;
		ArrayList<Integer> N = sample.N;
		ArrayList<Integer> S = sample.S;
		
		// this the number of samples
		int numSamples = N.size(); 

		// everything all right?
		assert (times.size() == N.size() + 1);
		assert (times.size() == S.size() + 1);

		// a and b
		ArrayList<ArrayList<BigDecimal>> a = new ArrayList<ArrayList<BigDecimal>> (numSamples + 1);
		ArrayList<ArrayList<BigDecimal>> b = new ArrayList<ArrayList<BigDecimal>> (numSamples + 1);
		
		// initialize for all timepoints
		for (int i = 0; i <= numSamples; i++) {
			ArrayList<BigDecimal> va = new ArrayList<BigDecimal>(maxN + 1);
			ArrayList<BigDecimal> vb = new ArrayList<BigDecimal>(maxN + 1);
			a.add(va);
			b.add(vb);
		}

		
		// get constants
		BigDecimal C = SelectionEigenfunktion.getC (alpha, beta, hetF, homF, matrixCutoff, maxM, mc);
		BigDecimal c0 = SelectionEigenfunktion.evaluateSquaredLength (alpha, beta, hetF, homF, matrixCutoff, maxM, 0, mc);
		// check C
		System.out.println ("# [see]\t" + C + "\t" + c0);
		
		BigDecimal[][] W = SelectionEigenfunktion.getEigenVectorMatrixCopy (alpha, beta, hetF, homF, matrixCutoff, maxN+1, maxM+1, mc);
		
		// transpose the matrix
		BigDecimal[][] tildeW = MatrixPower.transpose (W);
		for (int j = 0; j <= maxN; j++)	{
			BigDecimal cj = SelectionEigenfunktion.evaluateSquaredLength(alpha, beta, hetF, homF, matrixCutoff, maxM, j, mc);
			for (int i = 0; i <= maxM; i++)	{
				tildeW[i][j] = tildeW[i][j].divide(cj, mc);
			}
		}
		for (int i = 0; i <= maxM; i++)	{
			BigDecimal ciSquaredLength = JacobiPolynomials.evaluateSquaredLength (alpha, beta, i, mc);
			for (int j = 0; j <= maxN; j++)	{
				tildeW[i][j] = tildeW[i][j].multiply(ciSquaredLength, mc);
			}
		}
		
		
		// copy the initial values
		ArrayList<BigDecimal> vb0 = b.get(0);

		// and initialize them
		if (initCondition == InitialConditionEnum.MutationSelection) {
			// this is initial vector for mutation selection balance
			for (int n = 0; n <= maxN; n++)	{
				// don't forget the normalizing constant
				if (n == 0) vb0.add(C.divide(c0, mc));
				else vb0.add(BigDecimal.ZERO.setScale(mc.getPrecision()));
			}
		}
		else if (initCondition == InitialConditionEnum.InitialFrequency) {
			// this is initial vector for a given initial frequency
			for (int n = 0; n <= maxN; n++)	{
				BigDecimal BnX = SelectionEigenfunktion.evaluate (alpha, beta, hetF, homF, matrixCutoff, maxM, n, initFrequency, mc);
				BigDecimal cn = SelectionEigenfunktion.evaluateSquaredLength (alpha, beta, hetF, homF, matrixCutoff, maxM, n, mc);
				vb0.add (BnX.divide(cn, mc));
			}
		}
		else if (initCondition == InitialConditionEnum.MutationDrift) {
		// this is initial vector for a drift mutation balance
			//implicitly using heterozygous and homozygous fitnesses for hetF and homF
			BigDecimal mHetF = BigDecimal.ZERO.subtract(hetF, mc);
			BigDecimal mHomF = BigDecimal.ZERO.subtract(homF, mc);
			
			// get a copy of the eigenvector
			BigDecimal[] v = SelectionEigenfunktion.getEigenVectorCopy(alpha, beta, mHetF, mHomF, matrixCutoff, maxM, 0, mc);
			BigDecimal Cp = SelectionEigenfunktion.getC (alpha, beta, mHetF, mHomF, matrixCutoff, maxM, mc);
			BigDecimal thisBetaFunction = JacobiPolynomials.betaFunction(alpha, beta, mc);
			for (int i=0; i<v.length; i++) {
				v[i] = v[i].multiply(JacobiPolynomials.evaluateSquaredLength (alpha, beta, i, mc), mc);
			}

			// calculate entries
			for (int n = 0; n <= maxN; n++)	{
				BigDecimal value = BigDecimal.ZERO;
				for (int m=0; m<=maxM; m++) {
					BigDecimal increment = W[n][m].multiply(v[m], mc);
					value = value.add(increment, mc);
				}
				BigDecimal cn = SelectionEigenfunktion.evaluateSquaredLength (alpha, beta, hetF, homF, matrixCutoff, maxM, n, mc);
				value = value.divide(cn, mc).divide(thisBetaFunction, mc).divide(Cp, mc);
				vb0.add(value);
			}		
		}
		else {
			// no other initial condition allowed for now
			assert (false);
		}

		
		// analyze sample
		for (int k = 1; k <= numSamples; k++) {
			Integer Nk = N.get(k-1);
			Integer Sk = S.get(k-1);
			
			//update a's
			ArrayList<BigDecimal> va = a.get(k);
			for (int n = 0; n <= maxN; n++)	{
				BigDecimal LambdaN = SelectionEigenfunktion.getEigenValue (alpha, beta, hetF, homF, matrixCutoff, maxM, n, mc);
				BigDecimal val = BigDecimalMath.exp (BigDecimal.ZERO.setScale(mc.getPrecision()).subtract(times.get(k).subtract(times.get(k-1), mc).multiply(LambdaN , mc), mc), mc);
				val = val.multiply(b.get(k-1).get(n), mc);
				va.add(val);
			}

			//update b's
			ArrayList<BigDecimal> vb = b.get(k);
			
			// get the whole N
			BigDecimal[][] nMatrix = getNMatrix (alpha, beta, Nk, Sk, maxM, mc);
			
			BigDecimal[][] aMat = new BigDecimal[1][maxN+1];
			for (int i = 0; i <= maxN; i++)	{
				aMat[0][i] = a.get(k).get(i);
			}
			BigDecimal[][] ret = MatrixPower.multiplyMatricesBanded(aMat, W, maxN, maxM, mc);
			ret = MatrixPower.multiplyMatricesBanded(ret, nMatrix, maxM, Nk, mc);
			ret = MatrixPower.multiplyMatricesBanded(ret, tildeW, maxM, maxM, mc);
			for (int i = 0; i <= maxN; i++)	{
				if (! ignoreBinomials)	{
					// calculate the binomial
					BigDecimal binom = new BigDecimal (Combinatorics.choose(N.get(k-1), S.get(k-1)).toString());
					ret[0][i] = ret[0][i].multiply(binom, mc);
				}
				vb.add(ret[0][i]);
			}
			
		}
		
		// to get final answer, we see whether we want to condition
		BigDecimal ans = b.get(numSamples).get(0);
		
		// whats the conditioning
		if (!condOnLastSegregating) {
			// just normalize it with the right value
			ans = ans.multiply(c0, mc).divide(C, mc);
		}
		else {
			// if we condition on the last one segregating, we have to normalize our result
			// we have to start with the initial condition in b^(0) and evolve it all the way till the end
			
			// get the total time
			BigDecimal totalTime = times.get(numSamples).subtract (times.get (0), mc);
			
			// and evolve
			ArrayList<BigDecimal> totalA = new ArrayList<BigDecimal>(maxN + 1);
			for (int n = 0; n <= maxN; n++)	{
				BigDecimal LambdaN = SelectionEigenfunktion.getEigenValue (alpha, beta, hetF, homF, matrixCutoff, maxM, n, mc);
				BigDecimal val = BigDecimalMath.exp (BigDecimal.ZERO.subtract(totalTime, mc).multiply(LambdaN , mc), mc);
				val = val.multiply(b.get(0).get(n), mc);
				totalA.add(val);
			}

			//  then compute the probability of finding no derived allele at the end
			// we only need endB_0
			// thus, we only need some matrix product between some endA matrix
			BigDecimal[][] totalAMat = new BigDecimal[1][maxN+1];
			for (int i = 0; i <= maxN; i++)	{
				totalAMat[0][i] = totalA.get(i);
			}
			// end the corresponding sampling matrix
			int lastN = N.get(N.size()-1);
			BigDecimal[][] nMatrix = getNMatrix (alpha, beta, lastN, 0, maxM, mc);
			// with some W involved
			BigDecimal[][] ret = MatrixPower.multiplyMatricesBanded (totalAMat, W, maxN, maxM, mc);
			ret = MatrixPower.multiplyMatricesBanded(ret, nMatrix, maxM, lastN, mc);
			ret = MatrixPower.multiplyMatricesBanded(ret, tildeW, maxM, maxM, mc);

			// copy the last entry
			BigDecimal totalB_0 = ret[0][0];
			if (! ignoreBinomials)	{
				// calculate the binomial
				BigDecimal binom = new BigDecimal (Combinatorics.choose(lastN, 0).toString());
				totalB_0 = totalB_0.multiply(binom, mc);
			}

			// then adjust the final answer (divide by prob of segregating) [put the factor in there]
			ans = ans.divide (C.divide(c0, mc).subtract (totalB_0, mc), mc);

		}
		
		// and return the result
		return ans;
	}

	public static BigDecimal[][] getNMatrix(BigDecimal alpha, BigDecimal beta, Integer Nk, Integer Sk, int maxM, MathContext mc) {
		BigDecimal[][] L = JacobiPolynomials.coefficientMatrix(alpha, beta, Sk, maxM, mc, true);
		BigDecimal[][] M = JacobiPolynomials.coefficientMatrix(alpha, beta, Nk - Sk, maxM, mc, false);
		BigDecimal[][] ret = MatrixPower.multiplyMatricesBanded (L, M, Sk, Nk - Sk, mc);
		return ret;
	}

	@SuppressWarnings("unused")
	private static BigDecimal[][] getNMatrixRescaled(BigDecimal alpha, BigDecimal beta, Integer Nk, Integer Sk, int maxM, MathContext mc) {
		BigDecimal[][] L = JacobiPolynomials.coefficientMatrixRescaled(alpha, beta, Nk, Sk, maxM, mc, true, (Sk > Nk - Sk));
		BigDecimal[][] M = JacobiPolynomials.coefficientMatrixRescaled(alpha, beta, Nk, Nk - Sk, maxM, mc, false, (Sk <= Nk - Sk));
		BigDecimal[][] ret = MatrixPower.multiplyMatricesBanded (L, M, Sk, Nk - Sk, mc);
		return ret;
	}


//	@SuppressWarnings("unused")
//	private static BigDecimal[][] getEigenfunctions(BigDecimal alpha, BigDecimal beta, BigDecimal hetF, BigDecimal homF, int matrixCutoff, int maxM, int maxN, MathContext mc) throws EigenRefineException, EigenRefineError {
//		// make a grid
//		int resolution = 50;
//		BigDecimal bigResolution = new BigDecimal(resolution); 
//		BigDecimal[] theGrid = new BigDecimal[resolution + 1];
//		// beginning
//		theGrid[0] = BigDecimal.ZERO;
//		theGrid[0] = theGrid[0].setScale(mc.getPrecision());
//		// the middle
//		for (int i=1; i< resolution; i++) {
//			theGrid[i] = new BigDecimal(i).divide(bigResolution, mc);
//			theGrid[i] = theGrid[i].setScale(mc.getPrecision());
//		}
//		// end
//		theGrid[resolution] = BigDecimal.ONE;
//		theGrid[resolution] = theGrid[resolution].setScale(mc.getPrecision());
//
//		// we also need the weight function on the yGrid
//		BigDecimal[] weightFunction = new BigDecimal[theGrid.length];
//		for (int y=0; y<weightFunction.length; y++) {
////			 exp(\bar\sigma(y)) y^(alpha-1) (1-y)^(beta-1)
//			BigDecimal jacobiWeight = BigDecimalMath.pow(theGrid[y], alpha.subtract(BigDecimal.ONE, mc), mc).multiply(BigDecimalMath.pow (BigDecimal.ONE.subtract(theGrid[y], mc), beta.subtract(BigDecimal.ONE, mc), mc), mc);
//			weightFunction[y] = (BigDecimalMath.exp (SelectionEigenfunktion.meanFitness(hetF, homF, theGrid[y], mc), mc)).multiply (jacobiWeight, mc);
//		}
//		
//		// get the forward eigenfunctions on that grid
//		BigDecimal[][] results = new BigDecimal[maxN+1][theGrid.length];
//		for (int n=0; n<=maxN; n++) {
//			for (int y=0; y<theGrid.length; y++) {
//				results[n][y] = weightFunction[y].multiply(SelectionEigenfunktion.evaluate (alpha, beta, hetF, homF, matrixCutoff, maxM, n, theGrid[y], mc), mc);
//			}
//		}
//		return results;
//	}

	public static BigDecimal N(int j, int l, int nk, int sk, BigDecimal alpha, BigDecimal beta, MathContext mc) {
		
		// just to be sure
		assert (nk >= sk);

		// initialize return value
		BigDecimal ret = BigDecimal.ZERO;
		
		int offset = JacobiPolynomials.calculateOffset(alpha, beta);
		for (int mu = Math.max(0, Math.max(j - sk, l - (nk - sk))); mu <= Math.min(j + sk, l + (nk - sk)); mu++){
			ret = ret.add(JacobiPolynomials.L(alpha, beta, j, mu, sk, offset, mc).multiply(JacobiPolynomials.M(alpha, beta, mu, l, (nk-sk), offset, mc), mc), mc);
		}
		// multiply it with length
		BigDecimal clSquaredLength = JacobiPolynomials.evaluateSquaredLength (alpha, beta, l, mc);
		ret = ret.multiply(clSquaredLength, mc);
		
		// return the final thing
		return ret;
	}

	// debugging output of the vector
	@SuppressWarnings("unused")
	private static void printVector (String name, ArrayList<BigDecimal> vb0, int maxN, BigDecimal c0, BigDecimal C, MathContext mc) {
		System.out.println ("# " + name + "\t");
		for (BigDecimal value : vb0) {
			System.out.println (value + "\t");
		}
		System.out.println();
		assert (vb0.size() == maxN + 1);		
	}

	
	// wrapper to catch exceptions a bit better
	public static void main (String[] args) {
		try {
			realMain (args);
		} catch (FileNotFoundException e) {
			System.err.println ("Error: Invalid input file secified:\n\t" + e.getMessage());
		} catch (JSAPException e) {
			System.err.println ("Error while parsing command line arguments (--help for usage):\n\t" + e.getMessage());
		} catch (IOException e) {
			System.err.println ("I/O error:\n\t" + e.getMessage());
		}
	}

	
	public static void realMain (String[] args) throws JSAPException, FileNotFoundException, IOException {
		long startTime = System.currentTimeMillis();

		// build a parser for the input file, with the right parameters
		FileStringParser myFileParser = FileStringParser.getParser();
		myFileParser.setMustBeFile(true);
		myFileParser.setMustExist(true);
		
		// build the jsap object
		SimpleJSAP jsap = new SimpleJSAP( 
	            "spectralHMM", 
	            "Analyze temporal data using a spectral HMM method.",
	            new Parameter[] {
					new FlaggedOption( "inputFile", myFileParser, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'f', "inputFile", "Specify an input file. See documentation for formatting." ),
					new Switch ("mutSelBalance", 'j', "mutSelBalance", "Set the initial condition to be selection-drift balance."),
					new Switch ("mutDriftBalance", 'd', "mutDriftBalance", "Set the initial condition to be mutation-drift balance."),
					new FlaggedOption( "initFrequency", JSAP.BIGDECIMAL_PARSER, "-0.5", JSAP.NOT_REQUIRED, 'i', "initFrequency", "Set the initial frequency." ),
					new FlaggedOption( "initTime", JSAP.BIGDECIMAL_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 't', "initTime", "The initial time. If not given, the initial time is expected to be in the file with each dataset. Also, this only works with multiplexing." ),
					new FlaggedOption( "mutToBenef", JSAP.BIGDECIMAL_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'a', "mutToBenef", "The mutation rate from the wild type allele to the selected allele." ),
					new FlaggedOption( "mutFromBenef", JSAP.BIGDECIMAL_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'b', "mutFromBenef", "The mutation rate from the selected allele to the wild-type allele." ),
					new FlaggedOption( "selection", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 's', "selection", "Selection coefficient s (can be a range: [start:step:stop])." ),
					new FlaggedOption( "dominance", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'h', "dominance", "Dominance parameter h (can be a range: [start:step:stop])." ),
					new FlaggedOption( "hetF", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'w', "hetF", "Fitness of the heterozygote (can be a range: [start:step:stop])." ),
					new FlaggedOption( "homF", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'v', "homF", "Fitness of the homozygote (can be a range: [start:step:stop])." ),
					new FlaggedOption( "effPopSize", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'e', "effPopSize", "The effective population size (diploid)." ),
					new FlaggedOption( "yearsPerGen", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'y', "yearsPerGen", "Specify how many years a generation takes." ),
					new FlaggedOption( "matrixCutoff", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'c', "matrixCutoff", "The cutoff for the matrix whose eigenvectors yield the coefficients for the eigenfunctions." ),
					new FlaggedOption( "maxM", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'm', "maxM", "Specify how many summands to use in the infinite sum for each eigenfunction." ),
					new FlaggedOption( "maxN", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'n', "maxN", "Specify how many eigenfunctions/-values to use in the computations." ),
					new FlaggedOption( "precision", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, 'p', "precision", "Specify a precision to be used for the computations. This gives the number of significant digits used." ),
					new Switch ("ignoreBinomials", 'z', "ignoreBinomials", "If set, the likelihood is calculated without the binomial coefficients required to make it a probability."),
					new Switch ("multiplex", JSAP.NO_SHORTFLAG, "multiplex", "Analyze each of the multiple datasests given in the input file using all the selection parameters specified."),
					new Switch ("condOnLastSegregating", JSAP.NO_SHORTFLAG, "condOnLastSegregating", "When set, condition on the last sample being segregating."),
	            }
	        );

		JSAPResult config = jsap.parse(args);
		if (jsap.messagePrinted()) { System.exit(1); }
		
		// what about the binomials
		boolean ignoreBinomials = config.getBoolean ("ignoreBinomials");
		// what about multiplexing
		boolean multiplex = config.getBoolean ("multiplex");
		// what about conditioning
		boolean condOnLastSegregating = config.getBoolean ("condOnLastSegregating");
		
		// precision
		int precision = config.getInt("precision");
		int scale = precision;
		// set precision and some math context
		MathContext mc = new MathContext (precision, RoundingMode.HALF_EVEN);

		// what is the parametrization
		boolean paramByHS = false;
		String param1Name = null, param2Name = null;
		// selection parameterization
		if (config.contains("selection") && config.contains("dominance"))	{
			paramByHS = true;
			param1Name = "selection";
			param2Name = "dominance";
			if (config.contains("hetF") || config.contains("homF"))	{
				System.err.println("Specify exactly either the selection and dominance parameters, or the heterozygote and homozygote advantage parameters");
				System.exit(1);
			}
		}
		else if (config.contains("hetF") && config.contains("homF"))	{
			param1Name = "hetF";
			param2Name = "homF";
			if (config.contains("selection") || config.contains("dominance"))	{
				System.err.println("Specify exactly either the selection and dominance parameters, or the heterozygote and homozygote advantage parameters");
				System.exit(1);
			}
		}
		else {
			System.err.println("Specify exactly either the selection and dominance parameters, or the heterozygote and homozygote advantage parameters");
			System.exit(1);
		}
				
		
		// initial time?
		BigDecimal initTime = null;
		boolean initTimesInFile = !config.contains ("initTime");
		if (initTimesInFile) {
			if (!multiplex) {
				System.err.println ("The flag timesInFile only works together with multiplexing.");
				System.exit(-1);
			}
		}
		else {
			//we do have an initial time
			initTime = config.getBigDecimal("initTime");
		}
		
		// get the initial condition right
		BigDecimal initialFrequency = config.getBigDecimal("initFrequency");
		boolean initFreqSet = (initialFrequency.compareTo(BigDecimal.ZERO) >= 0d) && (initialFrequency.compareTo(BigDecimal.ONE) <= 0d);
		boolean mutSelSet = config.getBoolean ("mutSelBalance");
		boolean mutDriftSet = config.getBoolean ("mutDriftBalance");
		// not check whether exactly one is set, and which one it is
		InitialConditionEnum initialCondition = null;
		if (initFreqSet && !mutSelSet && !mutDriftSet) {
			// initial frequency given
			initialCondition = InitialConditionEnum.InitialFrequency;
		}
		else if (!initFreqSet && mutSelSet && !mutDriftSet) {
			// mutation selection balance requested
			initialCondition = InitialConditionEnum.MutationSelection;
		}
		else if (!initFreqSet && !mutSelSet && mutDriftSet) {
			// mutation drift balance requested
			initialCondition = InitialConditionEnum.MutationDrift;
		}
		else {
			System.err.println ("Must specify consistent initial conditions.");
			System.exit(1);
		}
		
		// get input file
		File inputFile = config.getFile("inputFile");
		
		// get raw mutation parameters
		BigDecimal Ne = (new BigDecimal(config.getString("effPopSize"))).setScale(scale);
		BigDecimal preAlpha = config.getBigDecimal("mutToBenef").setScale(scale);
		if (BigDecimal.ZERO.compareTo(preAlpha) >= 0) throw new IOException ("Zero mutation rate not implemented yet.");
		BigDecimal preBeta = config.getBigDecimal("mutFromBenef").setScale(scale);
		if (BigDecimal.ZERO.compareTo(preBeta) >= 0) throw new IOException ("Zero mutation rate not implemented yet.");
		// and rescale them to be population scaled
		BigDecimal alpha = (new BigDecimal("4")).multiply(Ne, mc).multiply(preAlpha, mc);
		BigDecimal beta = (new BigDecimal("4")).multiply(Ne, mc).multiply(preBeta, mc);
		
		BigDecimal yearsPerGen = (new BigDecimal(config.getString("yearsPerGen"))).setScale(scale);
		
		BigDecimal[] param1Range = null;
		BigDecimal[] param2Range = null;
		
		// some raw selection coefficients
		if (config.contains("selection"))	{
			assert(config.contains("dominance"));
			param1Range = parseRange (config.getString("selection"), scale, mc);
			param2Range = parseRange (config.getString("dominance"), scale, mc);
		}
		if (config.contains("hetF"))	{
			assert(config.contains("homF"));
			param1Range = parseRange (config.getString("hetF"), scale, mc);
			param2Range = parseRange (config.getString("homF"), scale, mc);
		}

		// fill the list of values
		ArrayList<BigDecimal> param1s = buildRange (param1Range[0], param1Range[1], param1Range[2], mc);

		// also we need the other parameters
		
		// make a list of other parameters
		ArrayList<BigDecimal> param2s = buildRange (param2Range[0], param2Range[1], param2Range[2], mc);
		
		// the remaining parameters
		int matrixCutoff = config.getInt("matrixCutoff");
		int maxM = config.getInt("maxM");
		int maxN = config.getInt("maxN");
		
		assert(maxM <= matrixCutoff);
		assert(maxN <= maxM);
		
		//print the command line arguments directly, and also the parameters nicely
//		PrintStream[] outStreams = new PrintStream[]{System.out, System.err};
		PrintStream[] outStreams = new PrintStream[]{System.out};
		for (PrintStream outStream: outStreams) {
			outStream.print("# Command-line arguments: ");
			for (String arg: args)	{
				outStream.print(arg + " ");
			}
			outStream.println();
			
			outStream.println("# Parameter values:");
			outStream.println("# mutToBenef = " + bigDecimalToString(preAlpha));
			outStream.println("# alpha = " + bigDecimalToString(alpha));
			outStream.println("# mutFromBenef = " + bigDecimalToString(preBeta));
			outStream.println("# beta = " + bigDecimalToString(beta));
			outStream.println("# selection param 1 (" + param1Name + ") = " + bigDecimalToString(param1Range[0]) + "\t" + bigDecimalToString(param1Range[1]) + "\t" + bigDecimalToString(param1Range[2]));
			outStream.println("# selection param 2 (" + param2Name + ") = " + bigDecimalToString(param2Range[0]) + "\t" + bigDecimalToString(param2Range[1]) + "\t" + bigDecimalToString(param2Range[2]));
			outStream.println("# precision = " + precision);
			outStream.println("# matrixCutoff = " + matrixCutoff);
			outStream.println("# maxM = " + maxM);
			outStream.println("# maxN = " + maxN);
			outStream.println("# effective Population size = " + bigDecimalToString(Ne) + " diploids");
			outStream.println("# Years per generation = " + bigDecimalToString(yearsPerGen));
			outStream.println("# ignore binomials = " + ignoreBinomials);
			outStream.println("# mulitplex = " + multiplex);
			outStream.println("# condition on last segregating = " + condOnLastSegregating);
			outStream.println ("# initial times in input file = " + initTimesInFile);
		}

		ArrayList<TimedSample> listOfSamples = new ArrayList<TimedSample>();
		// get the data
		if (!multiplex) {
			// only one dataset
			// read it
			ArrayList<BigDecimal> rawTimes = new ArrayList<BigDecimal>();
			ArrayList<Integer> N = new ArrayList<Integer>();
			ArrayList<Integer> S = new ArrayList<Integer>();
	
			// parse it from file
			parseInputSingleData (new FileReader(inputFile), rawTimes, N, S, scale);
			
			// add an additional time for the beginning
			assert (!rawTimes.isEmpty());
			rawTimes.add (0, initTime);
			// check here
			if (!ordered  (rawTimes)) {
				throw new IOException ("The times are not increasing.");
			}
	
			// adjust times in the sample for Ne and generation time
			BigDecimal divisor = (new BigDecimal("2")).multiply(Ne, mc).multiply(yearsPerGen, mc);
			ArrayList<BigDecimal> times = new ArrayList<BigDecimal>();
			for (BigDecimal time : rawTimes) {
				times.add (time.divide(divisor, mc));
			}
			
			// put it into a list
			listOfSamples.add (new TimedSample(times, N, S));
		}
		else {
			// a lot of datasets
			listOfSamples = parseInputMultipleSamples  (new FileReader(inputFile), initTimesInFile, initTime, Ne, yearsPerGen, mc);
		}
		
		if (condOnLastSegregating) {
			// see, what's the minimal sample that occurs at the last position
			int minLastPosition = scanLastPosForMin (listOfSamples);
			
			if (minLastPosition < 1) {
				System.err.println ("If we analyze samples conditional on segregation at last timepoint, then all your samples should be segregating.");
				System.exit(-1);
			}
		}
		
		// build some likelihood surface
		for (BigDecimal param1 : param1s) {
			for (BigDecimal param2 : param2s) {
				// some variables
				BigDecimal hetF, homF, h=null, s;
				
				//which parametrization?
				if (paramByHS)	{
					// homF = s
					homF = param1;
					// homF = s*h
					hetF = param2.multiply(param1, mc);
					// get the others for printing
					s = param1;
					h = param2;
				}
				else	{
					hetF = param1;
					homF = param2;
					// get the others to print nicely
					s = homF;
					if (homF.compareTo(BigDecimal.ZERO) != 0)	{
						h = hetF.divide(homF, mc);
					}
				}
				
				BigDecimal hetFScaled = (new BigDecimal("2")).multiply(Ne, mc).multiply(hetF, mc);
				BigDecimal homFScaled = (new BigDecimal("2")).multiply(Ne, mc).multiply(homF, mc);
				
				BigDecimal sigma  = (new BigDecimal("2")).setScale(scale).multiply(Ne, mc).multiply(s, mc);
				
				// print sigma into comment
				System.out.println("# " + param1Name + " = " + bigDecimalToString(hetF) + "\t" + param2Name + " = " + bigDecimalToString(homF));
				System.out.println("# hetF = " + bigDecimalToString(hetF) + "\thomF = " + bigDecimalToString(homF));
				System.out.println("# hetFScaled = " + bigDecimalToString(hetFScaled) + "\thomFScaled = " + bigDecimalToString(homFScaled));
				System.out.println("# s = " + bigDecimalToString(s) + "\th = " + (h == null ? "undefined" : bigDecimalToString(h)));
				System.out.println("# sigma = " + bigDecimalToString(sigma));

				// need a flag
				boolean successful = true;
				
				// now iterate over all samples
				for (int d=0; d<listOfSamples.size(); d++) {
					TimedSample thisSample = listOfSamples.get(d);
										
					// compute the likelihood
					BigDecimal ans;
					try {
						// have we been successful before
						if (successful) {
							// run it
							ans = run (alpha, beta, hetFScaled, homFScaled, initialCondition, initialFrequency, matrixCutoff, maxM, maxN, thisSample, ignoreBinomials, condOnLastSegregating, mc);
						}
						else {
							// if not, return the dummy
							ans = returnUnsuccessful;
						}
					} catch (EigenRefineException e) {
						// say something about not converged
						System.out.println ("# [EXCEPTION (maybe recoverable)] " + e.getMessage());
						// remember that it was not successful
						successful = false;
						// and give a dummy value
						ans = returnUnsuccessful;
					} catch (EigenRefineError e) {
						// say something about not converged
						System.out.println ("# [ERROR (not recoverable)] " + e.getMessage());
						// remember that it was not successful
						successful = false;
						// and give a dummy value
						ans = returnUnsuccessful;
					}
					
					// print sample
					printSample (thisSample, initialCondition, initialFrequency, System.out);
					
					// print result out
					String outPut = "";
					
					// add sample number if you want
					if (listOfSamples.size() > 1) {
						// plus one, cause we want a more natural numbering
						outPut += (d+1) + "\t";
					}
					
					// add first selection parameter
					outPut += bigDecimalToString(param1) + "\t";
					
					// add second selection parameter
					outPut += bigDecimalToString(param2) + "\t";
					
					// print it
					System.out.println (outPut + ans);
				}
				
				// clear caches
				// only selection, after all samples have been done
				SelectionEigenfunktion.clearCaches();
			}
		}
		
		long endTime = System.currentTimeMillis();
		
		System.out.printf("# Elapsed time = %d ms\n", endTime-startTime);
	}

	
	private static boolean ordered (ArrayList<BigDecimal> times) {
		for (int i=1; i<times.size(); i++) {
			if (times.get(i-1).compareTo(times.get(i)) >= 0) return false;
		}
		return true;
	}

	private static int scanLastPosForMin(ArrayList<TimedSample> listOfSamples) {
		// get the number of sampling times
		// WARNING: we allow this to vary
		
		// and start finding the minimum
		int thisNumSamplingTimes = listOfSamples.get(0).S.size();
		int min = listOfSamples.get(0).S.get(thisNumSamplingTimes-1);
		for (int i=1; i<listOfSamples.size(); i++) {
			thisNumSamplingTimes = listOfSamples.get(i).S.size();
			min = Math.min (min, listOfSamples.get(i).S.get(thisNumSamplingTimes-1));
		}

		// give it away
		return min;
	}

	private static ArrayList<TimedSample> parseInputMultipleSamples (FileReader fileReader, boolean initTimesInFile, BigDecimal initTime, BigDecimal Ne, BigDecimal yearsPerGen, MathContext mc) throws IOException {

		// the list with the results
		ArrayList<TimedSample> listOfSamples = new ArrayList<SelectionHMM.TimedSample>();
		
		// create an object to read lines
		LineNumberReader lineReader = new LineNumberReader (fileReader);
		
		// and go through lines
		String line;
		// each line should contain a sample
		while ((line = lineReader.readLine()) != null) {

			// ignore lines starting with '#'
			if (line.startsWith("#") || line.trim().isEmpty()) continue;
			
			// new containers
			ArrayList<BigDecimal> rawTimes = new ArrayList<BigDecimal>();
			if (!initTimesInFile) {
				// add the standard initTime
				rawTimes.add (initTime);
			}

			ArrayList<Integer> sampleSizes = new ArrayList<Integer>();
			ArrayList<Integer> numDerivedAlleles = new ArrayList<Integer>();
			
			// get the values
			readVectorFromLine (line, initTimesInFile, rawTimes, sampleSizes, numDerivedAlleles);
						
			// there should be some time
			if (rawTimes.isEmpty()) {
				throw new IOException("no times");
			}
			
			// adjust times in the sample for Ne and generation time
			BigDecimal divisor = (new BigDecimal("2")).multiply(Ne, mc).multiply(yearsPerGen, mc);
			ArrayList<BigDecimal> times = new ArrayList<BigDecimal>();
			for (BigDecimal time : rawTimes) {
				times.add (time.divide(divisor, mc));
			}
			
			// put it into a list
			listOfSamples.add (new TimedSample(times, sampleSizes, numDerivedAlleles));
		}
		
		// give away the list of samples read in
		return listOfSamples;
	}

	static void readVectorFromLine(String line, boolean initTimesInFile, ArrayList<BigDecimal> rawTimes, ArrayList<Integer> sampleSizes, ArrayList<Integer> numDerivedAlleles) throws IOException {
		// first clear the lists
		sampleSizes.clear();
		numDerivedAlleles.clear();
		
		// get the three-tuples from the vector
		// first split it into larger pieces
		String[] fields = line.split(";");
		
		int sampleIdx = 0;
		if (initTimesInFile) {
			// the first entry should only be the time
			rawTimes.add (new BigDecimal(fields[0].trim()));
			// set index to skip the first entry in the rest
			sampleIdx = 1;
		}
		
		// now, if we should have a time in the file, then it should be at the first position
		for (; sampleIdx<fields.length; sampleIdx++) {
			
			// then split the fields into hopefully three pieces
			String[] pieces = fields[sampleIdx].split("(\\))|(\\()|(\\,)");
			
			// we have to clean it a bit
			ArrayList<String> cleanPieces = new ArrayList<String>();
			for (String piece : pieces) {
				// throw away empty ones
				if (!piece.trim().isEmpty()) {
					cleanPieces.add(piece);
				}
			}
			
			// now it should be three pieces
			if (cleanPieces.size() != 3) {
				throw new IOException ("Component does not contain the right number of values: " + fields[sampleIdx] + "\nin line: " + line);
			}
			
			// now store them
			rawTimes.add (new BigDecimal(cleanPieces.get(0).trim()));
			sampleSizes.add (new Integer (cleanPieces.get(1).trim()));
			numDerivedAlleles.add (new Integer (cleanPieces.get(2).trim()));
		}
		
		// before we go, check the times
		if (!ordered  (rawTimes)) {
			throw new IOException ("The times are not increasing in at least one dataset.");
		}

	}

	static void printSample (TimedSample thisSample, InitialConditionEnum initialCondition, BigDecimal initialFrequency, PrintStream outStream) {
		// print initial condition
		if (initialCondition == InitialConditionEnum.InitialFrequency) {
			outStream.println("# initial condition: frequency (" + initialFrequency + ")");
		}
		else if (initialCondition == InitialConditionEnum.MutationSelection) {
			outStream.println("# initial condition: mutation selection balance");
		}
		else if (initialCondition == InitialConditionEnum.MutationDrift) {
			outStream.println("# initial condition: mutation drift balance");
			
		}
		else {
			assert (false);
		}
		outStream.println("# Sample");
		
		outStream.print("# ");
		for (BigDecimal time : thisSample.times)	{
			outStream.print(bigDecimalToString(time) + "\t");
		}
		outStream.println();
		
		outStream.print("# \t");
		for (Integer n : thisSample.N)	{
			outStream.print(n + "\t");
		}
		outStream.println();
		
		outStream.print("# \t");
		for (Integer k : thisSample.S)	{
			outStream.print (k + "\t");
		}
		outStream.println();
	}

	private static ArrayList<BigDecimal> buildRange (BigDecimal start, BigDecimal step, BigDecimal stop, MathContext mc) {
		// just be sure
		assert (start.compareTo(stop) < 1);
		assert (step.compareTo(BigDecimal.ZERO) > -1);
		
		ArrayList<BigDecimal> result = new ArrayList<BigDecimal>();
		// initial
		BigDecimal tmp = start;
		while (tmp.compareTo(stop) < 1){
			// add
			result.add(tmp);
			// increase
			tmp = tmp.add(step, mc);
		}

		return result;
	}

	
	private static BigDecimal[] parseRange(String rangeString, int bigDecimalScale, MathContext mc) throws JSAPException {
		if (rangeString.startsWith("[")) {
			// expect a real range
			if (!rangeString.endsWith("]")) { throw new JSAPException ("Range for selection has to be in format [start:step:stop]"); }
			
			// parse the range
			String[] fields = rangeString.substring(1, rangeString.length()-1).split(":");

			// get the array
			if (fields.length != 3) { throw new JSAPException ("Range for selection has to be in format [start:step:stop]"); }
			BigDecimal[] tmp = new BigDecimal[] { (new BigDecimal(fields[0])).setScale(bigDecimalScale), (new BigDecimal(fields[1])).setScale(bigDecimalScale), (new BigDecimal(fields[2])).setScale(bigDecimalScale) };
			
			// now check for the right stepping direction
			if (tmp[2].subtract(tmp[0], mc).multiply(tmp[1], mc).compareTo(BigDecimal.ZERO) < 0) { throw new JSAPException ("Range for selection has to be in format [start:step:stop]"); }
			// and change the order if necessary
			if (tmp[0].compareTo(tmp[2]) > 0) {
				BigDecimal swap = tmp[0];
				tmp[0] = tmp[2];
				tmp[1] = tmp[1].multiply(new BigDecimal("-1"), mc);
				tmp[2] = swap;
			}
			// special case, if boundaries are equal
			if (tmp[0].compareTo(tmp[2]) == 0 && tmp[1].compareTo(BigDecimal.ZERO) < 0) {
				// step is negative, so change it
				tmp[1] = tmp[1].multiply(BigDecimal.ZERO.subtract(BigDecimal.ONE, mc), mc);
			}
			
			// give it away now
			return tmp;
		}
		else {
			// should just be a single value
			BigDecimal val = new BigDecimal(rangeString);
			// give it like this
			return new BigDecimal[] {val, BigDecimal.ONE, val};
		}
	}

	private static void parseInputSingleData (Reader In, ArrayList<BigDecimal> times, ArrayList<Integer> sampleSizes, ArrayList<Integer> numDerivedAlleles, int bigDecimalScale) throws IOException {
		
		// clear containers
		times.clear();
		sampleSizes.clear();
		numDerivedAlleles.clear();
		
		// create an object to read lines
		LineNumberReader lineReader = new LineNumberReader (In);
		
		// and go through lines
		String line;
		while ((line = lineReader.readLine()) != null) {
			// ignore lines starting with '#'
			if (line.startsWith("#") || line.trim().isEmpty()) continue;
			
			// get the values
			String[] fields = line.split("\\s+");
			// should have 3 fields
			if (fields.length != 3) {
				throw new IOException ("Invalid input format (maybe missing multiplex).");
			}
			
			// add stuff
			times.add (new BigDecimal(fields[0]).setScale(bigDecimalScale));
			sampleSizes.add (new Integer(fields[1]));
			numDerivedAlleles.add (new Integer(fields[2]));
		}
		// done
	}
	
	
	public static String bigDecimalToString (BigDecimal x){
		StringBuilder sb = new StringBuilder();
		Formatter formatter = new Formatter(sb);
		
		formatter.format("%g", x);
		
		// close the formatter
		formatter.close();
		
		return sb.toString();
	}
	
	// class for a dataset
	public static class TimedSample {
		// members
		public ArrayList<BigDecimal> times;
		public ArrayList<Integer> N;
		public ArrayList<Integer> S;
		
		// easy constructor
		public TimedSample (ArrayList<BigDecimal> times, ArrayList<Integer> N, ArrayList<Integer> S) {
			// check some
			assert (times.size() == N.size()+1);
			assert (N.size() == S.size());
			
			// remember them
			this.times = times;
			this.N = N;
			this.S = S;
		}
	}
	
	public static int firstObservedIndex (ArrayList<BigDecimal> rawTimes, ArrayList<Integer> numDerivedAlleles) {
		// loop through
		int observationIndex = 0;
		while (numDerivedAlleles.get(observationIndex) == 0) {
			observationIndex++;
		}
		// and return the first time
		return observationIndex;
	}

	
	public static TimedSample buildSample (BigDecimal t0, ArrayList<BigDecimal> rawTimes, ArrayList<Integer> sampleSizes, ArrayList<Integer> numDerivedAlleles, BigDecimal Ne, BigDecimal yearsPerGen, MathContext mc) {
		// get the index of first observation
		int observationIndex = firstObservedIndex (rawTimes, numDerivedAlleles);
		// t0 should be right
		assert (t0.compareTo(rawTimes.get(observationIndex)) < 0);		
		
		// make the right times
		ArrayList<BigDecimal> diffusionTimes = new ArrayList<BigDecimal>();

		// get the divisor
		BigDecimal divisor = (new BigDecimal("2")).multiply(Ne).multiply(yearsPerGen);

		// put in the initial time
		BigDecimal diffusionT0 = t0.divide (divisor, mc);
		diffusionTimes.add (diffusionT0);
		// and rest
		for (int i=observationIndex; i<rawTimes.size(); i++) {
			diffusionTimes.add (rawTimes.get(i).divide(divisor, mc));
		}
		// now copy the sample size and derived number appropriately
		ArrayList<Integer> newSampleSizes = new ArrayList<Integer>();
		ArrayList<Integer> newNumDerivedAlleles = new ArrayList<Integer>();
		for (int i=observationIndex; i<sampleSizes.size(); i++) {
			newSampleSizes.add (sampleSizes.get(i));
			newNumDerivedAlleles.add (numDerivedAlleles.get(i));
		}
		
		// and return a new object
		return new TimedSample (diffusionTimes, newSampleSizes, newNumDerivedAlleles);
	}


	// enumeration for the initial condition
	public enum InitialConditionEnum { MutationDrift, MutationSelection, InitialFrequency};
	
	public static BigDecimal returnUnsuccessful = BigDecimal.ZERO.subtract(BigDecimal.ONE);
}
