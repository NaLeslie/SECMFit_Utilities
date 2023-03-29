/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package secm_grid_utilities;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Calendar;
import java.util.LinkedList;
import java.util.Scanner;

/**
 *
 * @author Nathaniel
 */
public class ThreeParameter {

    /**
     * Stand-in for auto-generated method
     * @param reactivity_data_filename
     * @param l The L parameter
     * @param xcoords The x coordinates of each sampled point
     * @param ycoords The y coordinates of each sampled point
     * @return 
     */
    static Model run(String reactivity_data_filename, double l, double[] xcoords, double[] ycoords){
        return new Model();
    }
    
    static Model runk(String reactivity_data_filename, double[] ks, double l, double xcoord, double ycoord){
        return new Model();
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main() throws NumberFormatException, FileNotFoundException {
        list_l = new LinkedList<Double>();
        list_logk = new LinkedList<Double>();
        list_data = new LinkedList<Double[]>();
        System.out.println(getDateStamp());
        try{
            readSECMInfo("Fitfile.csv");
            double amount = -2.0;
            for(int i = 0; i < 9; i++){
                int[][] newgrid = applyDilationErosion(amount);
                exportEditedGrid(i + "_" + amount + ".csv", newgrid);
                amount += 0.5;
            }
            //runFit("yep", 1.0, -5, false);
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    /*
    Fitting methods and fields
    */
    /**
     * Runs the fitting procedure using Levenberg-Marquardt
     * @param filename The name of the file containing the data to be fit to (for logging purposes). {@link #readSECMInfo(java.lang.String)} should be used before this method to load-in the necessary data.
     * @param firstl The initial guess for L
     * @param firstlogk The initial guess for logk. Use of {@link #findFirstLogK(double, boolean)} to obtain this value is recommended.
     * @param verbose flag to be ticked if extra data is desired (the current and derivatives at each iteration)
     * @return The status of the fitting procedure. 
     * <p>{@link #EXECUTED_OK} if the process converged,</p>
     * <p>{@link #MAX_ITERATIONS_REACHED} if there is no convergence after {@link #MAX_ITERATIONS} iterations,</p>
     * <p>or {@link #MAX_LAMBDA_REACHED} if lambda exceeds {@link #MAX_LAMBDA} when trying to compute the next iteration's parameters.</p>
     * @throws FileNotFoundException
     * @throws NumberFormatException
     * @throws IOException 
     */
    static int runFit(String filename, double firstl, double firstlogk, double firstr, boolean verbose) throws FileNotFoundException, NumberFormatException, IOException{
        double[] experimental = true_image;
        double lambda = 0;
        double ssr = 0;
        double last_ssr;
        double l = firstl;
        double logk = firstlogk;
        double r = firstr;
        double last_l;
        double last_logk;
        double last_r;
        //first iteration
        //compute the currents and residuals for the initial parameter guesses
        double[] curr = runModel(firstl, firstlogk, firstr);
        double[] residuals = subtract(experimental, curr);
        
        //compute the derivatives of the current with respect to the parameters
        double[] curr_dl = runModel(firstl + L_PERTURB, firstlogk, firstr);
        double[] dl = multiply(subtract(curr_dl, curr), 1.0/L_PERTURB);
        
        double[] curr_dk = runModel(firstl, firstlogk + LOGK_PERTURB, firstr);
        double[] dk = multiply(subtract(curr_dk, curr), 1.0/LOGK_PERTURB);
        
        double[] curr_dr = runModel(firstl, firstlogk, firstr + R_PERTURB);
        double[] dr = multiply(subtract(curr_dr, curr), 1.0/R_PERTURB);
        
        //construct the Jacobian, and perform an iteration of Levenberg-Marquardt
        //[dl dk dr]
        double[][] J = appendColumn(appendColumn(dl, dk), dr);
        double[][] JT = transpose(J);
        double[][] JTJ = multiply(JT, J);
        double[][] DTD = new double[JTJ.length][JTJ.length];
        for(int i = 0; i < DTD.length; i++){
            for(int j = 0; j < DTD.length; j++){
                if(i != j){
                    DTD[i][j] = 0;
                }
                else{
                    DTD[i][j] = JTJ[i][j];
                }
            }
        }
        double[][] lam_DTD = multiply(DTD, lambda);
        double[][] JTJinv = invert(add(JTJ, lam_DTD));
        double[] delta_c = multiply(multiply(JTJinv, JT), residuals);
        ssr = sumSquare(residuals);
        //initialize 'previous iteration' data
        last_ssr = ssr;
        last_l = l;
        last_logk = logk;
        last_r = r;
        l = last_l + round(delta_c[0], L_DECIMALS);
        logk = last_logk + round(delta_c[1], LOGK_DECIMALS);
        r = last_r + round(delta_c[2], R_DECIMALS);
        logInitialGuesses(filename, new String[]{"L", "log10k", "dilation/erosion"}, new double[]{l, logk, r}, ssr);
        boolean converged = false;
        int iterations = 1;
        if(verbose){
            writeIteration("Iteration_1.txt", physical_xs, physical_ys, curr, dl, dk, dr, residuals);
        }
        //subsequent iterations
        while(!converged && iterations <= MAX_ITERATIONS){
            iterations ++;
            System.out.println("Iteration " + iterations); 
            System.out.println("Lambda 0");
            //First lambda
            curr = runModel(l, logk, r);
            residuals = subtract(experimental, curr);
            ssr = sumSquare(residuals);
            boolean lambda_ok = ssr < last_ssr;
            logIteration(iterations, new double[]{DTD[0][0], DTD[1][1]}, lambda, new String[]{"L", "log10k", "dilation/erosion"}, new double[]{l, logk, r}, ssr, lambda_ok);

            //subsequent lambdas (if necessary)
            while(!lambda_ok && lambda <= MAX_LAMBDA){
                lambda = nextLambda(lambda);
                System.out.println("Lambda " + lambda);
                lam_DTD = multiply(DTD, lambda);
                JTJinv = invert(add(JTJ, lam_DTD));
                delta_c = multiply(multiply(JTJinv, JT), residuals);
                l = last_l + round(delta_c[0], L_DECIMALS);
                logk = last_logk + round(delta_c[1], LOGK_DECIMALS);
                r = last_r + round(delta_c[2], R_DECIMALS);
                curr = runModel(l, logk, r);
                residuals = subtract(experimental, curr);
                ssr = sumSquare(residuals);
                lambda_ok = ssr < last_ssr;
                logLambda(lambda, new String[]{"L", "log10k", "dilation/erosion"}, new double[]{l, logk, r}, ssr, lambda_ok);
                //if the prescribed changes to the parameters are small enough, declare convergence
                if(Math.abs(delta_c[0]) < 0.5*Math.pow(10, -L_DECIMALS) && Math.abs(delta_c[1]) < 0.5*Math.pow(10, -LOGK_DECIMALS) && Math.abs(delta_c[2]) < 0.5*Math.pow(10, -R_DECIMALS)){
                    System.out.println("DeltaC: " + delta_c[0] + "\t" + delta_c[1] + "\t" + delta_c[2]);
                    lambda_ok = true;
                    converged = true;
                    return EXECUTED_OK;
                }
            }
            //check if the loop ended due to lambda hitting its maximum
            if(!lambda_ok){
                return MAX_LAMBDA_REACHED;
            }
            //end the process if this was the last iteration.
            if(iterations >= MAX_ITERATIONS){
                return MAX_ITERATIONS_REACHED;
            }
            //Experimental: (use the set of parameters with the lowest residuals to-date)
            int lowest_sim = findLowestSSR();
            curr = convertToRegularDouble(list_data.get(lowest_sim));
            residuals = subtract(experimental, curr);
            ssr = sumSquare(residuals);
            l = list_l.get(lowest_sim);
            logk = list_logk.get(lowest_sim);
            //end of Experimental bit
            //update the previous iteration parameters
            last_ssr = ssr;
            last_l = l;
            last_logk = logk;
            int[][] de_grid = applyDilationErosion(r);
            int gridcount = getSum(de_grid);
            //compute or look up L derivative
            int query = checkList(l - L_PERTURB, logk, gridcount);
            if(query == -1){
                curr_dl = runModel(l + L_PERTURB, logk, r);
                dl = multiply(subtract(curr_dl, curr), 1.0/L_PERTURB);
            }
            else{
                curr_dl = convertToRegularDouble(list_data.get(query));
                //NOTE: this curr_dl is effectively simulated as being perturbed negatively
                dl = multiply(subtract(curr_dl, curr), -1.0/L_PERTURB);
            }
            
            //compute or look up logk derivative
            query = checkList(l, logk - LOGK_PERTURB, gridcount);
            if(query == -1){
                curr_dk = runModel(l, logk + LOGK_PERTURB, r);
                dk = multiply(subtract(curr_dk, curr), 1.0/LOGK_PERTURB);
            }
            else{
                curr_dk = convertToRegularDouble(list_data.get(query));
                //NOTE: this curr_dk is effectively simulated as being perturbed negatively
                dk = multiply(subtract(curr_dk, curr), -1.0/LOGK_PERTURB);
            }
            
            //compute or look up r derivative
            de_grid = applyDilationErosion(r - R_PERTURB);
            gridcount = getSum(de_grid);
            query = checkList(l, logk, gridcount);
            if(query == -1){
                curr_dr = runModel(l, logk, r + R_PERTURB);
                dr = multiply(subtract(curr_dr, curr), 1.0/R_PERTURB);
            }
            else{
                curr_dr = convertToRegularDouble(list_data.get(query));
                //NOTE: this curr_dr is effectively simulated as being perturbed negatively
                dr = multiply(subtract(curr_dr, curr), -1.0/R_PERTURB);
            }

            if(verbose){
                writeIteration("Iteration_" + iterations + ".txt", physical_xs, physical_ys, curr, dl, dk, dr, residuals);
            }
            //construct the jacobian
            J = appendColumn(appendColumn(dl, dk), dr);
            JT = transpose(J);
            JTJ = multiply(JT, J);
            for(int i = 0; i < DTD.length; i++){
                if(JTJ[i][i] > DTD[i][i]){
                    DTD[i][i] = JTJ[i][i];
                }
            }
            lambda = 0;
            lam_DTD = multiply(DTD, lambda);
            JTJinv = invert(add(JTJ, lam_DTD));
            delta_c = multiply(multiply(JTJinv, JT), residuals);
            l = last_l + round(delta_c[0], L_DECIMALS);
            logk = last_logk + round(delta_c[1], LOGK_DECIMALS);
            r = last_r + round(delta_c[2], R_DECIMALS);
            //if the prescribed changes to the parameters are small enough, declare convergence
            if(Math.abs(delta_c[0]) < 0.5*Math.pow(10, -L_DECIMALS) && Math.abs(delta_c[1]) < 0.5*Math.pow(10, -LOGK_DECIMALS) && Math.abs(delta_c[2]) < 0.5*Math.pow(10, -R_DECIMALS)){
                System.out.println("DeltaC: " + delta_c[0] + "\t" + delta_c[1] + "\t" + delta_c[2]);
                converged = true;
                return EXECUTED_OK;
            }
        }
        
        return MAX_ITERATIONS_REACHED;
    }
    
    /**
     * Checks if the parameters have already been used in a previous iteration, then either returns the previously calculated currents or simulates the currents.
     * @param l the L parameter
     * @param logk the base-ten logarithm of the k parameter
     * @return The currents of the microelectrode at the specified relative positions with respect to the reactive feature. 
     * The positions follow the same order that their currents are defined in the control file.
     * @throws FileNotFoundException 
     */
    static double[] runModel(double l, double logk, double r) throws FileNotFoundException{
        int[][] de_grid = applyDilationErosion(r);
        int gridcount = getSum(de_grid);
        int index = checkList(l, logk, gridcount);
        if(index == -1){
            System.out.println(getDateStamp() + ": simulating: L: " + l + "; logk: " + logk + "; erosion/dilation: " + r);
            Model temp = run(reactivity_mapfile, l, physical_xs, physical_ys);
            double[] data = readData();
            eraseDataFile();
            addToList(l, logk, gridcount, data);
            return data;
        }
        else{
            System.out.println(getDateStamp() + ":     lookup: L: " + l + "; logk: " + logk + "; erosion/dilation: " + r);
            return convertToRegularDouble(list_data.get(index));
        }
    }
    
    /**
     * Prints a progress bar to the console
     * @param prog the progress so far
     * @param total the maximum possible value for progress
     * @return The progress bar as a string. Note: the first character of this string is a '\r'.
     */
    static String pBar(int prog, int total){
	int len = 25;
	int filled = (len*prog) / total;
	int percent = (100*prog) / total;
	String bar = " ";
	for(int i = 1; i <= filled; i++){
	  bar = bar.concat("\u2588");
	}
	for(int i = filled + 1; i <= len; i++){
	  bar = bar.concat("\u2591");
	}
	bar = bar.concat("  " + percent + "% complete.");
	return "\r " + bar;
    }
    
    /**
     * Returns The current ISO8601 timestamp (minus the timezone information)
     * @return The current time and date as yyyy-mm-ddThh:mm:ss
     */
    static String getDateStamp(){
        Calendar cl = Calendar.getInstance();
        int year = cl.get(Calendar.YEAR);
        int month = cl.get(Calendar.MONTH) + 1;
        int day = cl.get(Calendar.DAY_OF_MONTH);
        int hour = cl.get(Calendar.HOUR_OF_DAY);
        int minute = cl.get(Calendar.MINUTE);
        int second = cl.get(Calendar.SECOND);
        String datetime = String.format("[%04d-%02d-%02dT%02d:%02d:%02d]",year, month, day, hour, minute, second);
        return datetime;
    }
    
    /**
     * Checks the list of previous simulations and returns where in the list a previous simulation with the same parameters occurred.
     * Or a -1 if no such simulation is in the list.
     * @param l The L parameter to find.
     * @param logk The logk to find.
     * @return -1 if there have been no previous simulations with both L and logk, 
     * or the index in {@link #list_data} that corresponds to the given L and logk.
     */
    static int checkList(double l, double logk, int gridsum){
        for(int i = 0; i < list_l.size(); i++){
            boolean leq = precisionEquals(l, list_l.get(i), L_DECIMALS);
            boolean logkeq = precisionEquals(l, list_logk.get(i), LOGK_DECIMALS);
            boolean gseq = list_gridsum.get(i).intValue() == gridsum;
            if(leq && logkeq && gseq){
                return i;
            }
        }
        return -1;
    }
    
    /**
     * Adds the l, logk, and data to the list of previous simulations.
     * @param l The l parameter.
     * @param logk The logk parameter.
     * @param data The simulated currents.
     */
    static void addToList(double l, double logk, int gridsum, double[] data){
        Double dl = l;
        Double dlogk = logk;
        Integer gsum = gridsum;
        Double[] ddata = convertToClassDouble(data);
        list_l.add(dl);
        list_logk.add(dlogk);
        list_gridsum.add(gsum);
        list_data.add(ddata);
    }
    
    /**
     * Method for converting double[] to Double[]
     * @param input The double[] to be converted.
     * @return The equivalent Double[].
     */
    static Double[] convertToClassDouble(double[] input){
        int len = input.length;
        Double[] output = new Double[len];
        for(int i = 0; i < len; i++){
            output[i] = input[i];
        }
        return output;
    }
    
    /**
     * Method for converting Double[] to double[]
     * @param input The Double[] to be converted.
     * @return The equivalent double[].
     */
    static double[] convertToRegularDouble(Double[] input){
        int len = input.length;
        double[] output = new double[len];
        for(int i = 0; i < len; i++){
            output[i] = (double)input[i];
        }
        return output;
    }
    
    /**
     * Determines if d1 and cd2 are equal to one-another to a certain number of decimal points.
     * @param d1 The double to be compared to cd2.
     * @param cd2 The Double to be compared to d1.
     * @param decimals The number of decimal points to which the equality test is done.
     * @return true if round(abs(d1-cd2)*(10**decimals)) == 0, false otherwise.
     */
    static boolean precisionEquals(double d1, Double cd2, int decimals){
        double d2 = (double)cd2;
        double factor = Math.pow(10, decimals);
        double difference = Math.abs(d1 - d2);
        long rounded_difference = Math.round(difference*factor);
        return rounded_difference == 0;
    }
    
    /**
     * Method for finding the best first guess for logk.
     * @param firstl the initial guess for the L parameter.
     * @param verbose flag for determining if extra data logging is desired.
     * @return The first guess for the logarithm of the k parameter.
     * @throws FileNotFoundException 
     */
    static double findFirstLogK(double firstl, boolean verbose) throws FileNotFoundException, IOException{
        //call the run method to produce a current when the electrode is at distance z and GridData.getCentre() relative to the reactive feature.
        writeReactivityFile(grid, 1.0);
        int[] centre = getCentre();
        Model model = runk(reactivity_mapfile, TEST_LOG_K, firstl, centre[0], centre[1]);
        double[] curr = readData();
        eraseDataFile();
        double max_derivative = curr[2] - curr[0];
        int max_derivative_index = 1;
        for(int i = 2; i < TEST_LOG_K.length -1; i++){
            double derivative = curr[i+1]-curr[i-1];
            if(derivative > max_derivative){
                max_derivative = derivative;
                max_derivative_index = i;
            }
        }
        if(verbose){
            writeKFit(KLOGFILE, curr);
        }
        return TEST_LOG_K[max_derivative_index];
    }
    
    /**
     * Method that controls the lambda-escalation policy for the Levenberg-Marquardt algorithm.
     * @param current_lambda The lambda that was just used.
     * @return The new lambda.
     */
    static double nextLambda(double current_lambda){
        if(current_lambda == 0.0){
            return 1E-4;
        }
        else if(current_lambda == 1E-4){
            return 1E-2;
        }
        else if(current_lambda == 1E-2){
            return 1E-1;
        }
        else if(current_lambda == 1E-1){
            return 0.3;
        }
        else if(current_lambda == 0.3){
            return 1;
        }
        else {
            return 10.0*current_lambda;
        }
    }
    
    /**
     * Computes the sum of squares for the residuals.
     * @param residuals The residuals to be squared and summed.
     * @return The sum of squares.
     */
    static double sumSquare(double[] residuals){
        double sum = 0;
        for(double r : residuals){
            sum += r*r;
        }
        return sum;
    }
    
    /**
     * Rounds value to decimals decimal places.
     * @param value The value to be rounded.
     * @param decimals The number of decimals to round to.
     * @return The rounded value.
     */
    static double round(double value, int decimals){
        return Math.rint(value*Math.pow(10, decimals))/Math.pow(10, decimals);
    }
    
    /**
     * Searches the simulation data for the parameters that have given the lowest sum of square residuals so far.
     * @return The index in the simulation data list that corresponds to the lowest SSR.
     */
    static int findLowestSSR(){
        if(list_l.isEmpty()){
            return -1;
        }
        double[] experimental = true_image;
        int len = list_l.size();
        int low_index = len - 1;
        double[] curr = convertToRegularDouble(list_data.get(low_index));
        double[] residuals = subtract(experimental, curr);
        double low_ssr = sumSquare(residuals);
        for(int i = 0; i < len; i++){
            curr = convertToRegularDouble(list_data.get(i));
            residuals = subtract(experimental, curr);
            double ssr = sumSquare(residuals);
            if(ssr < low_ssr){
                low_index = i;
                low_ssr = ssr;
            }
        }
        return low_index;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // CONSTANTS FOR FITTING ///////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    /**
     * The values for the base-ten logarithm of the k-parameter that are to be tested to find the most suitable initial guess.
     */
    static final double[] TEST_LOG_K = new double[]{-8, -7, -6, -5, -4, -3, -2, -1, 0, 1};
    
    /**
     * The maximum number of iterations to be taken by the fitting method.
     */
    static final int MAX_ITERATIONS = 15;
    
    /**
     * The maximum value of lambda that can be obtained before giving up.
     */
    static final double MAX_LAMBDA = 100.0;
    
    /**
     * The number of decimal points to which the L-parameter is to be considered.
     */
    static final int L_DECIMALS = 2;//in [um]
    
    /**
     * The perturbation that is made to the L-parameter
     */
    static final double L_PERTURB = 0.05;//in [um]
    
    /**
     * The number of decimal points to which the base-ten logarithm of k is to be considered.
     */
    static final int LOGK_DECIMALS = 3;//[k] = [m/s]
    
    /**
     * The perturbation that is to be introduced to the base-ten logarithm of k.
     */
    static final double LOGK_PERTURB = 0.005;//[k] = [m/s]
    
    /**
     * The number of decimal points to which the R-parameter is to be considered.
     */
    static final int R_DECIMALS = 1;//in grid pixels
    
    /**
     * The perturbation that is made to the R-parameter
     */
    static final double R_PERTURB = 0.6;//in grid pixels
    
    /*
    Lists for storing simulation history (to avoid redundant simulations)
    */
    
    /**
     * List for holding the l-parameters that have been simulated
     */
    static LinkedList<Double> list_l;
    
    /**
     * List for holding the logk parameters that have been simulated
     */
    static LinkedList<Double> list_logk;
    
    /**
     * List for holding the total number of grid points that are switched 'on'.
     * Since the possible grid configurations are controlled by erosion or dilation of the initial grid configuration,
     * The pixels will turn 'on' or 'off' in a very specific order.
     */
    static LinkedList<Integer> list_gridsum;
    
    /**
     * List for holding the current results for each secm image simulation.
     */
    static LinkedList<Double[]> list_data;
    
    /*
    Fitting status symbols
    */
    
    /**
     * Indicates that the fitting ended in convergence.
     */
    static final int EXECUTED_OK = 0;
    
    /**
     * Indicates that the fitting ended due to reaching the maximum number of iterations
     */
    static final int MAX_ITERATIONS_REACHED = 1;
    
    /**
     * Indicates that the fitting method ended due to the maximum value of lambda being reached prior to convergence.
     */
    static final int MAX_LAMBDA_REACHED = 2;
    
    
    /*
    Grid-handling methods and fields
    */
    
    /**
     * X-size of pixel {@link #grid} must be odd
     */
    static final int XSIZE = 201;
    
    /**
     * Y-size of pixel {@link #grid} must be odd
     */
    static final int YSIZE = 201;
    
    /**
     * Generates a list of boundary IDs that are to be switched on to position the {@link #grid} at (-x,-y) relative to the secm probe.
     * @param x The x-position of the secm probe relative to the {@link #grid}.
     * @param y The y-position of the secm probe relative to the {@link #grid}.
     * @return The list of boundaries to be used for the PS boundary condition.
     */
    public static int[] getPixels(int x, int y){
        //find the indices of grid_switches that correspond to (x,y) 
        int gridx = x - min_x;
        int gridy = y - min_y;
        
        int halfsizex = (XSIZE-1)/2;
        int halfsizey = (YSIZE-1)/2;
        
        int xmin = 0;
        if(halfsizex - gridx > 0){
            xmin = halfsizex - gridx;
        }
        int xmax = XSIZE;
        if(grid.length + halfsizex - gridx < XSIZE){
            xmax = grid.length + halfsizex - gridx;
        }
        
        int ymin = 0;
        if(halfsizey - gridy > 0){
            ymin = halfsizey - gridy;
        }
        int ymax = YSIZE;
        if(grid[0].length + halfsizey - gridy < YSIZE){
            ymax = grid[0].length + halfsizey - gridy;
        }
        
        LinkedList<Integer> x_list = new LinkedList<Integer>();
        LinkedList<Integer> y_list = new LinkedList<Integer>();
        
        for(int ix = xmin; ix < xmax; ix++){
            int xindex = ix - halfsizex + gridx;
            for(int iy = ymin; iy < ymax; iy++){
                int yindex = iy - halfsizey + gridy;
                if(grid[xindex][yindex] == 1){
                    x_list.add(ix);
                    y_list.add(iy);
                }
            }
        }
        
        ////////////////////////////////////////////////////////////////////////
        //CONVERTING (x,y) to BOUNDARIES MANY OF THESE LITERALS WILL CHANGE WHEN
        //CERTAIN CHANGES ARE MADE TO THE MODEL GEOMETRY!!
        ////////////////////////////////////////////////////////////////////////
        int[] boundaries = new int[x_list.size()];
        for(int i = 0; i < boundaries.length; i++){
            int bound = x_list.removeFirst()*YSIZE + y_list.removeFirst() + 6;
            if(bound > 16487){
                bound += 4;
            }
            if(bound > 18300){
                bound += 4;
            }
            if(bound > 20314){
                bound += 4;
            }
            boundaries[i] = bound;
        }
        
        return boundaries;
    }
    
    /**
     * Obtains the most central sampling point over the reactive feature for obtaining the initial logk guess.
     * This is achieved by scoring each point by the sum of the reciprocal distance squared to each of the 'on' pixels.
     * @return The position of the most central sampling point as an array {x, y}.
     */
    static int[] getCentre(){
        
        int bestx = 0;
        int besty = 0;
        double bestScore = 0;
        
        for(int sx = 0; sx < sample_xs.length; sx++){
            int ix = sample_xs[sx] - min_x;
            for(int sy = 0; sy < sample_ys.length; sy++){
                int iy = sample_ys[sy] - min_y;
                double score = 0;
                for(int x = 0; x < grid.length; x++){
                    for(int y = 0; y < grid[0].length; y++){
                        if((x != ix || y != iy) && grid[x][y] == 1){
                            double distsq = (ix-x)*(ix-x) + (iy-y)*(iy-y);
                            score += 1.0/distsq;
                        }
                    }
                }
                if(score > bestScore){
                    bestx = sample_xs[sx];
                    besty = sample_ys[sy];
                    bestScore = score;
                }
            }
        }
        return new int[]{bestx, besty};
    }
    
    /**
     * The experimental SECM currents
     */
    static double[] true_image;
    
    /**
     * The experimental x-coordinates in meters
     */
    static double[] physical_xs;
    
    /**
     * The experimental y-coordinates in meters
     */
    static double[] physical_ys;
    
    /**
     * The grid_switches x-indexes
     */
    static int[] sample_xs;
    
    /**
     * The grid_switches y-indexes
     */
    static int[] sample_ys;
    
    /**
     * The pixel grid_switches
     */
    static int[][] grid;
    
    /**
     * the minimum grid_switches x-index
     */
    static int min_x;
    
    /**
     * the minimum grid_switches y-index
     */
    static int min_y;
    
    /**
     * The filepath for the instruction file.
     */
    static String true_file;
    
    /*
    IO-handling methods and fields
    */
    
    /**
     * Reads-in the current data from the data file that is produced by the simulation.
     * @return The currents at each point in the same order as {@link #sample_xs}, {@link #sample_ys} and {@link #true_image}.
     * COMSOL prefixes comment lines with '%', so lines starting with '%' in the file are ignored.
     * @throws FileNotFoundException 
     */
    static double[] readData() throws FileNotFoundException{
        File f = new File("data.txt");
        Scanner s = new Scanner(f);
        int count = 0;
        while(s.hasNextLine()){
            String ln = s.nextLine();
            if(!ln.startsWith("%")){
                count ++;
            }
        }
        s.close();
        double[] data = new double[count];
        s = new Scanner(f);
        count = 0;
        while(s.hasNextLine()){
            String ln = s.nextLine();
            if(!ln.startsWith("%")){
                String[] linedata = ln.trim().split("\\s+");//split up the columns
                data[count] = Double.parseDouble(linedata[linedata.length - 1]);//take the data from the righmost column
                count ++;
            }
        }
        return data;
    }
    
    /**
     * Erases the data file produced by the simulation, writing a '%' to the file.
     * @throws FileNotFoundException 
     */
    static void eraseDataFile() throws FileNotFoundException{
        File f = new File("data.txt");
        PrintWriter pw = new PrintWriter(f);
            pw.print("%");
        pw.close();
    }
    
    /**
     * Logs iteration details into {@link #LOGFILE}.
     * @param iteration_num The number of the iteration.
     * @param diagonal The diagonal of the DTD matrix
     * @param lambda The lambda value
     * @param labels The labels for the parameters
     * @param params The values of the parameters
     * @param ssr The sum of square residuals for the parameters.
     * @param accepted Whether or not the parameters had an acceptable SSR.
     * @throws FileNotFoundException
     * @throws IOException 
     */
    static void logIteration(int iteration_num, double[] diagonal, double lambda, String[] labels, double[] params, double ssr, boolean accepted) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(new FileWriter(f, true));
            pw.append("\nIteration: " + iteration_num);
            pw.append("\n\t Diagonal: " + diagonal[0]);
            for(int i = 1; i < labels.length; i++){
                pw.append("," + diagonal[i]);
            }
        pw.close();
        logLambda(lambda, labels, params, ssr, accepted);
    }
    
    /**
     * Logs details regarding an attempt with a new lambda to {@link #LOGFILE}.
     * @param lambda The lambda value being used.
     * @param labels The parameter labels.
     * @param params The parameter values.
     * @param ssr The sum of square residuals.
     * @param accepted Whether or not the lambda attempt had an acceptable SSR.
     * @throws FileNotFoundException
     * @throws IOException 
     */
    static void logLambda(double lambda, String[] labels, double[] params, double ssr, boolean accepted) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(new FileWriter(f, true));
            pw.append("\n\tLAMBDA: " + lambda);
            for(int i = 0; i < labels.length; i++){
                pw.append("\n\t\t" + labels[i] + ": " + params[i]);
            }
            pw.append("\n\tSum. Square Residuals: " + ssr);
            if(accepted){
                pw.append(" (Accepted)");
            }
            else{
                pw.append(" (Rejected)");
            }
        pw.close();
    }
    
    /**
     * Logs the SSR for the initial parameters to {@link #LOGFILE}.
     * @param fname The filename containing the experimental data that is being fit. 
     * (Note that this parameter is not for the name of the file that this logging information is being written into)
     * @param labels The labels for the parameters.
     * @param params The values of the parameters.
     * @param ssr The sum of square residuals when using the parameters.
     * @throws FileNotFoundException
     * @throws IOException 
     */
    static void logInitialGuesses(String fname, String[] labels, double[] params, double ssr) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(f);
            pw.append("\nInitial guesses:");
            pw.append("\n\tFile: " + fname);
            for(int i = 0; i < labels.length; i++){
                pw.append("\n\t\t" + labels[i] + ": " + params[i]);
            }
            pw.append("\n\tSum. Square Residuals: " + ssr);
        pw.close();
    }
    
    /**
     * Logs the end result for the fitting to {@link #LOGFILE}.
     * @param condition The reason for the fit ending.
     * @throws FileNotFoundException
     * @throws IOException 
     * @see #EXECUTED_OK
     * @see #MAX_ITERATIONS_REACHED
     * @see #MAX_LAMBDA_REACHED
     */
    static void logEndCondition(int condition) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(new FileWriter(f, true));
            switch(condition){
                case EXECUTED_OK:
                    pw.append("\nPROCESS CONVERGED.");
                case MAX_ITERATIONS_REACHED:
                    pw.append("\nPROCESS STOPPED PREMATURELY AFTER " + MAX_ITERATIONS + " ITERATIONS.");
                case MAX_LAMBDA_REACHED:
                    pw.append("\nPROCESS STOPPED DUE TO MAXIMUM LAMBDA BEING REACHED.");
            }
        pw.close();
    }
    
    /**
     * Produces a copy of the file at original_filepath, replacing the currents at the sampled points with the currents in {@link #true_image}.
     * This copy will be found at new_filepath.
     * This method is intended to be used to easily generate simulated secm images for testing purposes.
     * @param original_filepath The path to a correctly formatted experimental control file.
     * @param new_filepath The path to the file that is to be created/written-to.
     * @throws FileNotFoundException
     * @throws IOException 
     */
    static void writeSECMInfo(String original_filepath, String new_filepath) throws FileNotFoundException, IOException{
        File originalfile = new File(original_filepath);
        File newfile = new File(new_filepath);
        if(original_filepath.equals(new_filepath)){
            throw new IOException("Cannot read and write to the same file");
        }
        Scanner s = new Scanner(originalfile);
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(newfile)));
        int index = 0;
        String sep = ",";
        
        if(s.hasNextLine()){
            String fl = s.nextLine();
            pw.print(fl);
            if(fl.equalsIgnoreCase("##ENCODING: csv")){
                sep = ",";
            }
            else if(fl.equalsIgnoreCase("##ENCODING: tsv")){
                sep = "\t";
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
        }
        else{
            throw new FileNotFoundException("File incorrectly formatted.");
        }
        
        while(s.hasNextLine()){
            String line = s.nextLine();
            if(!line.startsWith("#")){
                String[] linesplit = line.split(sep);
                int x = Integer.parseInt(linesplit[0]);
                int y = Integer.parseInt(linesplit[1]);
                if(x == sample_xs[index] && y == sample_ys[index]){
                    //copy-over all of the line save the current
                    pw.print("\n" + linesplit[0] + sep + linesplit[1] + sep + linesplit[2] + sep + linesplit[3] + sep + linesplit[4] + sep + true_image[index]);
                    index ++;
                    if(index >= sample_xs.length){
                        index = sample_xs.length - 1;
                    }
                }
                else{
                    pw.print("\n" + line);//copy-over non-sampled lines
                }
            }
            else{
                pw.print("\n" + line);//copy-over commented lines
            }
        }
        
        s.close();
        pw.close();
    }
    
    /**
     * Reads through the control file, pulling out information for the reactivity grid, the sampled points and the experimentally observed secm currents.
     * This method instantiates:
     * <ul>
     * <li>{@link #grid}</li>
     * <li>{@link #sample_xs}</li>
     * <li>{@link #sample_ys}</li>
     * <li>{@link #physical_xs}</li>
     * <li>{@link #physical_ys}</li>
     * <li>{@link #min_x}</li>
     * <li>{@link #min_y}</li>
     * </ul>
     * @param filepath The path to the control file.
     * @throws FileNotFoundException 
     */
    static void readSECMInfo(String filepath) throws FileNotFoundException{
        File f = new File(filepath);
        true_file = filepath;
        int minx = 0;
        int miny = 0;
        int max_x = 0;
        int max_y = 0;
        double[] trueimage;
        double[] physicalxs;
        double[] physicalys;
        int[] samplexs;
        int[] sampleys;
        String sep = ",";
        
        Scanner s = new Scanner(f);
            if(s.hasNextLine()){
                String fl = s.nextLine();
                if(fl.equalsIgnoreCase("##ENCODING: csv")){
                    sep = ",";
                }
                else if(fl.equalsIgnoreCase("##ENCODING: tsv")){
                    sep = "\t";
                }
                else{
                    throw new FileNotFoundException("File incorrectly formatted.");
                }
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            
            //X header
            int xstart = 0;
            int xstep = 0;
            int xnum = 0;
            if(s.hasNextLine()){
                String fl = s.nextLine();
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            if(s.hasNextLine()){
                String fl = s.nextLine();
                String[] tokens = fl.substring(1).trim().split(",");
                xstart = Integer.parseInt(tokens[0]);
                xstep = Integer.parseInt(tokens[1]);
                xnum = Integer.parseInt(tokens[2]);
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            
            //Y header
            int ystart = 0;
            int ystep = 0;
            int ynum = 0;
            if(s.hasNextLine()){
                String fl = s.nextLine();
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            if(s.hasNextLine()){
                String fl = s.nextLine();
                String[] tokens = fl.substring(1).trim().split(",");
                ystart = Integer.parseInt(tokens[0]);
                ystep = Integer.parseInt(tokens[1]);
                ynum = Integer.parseInt(tokens[2]);
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
            
            //read through the main data file to pull out trueimage, physicalxs, physicalys, samplexs, sampleys, minx, miny
            minx = xstart;
            max_x = xstart;
            miny = ystart;
            max_y = ystart;
            
            int asize = (xnum)*(ynum);
            int present_index = 0;
            
            trueimage = new double[asize];
            physicalxs = new double[asize];
            physicalys = new double[asize];
            samplexs = new int[asize];
            sampleys = new int[asize];
            
            while(s.hasNextLine()){
                String line = s.nextLine();
                if(!line.startsWith("#")){
                    String[] linesplit = line.split(sep);
                    int x = Integer.parseInt(linesplit[0]);
                    int y = Integer.parseInt(linesplit[1]);
                    double px = Double.parseDouble(linesplit[3]);
                    double py = Double.parseDouble(linesplit[4]);
                    double curr = Double.parseDouble(linesplit[5]);
                    
                    if(x < minx){
                        minx = x;
                    }
                    else if(x > max_x){
                        max_x = x;
                    }
                    
                    if(y < miny){
                        miny = y;
                    }
                    else if(y > max_y){
                        max_y = y;
                    }
                    
                    //check if x,y is one of the points to be sampled
                    boolean xvalid = (x >= xstart) && (x < xstart + xnum*xstep) && ((x-xstart)%xstep == 0);
                    boolean yvalid = (y >= ystart) && (y < ystart + ynum*ystep) && ((y-ystart)%ystep == 0);
                    if(xvalid && yvalid && present_index < asize){
                        trueimage[present_index] = curr;
                        physicalxs[present_index] = px;
                        physicalys[present_index] = py;
                        samplexs[present_index] = x;
                        sampleys[present_index] = y;
                        present_index ++;
                    }
                    
                }
            }
        s.close();
        //second read of the file to get grid_switches[][]
        int[][] grid_switches = new int[max_x - minx + 1][max_y - miny + 1];
        s = new Scanner(f);
            while(s.hasNextLine()){
                String line = s.nextLine();
                if(!line.startsWith("#")){
                    String[] linesplit = line.split(sep);
                    int x = Integer.parseInt(linesplit[0]);
                    int y = Integer.parseInt(linesplit[1]);
                    int rs = Integer.parseInt(linesplit[2]);
                    
                    grid_switches[x-minx][y-miny] = rs;
                    
                }
            }
        s.close();
        
        true_image = trueimage;
        physical_xs = physicalxs;
        physical_ys = physicalys;
        sample_xs = samplexs;
        sample_ys = sampleys;
        min_x = minx;
        min_y = miny;
        grid = grid_switches;
        
    }
    
    /**
     * Writes the current, derivatives and residuals to a file for plotting purposes.
     * @param fname The name of the file to be written to.
     * @param x The x-positions of the sampled points.
     * @param y The y-positions of the sampled points.
     * @param current The currents.
     * @param dl The Partial derivatives of current with respect to L.
     * @param dlogk The partial derivatives of current with respect to logk.
     * @param residual The residuals.
     * @throws FileNotFoundException 
     */
    static void writeIteration(String fname, double[] x, double[] y, double[] current, double[] dl, double[] dlogk, double[] dr, double[] residual) throws FileNotFoundException{
        File f = new File(fname);
        PrintWriter pw = new PrintWriter(f);
            pw.print("#x [m], y [m], i [A], dL [A], dlog10k [A], residual [A]");
            for(int i = 0; i < x.length; i ++){
                pw.print("\n" + x[i] + "," + y[i] + "," + current[i] + "," + dl[i] + "," + dlogk[i] + "," + residual[i]);
            }
        pw.close();
    }

    /**
     * Writes the logk-current curve data to a file.
     * @param fname The file to be written to.
     * @param current The currents from the logk simulations.
     * @throws FileNotFoundException 
     */
    static void writeKFit(String fname, double[] current) throws FileNotFoundException{
        File f = new File(fname);
        PrintWriter pw = new PrintWriter(f);
            pw.print("#log10(k/1[m/s]), i [A]");
            for(int i = 0; i < current.length; i ++){
                pw.print("\n" + TEST_LOG_K[i] + "," + current[i]);
            }
        pw.close();
    }
    
    static void writeReactivityFile(int[][] de_grid, double reactivity) throws IOException, FileNotFoundException{
        double[][] reactivitymap = new double[grid.length][grid[0].length];
        for(int x = 0; x < grid.length; x++){
            for(int y = 0; y < grid[0].length; y++){
                double pv = de_grid[x][y];
                reactivitymap[x][y] = reactivity*pv;
            }
        }
        File originalfile = new File(true_file);
        File newfile = new File(reactivity_mapfile);
        if(true_file.equals(reactivity_mapfile)){
            throw new IOException("Cannot read and write to the same file");
        }
        newfile.createNewFile();
        Scanner s = new Scanner(originalfile);
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(newfile)));
        String sep = ",";
        
        if(s.hasNextLine()){
            String fl = s.nextLine();
            pw.print(fl);
            if(fl.equalsIgnoreCase("##ENCODING: csv")){
                sep = ",";
            }
            else if(fl.equalsIgnoreCase("##ENCODING: tsv")){
                sep = "\t";
            }
            else{
                throw new FileNotFoundException("File incorrectly formatted.");
            }
        }
        else{
            throw new FileNotFoundException("File incorrectly formatted.");
        }
        
        boolean first = true;
        while(s.hasNextLine()){
            String line = s.nextLine();
            if(!line.startsWith("#")){
                String[] linesplit = line.split(sep);
                int x = Integer.parseInt(linesplit[0]);
                int y = Integer.parseInt(linesplit[1]);
                //copy-over the relevant parts of the file
                if(!first){
                    pw.print("\n" + linesplit[3] + "," + linesplit[4] + "," + reactivitymap[x][y]);
                }
                else{
                    pw.print(linesplit[3] + "," + linesplit[4] + "," + reactivitymap[x][y]);
                }
            }
        }
        
        s.close();
        pw.close();
    }
    
    /**
     * The file to which logging information will be written.
     */
    static final String LOGFILE = "fit.log";
    
    /**
     * The file to which logk-current curves will be written.
     */
    static final String KLOGFILE = "k-curve.csv";
    
    /**
     * The file to which the reactivity map will be written and read from.
     */
    static String reactivity_mapfile = "func.csv";
    
    /*
    Linear algebra-handling methods and fields
     */
    /**
     * Inverts the given matrix
     * @param a the matrix to be inverted
     * @return the inverse of a, a<sup>-1</sup>. a<sup>-1</sup>*a = a*a<sup>-1</sup>=I
     */
    static double[][] invert(double[][] a) throws NumberFormatException{
        double det = determinant(a);
        if(!Double.isFinite(det) || det == 0.0){
            throw new NumberFormatException("Matrix a is singular.");
        }
        else{
            double[][] cofactors = new double[a.length][a[0].length];
            for(int r = 0; r < a.length; r ++){
                for(int c = 0; c < a.length; c ++){
                    cofactors[r][c] = cofactor(a, r, c);
                }
            }
            return multiply(transpose(cofactors), 1.0/det);
        }
    }
    
    /**
     * Multiplies matrices a*b
     * @param a The left matrix. a[r][c] (r = row number; c = column number)
     * @param b The right matrix. b[r][c] (r = row number; c = column number)
     * @return The product matrix. p[r][c] (r = row number; c = column number)
     * @throws NumberFormatException 
     */
    static double[][] multiply(double[][] a, double[][] b) throws NumberFormatException{
        double[][] prod = new double[a.length][b[0].length];
        if(a[0].length != b.length){
            throw new NumberFormatException("columnspace of a, " + a[0].length + ", and rowspace of b, " + b.length + ", must be the same size.");
        }
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < b[0].length; c++){
                prod[r][c] = 0;
                for(int i = 0; i < b.length; i++){
                    prod[r][c] += a[r][i]*b[i][c];
                }
            }
        }
        return prod;
    }
    
    /**
     * Operates the matrix a on the vector v. a*v
     * @param a The left matrix. a[r][c] (r = row number; c = column number)
     * @param v The right vector. v[r] (r = row number)
     * @return The resulting product vector. p[r] (r = row number)
     * @throws NumberFormatException 
     */
    static double[] multiply(double[][]a, double[] v) throws NumberFormatException{
        if(a[0].length != v.length){
            throw new NumberFormatException("columnspace of a, " + a[0].length + ", and rowspace of v, " + v.length + ", must be the same size.");
        }
        double[] prod = new double[a.length];
        for(int r = 0; r < a.length; r++){
            prod[r] = 0;
            for(int c = 0; c < v.length; c++){
                prod[r] += a[r][c]*v[c];
            }
        }
        return prod;
    }
    
    /**
     * Multiplies a matrix a with a scalar s.
     * @param a The a matrix
     * @param s The s scalar
     * @return The product matrix
     */
    static double[][] multiply(double[][] a, double s){
        double[][] mul = new double[a.length][a[0].length];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                mul[r][c] = s*a[r][c];
            }
        }
        return mul;
    }
    
    /**
     * Multiplies a vector v with a scalar s.
     * @param v The a vector
     * @param s The s scalar
     * @return The product vector
     */
    static double[] multiply(double[] v, double s){
        double[] mul = new double[v.length];
        for(int i = 0; i < mul.length; i++){
            mul[i] = v[i]*s;
        }
        return mul;
    }
    
    /**
     * Adds two matrices a+b
     * @param a The left matrix. a[r][c] (r = row number; c = column number)
     * @param b The right matrix. b[r][c] (r = row number; c = column number)
     * @return The sum matrix. s[r][c] (r = row number; c = column number)
     * @throws NumberFormatException 
     */
    static double[][] add(double[][] a, double[][] b) throws NumberFormatException{
        if(a.length != b.length || a[0].length != b[0].length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[][] sum = new double[a.length][a[0].length];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                sum[r][c] = a[r][c] + b[r][c];
            }
        }
        return sum;
    }
    
    /**
     * Subtracts two matrices a-b
     * @param a The left matrix. a[r][c] (r = row number; c = column number)
     * @param b The right matrix. b[r][c] (r = row number; c = column number)
     * @return The difference matrix. d[r][c] (r = row number; c = column number)
     * @throws NumberFormatException 
     */
    static double[][] subtract(double[][] a, double[][] b) throws NumberFormatException{
        if(a.length != b.length || a[0].length != b[0].length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[][] dif = new double[a.length][a[0].length];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                dif[r][c] = a[r][c] - b[r][c];
            }
        }
        return dif;
    }
    
    /**
     * Adds two vectors a+b
     * @param a The left vector. a[i]
     * @param b The right vector. b[i]
     * @return The sum vector. s[i]
     * @throws NumberFormatException 
     */
    static double[] add(double[] a, double[] b) throws NumberFormatException{
        if(a.length != b.length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[] sum = new double[a.length];
        for(int r = 0; r < a.length; r++){
            sum[r] = a[r] + b[r];
        }
        return sum;
    }
    
    /**
     * Subtracts two vectors a-b
     * @param a The left vector. a[i]
     * @param b The right vector. b[i]
     * @return The difference vector. d[i]
     * @throws NumberFormatException 
     */
    static double[] subtract(double[] a, double[] b) throws NumberFormatException{
        if(a.length != b.length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[] dif = new double[a.length];
        for(int r = 0; r < a.length; r++){
            dif[r] = a[r] - b[r];
        }
        return dif;
    }
    
    /**
     * Swaps the column and rowspace of a.
     * @param a The matrix. a[r][c] (r = row number; c = column number) 
     * @return The transposed matrix. t[r][c] (r = row number; c = column number) 
     */
    static double[][] transpose(double[][] a){
        double[][] transp = new double[a[0].length][a.length];
        for(int r = 0; r < a.length; r ++){
            for(int c = 0; c < a[0].length; c ++){
                transp[c][r] = a[r][c];
            }
        }
        return transp;
    }
    
    /**
     * appends v to the right of a
     * @param a the matrix to which v is to be appended
     * @param v the vector to append
     * @return a copy of a with column v added to the right
     * @throws NumberFormatException 
     */
    static double[][] appendColumn(double[][] a, double[] v) throws NumberFormatException{
        if(a.length != v.length){
            throw new NumberFormatException("Rowspace of a, " + a.length + ", and v, " + v.length + ", must be the same");
        }
        double[][] augment = new double[a.length][a[0].length + 1];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                augment[r][c] = a[r][c];
            }
            augment[r][a[0].length] = v[r];
        }
        return augment;
    }
    
    /**
     * appends v to the right of a
     * @param a the vector to which v is to be appended
     * @param v the vector to append
     * @return a copy of a with column v added to the right
     * @throws NumberFormatException 
     */
    static double[][] appendColumn(double[] a, double[] v) throws NumberFormatException{
        if(a.length != v.length){
            throw new NumberFormatException("Rowspace of a, " + a.length + ", and v, " + v.length + ", must be the same");
        }
        double[][] augment = new double[a.length][2];
        for(int r = 0; r < a.length; r++){
            augment[r][0] = a[r];
            augment[r][1] = v[r];
        }
        return augment;
    }
    
    /**
     * Returns the nxn identity matrix, I.
     * @param n the dimension of the identity matrix
     * @return the nxn identity matrix
     */
    static double[][] identity(int n){
        double[][] ident = new double[n][n];
        for(int r = 0; r < n; r++){
            for(int c = 0; c < n; c++){
                if(r==c){
                    ident[r][c] = 1;
                }
                else{
                    ident[r][c] = 0;
                }
            }
        }
        return ident;
    }
    
    /**
     * Computes the determinant of a matrix, a.
     * @param a The matrix a[r][c] (r = row number; c = column number)
     * @return The determinant of a.
     * @throws NumberFormatException 
     * @see #inner_determinant(double[][]) 
     */
    static double determinant(double[][] a) throws NumberFormatException{
        if(a.length == a[0].length){
            return inner_determinant(a);
        }
        else{
            throw new NumberFormatException("Determinants may only be taken of square matrices.");
        }
    }
    
    /**
     * Please use this method's wrapper, {@link #determinant(double[][])} as this function doe not perform any input checking.
     * @param a
     * @return 
     */
    static double inner_determinant(double[][] a){
        if(a.length == 2){//2x2 matrix
            return a[0][0]*a[1][1] - a[0][1]*a[1][0];
        }
        else if(a.length == 1){//1x1 matrix
            return a[0][0];
        }
        else{
            double sum = 0;
            for(int c = 0; c < a[0].length; c++){
                if(c%2 == 0){
                    sum += a[0][c]*inner_determinant(minor(a,0,c));
                }
                else{
                    sum -= a[0][c]*inner_determinant(minor(a,0,c));
                }
            }
            return sum;
        }
    }
    
    /**
     * Generates the cofactor for the element of matrix a at r,c
     * @param a The matrix a[r][c] (r = row number; c = column number)
     * @param r The row index
     * @param c The column index
     * @return The cofactor for constructing the matrix of cofactors for a.
     * @throws NumberFormatException 
     */
    static double cofactor(double[][] a, int r, int c) throws NumberFormatException{
        double sign = -1.0;
        if ((r+c)%2 == 0){
            sign = 1.0;
        }
        return sign*determinant(minor(a, r, c));
    }
    
    /**
     * Returns a copy of a that omits row r and column c.
     * @param a the original matrix
     * @param r the row to be omitted
     * @param c the column to be omitted
     * @return
     */
    static double[][] minor(double[][] a, int r, int c){
        double[][] red = new double[a.length - 1][a[0].length - 1];
        
        for(int ri = 0; ri < a.length - 1; ri++){
            int ru = ri;//the r to be used
            if(ri >= r){
                ru ++;//if we are at or beyond the row to be omitted, increment the r to be read from a
            }
            for(int ci = 0; ci < a[0].length - 1; ci++){
                int cu = ci;//the c to be used
                if(ci >= c){
                    cu ++;//if we are at or beyond the column to be omitted, increment the c to be read from a
                }
                red[ri][ci] = a[ru][cu];
            }
        }
        
        return red;
    }
    
    /* 
    Erosion and dilation
    */
    /**
     * Applies a dilation or erosion by |amount|
     * If amount is negative, an erosion will be carried-out, otherwise, a dilation will.
     * @param amount the amount by which the grid is to be dilated or eroded in grid-space (i.e. the distance between adjacent grid elements is 1.
     */
    static int[][] applyDilationErosion(double amount){
        if(amount < 0){
            return erodeGrid(-amount);
        }
        else{
            return dilateGrid(amount);
        }
    }
    /**
     * dilates grid by amount, storing the result to edited_grid
     * @param amount the amount by which the grid is to be dilated in grid-space (i.e. the distance between adjacent grid elements is 1.
     */
    static int[][] dilateGrid(double amount){
        int upperbound = (int)Math.ceil(amount);
        double amountsq = (amount + 0.5) * (amount + 0.5);// adding 0.5 due to the dilation being from the center of the pixel
        //(re)-initialize the edited_grid
        int[][] edited_grid = new int[grid.length][grid[0].length];
        for(int x = 0; x < grid.length; x++){
            for(int y = 0; y < grid[0].length; y++){
                edited_grid[x][y] = grid[x][y];
            }
        }
        if(amount >= 0.2){ // don't bother with the dilation if the amount is too small to change anything
            //iterate through every point
            for(int x = 0; x < grid.length; x++){
                for(int y = 0; y < grid[0].length; y++){
                    if(grid[x][y] == 0){
                        //if the value of the grid point can be changed, iterate through its neighboring points
                        int[][] subdivisions = new int[5][5];
                        for(int xx = - upperbound; xx <= upperbound; xx++){
                            for(int yy = - upperbound; yy <= upperbound; yy++){
                                if(getGrid(x+xx, y+yy) == 1){
                                    /*if the point at (x+xx,y+yy) is different from 
                                    the one at (x,y), subdivide the point at (x,y)
                                    into a 5x5 grid and see how many of these grid
                                    points are within the erosion or dilation distance
                                    */
                                    for(int xxx = 0; xxx < 5; xxx ++){
                                        for(int yyy = 0; yyy < 5; yyy ++){
                                            double x_subdist = ((double)xxx)*0.2 - 0.4;
                                            double y_subdist = ((double)yyy)*0.2 - 0.4;
                                            double xdist = ((double)xx) + x_subdist;
                                            double ydist = ((double)yy) + y_subdist;
                                            double distancesq = xdist*xdist + ydist*ydist;
                                            if(distancesq < amountsq){
                                                subdivisions[xxx][yyy] = 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        /*compute the number of subdivisions at (x,y) that should change.
                        If they are in the majority, flip the pixel.
                        */
                        int sum = 0;
                        for(int i = 0; i < 5; i++){
                            for(int ii = 0; ii < 5; ii++){
                                sum += subdivisions[i][ii];
                            }
                        }
                        if(sum > 12){
                            edited_grid[x][y] = 1;
                        }
                    }
                }
            }
        }
        return edited_grid;
    }
    
    /**
     * erodes grid by amount, storing the result to edited_grid
     * @param amount the amount by which the grid is to be eroded in grid-space (i.e. the distance between adjacent grid elements is 1.
     */
    static int[][] erodeGrid(double amount){
        int upperbound = (int)Math.ceil(amount);
        double amountsq = (amount + 0.5) * (amount + 0.5);// adding 0.5 due to the erosion being from the center of the pixel
        //(re)-initialize the edited_grid
        int[][] edited_grid = new int[grid.length][grid[0].length];
        for(int x = 0; x < grid.length; x++){
            for(int y = 0; y < grid[0].length; y++){
                edited_grid[x][y] = grid[x][y];
            }
        }
        if(amount >= 0.2){// don't bother with the erosion if the amount is too small to change anything
            //iterate through every point
            for(int x = 0; x < grid.length; x++){
                for(int y = 0; y < grid[0].length; y++){
                    if(grid[x][y] == 1){
                        //if the value of the grid point can be changed, iterate through its neighboring points
                        int[][] subdivisions = new int[5][5];
                        for(int xx = - upperbound; xx <= upperbound; xx++){
                            for(int yy = - upperbound; yy <= upperbound; yy++){
                                if(getGrid(x+xx, y+yy) == 0){
                                    /*if the point at (x+xx,y+yy) is different from 
                                    the one at (x,y), subdivide the point at (x,y)
                                    into a 5x5 grid and see how many of these grid
                                    points are within the erosion or dilation distance
                                    */
                                    for(int xxx = 0; xxx < 5; xxx ++){
                                        for(int yyy = 0; yyy < 5; yyy ++){
                                            double x_subdist = ((double)xxx)*0.2 - 0.4;
                                            double y_subdist = ((double)yyy)*0.2 - 0.4;
                                            double xdist = ((double)xx) + x_subdist;
                                            double ydist = ((double)yy) + y_subdist;
                                            double distancesq = xdist*xdist + ydist*ydist;
                                            if(distancesq < amountsq){
                                                subdivisions[xxx][yyy] = 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        /*compute the number of subdivisions at (x,y) that should change.
                        If they are in the majority, flip the pixel.
                        */
                        int sum = 0;
                        for(int i = 0; i < 5; i++){
                            for(int ii = 0; ii < 5; ii++){
                                sum += subdivisions[i][ii];
                            }
                        }
                        if(sum > 12){
                            edited_grid[x][y] = 0;
                        }
                        
                    }
                }
            }
        }
        return edited_grid;
    }
    
    /**
     * Counts the number of pixels that are switched on in a given grid
     * @param de_grid A grid where '1' indicates the pixel is on and '0' indicates the pixel is off.
     * @return 
     */
    static int getSum(int[][] de_grid){
        int runningcount = 0;
        for(int x = 0; x < de_grid.length; x++){
            for(int y = 0; y < de_grid[0].length; y++){
                runningcount += de_grid[x][y];
            }
        }
        return runningcount;
    }
    
    /**
     * Gets the value of the {@link #grid} at (x,y). If (x,y) is outside of the grid, 0 will be returned instead.
     * @param x The x-index of interest.
     * @param y The y-index of interest.
     * @return The value of the grid[x][y] if x and y are in bounds. Zero otherwise.
     */
    static int getGrid(int x, int y){
        if(x > 0 && y > 0 && x < grid.length && y < grid[0].length){
            return grid[x][y];
        }
        else{
            return 0;
        }
    }
    
    /*
    TESTING AND DEBUGGING PURPOSES ONLY
    */
    static void exportEditedGrid(String filename, int[][] edited_grid) throws IOException{
        File f = new File(filename);
        f.createNewFile();
        PrintWriter pw = new PrintWriter(f);
        pw.print("#x,y,switch");
        for(int x = 0; x < edited_grid.length; x++){
            for(int y = 0; y < edited_grid[0].length; y++){
                pw.print("\n" + x + "," + y + "," + edited_grid[x][y]);
            }
        }
        pw.close();
    }
}
