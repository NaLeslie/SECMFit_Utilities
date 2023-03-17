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
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Scanner;

/**
 *
 * @author Nathaniel
 */
public class Testing {

    /**
     * Stand-in for auto-generated method
     * @param z
     * @param logk
     * @param testk
     * @return 
     */
    static Model run(double l, double logk, boolean testk){
        return new Model();
    }
    static Model run2(Model model, boolean testk){
        return new Model();
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws NumberFormatException, FileNotFoundException {
        list_l = new LinkedList<Double>();
        list_logk = new LinkedList<Double>();
        list_data = new LinkedList<Double[]>();
        try{
            readSECMInfo("Fitfile.csv");
            runFit("yep", 1.0, -5, false);
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
    /*
    Fitting methods and fields
    */
    
    static int runFit(String filename, double firstl, double firstlogk, boolean verbose) throws FileNotFoundException, NumberFormatException, IOException{
        double[] experimental = true_image;
        double lambda = 0;
        double ssr = 0;
        double last_ssr;
        double l = firstl;
        double logk = firstlogk;
        double last_l;
        double last_logk;
        //first iteration
        double[] curr = runModel(firstl, firstlogk);
        double[] residuals = subtract(experimental, curr);
        
        double[] curr_dl = runModel(firstl + L_PERTURB, firstlogk);
        double[] dl = multiply(subtract(curr_dl, curr), 1.0/L_PERTURB);
        
        double[] curr_dk = runModel(firstl, firstlogk + LOGK_PERTURB);
        double[] dk = multiply(subtract(curr_dk, curr), 1.0/LOGK_PERTURB);
        
        double[][] J = appendColumn(dl, dk);
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
        last_ssr = ssr;
        
        last_l = l;
        last_logk = logk;
        l = last_l + round(delta_c[0], L_DECIMALS);
        logk = last_logk + round(delta_c[1], LOGK_DECIMALS);
        boolean converged = false;
        int iterations = 1;
        
        //subsequent iterations
        while(!converged && iterations <= MAX_ITERATIONS){
            iterations ++;
            //First lambda
            curr = runModel(l, logk);
            residuals = subtract(experimental, curr);
            ssr = sumSquare(residuals);
            boolean lambda_ok = ssr < last_ssr;

            while(!lambda_ok && lambda <= MAX_LAMBDA){
                lambda = nextLambda(lambda);
                lam_DTD = multiply(DTD, lambda);
                JTJinv = invert(add(JTJ, lam_DTD));
                delta_c = multiply(multiply(JTJinv, JT), residuals);
                l = last_l + round(delta_c[0], L_DECIMALS);
                logk = last_logk + round(delta_c[1], LOGK_DECIMALS);
                curr = runModel(l, logk);
                residuals = subtract(experimental, curr);
                ssr = sumSquare(residuals);
                lambda_ok = ssr < last_ssr;
                logLambda(lambda, new String[]{"L", "log10k"}, new double[]{l, logk}, ssr, lambda_ok);
                if(delta_c[0] < Math.pow(10, L_DECIMALS) && delta_c[1] < Math.pow(10, LOGK_DECIMALS)){
                    lambda_ok = true;
                    converged = true;
                    return EXECUTED_OK;
                }
            }
            if(!lambda_ok){
                return MAX_LAMBDA_REACHED;
            }
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
            
            last_ssr = ssr;
            last_l = l;
            last_logk = logk;
            
            int query = checkList(l - L_PERTURB, logk);
            if(query == -1){
                curr_dl = runModel(l + L_PERTURB, logk);
                dl = multiply(subtract(curr_dl, curr), L_PERTURB);
            }
            else{
                curr_dl = convertToRegularDouble(list_data.get(query));
                //NOTE: this curr_dl is effectively simulated as being perturbed negatively
                dl = multiply(subtract(curr_dl, curr), -L_PERTURB);
            }
            

            query = checkList(l, logk - LOGK_PERTURB);
            if(query == -1){
                curr_dk = runModel(l, logk + LOGK_PERTURB);
                dk = multiply(subtract(curr_dk, curr), LOGK_PERTURB);
            }
            else{
                curr_dk = convertToRegularDouble(list_data.get(query));
                //NOTE: this curr_dk is effectively simulated as being perturbed negatively
                dk = multiply(subtract(curr_dk, curr), -LOGK_PERTURB);
            }
            

            J = appendColumn(dl, dk);
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
            if(delta_c[0] < Math.pow(10, L_DECIMALS) && delta_c[1] < Math.pow(10, LOGK_DECIMALS)){
                converged = true;
                return EXECUTED_OK;
            }
        }
        
        return MAX_ITERATIONS_REACHED;
    }
    
    static double[] runModel(double l, double logk) throws FileNotFoundException{
        int index = checkList(l, logk);
        if(index == -1){
            Model temp = run(l, logk, false);
            run2(temp, false);
            double[] data = readData();
            eraseDataFile();
            addToList(l, logk, data);
            return data;
        }
        else{
            return convertToRegularDouble(list_data.get(index));
        }
    }
    
    static int checkList(double l, double logk){
        for(int i = 0; i < list_l.size(); i++){
            boolean leq = precisionEquals(l, list_l.get(i), L_DECIMALS);
            boolean logkeq = precisionEquals(l, list_logk.get(i), LOGK_DECIMALS);
            if(leq && logkeq){
                return i;
            }
        }
        return -1;
    }
    
    static void addToList(double l, double logk, double[] data){
        Double dl = l;
        Double dlogk = logk;
        Double[] ddata = convertToClassDouble(data);
        list_l.add(dl);
        list_logk.add(dlogk);
        list_data.add(ddata);
    }
    
    static Double[] convertToClassDouble(double[] input){
        int len = input.length;
        Double[] output = new Double[len];
        for(int i = 0; i < len; i++){
            output[i] = input[i];
        }
        return output;
    }
    
    static double[] convertToRegularDouble(Double[] input){
        int len = input.length;
        double[] output = new double[len];
        for(int i = 0; i < len; i++){
            output[i] = (double)input[i];
        }
        return output;
    }
    
    static boolean precisionEquals(double d1, Double cd2, int decimals){
        double d2 = (double)cd2;
        double factor = Math.pow(10, decimals);
        double difference = Math.abs(d1 - d2);
        long rounded_difference = Math.round(difference*factor);
        return rounded_difference == 0;
    }
    
    static double findFirstLogK(double firstl, boolean verbose) throws FileNotFoundException{
        //call the run method to produce a current when the electrode is at distance z and GridData.getCentre() relative to the reactive feature.
        Model model = run(firstl, 1.0, true);
        run2(model, true);
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
    
    static double sumSquare(double[] residuals){
        double sum = 0;
        for(double r : residuals){
            sum += r*r;
        }
        return sum;
    }
    
    static double round(double value, int decimals){
        return Math.rint(value*Math.pow(10, decimals))/Math.pow(10, decimals);
    }
    
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
    
    static final double[] TEST_LOG_K = new double[]{-8, -7, -6, -5, -4, -3, -2, -1, 0, 1};
    
    static final int MAX_ITERATIONS = 15;
    
    static final double MAX_LAMBDA = 100.0;
    
    static final int L_DECIMALS = 2;//in [um]
    
    static final double L_PERTURB = 0.05;//in [um]
    
    static final int LOGK_DECIMALS = 3;//[k] = [m/s]
    
    static final double LOGK_PERTURB = 0.001;//[k] = [m/s]
    
    static final double[] DUMMY = new double[]{0};
    
    /*
    Lists for storing simulation history (to avoid redundant simulations)
    */
    
    static LinkedList<Double> list_l;
    static LinkedList<Double> list_logk;
    static LinkedList<Double[]> list_data;
    
    /*
    Fitting status symbols
    */
    
    static final int EXECUTED_OK = 0;
    static final int MAX_ITERATIONS_REACHED = 1;
    static final int MAX_LAMBDA_REACHED = 2;
    
    
    /*
    Grid-handling methods and fields
    */
    
    static final int XSIZE = 201;//Size of pixel grid_switches must be odd
    static final int YSIZE = 201;//Size of pixel grid_switches must be odd
    
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
     * The pixel grid that gets eroded or dilated relative to grid
     */
    static int[][] edited_grid;
    /**
     * the minimum grid_switches x-index
     */
    static int min_x;
    /**
     * the minimum grid_switches y-index
     */
    static int min_y;
    
    
    /*
    IO-handling methods and fields
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
                data[count] = Double.parseDouble(ln.trim());
                count ++;
            }
        }
        return data;
    }
    
    static void eraseDataFile() throws FileNotFoundException{
        File f = new File("data.txt");
        PrintWriter pw = new PrintWriter(f);
            pw.print("%");
        pw.close();
    }
    
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
    
    static void writeSECMInfo(String original_filepath, String new_filepath) throws FileNotFoundException, IOException{
        File originalfile = new File(original_filepath);
        File newfile = new File(new_filepath);
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
                    pw.print("\n" + linesplit[0] + sep + linesplit[1] + sep + linesplit[2] + sep + linesplit[3] + sep + linesplit[4] + sep + true_image[index]);
                    index ++;
                    if(index >= sample_xs.length){
                        index = sample_xs.length - 1;
                    }
                }
                else{
                    pw.print("\n" + line);
                }
            }
            else{
                pw.print("\n" + line);
            }
        }
        
        s.close();
        pw.close();
    }
    
    static void readSECMInfo(String filepath) throws FileNotFoundException{
        File f = new File(filepath);
        
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
    
    static void writeIteration(String fname, double[] x, double[] y, double[] current, double[] dl, double[] dlogk, double[] residual) throws FileNotFoundException{
        File f = new File(fname);
        PrintWriter pw = new PrintWriter(f);
            pw.print("#x [m], y [m], i [A], dL [A], dlog10k [A], residual [A]");
            for(int i = 0; i < x.length; i ++){
                pw.print("\n" + x[i] + "," + y[i] + "," + current[i] + "," + dl[i] + "," + dlogk[i] + "," + residual[i]);
            }
        pw.close();
    }

    static void writeKFit(String fname, double[] current) throws FileNotFoundException{
        File f = new File(fname);
        PrintWriter pw = new PrintWriter(f);
            pw.print("#log10(k/1[m/s]), i [A]");
            for(int i = 0; i < current.length; i ++){
                pw.print("\n" + TEST_LOG_K[i] + "," + current[i]);
            }
        pw.close();
    }
    
    static final String LOGFILE = "fit.log";
    
    static final String KLOGFILE = "k-curve.csv";
    
    
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
    
    static double[][] multiply(double[][] a, double[][] b) throws NumberFormatException{
        double[][] prod = new double[a.length][b[0].length];
        if(a[0].length != b.length){
            throw new NumberFormatException("columnspace of a, " + a[0].length + ", and rowspace of b, " + b.length + ", must be the same");
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
    
    static double[] multiply(double[][]a, double[] v) throws NumberFormatException{
        if(a[0].length != v.length){
            throw new NumberFormatException("columnspace of a, " + a[0].length + ", and rowspace of v, " + v.length + ", must be the same");
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
    
    static double[][] multiply(double[][] a, double s){
        double[][] mul = new double[a.length][a[0].length];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                mul[r][c] = s*a[r][c];
            }
        }
        return mul;
    }
    
    static double[] multiply(double[] v, double s){
        double[] mul = new double[v.length];
        for(int i = 0; i < mul.length; i++){
            mul[i] = v[i]*s;
        }
        return mul;
    }
    
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
    
    static double determinant(double[][] a) throws NumberFormatException{
        if(a.length == a[0].length){
            return inner_determinant(a);
        }
        else{
            throw new NumberFormatException("Determinants may only be taken of square matrices.");
        }
    }
    
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
     * dilates grid by amount, storing the result to edited_grid
     * @param amount 
     */
    static void dilateGrid(double amount){
        int upperbound = (int)Math.ceil(amount);
        double amountsq = (amount + 0.5) * (amount + 0.5);// adding 0.5 due to the dilation being from the center of the pixel
        for(int x = 0; x < grid.length; x++){
            for(int y = 0; y < grid[0].length; y++){
                edited_grid[x][y] = grid[x][y];
            }
        }
        if(amount >= 0.2){
            for(int x = 0; x < grid.length; x++){
                for(int y = 0; y < grid[0].length; y++){
                    if(grid[x][y] == 0){
                        int[][] subdivisions = new int[5][5];
                        for(int xx = - upperbound; xx <= upperbound; xx++){
                            for(int yy = - upperbound; yy <= upperbound; yy++){
                                if(getGrid(x+xx, y+yy) == 1){
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
                        int sum = 0;
                        for(int i = 0; i < 5; i++){
                            for(int ii = 0; ii < 5; ii++){
                                sum += subdivisions[i][ii];
                            }
                        }
                        if(sum > 12){
                            setGrid(x, y, 1);
                        }
                    }
                }
            }
        }
    }
    
    /**
     * dilates grid by amount, storing the result to edited_grid
     * @param amount 
     */
    static void erodeGrid(double amount){
        int upperbound = (int)Math.ceil(amount);
        double amountsq = (amount + 0.5) * (amount + 0.5);// adding 0.5 due to the dilation being from the center of the pixel
        for(int x = 0; x < grid.length; x++){
            for(int y = 0; y < grid[0].length; y++){
                edited_grid[x][y] = grid[x][y];
            }
        }
        if(amount >= 0.2){
            for(int x = 0; x < grid.length; x++){
                for(int y = 0; y < grid[0].length; y++){
                    if(grid[x][y] == 1){
                        int[][] subdivisions = new int[5][5];
                        for(int xx = - upperbound; xx <= upperbound; xx++){
                            for(int yy = - upperbound; yy <= upperbound; yy++){
                                if(getGrid(x+xx, y+yy) == 0){
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
                        int sum = 0;
                        for(int i = 0; i < 5; i++){
                            for(int ii = 0; ii < 5; ii++){
                                sum += subdivisions[i][ii];
                            }
                        }
                        if(sum > 12){
                            setGrid(x, y, 0);
                        }
                    }
                }
            }
        }
    }
    
    /**
     * 
     * @param x
     * @param y
     * @return 
     */
    static int getGrid(int x, int y){
        if(x > 0 && y > 0 && x < grid.length && y < grid[0].length){
            return grid[x][y];
        }
        else{
            return 0;
        }
    }
    
    /**
     * 
     * @param x
     * @param y
     * @param value 
     */
    static void setGrid(int x, int y, int value){
        if(x > 0 && y > 0 && x < edited_grid.length && y < edited_grid[0].length){
            edited_grid[x][y] = value;
        }
    }
    
}

class Model{
    public Model(){
        property = 0;
    }
    public int property;
}
