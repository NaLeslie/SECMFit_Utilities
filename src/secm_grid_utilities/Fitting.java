/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package secm_grid_utilities;

import java.io.FileNotFoundException;

/**
 * Class holding methods for fitting.
 * @author Nathaniel Leslie
 * Created:  2022-07-12
 * Modified: 2022-07-14
 */
public class Fitting {
    public static int runFit(GridData input, String filename, double firstl, double firstlogk, boolean verbose) throws FileNotFoundException, MatrixException{
        double[] experimental = input.getTrueImage();
        double lambda = 0;
        double ssr = 0;
        double last_ssr;
        double l = firstl;
        double logk = firstlogk;
        double last_l;
        double last_logk;
        //first iteration
        Dummy.run(input, firstl, firstlogk, false);
        double[] curr = IO.readData();
        double[] residuals = MatrixMath.subtract(experimental, curr);
        IO.eraseDataFile();
        
        Dummy.run(input, firstl + L_PERTURB, firstlogk, false);
        double[] curr_dl = IO.readData();
        double[] dl = MatrixMath.multiply(MatrixMath.subtract(curr_dl, curr), L_PERTURB);
        IO.eraseDataFile();
        
        Dummy.run(input, firstl, firstlogk + LOGK_PERTURB, false);
        double[] curr_dk = IO.readData();
        double[] dk = MatrixMath.multiply(MatrixMath.subtract(curr_dk, curr), LOGK_PERTURB);
        IO.eraseDataFile();
        
        double[][] J = MatrixMath.appendColumn(dl, dk);
        double[][] JT = MatrixMath.transpose(J);
        double[][] JTJ = MatrixMath.multiply(JT, J);
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
        double[][] lam_DTD = MatrixMath.multiply(DTD, lambda);
        double[][] JTJinv = MatrixMath.invert(MatrixMath.add(JTJ, lam_DTD));
        double[] delta_c = MatrixMath.multiply(MatrixMath.multiply(JTJinv, JT), residuals);
        last_ssr = ssr;
        ssr = sumSquare(residuals);
        IO.logInitialGuesses(filename, new String[]{"L", "log10k"}, new double[]{l,logk}, ssr);
        if(verbose){
            IO.writeIteration("Iteration_1.txt", input.getPhysicalXs(), input.getPhysicalYs(), curr, dl, dk, residuals);
        }
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
            Dummy.run(input, l, logk, false);
            curr = IO.readData();
            residuals = MatrixMath.subtract(experimental, curr);
            IO.eraseDataFile();
            ssr = sumSquare(residuals);
            boolean lambda_ok = ssr < last_ssr;
            IO.logIteration(iterations, new double[]{DTD[0][0], DTD[1][1]}, lambda, new String[]{"L", "log10k"}, new double[]{l, logk}, ssr, lambda_ok);
            if(verbose){
                IO.writeIteration("Iteration_" + iterations + ".txt", input.getPhysicalXs(), input.getPhysicalYs(), curr, dl, dk, residuals);
            }
            while(!lambda_ok && lambda <= MAX_LAMBDA){
                lambda = nextLambda(lambda);
                lam_DTD = MatrixMath.multiply(DTD, lambda);
                JTJinv = MatrixMath.invert(MatrixMath.add(JTJ, lam_DTD));
                delta_c = MatrixMath.multiply(MatrixMath.multiply(JTJinv, JT), residuals);
                l = last_l + round(delta_c[0], L_DECIMALS);
                logk = last_logk + round(delta_c[1], LOGK_DECIMALS);
                Dummy.run(input, l, logk, false);
                curr = IO.readData();
                residuals = MatrixMath.subtract(experimental, curr);
                IO.eraseDataFile();
                ssr = sumSquare(residuals);
                lambda_ok = ssr < last_ssr;
                IO.logLambda(lambda, new String[]{"L", "log10k"}, new double[]{l, logk}, ssr, lambda_ok);
                if(delta_c[0] < Math.pow(10, L_DECIMALS) && delta_c[1] < Math.pow(10, LOGK_DECIMALS)){
                    lambda_ok = true;
                    converged = true;
                    return EXECUTED_OK;
                }
            }
            if(!lambda_ok){
                return MAX_LAMBDA_REACHED;
            }
            last_ssr = ssr;
            last_l = l;
            last_logk = logk;
            
            Dummy.run(input, firstl + L_PERTURB, firstlogk, false);
            curr_dl = IO.readData();
            dl = MatrixMath.multiply(MatrixMath.subtract(curr_dl, curr), L_PERTURB);
            IO.eraseDataFile();

            Dummy.run(input, firstl, firstlogk + LOGK_PERTURB, false);
            curr_dk = IO.readData();
            dk = MatrixMath.multiply(MatrixMath.subtract(curr_dk, curr), LOGK_PERTURB);
            IO.eraseDataFile();

            J = MatrixMath.appendColumn(dl, dk);
            JT = MatrixMath.transpose(J);
            JTJ = MatrixMath.multiply(JT, J);
            for(int i = 0; i < DTD.length; i++){
                if(JTJ[i][i] > DTD[i][i]){
                    DTD[i][i] = JTJ[i][i];
                }
            }
            lambda = 0;
            lam_DTD = MatrixMath.multiply(DTD, lambda);
            JTJinv = MatrixMath.invert(MatrixMath.add(JTJ, lam_DTD));
            delta_c = MatrixMath.multiply(MatrixMath.multiply(JTJinv, JT), residuals);
            l = last_l + round(delta_c[0], L_DECIMALS);
            logk = last_logk + round(delta_c[1], LOGK_DECIMALS);
            if(delta_c[0] < Math.pow(10, L_DECIMALS) && delta_c[1] < Math.pow(10, LOGK_DECIMALS)){
                converged = true;
                return EXECUTED_OK;
            }
        }
        
        return MAX_ITERATIONS_REACHED;
    }
    
    public static double findFirstLogK(GridData input, double firstl) throws FileNotFoundException{
        //call the run method to produce a current when the electrode is at distance z and GridData.getCentre() relative to the reactive feature.
        Dummy.run(input, firstl, 1.0, true);
        double[] curr = IO.readData();
        IO.eraseDataFile();
        double max_derivative = curr[2] - curr[0];
        int max_derivative_index = 1;
        for(int i = 2; i < TEST_LOG_K.length -1; i++){
            double derivative = curr[i+1]-curr[i-1];
            if(derivative > max_derivative){
                max_derivative = derivative;
                max_derivative_index = i;
            }
        }
        return TEST_LOG_K[max_derivative_index];
    }
    
    public static double nextLambda(double current_lambda){
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
    
    private static double sumSquare(double[] residuals){
        double sum = 0;
        for(double r : residuals){
            sum += r*r;
        }
        return sum;
    }
    
    private static double round(double value, int decimals){
        return Math.rint(value*Math.pow(10, decimals))/Math.pow(10, decimals);
    }
    
    public static final double[] TEST_LOG_K = new double[]{-8, -7, -6, -5, -4, -5, -4, -3, -2, -1, 0, 1};
    
    public static final int MAX_ITERATIONS = 15;
    
    public static final double MAX_LAMBDA = 100.0;
    
    public static final int L_DECIMALS = 2;//in [um]
    
    public static final double L_PERTURB = 0.05;//in [um]
    
    public static final int LOGK_DECIMALS = 3;//[k] = [m/s]
    
    public static final double LOGK_PERTURB = 0.001;//[k] = [m/s]
    
    
    
    public static final int EXECUTED_OK = 0;
    public static final int MAX_ITERATIONS_REACHED = 1;
    public static final int MAX_LAMBDA_REACHED = 2;
}
