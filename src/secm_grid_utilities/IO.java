package secm_grid_utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Scanner;

/**
 * Class holding methods for data input and output.
 * @author Nathaniel Leslie
 * Created:  2022-07-12
 * Modified: 2022-07-12
 */
public class IO {
    public static double[] readData() throws FileNotFoundException{
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
    
    public static void eraseDataFile() throws FileNotFoundException{
        File f = new File("data.txt");
        try (PrintWriter pw = new PrintWriter(f)) {
            pw.write("%");
        }
    }
    
    public static void logIteration(int iteration_num, double[] diagonal, double lambda, String[] labels, double[] params, double ssr, boolean accepted) throws FileNotFoundException{
        File f = new File("fit.log");
        try (PrintWriter pw = new PrintWriter(f)) {
            pw.append("\nIteration: " + iteration_num);
            pw.append("\n\t Diagonal: " + diagonal[0]);
            for(int i = 1; i < labels.length; i++){
                pw.append("," + diagonal[i]);
            }
        }
        logLambda(lambda, labels, params, ssr, accepted);
    }
    
    public static void logLambda(double lambda, String[] labels, double[] params, double ssr, boolean accepted) throws FileNotFoundException{
        File f = new File("fit.log");
        try (PrintWriter pw = new PrintWriter(f)) {
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
        }
    }
    
    public static void logInitialGuesses(String fname, String[] labels, double[] params, double ssr) throws FileNotFoundException{
        File f = new File("fit.log");
        try (PrintWriter pw = new PrintWriter(f)) {
            pw.append("\nInitial guesses:");
            pw.append("\n\tFile: " + fname);
            for(int i = 0; i < labels.length; i++){
                pw.append("\n\t\t" + labels[i] + ": " + params[i]);
            }
            pw.append("\n\tSum. Square Residuals: " + ssr);
        }
    }
    
    public static double[] readTrue(String fname){
        
    }
    
    public static int[][] getGrid(String fname){
        
    }
    
    public static int[] getXs(String fname){
        
    }
    
    public static int[] getYs(String fname){
        
    }
    
    public static double[] getPhysicalXs(String fname){
        
    }
    
    public static double[] getPhysicalYs(String fname){
        
    }
}
