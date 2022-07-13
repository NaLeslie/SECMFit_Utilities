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
    
    public static GridData readSECMInfo(String filepath) throws FileNotFoundException{
        File f = new File(filepath);
        
        int min_x = 0;
        int min_y = 0;
        int max_x = 0;
        int max_y = 0;
        double[] trueimage;
        double[] physicalxs;
        double[] physicalys;
        int[] samplexs;
        int[] sampleys;
        String sep = ",";
        
        try (Scanner s = new Scanner(f)) {
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
            min_x = xstart;
            max_x = xstart;
            min_y = ystart;
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
                    
                    if(x < min_x){
                        min_x = x;
                    }
                    else if(x > max_x){
                        max_x = x;
                    }
                    
                    if(y < min_y){
                        min_y = y;
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
        }
        //second read of the file to get grid[][]
        int[][] grid = new int[max_x - min_x + 1][max_y - min_y + 1];
        try (Scanner s = new Scanner(f)) {
            while(s.hasNextLine()){
                String line = s.nextLine();
                if(!line.startsWith("#")){
                    String[] linesplit = line.split(sep);
                    int x = Integer.parseInt(linesplit[0]);
                    int y = Integer.parseInt(linesplit[1]);
                    int rs = Integer.parseInt(linesplit[2]);
                    
                    grid[x-min_x][y-min_y] = rs;
                    
                }
            }
        }
        
        return new GridData(trueimage, physicalxs, physicalys, samplexs, sampleys, min_x, min_y, grid);
        
    }
}
