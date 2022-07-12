/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package secm_grid_utilities;

/**
 * Class holding methods for fitting.
 * @author Nathaniel Leslie
 * Created:  2022-07-12
 * Modified: 2022-07-12
 */
public class Fitting {
    public static void runFit(){
        
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
    
    public static int[] getPixels(int[][] grid, int x, int y){
        
    }
    
    public static int bestScore(int[][] grid, int[] xs, int[] ys){
        
    }
    
}
