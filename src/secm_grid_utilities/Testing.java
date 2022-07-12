/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package secm_grid_utilities;

import java.util.Arrays;
import static secm_grid_utilities.MatrixMath.*;

/**
 *
 * @author Malak
 */
public class Testing {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws MatrixException {
        double[][] A = new double[][]{
            {1, 2, 3},
            {4, 9, 6},
            {9, 8, 7},
            {2, 6, 1}};
        double[][] B = new double[][]{
            {5, 4, 4, 9},
            {5, 6, 8, 8},
            {1, 5, 2, 7}};
        double[] v = new double[]{1,5,9};
        System.out.println("A01: " + A[0][1]);
        System.out.println("A21: " + A[2][1]);
        printmat(A);
        System.out.println("");
        System.out.println("A*B: ");
        printmat(multiply(A,B));
        System.out.println("");
        System.out.println("B*A: ");
        printmat(multiply(B,A));
        System.out.println("");
        System.out.println("A*v: ");
        System.out.println(Arrays.toString(multiply(A,v)));
    }
    
    private static void printmat(double[][] mat){
        for(int r = 0; r < mat.length; r++){
            for(int c = 0; c < mat[0].length; c++){
                System.out.print(mat[r][c] + "\t");
            }
            System.out.print("\n");
        }
    }
    
}
