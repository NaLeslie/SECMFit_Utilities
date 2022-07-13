/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package secm_grid_utilities;

import java.io.FileNotFoundException;
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
    public static void main(String[] args) throws MatrixException, FileNotFoundException {
        String fname = "C:\\Users\\Malak\\Documents\\NetBeansProjects\\SEM-SECM Align\\Ex\\test1.csv";
        GridData gd = IO.readSECMInfo(fname);
        int[] coords = gd.getCentre();
        System.out.println(Arrays.toString(coords));
        System.out.println("");
        gd.getPixels(coords[0], coords[1]);
    }
    
}
