/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package secm_grid_utilities;

import java.util.LinkedList;

/**
 * Holds relevant information for the switch grid.
 * @author Nathaniel Leslie
 * Created:  2022-07-13
 * Modified: 2022-07-13
 */
public class GridData {
    
    public GridData(double[] trueimage, double[] physicalxs, double[] physicalys, int[] samplexs, int[] sampleys, int minx, int miny, int[][] griddata){
        true_image = trueimage;
        physical_xs = physicalxs;
        physical_ys = physicalys;
        sample_xs = samplexs;
        sample_ys = sampleys;
        min_x = minx;
        min_y = miny;
        grid = griddata;
    }
    
    private static final int XSIZE = 201;//Size of pixel grid must be odd
    private static final int YSIZE = 201;//Size of pixel grid must be odd
        
    public int[] getPixels(int x, int y){
        //find the indices of grid that correspond to (x,y) 
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
        
        LinkedList<Integer> x_list = new LinkedList<>();
        LinkedList<Integer> y_list = new LinkedList<>();
        
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
            if(bound > 16286){
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
    
    public double[] getPhysicalXs(){
        return physical_xs;
    }
    
    public double[] getPhysicalYs(){
        return physical_ys;
    }
    
    public int[] getSampleXs(){
        return sample_xs;
    }
    
    public int[] getSampleYs(){
        return sample_ys;
    }
    
    public int[] getCentre(){
        
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
    
    public double[] getTrueImage(){
        return true_image;
    }
    
    private final double[] true_image;
    
    private final double[] physical_xs;
    private final double[] physical_ys;
    private final int[] sample_xs;
    private final int[] sample_ys;
    
    private final int[][] grid;
    private final int min_x;
    private final int min_y;
}
