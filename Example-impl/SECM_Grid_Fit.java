/*
 * SECM_Grid_Fit.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Scanner;

/** Model exported on Jul 15 2022, 11:07 by COMSOL 6.0.0.405. */
public class SECM_Grid_Fit {

  public static Model run(double l, double logk, boolean testk) {
    Model model = ModelUtil.create("Model");

    model.modelPath("C:\\Users\\Academic COMSOL\\Documents\\Nathaniel-COMSOL\\Raster\\RasterTest\\settings_3");

    model.label("SECM_HighMeshRasterTest_Model.mph");

    model.param().set("a", "UME_Diam / 2", "UME RADIUS");
    model.param().set("Rg", "1.87*a", "Normalized Glass RADIUS");
    model.param().set("UME_Diam", "4.7[um]", "UME DIAMETER");
    model.param().set("X_UME", "Max_XY/2", "X-position of UME");
    model.param().set("Y_UME", "Max_XY/2", "Y-position of UME");
    model.param().set("Z_UME", l + "*a", "Z-position of UME");//<MOD></MOD>
    model.param().set("Max_XY", "a*350", "Maximum x,y bound");
    model.param().set("Max_Z", "Max_XY *1.0", "Maximum z bound");
    model.param().set("L_UME", "Max_Z - Z_UME", "Length of UME");
    model.param().set("Cox_bulk", "0[mM]", "Bulk concentration of oxidized species");
    model.param().set("Cred_bulk", "1[mM]", "Bulk concentration of reduced species");
    model.param().set("D", "4.13E-11[m^2/s]", "Diffusion coefficient");
    model.param().set("k_UME", "1 [m/s]", "Rate constant at UME");
    model.param().set("k_PS", "10^logk [m/s]", "Rate constant at point source");
    model.param()
         .set("MaxMeshRadial_frac", "0.05", "Fraction of the radial features that will be the maximum mesh size");
    model.param().set("sl", "201", "grid side length");
    model.param().set("gridres", "0.1*a", "Grid resolution");
    model.param().set("logk", logk + "", "base-10 logarithm of k/1m/s");//<MOD></MOD>

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 3);

    model.component("comp1").curvedInterior(false);

    model.result().table().create("tbl1", "Table");

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").label("Simulation cell");
    model.component("comp1").geom("geom1").create("blk1", "Block");
    model.component("comp1").geom("geom1").feature("blk1").label("Bath");
    model.component("comp1").geom("geom1").feature("blk1").set("size", new String[]{"Max_XY", "Max_XY", "Max_Z"});
    model.component("comp1").geom("geom1").create("cyl1", "Cylinder");
    model.component("comp1").geom("geom1").feature("cyl1").label("UME");
    model.component("comp1").geom("geom1").feature("cyl1").set("pos", new String[]{"X_UME", "Y_UME", "Z_UME"});
    model.component("comp1").geom("geom1").feature("cyl1").set("r", "a");
    model.component("comp1").geom("geom1").feature("cyl1").set("h", "L_UME");
    model.component("comp1").geom("geom1").create("cyl3", "Cylinder");
    model.component("comp1").geom("geom1").feature("cyl3").label("UME_Ins");
    model.component("comp1").geom("geom1").feature("cyl3").set("pos", new String[]{"X_UME", "Y_UME", "Z_UME"});
    model.component("comp1").geom("geom1").feature("cyl3").set("r", "Rg");
    model.component("comp1").geom("geom1").feature("cyl3").set("h", "L_UME");
    model.component("comp1").geom("geom1").create("blk2", "Block");
    model.component("comp1").geom("geom1").feature("blk2").active(false);
    model.component("comp1").geom("geom1").feature("blk2").label("box00");
    model.component("comp1").geom("geom1").feature("blk2").set("pos", new String[]{"42.25*a", "43.75*a", "-1E-5"});
    model.component("comp1").geom("geom1").feature("blk2").set("size", new String[]{"0.5*a", "0.5*a", "1E-5"});
    model.component("comp1").geom("geom1").create("pol1", "Polygon");
    model.component("comp1").geom("geom1").feature("pol1").label("square");
    model.component("comp1").geom("geom1").feature("pol1").set("type", "closed");
    model.component("comp1").geom("geom1").feature("pol1").set("source", "table");
    model.component("comp1").geom("geom1").feature("pol1")
         .set("table", new String[][]{{"0.5*Max_XY-sl*0.5*gridres", "0.5*Max_XY-sl*0.5*gridres", "0"}, 
         {"0.5*Max_XY-sl*0.5*gridres + gridres", "0.5*Max_XY-sl*0.5*gridres", "0"}, 
         {"0.5*Max_XY-sl*0.5*gridres + gridres", "0.5*Max_XY-sl*0.5*gridres + gridres", "0"}, 
         {"0.5*Max_XY-sl*0.5*gridres", "0.5*Max_XY-sl*0.5*gridres + gridres", "0"}});
    model.component("comp1").geom("geom1").create("arr1", "Array");
    model.component("comp1").geom("geom1").feature("arr1").set("fullsize", new String[]{"sl", "sl", "1"});
    model.component("comp1").geom("geom1").feature("arr1").set("displ", new String[]{"gridres", "gridres", "0"});
    model.component("comp1").geom("geom1").feature("arr1").selection("input").set("pol1");
    model.component("comp1").geom("geom1").run();
    model.component("comp1").geom("geom1").run("fin");

    model.view().create("view2", 2);
    model.view().create("view3", 2);

    model.component("comp1").physics().create("tds", "DilutedSpecies", "geom1");
    model.component("comp1").physics("tds").field("concentration").field("Cox");
    model.component("comp1").physics("tds").field("concentration").component(new String[]{"Cox", "Cred"});
    model.component("comp1").physics("tds").selection().set(1);
    model.component("comp1").physics("tds").create("fl1", "FluxBoundary", 2);
    model.component("comp1").physics("tds").feature("fl1").selection().set(18303);
    model.component("comp1").physics("tds").create("fl2", "FluxBoundary", 2);
	
	//<MOD>
	int[] boundaries = getPixels(sample_xs[0], sample_ys[0]);
	if(testk){
		int[] coords = getCentre();
		boundaries = getPixels(coords[0], coords[1]);
		model.param().set("logk", TEST_LOG_K[0] + "", "base-10 logarithm of k/1m/s");
		System.out.println("Simulating k-i curve...");
		System.out.print(pBar(1, TEST_LOG_K.length + 1));
	}
	else{
		System.out.println("Simulating image...");
		System.out.print(pBar(1, sample_xs.length + 1));
	}
    model.component("comp1").physics("tds").feature("fl2").selection().set(boundaries);
	//</MOD>
	
    return model;
  }

  public static Model run2(Model model, boolean testk) {
    model.component("comp1").physics("tds").create("conc1", "Concentration", 2);
    model.component("comp1").physics("tds").feature("conc1").selection().set(1, 2, 4, 5, 40419);

    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");
    model.component("comp1").mesh("mesh1").feature("ftet1").selection().geom("geom1", 3);
    model.component("comp1").mesh("mesh1").feature("ftet1").selection().set(1);
    model.component("comp1").mesh("mesh1").feature("ftet1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").selection().geom("geom1", 2);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").selection().set(18303);

    model.component("comp1").probe().create("bnd1", "Boundary");
    model.component("comp1").probe("bnd1").selection().set(18303);

    model.result().table("tbl1").label("Probe Table 1");

    model.component("comp1").view("view1").set("transparency", true);
    model.view("view2").axis().set("xmin", 397.5);
    model.view("view2").axis().set("xmax", 865);
    model.view("view2").axis().set("ymin", 494.375);
    model.view("view2").axis().set("ymax", 755.625);
    model.view("view2").axis().set("viewscaletype", "automatic");
    model.view("view3").axis().set("xmin", 41.69999694824219);
    model.view("view3").axis().set("xmax", 59.30000305175781);
    model.view("view3").axis().set("ymin", 42.8861083984375);
    model.view("view3").axis().set("ymax", 57.1138916015625);
    model.view("view3").axis().set("viewscaletype", "automatic");

    model.component("comp1").physics("tds").prop("EquationForm").set("form", "Stationary");
    model.component("comp1").physics("tds").prop("TransportMechanism").set("Convection", false);
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_Cox", new String[][]{{"D"}, {"0"}, {"0"}, {"0"}, {"D"}, {"0"}, {"0"}, {"0"}, {"D"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_Cred", new String[][]{{"D"}, {"0"}, {"0"}, {"0"}, {"D"}, {"0"}, {"0"}, {"0"}, {"D"}});
    model.component("comp1").physics("tds").feature("init1")
         .set("initc", new String[][]{{"Cox_bulk"}, {"Cred_bulk"}});
    model.component("comp1").physics("tds").feature("fl1").set("species", new int[][]{{1}, {1}});
    model.component("comp1").physics("tds").feature("fl1")
         .set("J0", new String[][]{{"+k_UME*Cred"}, {"-k_UME*Cred"}});
    model.component("comp1").physics("tds").feature("fl1").label("Flux_UME");
    model.component("comp1").physics("tds").feature("fl2").set("species", new int[][]{{1}, {1}});
    model.component("comp1").physics("tds").feature("fl2").set("J0", new String[][]{{"-k_PS*Cox"}, {"+k_PS*Cox"}});
    model.component("comp1").physics("tds").feature("fl2").label("Flux_PS");
    model.component("comp1").physics("tds").feature("conc1").set("species", new int[][]{{1}, {1}});
    model.component("comp1").physics("tds").feature("conc1").set("c0", new String[][]{{"Cox_bulk"}, {"Cred_bulk"}});
    model.component("comp1").physics("tds").feature("conc1").label("Bulk");

    model.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("size").set("hmax", "80[um]");
    model.component("comp1").mesh("mesh1").feature("size").set("hmin", "10[um]");
    model.component("comp1").mesh("mesh1").feature("size").set("hgrad", 1.14);
    model.component("comp1").mesh("mesh1").feature("ftet1").set("optcurved", false);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").label("UME");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1")
         .set("hmax", "a*MaxMeshRadial_frac*0.18");
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("hmin", 1.8E-5);
    model.component("comp1").mesh("mesh1").feature("ftet1").feature("size1").set("hminactive", false);
    model.component("comp1").mesh("mesh1").run();

    model.component("comp1").probe("bnd1").label("UME_FLUX");
    model.component("comp1").probe("bnd1").set("type", "integral");
    model.component("comp1").probe("bnd1").set("probename", "I_UME");
    model.component("comp1").probe("bnd1").set("expr", "tds.ntflux_Cred*F_const");
    model.component("comp1").probe("bnd1").set("unit", "A");
    model.component("comp1").probe("bnd1").set("descractive", true);
    model.component("comp1").probe("bnd1").set("descr", "Current at UME");
    model.component("comp1").probe("bnd1").set("table", "tbl1");
    model.component("comp1").probe("bnd1").set("window", "window1");

    model.study().create("std1");
    model.study("std1").create("stat", "Stationary");

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").attach("std1");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("i1", "Iterative");
    model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").create("sl1", "SORLine");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature().remove("fcDef");

    model.result().dataset().create("dset3", "Solution");
    model.result().dataset().create("int1", "Integral");
    model.result().dataset("dset3").set("probetag", "bnd1");
    model.result().dataset("int1").set("probetag", "bnd1");
    model.result().dataset("int1").set("data", "dset3");
    model.result().dataset("int1").selection().geom("geom1", 2);
    model.result().dataset("int1").selection().set(18303);
    model.result().numerical().create("pev1", "EvalPoint");
    model.result().numerical("pev1").set("probetag", "bnd1");
    model.result().create("pg6", "PlotGroup2D");
    model.result().create("pg7", "PlotGroup2D");
    model.result().create("pg3", "PlotGroup1D");
    model.result("pg6").create("tbls1", "TableSurface");
    model.result("pg6").create("surf1", "Surface");
    model.result("pg6").feature("surf1").set("expr", "comp1.Cox");
    model.result("pg7").create("tbls1", "TableSurface");
    model.result("pg3").set("probetag", "window1");
    model.result("pg3").create("tblp1", "Table");
    model.result("pg3").feature("tblp1").set("probetag", "bnd1");
    model.result().export().create("tbl1", "Table");

    model.component("comp1").probe("bnd1").genResult(null);

    model.result("pg8").tag("pg3");

    model.sol("sol1").attach("std1");
    model.sol("sol1").feature("st1").label("Compile Equations: Stationary");
    model.sol("sol1").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol1").feature("s1").label("Stationary Solver 1.1");
    model.sol("sol1").feature("s1").set("stol", "5e-12");
    model.sol("sol1").feature("s1").feature("dDef").label("Direct 1");
    model.sol("sol1").feature("s1").feature("dDef").set("thresh", 0.1);
    model.sol("sol1").feature("s1").feature("aDef").label("Advanced 1");
    model.sol("sol1").feature("s1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol1").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol1").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol1").feature("s1").feature("i1").label("Iterative 1.1");
    model.sol("sol1").feature("s1").feature("i1").set("nlinnormuse", true);
    model.sol("sol1").feature("s1").feature("i1").set("maxlinit", 400);
    model.sol("sol1").feature("s1").feature("i1").set("rhob", 40);
    model.sol("sol1").feature("s1").feature("i1").feature("ilDef").label("Incomplete LU 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").label("Multigrid 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").label("Presmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("linerelax", 0.2);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("pr").feature("sl1").set("relax", 0.4);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").label("Postsmoother 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("soDef").label("SOR 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").label("SOR Line 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("linerelax", 0.2);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("seconditer", 2);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("po").feature("sl1").set("relax", 0.4);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").label("Coarse Solver 1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("dDef").label("Direct 2");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("dDef").set("thresh", 0.1);
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1").label("Direct 1.1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("linsolver", "pardiso");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").feature("cs").feature("d1")
         .set("pivotperturb", 1.0E-13);
    model.sol("sol1").runAll();

    model.result().dataset("dset3").label("Probe Solution 3");
    model.result("pg6").set("view", "view2");
    model.result("pg6").feature("tbls1").set("dataformat", "filledtable");
    model.result("pg6").feature("surf1").set("resolution", "normal");
    model.result("pg7").set("view", "view3");
    model.result("pg7").set("xlabel", "X_Norm");
    model.result("pg7").set("ylabel", "Y_Norm");
    model.result("pg7").set("xlabelactive", false);
    model.result("pg7").set("ylabelactive", false);
    model.result("pg7").feature("tbls1").set("dataformat", "filledtable");
    model.result("pg3").label("Probe Plot Group 3");
    model.result("pg3").set("xlabel", "Current at UME (A), UME_FLUX");
    model.result("pg3").set("ylabel", "Current at UME (A), UME_FLUX");
    model.result("pg3").set("xlabelactive", false);
    model.result("pg3").set("ylabelactive", false);
    model.result().export("tbl1").set("filename", "data.txt");
    model.result().export("tbl1").set("ifexists", "append");
	model.result().export("tbl1").run();
	
	//<MOD>
	if(!testk){
		for(int i = 1; i < sample_xs.length; i++){
			System.out.print(pBar(i + 1, sample_xs.length + 1));
			model.component("comp1").physics("tds").feature("fl2").selection().set(getPixels(sample_xs[i], sample_ys[i]));
			model.sol("sol1").runAll();
			model.result().export("tbl1").run();
		}
		System.out.print(pBar(1, 1) + "\n");
	}
	else{
	  for(int i = 1; i < TEST_LOG_K.length; i++){
		System.out.print(pBar(i + 1, TEST_LOG_K.length + 1));
		model.param().set("logk", TEST_LOG_K[i] + "", "base-10 logarithm of k/1m/s");
		model.sol("sol1").runAll();
		model.result().export("tbl1").run();
	  }
	}
	//</MOD>
	
    return model;
  }

  public static void main(String[] args) {
	  System.out.println("start");
    try{
//		eraseDataFile();
//		String filename = "Fitfile.csv";
	
//		double L = 1.0;//Normalized tip-to-substrate distance
//		readSECMInfo(filename);
//		System.out.println("Grid data read");
//		double k = -5;//findFirstLogK(L);
//		System.out.println("initial log10k: " + k);
//		int result = runFit(filename, L, k, true);
//		logEndCondition(result);
		
		
		eraseDataFile();
		String filename = "sm_original.csv";
		String newfile = "true-6.csv";
		Model model = run(1.0, -6, false);
        run2(model, false);
		true_image = readData();
		writeSECMInfo(filename, newfile);
		eraseDataFile();
		
		newfile = "true-5.csv";
		model = run(1.0, -5, false);
        run2(model, false);
		true_image = readData();
		writeSECMInfo(filename, newfile);
		eraseDataFile();
		
		newfile = "true-4.csv";
		model = run(1.0, -4, false);
        run2(model, false);
		true_image = readData();
		writeSECMInfo(filename, newfile);
		eraseDataFile();
		
		newfile = "true-3.csv";
		model = run(1.0, -3, false);
        run2(model, false);
		true_image = readData();
		writeSECMInfo(filename, newfile);
		eraseDataFile();
		
		newfile = "true0.csv";
		model = run(1.0, 0, false);
        run2(model, false);
		true_image = readData();
		writeSECMInfo(filename, newfile);
		eraseDataFile();
		
		newfile = "true-3554.csv";
		model = run(1.0, -3.554, false);
        run2(model, false);
		true_image = readData();
		writeSECMInfo(filename, newfile);
		eraseDataFile();
		
		newfile = "true-5569.csv";
		model = run(1.0, -5.569, false);
        run2(model, false);
		true_image = readData();
		writeSECMInfo(filename, newfile);
		eraseDataFile();
		
		newfile = "true-5119.csv";
		model = run(1.0, -5.119, false);
        run2(model, false);
		true_image = readData();
		writeSECMInfo(filename, newfile);
		eraseDataFile();
	}
	catch(Exception e){
		e.printStackTrace();
	}
  }
  
  public static String pBar(int prog, int total){
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
  
  public static int[] fullsurface(){
    int[] sel = new int[40402]; 
	sel[0] = 3;
	for(int i = 6; i <= 16487; i++){
		sel[i-5] = i;
	}
	for(int i = 16492; i <= 18300; i++){
		sel[i-9] = i;
	}
	for(int i = 18305; i <= 20314; i++){
		sel[i-13] = i;
	}
	for(int i = 20319; i <= 40418; i++){
		sel[i-17] = i;
	}
	return sel;
  }

/*
    Fitting methods and fields
    */
    
    public static int runFit(String filename, double firstl, double firstlogk, boolean verbose) throws FileNotFoundException, IOException, NumberFormatException{
        double[] experimental = true_image;
        double lambda = 0;
        double ssr = 0;
        double last_ssr;
        double l = firstl;
        double logk = firstlogk;
        double last_l;
        double last_logk;
        //first iteration
        Model model = run(firstl, firstlogk, false);
        run2(model, false);
        double[] curr = readData();
        double[] residuals = subtract(experimental, curr);
        eraseDataFile();
        
        model = run(firstl + L_PERTURB, firstlogk, false);
        run2(model, false);
        double[] curr_dl = readData();
        double[] dl = multiply(subtract(curr_dl, curr), 1.0/L_PERTURB);
        eraseDataFile();
        
        model = run(firstl, firstlogk + LOGK_PERTURB, false);
        run2(model, false);
        double[] curr_dk = readData();
        double[] dk = multiply(subtract(curr_dk, curr), 1.0/LOGK_PERTURB);
        eraseDataFile();
		
		System.out.println("About to compute SSR");
        last_ssr = ssr;
        ssr = sumSquare(residuals);
        logInitialGuesses(filename, new String[]{"L", "log10k"}, new double[]{l,logk}, ssr);
		
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
		
        System.out.println("About to export i,dl, dlogk, resid");
        if(verbose){
            writeIteration("Iteration_1.txt", physical_xs, physical_ys, curr, dl, dk, residuals);
        }
        last_l = l;
        last_logk = logk;
        l = last_l + round(delta_c[0], L_DECIMALS);
        logk = last_logk + round(delta_c[1], LOGK_DECIMALS);
        boolean converged = false;
        int iterations = 1;
        System.out.println("About to start loop");
		last_ssr = ssr;
        //subsequent iterations
        while(!converged && iterations <= MAX_ITERATIONS){
            iterations ++;
			System.out.println("Iteration " + iterations); 
			System.out.println("Lambda 0");
            //First lambda
            model = run(l, logk, false);
            run2(model, false);
            curr = readData();
            residuals = subtract(experimental, curr);
            eraseDataFile();
            ssr = sumSquare(residuals);
            boolean lambda_ok = ssr < last_ssr;
            logIteration(iterations, new double[]{DTD[0][0], DTD[1][1]}, lambda, new String[]{"L", "log10k"}, new double[]{l, logk}, ssr, lambda_ok);
            
            while(!lambda_ok && lambda <= MAX_LAMBDA){
                lambda = nextLambda(lambda);
				System.out.println("Lambda " + lambda);
                lam_DTD = multiply(DTD, lambda);
                JTJinv = invert(add(JTJ, lam_DTD));
                delta_c = multiply(multiply(JTJinv, JT), residuals);
                l = last_l + round(delta_c[0], L_DECIMALS);
                logk = last_logk + round(delta_c[1], LOGK_DECIMALS);
                model = run(l, logk, false);
                run2(model, false);
                curr = readData();
                residuals = subtract(experimental, curr);
                eraseDataFile();
                ssr = sumSquare(residuals);
                lambda_ok = ssr < last_ssr;
                logLambda(lambda, new String[]{"L", "log10k"}, new double[]{l, logk}, ssr, lambda_ok);
                if(Math.abs(delta_c[0]) < 0.5*Math.pow(10, -L_DECIMALS) && Math.abs(delta_c[1]) < 0.5*Math.pow(10, -LOGK_DECIMALS)){
					System.out.println("DeltaC: " + delta_c[0] + "\t" + delta_c[1]);
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
            
            model = run(firstl + L_PERTURB, firstlogk, false);
            run2(model, false);
            curr_dl = readData();
            dl = multiply(subtract(curr_dl, curr), 1.0/L_PERTURB);
            eraseDataFile();

            model = run(firstl, firstlogk + LOGK_PERTURB, false);
            run2(model, false);
            curr_dk = readData();
            dk = multiply(subtract(curr_dk, curr), 1.0/LOGK_PERTURB);
            eraseDataFile();

			if(verbose){
                writeIteration("Iteration_" + iterations + ".txt", physical_xs, physical_ys, curr, dl, dk, residuals);
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
            if(Math.abs(delta_c[0]) < 0.5*Math.pow(10, -L_DECIMALS) && Math.abs(delta_c[1]) < 0.5*Math.pow(10, -LOGK_DECIMALS)){
                System.out.println("DeltaC: " + delta_c[0] + "\t" + delta_c[1]);
				converged = true;
                return EXECUTED_OK;
            }
        }
        
        return MAX_ITERATIONS_REACHED;
    }
    
    public static double findFirstLogK(double firstl) throws FileNotFoundException{
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
    
    public static final double[] TEST_LOG_K = new double[]{-8, -7, -6, -5, -4, -3, -2, -1, 0, 1};
    
    public static final int MAX_ITERATIONS = 8;
    
    public static final double MAX_LAMBDA = 100.0;
    
    public static final int L_DECIMALS = 2;//in [um]
    
    public static final double L_PERTURB = 0.05;//in [um]
    
    public static final int LOGK_DECIMALS = 3;//[k] = [m/s]
    
    public static final double LOGK_PERTURB = 0.001;//[k] = [m/s]
    
    public static final int EXECUTED_OK = 0;
    public static final int MAX_ITERATIONS_REACHED = 1;
    public static final int MAX_LAMBDA_REACHED = 2;
    
    
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
        PrintWriter pw = new PrintWriter(f);
            pw.print("%");
        pw.close();
    }
    
    public static void logIteration(int iteration_num, double[] diagonal, double lambda, String[] labels, double[] params, double ssr, boolean accepted) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(new FileWriter(f, true));
            pw.append("\nIteration: " + iteration_num);
            pw.append("\n\tDiagonal: " + diagonal[0]);
            for(int i = 1; i < labels.length; i++){
                pw.append("," + diagonal[i]);
            }
        pw.close();
        logLambda(lambda, labels, params, ssr, accepted);
    }
    
    public static void logLambda(double lambda, String[] labels, double[] params, double ssr, boolean accepted) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(new FileWriter(f, true));
            pw.append("\n\tLAMBDA: " + lambda);
            for(int i = 0; i < labels.length; i++){
                pw.append("\n\t\t" + labels[i] + ": " + params[i]);
            }
            pw.append("\n\t\tSum. Square Residuals: " + ssr);
            if(accepted){
                pw.append(" (Accepted)");
            }
            else{
                pw.append(" (Rejected)");
            }
        pw.close();
    }
    
    public static void logInitialGuesses(String fname, String[] labels, double[] params, double ssr) throws FileNotFoundException{
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
    
    public static void logEndCondition(int condition) throws FileNotFoundException, IOException{
        File f = new File(LOGFILE);
        PrintWriter pw = new PrintWriter(new FileWriter(f, true));
            switch(condition){
                case EXECUTED_OK:
                    pw.append("\nPROCESS CONVERGED.");
					break;
                case MAX_ITERATIONS_REACHED:
                    pw.append("\nPROCESS STOPPED PREMATURELY AFTER " + MAX_ITERATIONS + " ITERATIONS.");
					break;
                case MAX_LAMBDA_REACHED:
                    pw.append("\nPROCESS STOPPED DUE TO MAXIMUM LAMBDA BEING REACHED.");
					break;
            }
        pw.close();
    }
    
	public static void writeSECMInfo(String original_filepath, String new_filepath) throws FileNotFoundException, IOException{
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
	
    public static void readSECMInfo(String filepath) throws FileNotFoundException{
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
    
    public static void writeIteration(String fname, double[] x, double[] y, double[] current, double[] dl, double[] dlogk, double[] residual) throws FileNotFoundException{
        File f = new File(fname);
        PrintWriter pw = new PrintWriter(f);
            pw.print("#x [m], y [m], i [A], dL [A], dlog10k [A], residual [A]");
            for(int i = 0; i < x.length; i ++){
                pw.print("\n" + x[i] + "," + y[i] + "," + current[i] + "," + dl[i] + "," + dlogk[i] + "," + residual[i]);
            }
        pw.close();
    }

    public static final String LOGFILE = "fit.log";
    
    
    /*
    Linear algebra-handling methods and fields
     */
    /**
     * Inverts the given matrix
     * @param a the matrix to be inverted
     * @return the inverse of a, a<sup>-1</sup>. a<sup>-1</sup>*a = a*a<sup>-1</sup>=I
     * @throws MatrixException if a is not square or determinant(a)=0
     */
    public static double[][] invert(double[][] a) throws NumberFormatException{
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
    
    public static double[][] multiply(double[][] a, double[][] b) throws NumberFormatException{
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
    
    public static double[] multiply(double[][]a, double[] v) throws NumberFormatException{
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
    
    public static double[][] multiply(double[][] a, double s){
        double[][] mul = new double[a.length][a[0].length];
        for(int r = 0; r < a.length; r++){
            for(int c = 0; c < a[0].length; c++){
                mul[r][c] = s*a[r][c];
            }
        }
        return mul;
    }
    
    public static double[] multiply(double[] v, double s){
        double[] mul = new double[v.length];
        for(int i = 0; i < mul.length; i++){
            mul[i] = v[i]*s;
        }
        return mul;
    }
    
    public static double[][] add(double[][] a, double[][] b) throws NumberFormatException{
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
    
    public static double[][] subtract(double[][] a, double[][] b) throws NumberFormatException{
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
    
    public static double[] add(double[] a, double[] b) throws NumberFormatException{
        if(a.length != b.length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[] sum = new double[a.length];
        for(int r = 0; r < a.length; r++){
            sum[r] = a[r] + b[r];
        }
        return sum;
    }
    
    public static double[] subtract(double[] a, double[] b) throws NumberFormatException{
        if(a.length != b.length){
            throw new NumberFormatException("a and b must be the same size!");
        }
        double[] dif = new double[a.length];
        for(int r = 0; r < a.length; r++){
            dif[r] = a[r] - b[r];
        }
        return dif;
    }
    
    public static double[][] transpose(double[][] a){
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
     * @throws MatrixException 
     */
    public static double[][] appendColumn(double[][] a, double[] v) throws NumberFormatException{
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
     * @throws MatrixException 
     */
    public static double[][] appendColumn(double[] a, double[] v) throws NumberFormatException{
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
    public static double[][] identity(int n){
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
    
    public static double determinant(double[][] a) throws NumberFormatException{
        if(a.length == a[0].length){
            return inner_determinant(a);
        }
        else{
            throw new NumberFormatException("Determinants may only be taken of square matrices.");
        }
    }
    
    private static double inner_determinant(double[][] a){
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
    
    public static double cofactor(double[][] a, int r, int c) throws NumberFormatException{
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
     * @throws MatrixException 
     */
    private static double[][] minor(double[][] a, int r, int c){
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

}
