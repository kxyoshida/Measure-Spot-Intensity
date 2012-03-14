import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.measure.ResultsTable;
import ij.text.TextWindow;
import ij.text.TextPanel;
import java.util.*;
import java.lang.Math;



public class Measure_Spot_Intensity implements PlugInFilter {
    ImagePlus imp;
    int w,h,fr;
    TextWindow mtw;
    static int r=4;
    static double br=0.25;   //25th percentile
    static int avant=40;
    static int apres=40;
    static Boolean fission_at_start=false;
    static Boolean dump_pixel_value=false;
    static double fps=2.0;

    public int setup(String arg, ImagePlus imp) {
	this.imp = imp;
	return DOES_16+DOES_32;
    }
    
    public void run(ImageProcessor ip) {
	w=imp.getWidth()/2;
	h=imp.getHeight();
	
	GenericDialog gd = new GenericDialog("Mask_Spot", IJ.getInstance());
	gd.addNumericField("Target window radius:", r, 0);
	gd.addNumericField("Background rank", br, 2);
	gd.addNumericField("Frame per second (Hz)", fps, 2);
	gd.addNumericField("Pre-fission frames", avant, 0);
	gd.addNumericField("Post-fission frames", apres, 0);
	gd.addCheckbox("Display as Pre-fission + fission",fission_at_start);
	gd.addCheckbox("Dump Pixel Values For Background Estimation",dump_pixel_value);
	gd.showDialog();
	if (gd.wasCanceled()) 
	    return;
	
	r=(int)gd.getNextNumber();
	br=(double)gd.getNextNumber();
	fps=(double)gd.getNextNumber();
	avant=(int)gd.getNextNumber();
	apres=(int)gd.getNextNumber();
	fission_at_start=gd.getNextBoolean();
	dump_pixel_value=gd.getNextBoolean();
	
	int rw=2*r+1;
	int rx=(w-rw)/2;
	
	int nSlices = imp.getStackSize();

	double[] time=new double[nSlices];
	double[] gf=new double[nSlices];
	double[] rf=new double[nSlices];

	double rb=getPercentileInOvalRoiStack(w,0,w,h,br,"red",avant,apres);
	double gb=getPercentileInOvalRoiStack(0,0,w,h,br,"grn",avant,apres);

	if (dump_pixel_value) {
	    dumpPixelValuesInOvalRoiStack(w,0,w,h,"red");
	    dumpPixelValuesInOvalRoiStack(0,0,w,h,"grn");
	}
	
	for (int slice=1;slice<=nSlices;slice++) {
	    IJ.showStatus("Processing "+slice+"/"+nSlices+"");
	    IJ.showProgress(slice, nSlices);
	    ip=imp.getStack().getProcessor(slice).convertToFloat();
	    double gwm=getMeanInOvalRoi(ip, rx, rx, rw, rw);
	    double rwm=getMeanInOvalRoi(ip, w+rx, rx, rw, rw);
	    gf[slice-1]=gwm-gb;
	    rf[slice-1]=rwm-rb;

	    if (fission_at_start) {
		fr=slice-avant;
	    } else {
		fr=slice-(nSlices-apres);
	    }
	    writeResults(fr,gwm,rwm,gb,rb);
	    time[slice-1]=fr/fps;
	}

        Plot plot = new Plot(imp.getShortTitle()+"_plot","frame number","fluorescence", time, gf);
	    if (fission_at_start) {
		plot.setLimits(-avant/fps, 100/fps, 0, 150);
	    } else {
		plot.setLimits(-100/fps, apres/fps, 0, 150);
	    }

        // add a second curve
        plot.setColor(Color.red);
        plot.addPoints(time,rf,Plot.LINE);
        plot.setColor(Color.green);
        plot.draw();
	plot.show();
    }

    double getMeanInOvalRoi(ImageProcessor fip, int rx, int ry, int rw, int rh) {
	OvalRoi or=new OvalRoi(rx,ry,rw,rh);	
	double sum=0;
	int count=0;
	for (int y=ry;y<ry+rh;y++) {
	    for (int x=rx;x<rx+rw;x++) {
		if (or.contains(x,y)) {
		    sum+=fip.getPixelValue(x,y);
		    count++;
		}
	    }
	}
	return (sum/count);
    }

    double getPercentileInOvalRoiStack(int rx, int ry, int rw, int rh, double pp, String label, int avant, int apres) {
	int nSlices=imp.getStackSize();
	OvalRoi or=new OvalRoi(rx,ry,rw,rh);
	int count=0;
	for (int y=ry;y<ry+rh;y++) {
	    for (int x=rx;x<rx+rw;x++) {
		if (or.contains(x,y)) {
		    count++;
		}
	    }
	}

	float[] a=new float[count*(nSlices-avant-apres)];
	float[] cdf=new float[count*(nSlices-avant-apres)];
	int i=0;
	float suma=0;
	for (int slice=avant+1;slice<=nSlices-apres;slice++) {
	    ImageProcessor ip=imp.getStack().getProcessor(slice).convertToFloat();
	    for (int y=ry;y<ry+rh;y++) {
		for (int x=rx;x<rx+rw;x++) {
		    if (or.contains(x,y)) {
			suma+=a[i++]=ip.getPixelValue(x,y);
		    }
		}
	    }
	}

	Arrays.sort(a);                          // sort the data

	// plot cdf
	
	float csum=0;
	for (i=0;i<a.length;i++) {
	    csum+=a[i]/suma;
	    cdf[i]=csum;
	}
        Plot cdfplot = new Plot(imp.getShortTitle()+"_"+label+"_cdf","I","Fn(I)", a, cdf);
	cdfplot.setLimits(a[0], a[a.length-1], 0.0, 1.0);

	// calculate the percentile
	double pt=0;
	int pos=(int)(count*(nSlices-avant-apres)*pp);
	pt=a[pos];

	float[] ptx={a[pos]};
	float[] pty={cdf[pos]};
	cdfplot.setColor(Color.red);
	cdfplot.addPoints(ptx,pty,Plot.X);
	cdfplot.setColor(Color.blue);
        cdfplot.draw();
	cdfplot.show();

	return pt;
    }

    void dumpPixelValuesInOvalRoiStack(int rx, int ry, int rw, int rh, String label) {
	int nSlices=imp.getStackSize();
	OvalRoi or=new OvalRoi(rx,ry,rw,rh);
	String aLine="#PixelValuesInOvalRoiStack";
	TextWindow dpvtw=new TextWindow("PixelValuesInOval_"+rw+"_"+label, "",aLine,300,180);
	for (int slice=1;slice<=nSlices;slice++) {
	    ImageProcessor ip=imp.getStack().getProcessor(slice).convertToFloat();
	    for (int y=ry;y<ry+rh;y++) {
		for (int x=rx;x<rx+rw;x++) {
		    if (or.contains(x,y)) {
			double a=ip.getPixelValue(x,y);
			dpvtw.append(IJ.d2s(a,4));
		    }
		}
	    }
	}
    }

    void writeResults(int index, double gf, double rf, double gb, double rb){
	String title,headings,aLine ;
	aLine = index+"\t"+IJ.d2s(gf,4)+"\t"+IJ.d2s(rf,4);
	if (mtw==null) {
	    title = imp.getShortTitle();
	    //	    headings = "#Index\tFrame\tCenter_X\tCenter_Y\tArea";
	    //	    mtw = new TextWindow(title, headings, aLine, 550, 180);
	    mtw = new TextWindow(title, "#background=\t"+IJ.d2s(gb,4)+"\t"+IJ.d2s(rb,4), aLine, 550, 180);
	} else
	    mtw.append(aLine);
    }

}

